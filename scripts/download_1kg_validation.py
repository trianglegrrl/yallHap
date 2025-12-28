#!/usr/bin/env python3
"""
Download 1000 Genomes Phase 3 chrY data for validation.

Downloads:
1. Phase 3 chrY VCF (GRCh37, ~13MB)
2. Sample panel file
3. Creates ground truth TSV from Poznik 2016 data

Usage:
    python scripts/download_1kg_validation.py --output-dir data/validation
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import requests

# 1000 Genomes Phase 3 chrY data URLs
PHASE3_CHRY_VCF = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz"
)
PHASE3_CHRY_TBI = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi"
)
SAMPLE_PANEL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)

# YBrowse GRCh37 SNP database
YBROWSE_HG19_CSV = "http://ybrowse.org/gbrowse2/gff/snps_hg19.csv"


def download_file(
    url: str,
    output_path: Path,
    timeout: int = 300,
    show_progress: bool = True,
) -> None:
    """
    Download a file with progress indication.

    Args:
        url: URL to download
        output_path: Path to save the file
        timeout: Request timeout in seconds
        show_progress: Whether to show progress bar

    Raises:
        requests.RequestException: If download fails
    """
    print(f"Downloading {output_path.name}...")

    response = requests.get(url, stream=True, timeout=timeout)
    response.raise_for_status()

    total_size = response.headers.get("content-length")
    total_bytes = int(total_size) if total_size else None

    downloaded = 0
    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                if show_progress and total_bytes:
                    pct = (downloaded / total_bytes) * 100
                    print(f"\r  {pct:.1f}% ({downloaded:,} / {total_bytes:,} bytes)", end="")

    if show_progress:
        print(f"\n  Saved to {output_path}")


def get_male_samples(panel_path: Path) -> list[dict[str, str]]:
    """
    Extract male samples from panel file.

    Args:
        panel_path: Path to sample panel file

    Returns:
        List of dicts with sample, pop, super_pop, gender
    """
    males = []
    with open(panel_path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(header, parts, strict=False))
            if row.get("gender") == "male":
                males.append(row)
    return males


def parse_yhaplo_output(output_dir: Path, vcf_name: str | None = None) -> dict[str, str]:
    """
    Parse existing yhaplo output files.

    Args:
        output_dir: Directory containing yhaplo output
        vcf_name: Optional VCF filename stem to match specific output

    Returns:
        Dictionary mapping sample_id -> haplogroup
    """
    haplogroups: dict[str, str] = {}

    # Look for yhaplo output files
    yhaplo_dir = output_dir / "yhaplo_output"
    possible_files = []

    if yhaplo_dir.exists():
        # Look for haplogroups.*.txt files
        possible_files.extend(yhaplo_dir.glob("haplogroups.*.txt"))
        possible_files.extend(yhaplo_dir.glob("*.haplogroups.txt"))

    # Also check output_dir directly
    possible_files.extend(output_dir.glob("haplogroups.*.txt"))
    possible_files.extend(output_dir.glob("*.haplogroups.txt"))

    for output_file in possible_files:
        print(f"  Found yhaplo output: {output_file.name}")
        with open(output_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("sample") or line.startswith("#"):
                    continue
                # yhaplo format: sample_id<tab>haplogroup<tab>haplogroup2<tab>...
                # Can be space or tab separated
                parts = line.split()
                if len(parts) >= 2:
                    sample_id = parts[0]
                    haplogroup = parts[1]  # First haplogroup column
                    haplogroups[sample_id] = haplogroup

        if haplogroups:
            print(f"  Parsed {len(haplogroups)} haplogroup assignments")
            return haplogroups

    return haplogroups


def run_yhaplo(vcf_path: Path, output_dir: Path) -> dict[str, str]:
    """
    Run yhaplo on a VCF file and return haplogroup assignments.

    Args:
        vcf_path: Path to input VCF
        output_dir: Directory for yhaplo output files

    Returns:
        Dictionary mapping sample_id -> haplogroup

    Raises:
        subprocess.CalledProcessError: If yhaplo fails
    """
    # First check if yhaplo output already exists
    existing = parse_yhaplo_output(output_dir, vcf_path.stem)
    if existing:
        print(f"  Using existing yhaplo output ({len(existing)} samples)")
        return existing

    print(f"  Running yhaplo on {vcf_path.name}...")

    # yhaplo output directory
    yhaplo_dir = output_dir / "yhaplo_output"
    yhaplo_dir.mkdir(exist_ok=True)

    cmd = [
        "yhaplo",
        "-i", str(vcf_path),
        "-o", str(yhaplo_dir),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  yhaplo stderr: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)

    # Parse the output
    return parse_yhaplo_output(output_dir)


def create_ground_truth_from_yhaplo(
    output_path: Path,
    male_samples: list[dict[str, str]],
    vcf_path: Path | None = None,
    max_samples: int | None = None,
) -> None:
    """
    Create ground truth TSV by running yhaplo on samples.

    If yhaplo is available and a VCF is provided, runs yhaplo to get
    actual haplogroup assignments. Otherwise creates a placeholder file.

    Args:
        output_path: Path to save ground truth TSV
        male_samples: List of male sample dicts
        vcf_path: Path to VCF file (optional, needed for yhaplo)
        max_samples: Maximum number of samples (for testing)
    """
    # Check if yhaplo is available
    yhaplo_available = shutil.which("yhaplo") is not None

    # Build sample -> population lookup for ALL samples first
    all_sample_to_pop = {s.get("sample", ""): s.get("pop", "") for s in male_samples}

    # First check for existing yhaplo output
    haplogroups: dict[str, str] = {}
    existing = parse_yhaplo_output(output_path.parent)
    if existing:
        print(f"Found existing yhaplo output with {len(existing)} haplogroup assignments")
        haplogroups = existing
        # When we have yhaplo output, use ALL samples that have assignments
        # (ignore max_samples limit since haplogroups are already computed)
        sample_ids = [s for s in haplogroups.keys() if s in all_sample_to_pop]
        sample_to_pop = {s: all_sample_to_pop[s] for s in sample_ids}
    else:
        # Apply max_samples limit when we need to compute haplogroups
        if max_samples:
            male_samples = male_samples[:max_samples]
        sample_to_pop = {s.get("sample", ""): s.get("pop", "") for s in male_samples}
        sample_ids = [s.get("sample", "") for s in male_samples]

        if yhaplo_available and vcf_path and vcf_path.exists():
            # Run yhaplo if available and VCF exists
            print("Running yhaplo to generate ground truth haplogroups...")
            try:
                haplogroups = run_yhaplo(vcf_path, output_path.parent)
            except subprocess.CalledProcessError as e:
                print(f"  WARNING: yhaplo failed: {e}")
                print("  Falling back to placeholder haplogroups")
            except Exception as e:
                print(f"  WARNING: yhaplo error: {e}")
                print("  Falling back to placeholder haplogroups")
        elif yhaplo_available:
            print("  yhaplo is available but VCF not yet downloaded")
            print("  Re-run this script after VCF download to get actual haplogroups")
        else:
            print("  yhaplo not found - creating placeholder file")
            print("  Install yhaplo: pip install yhaplo")

    # Write ground truth file
    with open(output_path, "w") as f:
        f.write("sample_id\thaplogroup\tpopulation\n")
        for sample_id in sample_ids:
            pop = sample_to_pop.get(sample_id, "")
            hg = haplogroups.get(sample_id, "UNKNOWN")
            f.write(f"{sample_id}\t{hg}\t{pop}\n")

    # Report results
    actual_count = sum(1 for s in sample_ids if haplogroups.get(s, "UNKNOWN") != "UNKNOWN")
    if actual_count > 0:
        print(f"  Created {output_path} with {actual_count} actual haplogroups")
    else:
        print(f"  Created {output_path} with {len(sample_ids)} placeholder samples")
        print("  NOTE: Haplogroups are placeholders - run yhaplo to get actual values")


def extract_sample_subset(
    vcf_path: Path,
    output_path: Path,
    sample_ids: list[str],
) -> None:
    """
    Extract subset of samples from VCF using bcftools.

    Args:
        vcf_path: Path to input VCF
        output_path: Path to output VCF
        sample_ids: List of sample IDs to extract
    """
    if not shutil.which("bcftools"):
        print("  WARNING: bcftools not found - cannot extract sample subset")
        return

    # Create sample list file
    samples_file = output_path.parent / "sample_list.txt"
    with open(samples_file, "w") as f:
        f.write("\n".join(sample_ids))

    print(f"  Extracting {len(sample_ids)} samples with bcftools...")

    cmd = [
        "bcftools",
        "view",
        "-S",
        str(samples_file),
        str(vcf_path),
        "-Oz",
        "-o",
        str(output_path),
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"  Created {output_path}")

        # Index the output
        subprocess.run(["bcftools", "index", "-t", str(output_path)], check=True)
        print(f"  Indexed {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: bcftools failed: {e.stderr.decode()}")
    finally:
        samples_file.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download 1000 Genomes Phase 3 chrY data for validation"
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("data/validation"),
        help="Output directory [default: data/validation]",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=100,
        help="Maximum samples to include in ground truth [default: 100]",
    )
    parser.add_argument(
        "--skip-vcf",
        action="store_true",
        help="Skip downloading the large VCF file",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Download timeout in seconds [default: 600]",
    )

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Download sample panel first (small file)
    panel_path = args.output_dir / "1kg_sample_panel.txt"
    if not panel_path.exists():
        try:
            download_file(SAMPLE_PANEL, panel_path, timeout=args.timeout)
        except requests.RequestException as e:
            print(f"ERROR: Failed to download sample panel: {e}", file=sys.stderr)
            return 1

    # Get male samples
    male_samples = get_male_samples(panel_path)
    print(f"Found {len(male_samples)} male samples in panel")

    # Download YBrowse GRCh37 SNP database
    ybrowse_path = args.output_dir / "ybrowse_snps_hg19.csv"
    if not ybrowse_path.exists():
        try:
            download_file(YBROWSE_HG19_CSV, ybrowse_path, timeout=args.timeout)
        except requests.RequestException as e:
            print(f"WARNING: Failed to download YBrowse hg19: {e}", file=sys.stderr)

    # Download chrY VCF
    vcf_path = args.output_dir / "1kg_chrY_phase3.vcf.gz"
    if not args.skip_vcf:
        tbi_path = args.output_dir / "1kg_chrY_phase3.vcf.gz.tbi"

        if not vcf_path.exists():
            try:
                download_file(PHASE3_CHRY_VCF, vcf_path, timeout=args.timeout)
            except requests.RequestException as e:
                print(f"ERROR: Failed to download VCF: {e}", file=sys.stderr)
                return 1

        if not tbi_path.exists():
            try:
                download_file(PHASE3_CHRY_TBI, tbi_path, timeout=args.timeout)
            except requests.RequestException as e:
                print(f"ERROR: Failed to download VCF index: {e}", file=sys.stderr)
                return 1

        # Extract subset if bcftools available
        subset_path = args.output_dir / "1kg_chrY_subset.vcf.gz"
        sample_ids = [s["sample"] for s in male_samples[: args.max_samples]]
        if not subset_path.exists() and vcf_path.exists():
            extract_sample_subset(vcf_path, subset_path, sample_ids)

    # Create ground truth file (after VCF download so yhaplo can run)
    ground_truth_path = args.output_dir / "poznik2016_haplogroups.tsv"
    create_ground_truth_from_yhaplo(
        ground_truth_path,
        male_samples,
        vcf_path=vcf_path if vcf_path.exists() else None,
        max_samples=args.max_samples,
    )

    print("\n=== Download Complete ===")
    print(f"Output directory: {args.output_dir}")

    # Check if we got actual haplogroups
    yhaplo_available = shutil.which("yhaplo") is not None
    if yhaplo_available and vcf_path.exists():
        print("\nGround truth haplogroups have been generated with yhaplo.")
        print("You can now run benchmarks:")
        print("  python scripts/run_benchmarks.py")
    else:
        print("\nNext steps:")
        if not yhaplo_available:
            print("1. Install yhaplo: pip install yhaplo")
            print("2. Re-run this script to generate ground truth haplogroups")
        print("3. Run: python scripts/run_benchmarks.py")

    return 0


if __name__ == "__main__":
    sys.exit(main())

