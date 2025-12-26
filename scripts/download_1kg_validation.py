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


def create_ground_truth_from_yhaplo(
    output_path: Path,
    male_samples: list[dict[str, str]],
    max_samples: int | None = None,
) -> None:
    """
    Create ground truth TSV by running yhaplo on samples.

    Note: This requires yhaplo to be installed. For now, we create a
    placeholder file that can be populated with actual yhaplo output.

    Args:
        output_path: Path to save ground truth TSV
        male_samples: List of male sample dicts
        max_samples: Maximum number of samples (for testing)
    """
    print("Creating ground truth placeholder...")

    # Check if yhaplo is available
    yhaplo_available = shutil.which("yhaplo") is not None

    if max_samples:
        male_samples = male_samples[:max_samples]

    with open(output_path, "w") as f:
        f.write("sample_id\thaplogroup\tpopulation\n")

        if yhaplo_available:
            print("  yhaplo is available - you can run it to get ground truth")
        else:
            print("  yhaplo not found - creating placeholder file")
            print("  Install yhaplo: pip install yhaplo")
            print("  Then run: yhaplo --example_text")

        # Write placeholder entries (actual haplogroups to be filled in)
        for sample in male_samples:
            sample_id = sample.get("sample", "")
            pop = sample.get("pop", "")
            # Placeholder - actual haplogroup would come from yhaplo
            f.write(f"{sample_id}\tUNKNOWN\t{pop}\n")

    print(f"  Created {output_path} with {len(male_samples)} samples")
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

    # Create ground truth file
    ground_truth_path = args.output_dir / "poznik2016_haplogroups.tsv"
    create_ground_truth_from_yhaplo(
        ground_truth_path,
        male_samples,
        max_samples=args.max_samples,
    )

    # Download YBrowse GRCh37 SNP database
    ybrowse_path = args.output_dir / "ybrowse_snps_hg19.csv"
    if not ybrowse_path.exists():
        try:
            download_file(YBROWSE_HG19_CSV, ybrowse_path, timeout=args.timeout)
        except requests.RequestException as e:
            print(f"WARNING: Failed to download YBrowse hg19: {e}", file=sys.stderr)

    # Download chrY VCF
    if not args.skip_vcf:
        vcf_path = args.output_dir / "1kg_chrY_phase3.vcf.gz"
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

    print("\n=== Download Complete ===")
    print(f"Output directory: {args.output_dir}")
    print("\nNext steps:")
    print("1. Install yhaplo: pip install yhaplo")
    print("2. Run yhaplo on the VCF to get ground truth haplogroups")
    print("3. Update poznik2016_haplogroups.tsv with actual haplogroup calls")
    print("4. Run: python scripts/validate_1kg.py")

    return 0


if __name__ == "__main__":
    sys.exit(main())

