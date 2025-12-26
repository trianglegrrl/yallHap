#!/usr/bin/env python3
"""
Download gnomAD v3.1.2 HGDP/1000 Genomes high-coverage chrY VCF.

This dataset contains:
- 4,151 samples with high-coverage whole genome sequencing
- HGDP (Human Genome Diversity Project) + 1000 Genomes samples
- GRCh38 reference
- AD (allelic depth) fields for Bayesian classification

Downloads:
1. gnomAD HGDP/1000G chrY VCF (~9GB compressed)
2. Rebuilds index locally (source index may be outdated)

Usage:
    python scripts/download_gnomad_highcov.py --output-dir data/validation_highcov

Note: This is a large download (~9GB). Use --sample-only to download and extract
a subset of samples for testing.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# gnomAD v3.1.2 HGDP/1000G chrY VCF (GRCh38)
GNOMAD_HGDP_1KG_CHRY = (
    "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/"
    "gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz"
)

# Expected samples with published haplogroup calls (from 1000 Genomes high-coverage)
# These are GBR samples with well-characterized R and I haplogroups
VALIDATION_SAMPLES = {
    "HG00096": "R-CTS3057",
    "HG00097": "R-CTS3087",
    "HG00099": "R-Z278",
    "HG00100": "R-L21",
    "HG00101": "I-BY266",
    "HG00102": "R-Z278",
    "HG00103": "R-CTS3087",
    "HG00105": "R-BY157",
    "HG00106": "R-Z302",
    "HG00107": "R-Z2571",
    "HG00108": "R-Z372",
    "HG00109": "R-Z287",
    "HG00110": "R-L2",
    "HG00111": "R-CTS241",
    "HG00112": "R-Z198",
    "HG00113": "R-Z87",
    "HG00114": "R-Z87",
    "HG00115": "R-CTS241",
    "HG00116": "R-Z16539",
    "HG00117": "I-M253",
    "HG00118": "R-Z302",
    "HG00119": "R-CTS241",
}


def download_with_curl(url: str, output_path: Path, timeout: int = 7200) -> bool:
    """
    Download a file using curl with progress bar.

    Args:
        url: URL to download
        output_path: Path to save the file
        timeout: Request timeout in seconds

    Returns:
        True if successful, False otherwise
    """
    print(f"Downloading {output_path.name}...")
    print(f"  Source: {url}")
    print(f"  This is a large file (~9GB), please be patient...")

    cmd = [
        "curl",
        "-L",
        "--progress-bar",
        "--max-time",
        str(timeout),
        "-o",
        str(output_path),
        url,
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"  Downloaded to {output_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: Download failed with code {e.returncode}")
        return False
    except FileNotFoundError:
        print("  ERROR: curl not found. Please install curl.")
        return False


def index_vcf(vcf_path: Path) -> bool:
    """
    Index a VCF file using tabix.

    Args:
        vcf_path: Path to VCF file

    Returns:
        True if successful, False otherwise
    """
    if not shutil.which("tabix"):
        print("  ERROR: tabix not found. Please install htslib/tabix.")
        return False

    print(f"Indexing {vcf_path.name}...")
    cmd = ["tabix", "-f", "-p", "vcf", str(vcf_path)]

    try:
        subprocess.run(cmd, check=True)
        print(f"  Created {vcf_path}.tbi")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: Indexing failed: {e}")
        return False


def extract_sample_subset(
    vcf_path: Path,
    output_path: Path,
    sample_ids: list[str],
) -> bool:
    """
    Extract subset of samples from VCF using bcftools.

    Args:
        vcf_path: Path to input VCF
        output_path: Path to output VCF
        sample_ids: List of sample IDs to extract

    Returns:
        True if successful, False otherwise
    """
    if not shutil.which("bcftools"):
        print("  ERROR: bcftools not found. Please install bcftools.")
        return False

    # Create sample list file
    samples_file = output_path.parent / "gnomad_sample_list.txt"
    with open(samples_file, "w") as f:
        f.write("\n".join(sample_ids))

    print(f"  Extracting {len(sample_ids)} samples...")

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
        samples_file.unlink(missing_ok=True)

        # Index the output
        return index_vcf(output_path)
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: bcftools failed: {e.stderr.decode()}")
        samples_file.unlink(missing_ok=True)
        return False


def create_ground_truth(output_path: Path, samples: dict[str, str]) -> None:
    """
    Create ground truth TSV file.

    Args:
        output_path: Path to save ground truth TSV
        samples: Dict of sample_id -> haplogroup
    """
    print(f"Creating ground truth file...")

    with open(output_path, "w") as f:
        f.write("sample_id\thaplogroup\tsource\n")
        for sample_id, haplogroup in sorted(samples.items()):
            f.write(f"{sample_id}\t{haplogroup}\t1000G_highcov\n")

    print(f"  Created {output_path} with {len(samples)} samples")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download gnomAD v3.1.2 HGDP/1000G high-coverage chrY VCF"
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("data/validation_highcov"),
        help="Output directory [default: data/validation_highcov]",
    )
    parser.add_argument(
        "--sample-only",
        action="store_true",
        help="Only download/extract validation sample subset (faster for testing)",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip download, only create ground truth and index",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=7200,
        help="Download timeout in seconds [default: 7200 (2 hours)]",
    )

    args = parser.parse_args()

    # Create output directories
    vcf_dir = args.output_dir / "vcf"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    # Full VCF path
    full_vcf = vcf_dir / "gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz"
    subset_vcf = vcf_dir / "gnomad_validation_subset.vcf.gz"
    ground_truth_path = args.output_dir / "gnomad_highcov_ground_truth.tsv"

    # Create ground truth file
    create_ground_truth(ground_truth_path, VALIDATION_SAMPLES)

    if not args.skip_download:
        if args.sample_only:
            # For sample-only mode, we still need the full file first
            # unless it already exists
            if not full_vcf.exists():
                print("\nNote: Full VCF download required even for subset extraction.")
                print("Consider using --skip-download if you already have the file.\n")
                if not download_with_curl(GNOMAD_HGDP_1KG_CHRY, full_vcf, args.timeout):
                    return 1

            # Index the full VCF if needed
            if not full_vcf.with_suffix(".bgz.tbi").exists():
                if not index_vcf(full_vcf):
                    return 1

            # Extract subset
            if not extract_sample_subset(
                full_vcf, subset_vcf, list(VALIDATION_SAMPLES.keys())
            ):
                return 1
        else:
            # Download full VCF
            if not full_vcf.exists():
                if not download_with_curl(GNOMAD_HGDP_1KG_CHRY, full_vcf, args.timeout):
                    return 1

            # Index the full VCF
            if not index_vcf(full_vcf):
                return 1

    print("\n=== Download Complete ===")
    print(f"Output directory: {args.output_dir}")
    print(f"\nFiles created:")
    print(f"  Ground truth: {ground_truth_path}")
    if full_vcf.exists():
        print(f"  Full VCF: {full_vcf}")
    if subset_vcf.exists():
        print(f"  Subset VCF: {subset_vcf}")

    print("\nNext steps:")
    print("  1. Run validation: python scripts/run_benchmarks.py")
    print("  2. Or quick test:")
    print("     yallhap classify --vcf data/validation_highcov/vcf/gnomad_*.vcf.* \\")
    print("       --sample HG00096 --bayesian")

    return 0


if __name__ == "__main__":
    sys.exit(main())

