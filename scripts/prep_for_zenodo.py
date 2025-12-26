#!/usr/bin/env python3
"""
Prepare validation data for Zenodo upload.

This script bundles all validation datasets into a tarball with checksums
for upload to Zenodo. The resulting archive can be downloaded by users
to reproduce yallHap validation tests.

Usage:
    python scripts/prep_for_zenodo.py

Output:
    - yallhap-validation-v1.tar.gz  (~165 MB)
    - yallhap-validation-v1.md5

Files included:
    - validation/1kg_chrY_phase3.vcf.gz (1000 Genomes VCF)
    - validation/1kg_chrY_phase3.vcf.gz.tbi
    - validation/poznik2016_haplogroups.tsv (ground truth)
    - validation/ybrowse_snps_hg19.csv (GRCh37 SNP positions)
    - ancient/aadr_chrY_v2.vcf.gz (ancient DNA VCF)
    - ancient/aadr_chrY_v2.vcf.gz.tbi
    - ancient/aadr_1240k_ground_truth.tsv (ancient ground truth)
    - liftover/grch38-chm13v2.chain (T2T liftover)
    - liftover/hg19-chm13v2.chain (T2T liftover)
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
import tarfile
from pathlib import Path

# Version for the validation bundle
VERSION = "v1"

# Files to include in the bundle (relative to data/)
VALIDATION_FILES = [
    # 1000 Genomes validation
    "validation/1kg_chrY_phase3.vcf.gz",
    "validation/1kg_chrY_phase3.vcf.gz.tbi",
    "validation/poznik2016_haplogroups.tsv",
    "validation/ybrowse_snps_hg19.csv",
    # Ancient DNA validation
    "ancient/aadr_chrY_v2.vcf.gz",
    "ancient/aadr_chrY_v2.vcf.gz.tbi",
    "ancient/aadr_1240k_ground_truth.tsv",
    # T2T liftover chains
    "liftover/grch38-chm13v2.chain",
    "liftover/hg19-chm13v2.chain",
]


def get_file_size_mb(path: Path) -> float:
    """Get file size in megabytes."""
    return path.stat().st_size / (1024 * 1024)


def compute_md5(path: Path) -> str:
    """Compute MD5 checksum of a file."""
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def main() -> int:
    """Create validation data bundle for Zenodo."""
    # Determine paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    data_dir = project_root / "data"

    output_name = f"yallhap-validation-{VERSION}"
    tarball_path = project_root / f"{output_name}.tar.gz"
    md5_path = project_root / f"{output_name}.md5"

    print("=" * 60)
    print("yallHap Validation Data Bundler")
    print("=" * 60)
    print()

    # Check all required files exist
    print("Checking required files...")
    missing_files = []
    total_size = 0.0

    for rel_path in VALIDATION_FILES:
        full_path = data_dir / rel_path
        if full_path.exists():
            size_mb = get_file_size_mb(full_path)
            total_size += size_mb
            print(f"  ✓ {rel_path} ({size_mb:.1f} MB)")
        else:
            missing_files.append(rel_path)
            print(f"  ✗ {rel_path} (MISSING)")

    print()

    if missing_files:
        print("ERROR: Missing files!")
        print()
        print("To download missing files, run:")
        print()
        if any("1kg" in f or "poznik" in f for f in missing_files):
            print("  python scripts/download_1kg_validation.py")
        if any("aadr" in f or "ancient" in f for f in missing_files):
            print("  python scripts/download_ancient_test_data.py")
        if any("liftover" in f or "chain" in f for f in missing_files):
            print("  python scripts/download_liftover_chains.py")
        print()
        return 1

    print(f"Total size: {total_size:.1f} MB")
    print()

    # Create tarball
    print(f"Creating {tarball_path.name}...")

    with tarfile.open(tarball_path, "w:gz") as tar:
        for rel_path in VALIDATION_FILES:
            full_path = data_dir / rel_path
            # Archive with path relative to data/ so it extracts correctly
            arcname = rel_path
            tar.add(full_path, arcname=arcname)
            print(f"  Added: {rel_path}")

    tarball_size = get_file_size_mb(tarball_path)
    print()
    print(f"Created: {tarball_path.name} ({tarball_size:.1f} MB)")

    # Compute and save MD5
    print()
    print("Computing MD5 checksum...")
    md5_hash = compute_md5(tarball_path)

    with open(md5_path, "w") as f:
        f.write(f"{md5_hash}  {tarball_path.name}\n")

    print(f"MD5: {md5_hash}")
    print(f"Saved: {md5_path.name}")

    # Print upload instructions
    print()
    print("=" * 60)
    print("Zenodo Upload Instructions")
    print("=" * 60)
    print()
    print("1. Go to https://zenodo.org/deposit/new")
    print()
    print("2. Upload these files:")
    print(f"   - {tarball_path.name}")
    print(f"   - {md5_path.name}")
    print()
    print("3. Fill in metadata:")
    print("   Title: yallHap Validation Data")
    print(f"   Version: {VERSION}")
    print("   Description: Validation datasets for yallHap Y-chromosome")
    print("                haplogroup classifier. Includes 1000 Genomes")
    print("                Phase 3 and AADR ancient DNA samples.")
    print("   Keywords: Y-chromosome, haplogroup, validation, 1000 Genomes, AADR")
    print("   License: CC-BY-4.0 (for derived data)")
    print()
    print("4. Set access rights: Open Access")
    print()
    print("5. Publish and note the DOI")
    print()
    print("6. Update VALIDATION_TESTING.md with the DOI:")
    print("   Replace 'XXXXXX' with the Zenodo record number")
    print()
    print("Done!")

    return 0


if __name__ == "__main__":
    sys.exit(main())

