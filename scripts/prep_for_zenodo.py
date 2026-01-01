#!/usr/bin/env python3
"""
Prepare validation data for Zenodo upload.

This script bundles all validation datasets into a tarball with checksums
for upload to Zenodo. The resulting archive can be downloaded by users
to reproduce yallHap validation tests.

Usage:
    python scripts/prep_for_zenodo.py

Output:
    - yallhap-validation-v4.tar.gz  (~650 MB)
    - yallhap-validation-v4.md5

Files included:
    Core reference data:
    - yfull_tree.json (YFull phylogeny, ~14 MB)
    - ybrowse_snps.csv (GRCh38 SNP database, ~400 MB)
    - isogg_snps_grch38.txt (ISOGG SNP database)

    1000 Genomes validation:
    - validation/1kg_chrY_phase3.vcf.gz (~5 MB)
    - validation/1kg_chrY_phase3.vcf.gz.tbi
    - validation/poznik2016_haplogroups.tsv (ground truth)
    - validation/1kg_sample_panel.txt (population info)
    - validation/ybrowse_snps_hg19.csv (GRCh37 SNPs, ~48 MB)

    AADR Ancient DNA validation:
    - ancient/aadr_chrY_v2.vcf.gz (~97 MB)
    - ancient/aadr_chrY_v2.vcf.gz.tbi
    - ancient/aadr_1240k_ground_truth.tsv

    Ancient genomes for tool comparison (BAM + VCF, ~85 MB):
    - ancient_genomes/*.bam, *.bam.bai, *.vcf.gz, *.vcf.gz.tbi
    - 9 samples: I0231, I0443, Kennewick, SB524A, SB524A2, VK287, VK292, VK296, VK582

    gnomAD High-Coverage validation (pre-extracted subset):
    - validation_highcov/vcf/gnomad_1kg_shared_diagnostic.vcf.gz
    - validation_highcov/vcf/gnomad_1kg_shared_diagnostic.vcf.gz.tbi
    - validation_highcov/vcf/diagnostic_positions.tsv
    (Full 9GB VCF not included - download with scripts/download_gnomad_highcov.py)

    T2T liftover chains:
    - liftover/grch38-chm13v2.chain
    - liftover/hg19-chm13v2.chain

    Benchmark results:
    - results/benchmark_results.json (benchmark results JSON)
    - results/*.md (all markdown reports in results/ directory)

    Validation scripts (for reproducibility):
    - scripts/download_*.py (data download helpers)
    - scripts/run_benchmarks.py (main benchmark runner)
    - scripts/validate_1kg.py (1000 Genomes validation)
    - scripts/validate_aadr_stratified.py (stratified ancient DNA validation)

    Documentation:
    - ZENODO_MANIFEST.md (this file manifest with descriptions)
"""

from __future__ import annotations

import hashlib
import sys
import tarfile
from pathlib import Path

# Version for the validation bundle
VERSION = "v4"

# Ancient genome samples for tool comparison
ANCIENT_SAMPLES = [
    "I0231.390k.chrY",
    "I0443.390k.chrY",
    "Kennewick_defaultMap1extr.realign.md.head.rmdup.chrY",
    "SB524A_lib.merged.markdup.chrY",
    "SB524A2_lib.merged.markdup.chrY",
    "VK287.final.chrY",
    "VK292.final.chrY",
    "VK296.final.chrY",
    "VK582.final.chrY",
]

# Files to include in the bundle (relative to data/)
VALIDATION_FILES = [
    # Core reference data
    "yfull_tree.json",
    "ybrowse_snps.csv",
    "isogg_snps_grch38.txt",
    # 1000 Genomes validation
    "validation/1kg_chrY_phase3.vcf.gz",
    "validation/1kg_chrY_phase3.vcf.gz.tbi",
    "validation/poznik2016_haplogroups.tsv",
    "validation/1kg_sample_panel.txt",
    "validation/ybrowse_snps_hg19.csv",
    # Ancient DNA validation (AADR multi-sample VCF)
    "ancient/aadr_chrY_v2.vcf.gz",
    "ancient/aadr_chrY_v2.vcf.gz.tbi",
    "ancient/aadr_1240k_ground_truth.tsv",
    # gnomAD High-Coverage validation (pre-extracted subset, not full 9GB VCF)
    "validation_highcov/vcf/gnomad_1kg_shared_diagnostic.vcf.gz",
    "validation_highcov/vcf/gnomad_1kg_shared_diagnostic.vcf.gz.tbi",
    "validation_highcov/vcf/diagnostic_positions.tsv",
    # T2T liftover chains
    "liftover/grch38-chm13v2.chain",
    "liftover/hg19-chm13v2.chain",
]

# Ancient genome files (BAM + VCF for each sample) - built dynamically
ANCIENT_GENOME_FILES: list[str] = []
for sample in ANCIENT_SAMPLES:
    ANCIENT_GENOME_FILES.extend([
        f"ancient_genomes/{sample}.bam",
        f"ancient_genomes/{sample}.bam.bai",
        f"ancient_genomes/{sample}.vcf.gz",
        f"ancient_genomes/{sample}.vcf.gz.tbi",
    ])

# Results files to include (relative to project root)
# Note: All results/*.md files are included dynamically
# Note: results/validation_temp/ is excluded (temporary tool outputs)
RESULTS_FILES = [
    "results/benchmark_results.json",
    "results/coalescent_prior_comparison.json",
    "results/isogg_mismatch_analysis.json",
    "results/isogg_mismatches_detailed.json",
    "results/isogg_validation_full.json",
    "results/isogg_validation.json",
    "results/power_analysis_ancient.json",
    "results/significance_tests.json",
]

# Validation scripts to include for reproducibility (relative to project root)
# These scripts allow users to reproduce validation from scratch
VALIDATION_SCRIPTS = [
    # Data download scripts
    "scripts/download_1kg_validation.py",
    "scripts/download_ancient_test_data.py",
    "scripts/download_gnomad_highcov.py",
    "scripts/download_liftover_chains.py",
    "scripts/download_yfull_tree.py",
    # Core validation runners
    "scripts/run_benchmarks.py",
    "scripts/validate_1kg.py",
    "scripts/validate_aadr_stratified.py",
    "scripts/gather_validation_and_comparative_data.py",
    # ISOGG validation
    "scripts/validate_isogg.py",
    "scripts/analyze_isogg_mismatches.py",
    "scripts/generate_isogg_mismatch_table.py",
    # Statistical analysis
    "scripts/power_analysis_ancient.py",
    "scripts/calculate_significance.py",
    "scripts/test_coalescent_priors.py",
    # Asset generation
    "scripts/generate_assets.py",
    # Bundle creation (this script)
    "scripts/prep_for_zenodo.py",
]

# Documentation files (relative to project root)
DOCUMENTATION_FILES = [
    "ZENODO_MANIFEST.md",
]

# Optional files (include if present, don't fail if missing)
OPTIONAL_FILES: list[str] = []


def get_file_size_mb(path: Path) -> float:
    """Get file size in megabytes."""
    return path.stat().st_size / (1024 * 1024)


def get_file_size_gb(path: Path) -> float:
    """Get file size in gigabytes."""
    return path.stat().st_size / (1024 * 1024 * 1024)


def compute_md5(path: Path) -> str:
    """Compute MD5 checksum of a file."""
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def format_size(size_mb: float) -> str:
    """Format size as MB or GB as appropriate."""
    if size_mb >= 1024:
        return f"{size_mb / 1024:.2f} GB"
    return f"{size_mb:.1f} MB"


def main() -> int:
    """Create validation data bundle for Zenodo."""
    # Determine paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    data_dir = project_root / "data"

    output_name = f"yallhap-validation-{VERSION}"
    tarball_path = project_root / f"{output_name}.tar.gz"
    md5_path = project_root / f"{output_name}.md5"

    print("=" * 70)
    print("yallHap Validation Data Bundler")
    print("=" * 70)
    print()
    print(f"Version: {VERSION}")
    print(f"Output: {tarball_path.name}")
    print()

    # Check all required files exist
    print("Checking validation data files...")
    missing_files = []
    total_size = 0.0
    files_to_add: list[tuple[Path, str]] = []  # (full_path, archive_name)

    for rel_path in VALIDATION_FILES:
        full_path = data_dir / rel_path
        if full_path.exists():
            size_mb = get_file_size_mb(full_path)
            total_size += size_mb
            size_str = format_size(size_mb)
            print(f"  ✓ {rel_path} ({size_str})")
            files_to_add.append((full_path, rel_path))
        else:
            missing_files.append(rel_path)
            print(f"  ✗ {rel_path} (MISSING)")

    print("\nChecking ancient genome files (for tool comparison)...")
    for rel_path in ANCIENT_GENOME_FILES:
        full_path = data_dir / rel_path
        if full_path.exists():
            size_mb = get_file_size_mb(full_path)
            total_size += size_mb
            size_str = format_size(size_mb)
            print(f"  ✓ {rel_path} ({size_str})")
            files_to_add.append((full_path, rel_path))
        else:
            missing_files.append(rel_path)
            print(f"  ✗ {rel_path} (MISSING)")

    print("\nChecking benchmark results files...")
    for rel_path in RESULTS_FILES:
        full_path = project_root / rel_path
        if full_path.exists():
            size_mb = get_file_size_mb(full_path)
            total_size += size_mb
            print(f"  ✓ {rel_path} ({size_mb:.2f} MB)")
            # Keep the path as-is since it's already in results/
            files_to_add.append((full_path, rel_path))
        else:
            missing_files.append(rel_path)
            print(f"  ✗ {rel_path} (MISSING)")

    print("\nChecking results/*.md files...")
    results_dir = project_root / "results"
    if results_dir.exists():
        for md_file in sorted(results_dir.glob("*.md")):
            size_mb = get_file_size_mb(md_file)
            total_size += size_mb
            print(f"  ✓ results/{md_file.name} ({size_mb:.2f} MB)")
            files_to_add.append((md_file, f"results/{md_file.name}"))
    else:
        print("  - results/ directory not found")

    print("\nChecking validation scripts...")
    for rel_path in VALIDATION_SCRIPTS:
        full_path = project_root / rel_path
        if full_path.exists():
            size_mb = get_file_size_mb(full_path)
            total_size += size_mb
            print(f"  ✓ {rel_path} ({size_mb:.2f} MB)")
            files_to_add.append((full_path, rel_path))
        else:
            missing_files.append(rel_path)
            print(f"  ✗ {rel_path} (MISSING)")

    print("\nChecking documentation files...")
    for rel_path in DOCUMENTATION_FILES:
        full_path = project_root / rel_path
        if full_path.exists():
            size_mb = get_file_size_mb(full_path)
            total_size += size_mb
            print(f"  ✓ {rel_path} ({size_mb:.2f} MB)")
            # Place at root of archive (no subdirectory)
            files_to_add.append((full_path, rel_path))
        else:
            missing_files.append(rel_path)
            print(f"  ✗ {rel_path} (MISSING)")

    print("\nChecking optional files...")
    for rel_path in OPTIONAL_FILES:
        full_path = project_root / rel_path
        if full_path.exists():
            size_mb = get_file_size_mb(full_path)
            total_size += size_mb
            print(f"  ✓ {rel_path} ({size_mb:.2f} MB)")
            files_to_add.append((full_path, f"results/{rel_path}"))
        else:
            print(f"  - {rel_path} (not present, skipping)")

    print()

    if missing_files:
        print("ERROR: Missing required files!")
        print()
        print("To download missing files, run:")
        print()
        if any("yfull" in f or "ybrowse_snps.csv" in f for f in missing_files):
            print("  yallhap download --output-dir data/")
        if any("isogg" in f for f in missing_files):
            print("  # isogg_snps_grch38.txt - download from ISOGG website")
        if any("1kg" in f or "poznik" in f for f in missing_files):
            print("  python scripts/download_1kg_validation.py")
        if any("aadr" in f for f in missing_files):
            print("  python scripts/download_ancient_test_data.py")
        if any("ancient_genomes" in f for f in missing_files):
            print("  # Ancient genome BAMs/VCFs - see scripts/download_ancient_genomes.py")
        if any("gnomad" in f or "validation_highcov" in f or "diagnostic" in f for f in missing_files):
            print("  python scripts/download_gnomad_highcov.py")
            print("  # Then run: python scripts/run_benchmarks.py --threads 16")
            print("  # (creates gnomad_1kg_shared_diagnostic.vcf.gz)")
        if any("liftover" in f or "chain" in f for f in missing_files):
            print("  python scripts/download_liftover_chains.py")
        if any("benchmark" in f or "results/" in f for f in missing_files):
            print("  python scripts/run_benchmarks.py --threads 16")
            print("  python scripts/validate_aadr_stratified.py -o results/aadr_density_stratified.md")
        print()
        return 1

    print(f"Total size: {format_size(total_size)}")
    print()

    # Warn about large files
    if total_size > 2000:  # > 2 GB
        print("⚠️  WARNING: Bundle is larger than expected.")
        print("   Check if any unexpected large files are included.")
        print()

    # Confirm before creating large archive
    if total_size > 1000:  # > 1 GB
        response = input(f"Create {format_size(total_size)} archive? [y/N] ")
        if response.lower() != "y":
            print("Aborted.")
            return 0

    # Create tarball
    print(f"Creating {tarball_path.name}...")
    print("(This may take several minutes for large files)")
    print()

    with tarfile.open(tarball_path, "w:gz") as tar:
        for full_path, arcname in files_to_add:
            size_str = format_size(get_file_size_mb(full_path))
            print(f"  Adding: {arcname} ({size_str})...")
            tar.add(full_path, arcname=arcname)

    tarball_size = get_file_size_mb(tarball_path)
    print()
    print(f"Created: {tarball_path.name} ({format_size(tarball_size)})")

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
    print("=" * 70)
    print("Zenodo Upload Instructions")
    print("=" * 70)
    print()
    print("1. Go to https://zenodo.org/deposit/new")
    print()
    print("2. Upload these files:")
    print(f"   - {tarball_path.name} ({format_size(tarball_size)})")
    print(f"   - {md5_path.name}")
    print()
    print("   Note: ZENODO_MANIFEST.md is included IN the tarball for users")
    print()
    print("3. Fill in metadata:")
    print("   Title: yallHap Validation Data Bundle")
    print(f"   Version: {VERSION}")
    print("   Description:")
    print("     Complete validation dataset bundle for yallHap Y-chromosome")
    print("     haplogroup classifier. Includes:")
    print("     - 1000 Genomes Phase 3 (1,233 samples, GRCh37)")
    print("     - gnomAD HGDP/1KG high-coverage subset (GRCh38)")
    print("     - AADR v54 ancient DNA (7,333 samples, GRCh37)")
    print("     - 9 ancient genome BAM/VCF files for tool comparison")
    print("     - YFull tree, YBrowse SNPs, and ISOGG SNP database")
    print("     - T2T liftover chain files")
    print("     - Benchmark results and validation reports")
    print()
    print("   Keywords: Y-chromosome, haplogroup, validation, 1000 Genomes,")
    print("             AADR, gnomAD, ancient DNA, phylogenetics")
    print()
    print("   License:")
    print("     - yallHap software: PolyForm Noncommercial License 1.0.0")
    print("     - Derived data/results: CC-BY-4.0")
    print("     - Original consortium data: see respective provider terms")
    print()
    print("4. Set access rights: Open Access")
    print()
    print("5. Publish and note the DOI")
    print()
    print("6. Update VALIDATION_TESTING.md with the DOI:")
    print("   Replace 'XXXXXX' with the Zenodo record number")
    print()
    print("Done!")

    # Print file manifest
    print()
    print("=" * 70)
    print("File Manifest")
    print("=" * 70)
    print()
    for full_path, arcname in files_to_add:
        size_str = format_size(get_file_size_mb(full_path))
        print(f"  {arcname}: {size_str}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
