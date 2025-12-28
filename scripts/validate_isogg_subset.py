#!/usr/bin/env python3
"""
Validate ISOGG output on a subset of well-known haplogroups.

Tests the ISOGG mapping functionality by:
1. Selecting samples with unambiguous major haplogroups
2. Running yallHap with --isogg output
3. Checking that ISOGG names are reasonable for the major clade

Usage:
    python scripts/validate_isogg_subset.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from yallhap.classifier import HaplogroupClassifier
from yallhap.isogg import ISOGGDatabase, ISOGGMapper
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


# Expected ISOGG prefix patterns for major haplogroups
# Format: YFull major clade letter -> list of valid ISOGG prefixes
EXPECTED_ISOGG_PREFIXES = {
    "A": ["A", "A0", "A00", "A0-T"],
    "B": ["B", "B2"],
    "C": ["C", "C1", "C2"],
    "D": ["D", "D1", "D2"],
    "E": ["E", "E1", "E2"],
    "G": ["G", "G1", "G2"],
    "H": ["H", "H1", "H2"],
    "I": ["I", "I1", "I2"],
    "J": ["J", "J1", "J2"],
    "L": ["L", "L1"],
    "N": ["N", "N1"],
    "O": ["O", "O1", "O2"],
    "Q": ["Q", "Q1"],
    "R": ["R", "R1", "R2"],
    "T": ["T", "T1"],
}


def load_ground_truth(tsv_path: Path) -> dict[str, str]:
    """Load ground truth haplogroups from TSV."""
    ground_truth = {}
    with open(tsv_path) as f:
        header = True
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            if header:
                header = False
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sample_id = parts[0]
                haplogroup = parts[1]
                ground_truth[sample_id] = haplogroup
    return ground_truth


def get_major_clade(haplogroup: str) -> str:
    """Extract major clade letter from haplogroup."""
    if not haplogroup:
        return ""
    hg = haplogroup.upper()
    if hg.startswith("ROOT"):
        return "ROOT"
    # Return first letter
    return hg[0] if hg and hg[0].isalpha() else ""


def validate_isogg_mapping(
    yfull_hg: str,
    isogg_hg: str,
    expected_major: str,
) -> tuple[bool, str]:
    """
    Validate that ISOGG mapping is reasonable for the major clade.

    Returns:
        (is_valid, reason)
    """
    if not isogg_hg:
        return False, "No ISOGG mapping"

    # Check if ISOGG starts with expected prefix
    expected_prefixes = EXPECTED_ISOGG_PREFIXES.get(expected_major, [expected_major])

    for prefix in expected_prefixes:
        if isogg_hg.upper().startswith(prefix.upper()):
            return True, f"Valid prefix: {prefix}"

    # Check if it's just the major clade letter
    if isogg_hg.upper() == expected_major:
        return True, "Major clade match"

    return False, f"Expected {expected_prefixes}, got {isogg_hg}"


def main() -> None:
    """Run ISOGG validation on subset."""
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"
    validation_dir = data_dir / "validation"
    results_dir = base_dir / "results"

    # File paths
    vcf_path = validation_dir / "1kg_chrY_phase3.vcf.gz"
    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = validation_dir / "ybrowse_snps_hg19.csv"
    isogg_path = data_dir / "isogg_snps_grch38.txt"
    ground_truth_path = validation_dir / "poznik2016_haplogroups.tsv"
    output_path = results_dir / "isogg_validation.json"

    print("=" * 60)
    print("ISOGG Validation Test")
    print("=" * 60)

    # Check files exist
    for path, name in [
        (vcf_path, "VCF"),
        (tree_path, "Tree"),
        (snp_db_path, "SNP DB"),
        (isogg_path, "ISOGG DB"),
        (ground_truth_path, "Ground truth"),
    ]:
        if not path.exists():
            print(f"Error: {name} not found at {path}")
            sys.exit(1)

    # Load resources
    print("\nLoading resources...")
    tree = Tree.from_json(tree_path)
    print(f"  Tree: {len(tree)} nodes")

    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)
    print(f"  SNP DB: {len(snp_db)} SNPs")

    isogg_db = ISOGGDatabase.from_file(isogg_path)
    print(f"  ISOGG DB: {len(isogg_db)} SNPs")

    isogg_mapper = ISOGGMapper(tree, isogg_db)

    ground_truth = load_ground_truth(ground_truth_path)
    print(f"  Ground truth: {len(ground_truth)} samples")

    # Create classifier
    classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference="grch37",
    )

    # Select samples covering major haplogroups
    # Group ground truth by major clade
    samples_by_clade: dict[str, list[str]] = {}
    for sample, hg in ground_truth.items():
        major = get_major_clade(hg)
        if major and major in EXPECTED_ISOGG_PREFIXES:
            if major not in samples_by_clade:
                samples_by_clade[major] = []
            samples_by_clade[major].append(sample)

    # Select 2-3 samples per major clade
    test_samples = []
    for clade, samples in samples_by_clade.items():
        test_samples.extend(samples[:3])

    print(f"\nSelected {len(test_samples)} samples across {len(samples_by_clade)} major clades")

    # Classify samples
    print("Classifying samples...")
    results = classifier.classify_batch(vcf_path, test_samples, threads=8)

    # Validate ISOGG mappings
    print("\nValidating ISOGG mappings...")

    valid_count = 0
    invalid_count = 0
    validation_details = []

    for result in results:
        gt = ground_truth.get(result.sample, "")
        gt_major = get_major_clade(gt)

        # Get ISOGG mapping
        isogg_hg = isogg_mapper.to_isogg(result.haplogroup)

        # Validate
        is_valid, reason = validate_isogg_mapping(result.haplogroup, isogg_hg, gt_major)

        if is_valid:
            valid_count += 1
        else:
            invalid_count += 1

        validation_details.append({
            "sample": result.sample,
            "ground_truth": gt,
            "yfull": result.haplogroup,
            "isogg": isogg_hg,
            "valid": is_valid,
            "reason": reason,
        })

    # Print results
    print("\n" + "=" * 60)
    print("VALIDATION RESULTS")
    print("=" * 60)

    total = valid_count + invalid_count
    accuracy = valid_count / total * 100 if total > 0 else 0

    print(f"\nTotal samples: {total}")
    print(f"Valid mappings: {valid_count} ({accuracy:.1f}%)")
    print(f"Invalid mappings: {invalid_count}")

    # Show some examples
    print("\nSample mappings by clade:")
    shown_clades: set[str] = set()
    for d in validation_details:
        clade = get_major_clade(d["ground_truth"])
        if clade not in shown_clades and clade:
            shown_clades.add(clade)
            status = "✓" if d["valid"] else "✗"
            print(f"  {clade}: {d['yfull']} -> {d['isogg']} {status}")
            if len(shown_clades) >= 10:
                break

    if invalid_count > 0:
        print("\nInvalid mappings:")
        for d in validation_details:
            if not d["valid"]:
                print(f"  {d['sample']}: {d['yfull']} -> {d['isogg']} ({d['reason']})")
                if len([x for x in validation_details if not x["valid"]]) > 5:
                    break

    # Save results
    results_dir.mkdir(exist_ok=True)
    output_data = {
        "total_samples": total,
        "valid_count": valid_count,
        "invalid_count": invalid_count,
        "accuracy_pct": round(accuracy, 2),
        "details": validation_details[:20],  # First 20 for review
    }

    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()

