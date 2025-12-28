#!/usr/bin/env python3
"""
Validate yallHap against 1000 Genomes Phase 3 ground truth.

Runs yallHap classification on samples and compares against
Poznik 2016 haplogroup assignments.

Usage:
    python scripts/validate_1kg.py --vcf data/validation/1kg_chrY_subset.vcf.gz
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from yallhap.snps import SNPDatabase
from yallhap.tree import Tree
from yallhap.validation import (
    ValidationRunner,
    compute_metrics,
    load_ground_truth,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate yallHap against 1000 Genomes Phase 3"
    )
    parser.add_argument(
        "--vcf",
        type=Path,
        required=True,
        help="Path to chrY VCF file",
    )
    parser.add_argument(
        "--tree",
        type=Path,
        default=Path("data/yfull_tree.json"),
        help="Path to YFull tree JSON [default: data/yfull_tree.json]",
    )
    parser.add_argument(
        "--snp-db",
        type=Path,
        default=Path("data/validation/ybrowse_snps_hg19.csv"),
        help="Path to YBrowse SNP CSV [default: data/validation/ybrowse_snps_hg19.csv]",
    )
    parser.add_argument(
        "--ground-truth",
        type=Path,
        default=Path("data/validation/poznik2016_haplogroups.tsv"),
        help="Path to ground truth TSV [default: data/validation/poznik2016_haplogroups.tsv]",
    )
    parser.add_argument(
        "--reference",
        type=str,
        choices=["grch37", "grch38", "t2t"],
        default="grch37",
        help="Reference genome [default: grch37]",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Maximum samples to validate [default: all]",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON file for results [default: stdout]",
    )

    args = parser.parse_args()

    # Validate input files exist
    if not args.vcf.exists():
        print(f"ERROR: VCF file not found: {args.vcf}", file=sys.stderr)
        return 1
    if not args.tree.exists():
        print(f"ERROR: Tree file not found: {args.tree}", file=sys.stderr)
        return 1
    if not args.snp_db.exists():
        print(f"ERROR: SNP database not found: {args.snp_db}", file=sys.stderr)
        return 1
    if not args.ground_truth.exists():
        print(f"ERROR: Ground truth file not found: {args.ground_truth}", file=sys.stderr)
        return 1

    # Load resources
    print("Loading YFull tree...", file=sys.stderr)
    tree = Tree.from_json(args.tree)
    print(f"  Loaded {len(tree)} nodes", file=sys.stderr)

    print("Loading SNP database...", file=sys.stderr)
    # Detect format based on file content
    snp_db = SNPDatabase.from_ybrowse_gff_csv(args.snp_db)
    print(f"  Loaded {len(snp_db)} SNPs", file=sys.stderr)

    print("Loading ground truth...", file=sys.stderr)
    ground_truth = load_ground_truth(args.ground_truth)
    print(f"  Loaded {len(ground_truth)} samples", file=sys.stderr)

    # Filter to only samples with valid haplogroups (not UNKNOWN)
    valid_ground_truth = {
        k: v for k, v in ground_truth.items() if v and v != "UNKNOWN"
    }
    if len(valid_ground_truth) < len(ground_truth):
        print(
            f"  Filtered to {len(valid_ground_truth)} samples with valid haplogroups",
            file=sys.stderr,
        )

    if not valid_ground_truth:
        print(
            "ERROR: No valid ground truth haplogroups found. "
            "Run yhaplo first to populate haplogroups.",
            file=sys.stderr,
        )
        return 1

    # Limit samples if requested
    sample_ids = list(valid_ground_truth.keys())
    if args.max_samples:
        sample_ids = sample_ids[: args.max_samples]

    print(f"\nValidating {len(sample_ids)} samples...", file=sys.stderr)

    # Create runner
    runner = ValidationRunner(
        tree=tree,
        snp_db=snp_db,
        ground_truth=valid_ground_truth,
        reference=args.reference,  # type: ignore
    )

    # Run validation
    results = []
    for i, sample_id in enumerate(sample_ids, 1):
        try:
            result = runner.validate_sample(sample_id, args.vcf)
            results.append(result)

            status = "✓" if result.comparison.major_match else "✗"
            print(
                f"  [{i}/{len(sample_ids)}] {sample_id}: "
                f"{result.expected_haplogroup} → {result.called_haplogroup} "
                f"({result.confidence:.2f}) {status}",
                file=sys.stderr,
            )
        except Exception as e:
            print(f"  [{i}/{len(sample_ids)}] {sample_id}: ERROR - {e}", file=sys.stderr)

    # Compute metrics
    metrics = compute_metrics(results)

    print("\n=== Validation Results ===", file=sys.stderr)
    print(f"Total samples: {metrics.total_samples}", file=sys.stderr)
    print(f"Exact match rate: {metrics.exact_match_rate:.1%}", file=sys.stderr)
    print(f"Major match rate: {metrics.major_match_rate:.1%}", file=sys.stderr)
    print(f"Mean confidence: {metrics.mean_confidence:.2f}", file=sys.stderr)

    # Prepare output
    output_data = {
        "metrics": {
            "total_samples": metrics.total_samples,
            "exact_match_rate": metrics.exact_match_rate,
            "major_match_rate": metrics.major_match_rate,
            "mean_confidence": metrics.mean_confidence,
            "mean_depth_difference": metrics.mean_depth_difference,
        },
        "results": [
            {
                "sample_id": r.sample_id,
                "expected": r.expected_haplogroup,
                "called": r.called_haplogroup,
                "confidence": r.confidence,
                "exact_match": r.comparison.exact_match,
                "major_match": r.comparison.major_match,
            }
            for r in results
        ],
    }

    # Write output
    if args.output:
        with open(args.output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults written to {args.output}", file=sys.stderr)
    else:
        print(json.dumps(output_data, indent=2))

    # Return exit code based on success
    if metrics.major_match_rate >= 0.9:
        return 0  # Success: >90% major match
    elif metrics.major_match_rate >= 0.7:
        return 1  # Warning: 70-90% match
    else:
        return 2  # Failure: <70% match


if __name__ == "__main__":
    sys.exit(main())
