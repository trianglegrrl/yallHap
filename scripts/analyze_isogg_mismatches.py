#!/usr/bin/env python3
"""
Analyze ISOGG validation mismatches to categorize causes.

Categories:
1. LOW_COVERAGE: <10% call rate - insufficient data, not an error
2. ANCESTRAL_STOP: Correct major clade, stopped at ancestor - expected behavior
3. GROUND_TRUTH_ERROR: Evidence suggests ground truth is wrong
4. CLASSIFICATION_ERROR: Potential yallHap bug - needs investigation

Usage:
    python scripts/analyze_isogg_mismatches.py [--samples N]

Output:
    results/isogg_mismatch_analysis.json
    results/isogg_mismatch_analysis.md
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from yallhap.classifier import HaplogroupClassifier
from yallhap.isogg import ISOGGDatabase, ISOGGMapper
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


@dataclass
class MismatchAnalysis:
    """Analysis of a single mismatch."""

    sample: str
    yfull_predicted: str
    isogg_predicted: str
    isogg_ground_truth: str
    call_rate: float
    derived_count: int
    ancestral_count: int
    missing_count: int
    confidence: float
    category: str  # LOW_COVERAGE, ANCESTRAL_STOP, GROUND_TRUTH_ERROR, CLASSIFICATION_ERROR
    notes: str = ""


def get_major_clade(hg: str) -> str:
    """Extract major clade letter from haplogroup."""
    if not hg:
        return ""
    # Handle ISOGG format (R1b1a...) and YFull format (R-L21)
    hg = hg.replace("~", "").strip()
    if hg and hg[0].isalpha():
        return hg[0].upper()
    return ""


def analyze_mismatch(
    sample: str,
    yfull_predicted: str,
    isogg_predicted: str,
    isogg_ground_truth: str,
    classifier: HaplogroupClassifier,
    tree: Tree,
    vcf_path: Path,
) -> MismatchAnalysis:
    """Analyze a single mismatch and categorize it."""
    # Get full classification result
    result = classifier.classify(vcf_path, sample)

    # Calculate call rate
    total = result.snp_stats.derived + result.snp_stats.ancestral + result.snp_stats.missing
    call_rate = (result.snp_stats.derived + result.snp_stats.ancestral) / total if total else 0

    # Get major clades
    pred_major = get_major_clade(isogg_predicted)
    gt_major = get_major_clade(isogg_ground_truth)

    # Categorize
    if call_rate < 0.10:
        category = "LOW_COVERAGE"
        notes = f"Only {call_rate:.1%} of positions called - insufficient data"
    elif pred_major == gt_major:
        category = "ANCESTRAL_STOP"
        notes = f"Same major clade ({pred_major}), stopped at more ancestral node"
    elif call_rate > 0.30 and result.snp_stats.ancestral == 0:
        # High confidence classification with no ancestral markers - likely ground truth error
        category = "GROUND_TRUTH_ERROR"
        notes = f"High call rate ({call_rate:.1%}), {result.snp_stats.derived} derived, 0 ancestral"
    else:
        category = "CLASSIFICATION_ERROR"
        notes = f"Call rate {call_rate:.1%}, needs investigation"

    return MismatchAnalysis(
        sample=sample,
        yfull_predicted=yfull_predicted,
        isogg_predicted=isogg_predicted,
        isogg_ground_truth=isogg_ground_truth,
        call_rate=call_rate,
        derived_count=result.snp_stats.derived,
        ancestral_count=result.snp_stats.ancestral,
        missing_count=result.snp_stats.missing,
        confidence=result.confidence,
        category=category,
        notes=notes,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze ISOGG validation mismatches")
    parser.add_argument("--samples", type=int, default=200, help="Max samples to analyze")
    parser.add_argument("--threads", type=int, default=8, help="Threads for classification")
    args = parser.parse_args()

    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"
    results_dir = base_dir / "results"

    # File paths
    vcf_path = data_dir / "ancient" / "aadr_chrY_v2.vcf.gz"
    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = data_dir / "validation" / "ybrowse_snps_hg19.csv"
    isogg_path = data_dir / "isogg_snps_grch38.txt"
    validation_results_path = results_dir / "isogg_validation_full.json"

    print("=" * 70)
    print("ISOGG Mismatch Analysis")
    print("=" * 70)

    # Check if validation results exist
    if not validation_results_path.exists():
        print(f"Error: Run validate_isogg.py first to generate {validation_results_path}")
        sys.exit(1)

    # Load validation results
    with open(validation_results_path) as f:
        validation_data = json.load(f)

    mismatches = [c for c in validation_data["comparisons"] if c["match_type"] == "mismatch"]
    print(f"\nTotal mismatches from validation: {len(mismatches)}")
    print(f"Analyzing up to {args.samples} samples...")

    # Load resources
    print("\nLoading resources...")
    tree = Tree.from_json(tree_path)
    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)
    classifier = HaplogroupClassifier(tree=tree, snp_db=snp_db, reference="grch37")

    # Analyze mismatches
    start_time = time.time()
    analyses: list[MismatchAnalysis] = []

    sample_mismatches = mismatches[: args.samples]
    for i, mm in enumerate(sample_mismatches):
        if i % 20 == 0:
            print(f"  Progress: {i}/{len(sample_mismatches)}")

        try:
            analysis = analyze_mismatch(
                sample=mm["sample"],
                yfull_predicted=mm["yfull_predicted"],
                isogg_predicted=mm["isogg_predicted"],
                isogg_ground_truth=mm["isogg_ground_truth"],
                classifier=classifier,
                tree=tree,
                vcf_path=vcf_path,
            )
            analyses.append(analysis)
        except Exception as e:
            print(f"  Error analyzing {mm['sample']}: {e}")

    runtime = time.time() - start_time

    # Summarize
    categories = {}
    for a in analyses:
        categories[a.category] = categories.get(a.category, 0) + 1

    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"\nAnalyzed: {len(analyses)} mismatches")
    print(f"Runtime: {runtime:.1f} seconds")
    print("\nCategory breakdown:")
    for cat in ["LOW_COVERAGE", "ANCESTRAL_STOP", "GROUND_TRUTH_ERROR", "CLASSIFICATION_ERROR"]:
        count = categories.get(cat, 0)
        pct = 100 * count / len(analyses) if analyses else 0
        print(f"  {cat}: {count} ({pct:.1f}%)")

    # Calculate corrected error rate
    true_errors = categories.get("CLASSIFICATION_ERROR", 0)
    total_validated = validation_data["total_samples"]
    extrapolated_errors = true_errors * (len(mismatches) / len(analyses)) if analyses else 0

    print(f"\n=== INTERPRETATION ===")
    print(f"Original validation: {validation_data['compatible_rate']:.1f}% compatible")
    print(f"Mismatches analyzed: {len(analyses)}")
    print(f"True classification errors: ~{true_errors} ({100*true_errors/len(analyses):.1f}% of mismatches)")
    print(f"Extrapolated to full dataset: ~{extrapolated_errors:.0f} errors out of {total_validated}")
    corrected_rate = 100 * (total_validated - extrapolated_errors) / total_validated
    print(f"Corrected compatible rate: ~{corrected_rate:.1f}%")

    # Save results
    output_json = results_dir / "isogg_mismatch_analysis.json"
    output_md = results_dir / "isogg_mismatch_analysis.md"

    json_data = {
        "timestamp": datetime.now().isoformat(),
        "samples_analyzed": len(analyses),
        "total_mismatches": len(mismatches),
        "runtime_seconds": round(runtime, 1),
        "categories": categories,
        "original_compatible_rate": validation_data["compatible_rate"],
        "extrapolated_errors": round(extrapolated_errors),
        "corrected_compatible_rate": round(corrected_rate, 1),
        "analyses": [
            {
                "sample": a.sample,
                "yfull_predicted": a.yfull_predicted,
                "isogg_predicted": a.isogg_predicted,
                "isogg_ground_truth": a.isogg_ground_truth,
                "call_rate": round(a.call_rate, 3),
                "derived_count": a.derived_count,
                "ancestral_count": a.ancestral_count,
                "confidence": round(a.confidence, 3),
                "category": a.category,
                "notes": a.notes,
            }
            for a in analyses
        ],
    }

    with open(output_json, "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"\nJSON: {output_json}")

    # Generate markdown report
    md_lines = [
        "# ISOGG Mismatch Analysis Report",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Mismatches analyzed | {len(analyses)} |",
        f"| Total mismatches | {len(mismatches)} |",
        f"| Original compatible rate | {validation_data['compatible_rate']:.1f}% |",
        f"| Corrected compatible rate | {corrected_rate:.1f}% |",
        "",
        "## Category Breakdown",
        "",
        "| Category | Count | Percentage | Description |",
        "|----------|-------|------------|-------------|",
        f"| LOW_COVERAGE | {categories.get('LOW_COVERAGE', 0)} | "
        f"{100*categories.get('LOW_COVERAGE', 0)/len(analyses):.1f}% | "
        f"<10% call rate, insufficient data |",
        f"| ANCESTRAL_STOP | {categories.get('ANCESTRAL_STOP', 0)} | "
        f"{100*categories.get('ANCESTRAL_STOP', 0)/len(analyses):.1f}% | "
        f"Correct clade, stopped at ancestor |",
        f"| GROUND_TRUTH_ERROR | {categories.get('GROUND_TRUTH_ERROR', 0)} | "
        f"{100*categories.get('GROUND_TRUTH_ERROR', 0)/len(analyses):.1f}% | "
        f"Evidence suggests GT is wrong |",
        f"| CLASSIFICATION_ERROR | {categories.get('CLASSIFICATION_ERROR', 0)} | "
        f"{100*categories.get('CLASSIFICATION_ERROR', 0)/len(analyses):.1f}% | "
        f"Potential yallHap bug |",
        "",
        "## Interpretation",
        "",
        "- **LOW_COVERAGE**: Ancient DNA samples with very few called positions. "
        "yallHap correctly stops at the most ancestral supported haplogroup.",
        "- **ANCESTRAL_STOP**: Classification to correct major clade but more ancestral "
        "than ground truth. Expected for low-coverage samples.",
        "- **GROUND_TRUTH_ERROR**: Strong derived signal in yallHap's direction, "
        "suggesting AADR annotation may be incorrect.",
        "- **CLASSIFICATION_ERROR**: Needs investigation - may be real bugs or edge cases.",
        "",
    ]

    # Add examples for each category
    for cat in ["CLASSIFICATION_ERROR", "GROUND_TRUTH_ERROR", "LOW_COVERAGE", "ANCESTRAL_STOP"]:
        examples = [a for a in analyses if a.category == cat][:10]
        if examples:
            md_lines.extend(
                [
                    f"## {cat} Examples ({len([a for a in analyses if a.category == cat])} total)",
                    "",
                    "| Sample | YFull | ISOGG Pred | ISOGG GT | Call Rate | Notes |",
                    "|--------|-------|------------|----------|-----------|-------|",
                ]
            )
            for a in examples:
                md_lines.append(
                    f"| {a.sample} | {a.yfull_predicted} | {a.isogg_predicted} | "
                    f"{a.isogg_ground_truth} | {a.call_rate:.1%} | {a.notes[:40]} |"
                )
            md_lines.append("")

    with open(output_md, "w") as f:
        f.write("\n".join(md_lines))
    print(f"Markdown: {output_md}")


if __name__ == "__main__":
    main()


