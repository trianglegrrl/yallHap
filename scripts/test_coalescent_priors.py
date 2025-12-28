#!/usr/bin/env python3
"""
Test coalescent vs uniform priors for Bayesian haplogroup classification.

Compares accuracy and confidence distributions between:
- Uniform priors: Equal probability for all haplogroups
- Coalescent priors: Weighted by branch length (SNP count)

Usage:
    python scripts/test_coalescent_priors.py [--samples N]
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import pysam

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from yallhap.classifier import HaplogroupClassifier
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


@dataclass
class ClassificationResult:
    """Result of a single classification."""

    sample: str
    ground_truth: str
    predicted_uniform: str
    confidence_uniform: float
    predicted_coalescent: str
    confidence_coalescent: float


def load_ground_truth(tsv_path: Path) -> dict[str, str]:
    """Load ground truth haplogroups from Poznik TSV."""
    ground_truth = {}
    with open(tsv_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sample_id = parts[0]
                haplogroup = parts[1]
                ground_truth[sample_id] = haplogroup
    return ground_truth


def get_major_clade(haplogroup: str) -> str:
    """Extract major clade (first letter) from haplogroup."""
    if not haplogroup:
        return ""
    # Handle special cases
    hg = haplogroup.upper()
    if hg.startswith("ROOT"):
        return "ROOT"
    # First letter or letter combo (e.g., "A0", "CT", "BT")
    for i, c in enumerate(hg):
        if c.isdigit() or c == "-":
            return hg[:i] if i > 0 else hg[0]
    return hg[0] if hg else ""


def classify_samples(
    vcf_path: Path,
    tree: Tree,
    snp_db: SNPDatabase,
    ground_truth: dict[str, str],
    sample_ids: list[str],
    reference: str = "grch37",
    threads: int = 8,
) -> list[ClassificationResult]:
    """
    Classify samples with both uniform and coalescent priors.
    """
    results = []

    # Create classifiers with different priors
    print("  Creating uniform prior classifier...")
    classifier_uniform = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference=reference,  # type: ignore
        bayesian=True,
        prior_type="uniform",
    )

    print("  Creating coalescent prior classifier...")
    classifier_coalescent = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference=reference,  # type: ignore
        bayesian=True,
        prior_type="coalescent",
    )

    # Filter to valid samples
    valid_samples = [s for s in sample_ids if s in ground_truth]

    # Classify with uniform prior
    print(f"  Classifying {len(valid_samples)} samples with uniform prior...")
    uniform_results = classifier_uniform.classify_batch(
        vcf_path, valid_samples, threads=threads
    )
    uniform_map = {r.sample: r for r in uniform_results}

    # Classify with coalescent prior
    print(f"  Classifying {len(valid_samples)} samples with coalescent prior...")
    coalescent_results = classifier_coalescent.classify_batch(
        vcf_path, valid_samples, threads=threads
    )
    coalescent_map = {r.sample: r for r in coalescent_results}

    # Combine results
    for sample_id in valid_samples:
        if sample_id not in uniform_map or sample_id not in coalescent_map:
            continue

        u_result = uniform_map[sample_id]
        c_result = coalescent_map[sample_id]

        results.append(
            ClassificationResult(
                sample=sample_id,
                ground_truth=ground_truth[sample_id],
                predicted_uniform=u_result.haplogroup,
                confidence_uniform=u_result.confidence,
                predicted_coalescent=c_result.haplogroup,
                confidence_coalescent=c_result.confidence,
            )
        )

    return results


def analyze_results(results: list[ClassificationResult]) -> dict:
    """Analyze and compare results from both prior types."""
    uniform_correct = 0
    coalescent_correct = 0
    uniform_confidences = []
    coalescent_confidences = []

    # Track disagreements
    disagreements = []

    for r in results:
        gt_major = get_major_clade(r.ground_truth)
        uniform_major = get_major_clade(r.predicted_uniform)
        coalescent_major = get_major_clade(r.predicted_coalescent)

        uniform_match = gt_major == uniform_major
        coalescent_match = gt_major == coalescent_major

        if uniform_match:
            uniform_correct += 1
        if coalescent_match:
            coalescent_correct += 1

        uniform_confidences.append(r.confidence_uniform)
        coalescent_confidences.append(r.confidence_coalescent)

        # Track where they disagree
        if uniform_major != coalescent_major:
            disagreements.append(
                {
                    "sample": r.sample,
                    "ground_truth": r.ground_truth,
                    "uniform": r.predicted_uniform,
                    "coalescent": r.predicted_coalescent,
                    "uniform_correct": uniform_match,
                    "coalescent_correct": coalescent_match,
                }
            )

    n = len(results)
    uniform_acc = uniform_correct / n * 100 if n > 0 else 0
    coalescent_acc = coalescent_correct / n * 100 if n > 0 else 0
    avg_conf_uniform = sum(uniform_confidences) / n if n > 0 else 0
    avg_conf_coalescent = sum(coalescent_confidences) / n if n > 0 else 0

    return {
        "total_samples": n,
        "uniform_prior": {
            "correct": uniform_correct,
            "accuracy_pct": round(uniform_acc, 2),
            "mean_confidence": round(avg_conf_uniform, 4),
        },
        "coalescent_prior": {
            "correct": coalescent_correct,
            "accuracy_pct": round(coalescent_acc, 2),
            "mean_confidence": round(avg_conf_coalescent, 4),
        },
        "difference_pp": round(coalescent_acc - uniform_acc, 2),
        "disagreements": len(disagreements),
        "disagreement_details": disagreements[:10],  # First 10 for review
    }


def main() -> None:
    """Run coalescent prior comparison test."""
    parser = argparse.ArgumentParser(description="Test coalescent vs uniform priors")
    parser.add_argument(
        "--samples",
        type=int,
        default=200,
        help="Number of samples to test (default: 200)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads for parallel processing (default: 8)",
    )
    args = parser.parse_args()

    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"

    # File paths
    vcf_path = data_dir / "validation" / "1kg_chrY_phase3.vcf.gz"
    tree_path = data_dir / "yfull_tree.json"
    snp_path = data_dir / "validation" / "ybrowse_snps_hg19.csv"
    ground_truth_path = data_dir / "validation" / "poznik2016_haplogroups.tsv"
    output_path = base_dir / "results" / "coalescent_prior_comparison.json"

    print("=" * 60)
    print("Coalescent Prior Comparison Test")
    print("=" * 60)

    # Check files exist
    for path, name in [
        (vcf_path, "VCF"),
        (tree_path, "Tree"),
        (snp_path, "SNP DB"),
        (ground_truth_path, "Ground truth"),
    ]:
        if not path.exists():
            print(f"Error: {name} not found at {path}")
            sys.exit(1)

    # Load resources
    print("\nLoading resources...")
    tree = Tree.from_json(tree_path)
    print(f"  Tree: {len(tree)} nodes")

    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_path)
    print(f"  SNP DB: {len(snp_db)} SNPs")

    ground_truth = load_ground_truth(ground_truth_path)
    print(f"  Ground truth: {len(ground_truth)} samples")

    # Get sample list from VCF
    vcf = pysam.VariantFile(str(vcf_path))
    all_samples = list(vcf.header.samples)
    vcf.close()

    # Filter to samples with ground truth
    valid_samples = [s for s in all_samples if s in ground_truth]
    print(f"  Valid samples: {len(valid_samples)}")

    # Random subsample
    random.seed(args.seed)
    test_samples = random.sample(valid_samples, min(args.samples, len(valid_samples)))
    print(f"  Testing: {len(test_samples)} samples")

    # Run classification
    print("\nClassifying samples with both priors...")
    results = classify_samples(
        vcf_path=vcf_path,
        tree=tree,
        snp_db=snp_db,
        ground_truth=ground_truth,
        sample_ids=test_samples,
        reference="grch37",
        threads=args.threads,
    )

    print(f"  Completed: {len(results)} samples")

    # Analyze results
    print("\nAnalyzing results...")
    analysis = analyze_results(results)

    # Print summary
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"\nSamples tested: {analysis['total_samples']}")
    print(f"\nUniform Prior:")
    print(f"  Accuracy: {analysis['uniform_prior']['accuracy_pct']:.1f}%")
    print(f"  Mean confidence: {analysis['uniform_prior']['mean_confidence']:.4f}")
    print(f"\nCoalescent Prior:")
    print(f"  Accuracy: {analysis['coalescent_prior']['accuracy_pct']:.1f}%")
    print(f"  Mean confidence: {analysis['coalescent_prior']['mean_confidence']:.4f}")
    print(f"\nDifference: {analysis['difference_pp']:+.2f} pp")
    print(f"Disagreements: {analysis['disagreements']}")

    if analysis["disagreement_details"]:
        print("\nSample disagreements:")
        for d in analysis["disagreement_details"][:5]:
            u_mark = "✓" if d["uniform_correct"] else "✗"
            c_mark = "✓" if d["coalescent_correct"] else "✗"
            print(
                f"  {d['sample']}: GT={d['ground_truth']}, "
                f"Uniform={d['uniform']} {u_mark}, "
                f"Coalescent={d['coalescent']} {c_mark}"
            )

    # Interpretation
    print("\n" + "=" * 60)
    print("INTERPRETATION")
    print("=" * 60)
    diff = analysis["difference_pp"]
    if abs(diff) < 1.0:
        print("No significant difference between priors on this dataset.")
        print("Uniform priors recommended as the simpler, unbiased default.")
    elif diff > 0:
        print(f"Coalescent priors improved accuracy by {diff:.1f} pp.")
        print("Consider offering as an option for population-matched datasets.")
    else:
        print(f"Coalescent priors decreased accuracy by {abs(diff):.1f} pp.")
        print("This may indicate population bias against rare haplogroups.")
        print("Uniform priors recommended as the unbiased default.")

    # Save results
    output_path.parent.mkdir(exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(analysis, f, indent=2)
    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
