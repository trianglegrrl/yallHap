#!/usr/bin/env python3
"""
Power analysis for ancient DNA classification.

Answers: At what confidence threshold does classification become meaningful
at each variant density level?

For each density bin, calculates:
- Accuracy at confidence thresholds (0.5, 0.7, 0.8, 0.9, 0.95)
- False discovery rate at each threshold
- Proportion of samples passing each threshold

Usage:
    python scripts/power_analysis_ancient.py [--threads N]
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from yallhap.classifier import HaplogroupClassifier
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


# Density bins matching the paper
DENSITY_BINS = {
    "<1%": (0.0, 1.0),
    "1-4%": (1.0, 4.0),
    "4-10%": (4.0, 10.0),
    "10-50%": (10.0, 50.0),
    ">=50%": (50.0, 100.1),
}

# Confidence thresholds to test
CONFIDENCE_THRESHOLDS = [0.5, 0.7, 0.8, 0.9, 0.95]


def get_density_bin(density_pct: float) -> str | None:
    """Assign variant density percentage to bin."""
    for bin_name, (low, high) in DENSITY_BINS.items():
        if low <= density_pct < high:
            return bin_name
    return None


def load_ground_truth(tsv_path: Path) -> dict[str, str]:
    """Load ground truth haplogroups from AADR TSV."""
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


def load_variant_density(cache_path: Path) -> dict[str, float]:
    """Load cached variant density from JSON."""
    if cache_path.exists():
        with open(cache_path) as f:
            data = json.load(f)
        # Handle both formats: flat dict or nested under "densities"
        if isinstance(data, dict):
            if "densities" in data:
                return data["densities"]
            # Flat format - sample -> density
            return data
    return {}


def get_major_clade(haplogroup: str) -> str:
    """Extract major clade from haplogroup."""
    if not haplogroup:
        return ""
    hg = haplogroup.upper()
    if hg.startswith("ROOT"):
        return "ROOT"
    for i, c in enumerate(hg):
        if c.isdigit() or c == "-":
            return hg[:i] if i > 0 else hg[0]
    return hg[0] if hg else ""


@dataclass
class ThresholdStats:
    """Statistics for a single confidence threshold."""

    threshold: float
    samples_passing: int
    samples_total: int
    correct: int
    incorrect: int
    accuracy_pct: float
    false_discovery_rate: float
    retention_rate: float


@dataclass
class BinAnalysis:
    """Analysis for a single density bin."""

    bin_name: str
    total_samples: int
    thresholds: list[ThresholdStats] = field(default_factory=list)


def analyze_threshold(
    results: list[tuple[str, str, float, bool]],  # (sample, predicted, conf, correct)
    threshold: float,
) -> ThresholdStats:
    """Analyze accuracy at a specific confidence threshold."""
    passing = [(s, p, c, correct) for s, p, c, correct in results if c >= threshold]
    total = len(results)
    n_passing = len(passing)

    correct = sum(1 for _, _, _, c in passing if c)
    incorrect = n_passing - correct

    accuracy = correct / n_passing * 100 if n_passing > 0 else 0.0
    fdr = incorrect / n_passing if n_passing > 0 else 0.0
    retention = n_passing / total * 100 if total > 0 else 0.0

    return ThresholdStats(
        threshold=threshold,
        samples_passing=n_passing,
        samples_total=total,
        correct=correct,
        incorrect=incorrect,
        accuracy_pct=round(accuracy, 2),
        false_discovery_rate=round(fdr, 4),
        retention_rate=round(retention, 2),
    )


def main() -> None:
    """Run power analysis on ancient DNA samples."""
    parser = argparse.ArgumentParser(description="Power analysis for ancient DNA")
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads (default: 16)",
    )
    args = parser.parse_args()

    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"
    ancient_dir = data_dir / "ancient"
    validation_dir = data_dir / "validation"
    results_dir = base_dir / "results"

    # File paths
    vcf_path = ancient_dir / "aadr_chrY_v2.vcf.gz"
    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = validation_dir / "ybrowse_snps_hg19.csv"
    ground_truth_path = ancient_dir / "aadr_1240k_ground_truth.tsv"
    density_cache_path = ancient_dir / "aadr_variant_density.json"
    output_path = results_dir / "power_analysis_ancient.json"

    print("=" * 70)
    print("Power Analysis for Ancient DNA Classification")
    print("=" * 70)

    # Check files exist
    for path, name in [
        (vcf_path, "VCF"),
        (tree_path, "Tree"),
        (snp_db_path, "SNP DB"),
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

    ground_truth = load_ground_truth(ground_truth_path)
    print(f"  Ground truth: {len(ground_truth)} samples")

    # Load variant density cache
    density_cache = load_variant_density(density_cache_path)
    print(f"  Density cache: {len(density_cache)} samples")

    if not density_cache:
        print("Error: No density cache found. Run validate_aadr_stratified.py first.")
        sys.exit(1)

    # Create classifier (Bayesian ancient mode)
    print("\nCreating classifier...")
    classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference="grch37",
        bayesian=True,
        ancient_mode=True,
        min_depth=1,
        min_quality=0,
    )

    # Group samples by density bin
    samples_by_bin: dict[str, list[str]] = {bin_name: [] for bin_name in DENSITY_BINS}
    for sample, density in density_cache.items():
        if sample in ground_truth:
            bin_name = get_density_bin(density)
            if bin_name:
                samples_by_bin[bin_name].append(sample)

    print("\nSamples by density bin:")
    for bin_name, samples in samples_by_bin.items():
        print(f"  {bin_name}: {len(samples)}")

    # Classify all samples
    all_samples = [s for samples in samples_by_bin.values() for s in samples]
    print(f"\nClassifying {len(all_samples)} samples...")

    results = classifier.classify_batch(vcf_path, all_samples, threads=args.threads)
    result_map = {r.sample: r for r in results}

    # Analyze each bin
    bin_analyses: list[BinAnalysis] = []

    for bin_name, samples in samples_by_bin.items():
        print(f"\nAnalyzing {bin_name}...")

        # Collect results for this bin
        bin_results: list[tuple[str, str, float, bool]] = []
        for sample in samples:
            if sample not in result_map:
                continue
            result = result_map[sample]
            gt = ground_truth.get(sample, "")

            gt_major = get_major_clade(gt)
            pred_major = get_major_clade(result.haplogroup)
            correct = gt_major == pred_major

            bin_results.append((sample, result.haplogroup, result.confidence, correct))

        # Analyze at each threshold
        analysis = BinAnalysis(bin_name=bin_name, total_samples=len(bin_results))

        for threshold in CONFIDENCE_THRESHOLDS:
            stats = analyze_threshold(bin_results, threshold)
            analysis.thresholds.append(stats)

        bin_analyses.append(analysis)

    # Print summary
    print("\n" + "=" * 70)
    print("POWER ANALYSIS RESULTS")
    print("=" * 70)

    for analysis in bin_analyses:
        print(f"\n{analysis.bin_name} (n={analysis.total_samples}):")
        print(
            f"  {'Threshold':<10} {'Passing':<10} {'Accuracy':<12} {'FDR':<10} {'Retention'}"
        )
        print(f"  {'-' * 55}")
        for t in analysis.thresholds:
            print(
                f"  {t.threshold:<10.2f} {t.samples_passing:<10} {t.accuracy_pct:>6.1f}%"
                f"      {t.false_discovery_rate:>6.2%}    {t.retention_rate:>5.1f}%"
            )

    # Find minimum density for reliable classification
    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    for analysis in bin_analyses:
        # Find threshold where accuracy > 80%
        best_threshold = None
        for t in analysis.thresholds:
            if t.accuracy_pct >= 80 and t.samples_passing >= 10:
                best_threshold = t
                break

        if best_threshold:
            print(
                f"\n{analysis.bin_name}: At threshold {best_threshold.threshold}, "
                f"accuracy = {best_threshold.accuracy_pct}% (n={best_threshold.samples_passing})"
            )
        else:
            # Check if any threshold achieves >70%
            for t in analysis.thresholds:
                if t.accuracy_pct >= 70 and t.samples_passing >= 5:
                    print(
                        f"\n{analysis.bin_name}: At threshold {t.threshold}, "
                        f"accuracy = {t.accuracy_pct}% (n={t.samples_passing}, best available)"
                    )
                    break
            else:
                print(
                    f"\n{analysis.bin_name}: No threshold achieves >70% accuracy "
                    f"with sufficient samples"
                )

    # Save results
    output_data = {
        "density_bins": [
            {
                "bin": a.bin_name,
                "total_samples": a.total_samples,
                "thresholds": [
                    {
                        "threshold": t.threshold,
                        "samples_passing": t.samples_passing,
                        "accuracy_pct": t.accuracy_pct,
                        "false_discovery_rate": t.false_discovery_rate,
                        "retention_rate": t.retention_rate,
                    }
                    for t in a.thresholds
                ],
            }
            for a in bin_analyses
        ],
    }

    results_dir.mkdir(exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()

