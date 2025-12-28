#!/usr/bin/env python3
"""
AADR Stratified Validation for yallHap Ancient DNA Modes.

Stratifies samples by variant density (called variants / total variants in Y-VCF)
and tests all ancient DNA modes:
- Heuristic (transversions-only)
- Bayesian + ancient mode
- Bayesian + transversions-only

Reports per-bin accuracy with variant density.

Usage:
    python scripts/validate_aadr_stratified.py -o results/aadr_stratified.md --samples-per-bin 100
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import pysam

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


@dataclass
class AADRSample:
    """AADR sample with variant density."""

    genetic_id: str
    variant_density: float  # Percentage (0-100) of called variants
    ground_truth_terminal: str
    ground_truth_isogg: str


@dataclass
class ClassificationResult:
    """Result from a classification run."""

    sample: str
    haplogroup: str
    confidence: float | None = None
    ground_truth: str = ""
    correct_major: bool = False
    variant_density: float = 0.0
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class BinResult:
    """Results for a variant density bin."""

    bin_name: str
    density_range: tuple[float, float]
    samples: int
    correct: int
    mean_confidence: float | None
    mean_density: float
    results: list[ClassificationResult] = field(default_factory=list)

    @property
    def accuracy(self) -> float:
        return 100 * self.correct / self.samples if self.samples > 0 else 0.0


def calculate_variant_density(
    vcf_path: Path,
    sample_ids: list[str],
    cache_path: Path | None = None,
) -> dict[str, float]:
    """
    Calculate variant density for each sample from the VCF.

    Variant density = (called variants / total variants) * 100

    Args:
        vcf_path: Path to VCF file
        sample_ids: List of sample IDs to calculate density for
        cache_path: Optional path to cache results

    Returns:
        Dictionary mapping sample_id -> density percentage (0-100)
    """
    # Check cache
    if cache_path and cache_path.exists():
        vcf_mtime = vcf_path.stat().st_mtime
        cache_mtime = cache_path.stat().st_mtime
        if cache_mtime > vcf_mtime:
            print(f"  Loading cached variant density from {cache_path.name}...")
            with open(cache_path) as f:
                cached = json.load(f)
            # Filter to requested samples
            return {s: cached[s] for s in sample_ids if s in cached}

    print("  Calculating variant density from VCF (this may take a few minutes)...")

    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)

    # Filter to samples in VCF
    target_samples = [s for s in sample_ids if s in vcf_samples]

    # Initialize counts
    called_counts: dict[str, int] = {s: 0 for s in target_samples}
    total_variants = 0

    # Count called variants per sample
    for i, rec in enumerate(vcf.fetch()):
        total_variants += 1
        if i % 5000 == 0 and i > 0:
            print(f"    Processed {i} variants...")
        for sample in target_samples:
            gt = rec.samples[sample]["GT"]
            if gt is not None and gt != (None,) and gt != (None, None):
                called_counts[sample] += 1

    vcf.close()
    print(f"    Processed {total_variants} total variants")

    # Calculate density percentages
    density_map: dict[str, float] = {}
    for sample, called in called_counts.items():
        density_map[sample] = (called / total_variants * 100) if total_variants > 0 else 0.0

    # Cache results (cache all samples for future use)
    if cache_path:
        print(f"  Caching variant density to {cache_path.name}...")
        # Load existing cache and merge
        all_densities = {}
        if cache_path.exists():
            with open(cache_path) as f:
                all_densities = json.load(f)
        all_densities.update(density_map)
        with open(cache_path, "w") as f:
            json.dump(all_densities, f)

    return density_map


def load_aadr_samples(
    vcf_path: Path,
    ground_truth_path: Path,
    cache_path: Path | None = None,
) -> list[AADRSample]:
    """Load AADR samples with variant density and ground truth."""
    # Load ground truth
    ground_truth: dict[str, tuple[str, str]] = {}
    with open(ground_truth_path) as f:
        header = f.readline().strip().split("\t")
        sample_idx = header.index("sample_id") if "sample_id" in header else 0
        terminal_idx = (
            header.index("haplogroup_terminal") if "haplogroup_terminal" in header else 1
        )
        isogg_idx = header.index("haplogroup_isogg") if "haplogroup_isogg" in header else -1

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > terminal_idx:
                sample_id = parts[sample_idx]
                terminal = parts[terminal_idx] if len(parts) > terminal_idx else ""
                isogg = parts[isogg_idx] if isogg_idx >= 0 and len(parts) > isogg_idx else ""
                if terminal and "-" in terminal and terminal[0].isalpha():
                    ground_truth[sample_id] = (terminal, isogg)

    # Get samples in VCF
    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)
    vcf.close()

    # Filter to samples in VCF
    valid_samples = [s for s in ground_truth if s in vcf_samples]
    print(f"  {len(valid_samples)} samples with ground truth in VCF")

    # Calculate variant density
    density_map = calculate_variant_density(vcf_path, valid_samples, cache_path)

    # Build sample list
    samples = []
    for sample_id in valid_samples:
        if sample_id not in density_map:
            continue
        terminal, isogg = ground_truth[sample_id]
        samples.append(
            AADRSample(
                genetic_id=sample_id,
                variant_density=density_map[sample_id],
                ground_truth_terminal=terminal,
                ground_truth_isogg=isogg,
            )
        )

    return samples


def stratify_samples(
    samples: list[AADRSample],
    bins: list[tuple[str, float, float]],
    samples_per_bin: int,
    seed: int = 42,
) -> dict[str, list[AADRSample]]:
    """Stratify samples into variant density bins."""
    random.seed(seed)
    stratified: dict[str, list[AADRSample]] = {name: [] for name, _, _ in bins}

    for sample in samples:
        density = sample.variant_density
        for name, low, high in bins:
            if low <= density < high:
                stratified[name].append(sample)
                break

    # Subsample each bin
    for name in stratified:
        if len(stratified[name]) > samples_per_bin:
            stratified[name] = random.sample(stratified[name], samples_per_bin)

    return stratified


def run_yallhap_on_samples(
    samples: list[AADRSample],
    vcf_path: Path,
    tree_path: Path,
    snp_db_path: Path,
    mode: str,
    threads: int = 1,
) -> list[ClassificationResult]:
    """Run yallHap on samples with specified mode."""
    from yallhap.classifier import HaplogroupClassifier
    from yallhap.snps import SNPDatabase
    from yallhap.tree import Tree

    tree = Tree.from_json(tree_path)
    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)

    # Configure based on mode
    bayesian = mode in ("bayesian", "bayesian_tv")
    ancient_mode = mode in ("bayesian", "bayesian_tv", "transversions")
    transversions_only = mode in ("transversions", "bayesian_tv")

    # For ancient+bayesian, use relaxed filtering
    min_depth = 1 if (ancient_mode and bayesian) else 1
    min_quality = 0 if (ancient_mode and bayesian) else 10

    classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference="grch37",
        bayesian=bayesian,
        ancient_mode=ancient_mode,
        transversions_only=transversions_only,
        min_depth=min_depth,
        min_quality=min_quality,
    )

    sample_ids = [s.genetic_id for s in samples]
    sample_map = {s.genetic_id: s for s in samples}

    results = classifier.classify_batch(vcf_path, sample_ids, threads=threads)

    classification_results = []
    for result in results:
        sample = sample_map[result.sample]
        correct = compare_major_haplogroup(
            sample.ground_truth_terminal, result.haplogroup
        )
        classification_results.append(
            ClassificationResult(
                sample=result.sample,
                haplogroup=result.haplogroup,
                confidence=result.confidence,
                ground_truth=sample.ground_truth_terminal,
                correct_major=correct,
                variant_density=sample.variant_density,
            )
        )

    return classification_results


def compare_major_haplogroup(expected: str, predicted: str) -> bool:
    """Compare major haplogroup letters."""
    if not expected or not predicted:
        return False
    if predicted in ("NA", "ERROR", "TIMEOUT", "ROOT"):
        return False
    return expected[0].upper() == predicted[0].upper()


def generate_report(
    mode_results: dict[str, dict[str, BinResult]],
    output_path: Path,
    bins: list[tuple[str, float, float]],
) -> None:
    """Generate markdown report."""
    with open(output_path, "w") as f:
        f.write("# AADR Stratified Validation Report\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Summary table
        f.write("## Summary by Variant Density Bin\n\n")

        # Header
        modes = list(mode_results.keys())
        f.write("| Density Bin | Mean Density | Samples |")
        for mode in modes:
            f.write(f" {mode} Acc. | {mode} Conf. |")
        f.write("\n")
        f.write("|-------------|--------------|---------|")
        for _ in modes:
            f.write("------------|------------|")
        f.write("\n")

        # Rows by bin
        for bin_name, low, high in bins:
            # Get first mode's bin for sample count and density
            first_mode = modes[0]
            if bin_name not in mode_results[first_mode]:
                continue

            bin_data = mode_results[first_mode][bin_name]
            f.write(f"| {bin_name} | {bin_data.mean_density:.1f}% | {bin_data.samples} |")

            for mode in modes:
                if bin_name in mode_results[mode]:
                    br = mode_results[mode][bin_name]
                    conf_str = f"{br.mean_confidence:.3f}" if br.mean_confidence else "N/A"
                    f.write(f" {br.accuracy:.1f}% | {conf_str} |")
                else:
                    f.write(" N/A | N/A |")
            f.write("\n")

        # Overall row
        f.write("| **Overall** | - | ")
        total_samples = sum(
            br.samples for br in mode_results[modes[0]].values()
        )
        f.write(f"{total_samples} |")
        for mode in modes:
            total_correct = sum(br.correct for br in mode_results[mode].values())
            total = sum(br.samples for br in mode_results[mode].values())
            acc = 100 * total_correct / total if total > 0 else 0
            confs = [
                br.mean_confidence
                for br in mode_results[mode].values()
                if br.mean_confidence
            ]
            mean_conf = sum(confs) / len(confs) if confs else None
            conf_str = f"{mean_conf:.3f}" if mean_conf else "N/A"
            f.write(f" {acc:.1f}% | {conf_str} |")
        f.write("\n\n")

        # Per-mode details
        for mode, bins_data in mode_results.items():
            f.write(f"## {mode} Mode Details\n\n")
            f.write("| Bin | Samples | Correct | Accuracy | Mean Conf | Mean Density |\n")
            f.write("|-----|---------|---------|----------|-----------|-------------|\n")
            for bin_name, _, _ in bins:
                if bin_name not in bins_data:
                    continue
                br = bins_data[bin_name]
                conf_str = f"{br.mean_confidence:.3f}" if br.mean_confidence else "N/A"
                f.write(
                    f"| {bin_name} | {br.samples} | {br.correct} | "
                    f"{br.accuracy:.1f}% | {conf_str} | {br.mean_density:.1f}% |\n"
                )
            f.write("\n")

        # Methods
        f.write("## Methods\n\n")
        f.write("### yallHap Modes Tested\n\n")
        f.write("- **heuristic_tv**: Heuristic classifier, transversions-only\n")
        f.write("- **bayesian_ancient**: Bayesian classifier with ancient mode ")
        f.write("(min-depth=1, min-quality=0)\n")
        f.write("- **bayesian_tv**: Bayesian classifier with transversions-only\n\n")
        f.write("### Variant Density\n\n")
        f.write("Variant density = (called variants / total variants in Y-VCF) × 100.\n")
        f.write("This directly measures data completeness for each sample.\n\n")
        f.write("### Density Bins\n\n")
        for name, low, high in bins:
            f.write(f"- **{name}**: {low}% to {high}% of variants called\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="AADR stratified validation by variant density")
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--samples-per-bin", type=int, default=100)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Don't use cached variant density calculations",
    )
    parser.add_argument(
        "--deciles",
        action="store_true",
        help="Use 10 equal-width bins (0-10%, 10-20%, ..., 90-100%) instead of default bins",
    )
    parser.add_argument(
        "--range",
        type=str,
        help="Custom range to analyze with 1%% bins, e.g., '1-10' for 1-10%% in 1%% increments",
    )
    args = parser.parse_args()

    project_root = Path(__file__).parent.parent

    # Paths
    vcf_path = project_root / "data" / "ancient" / "aadr_chrY_v2.vcf.gz"
    ground_truth_path = project_root / "data" / "ancient" / "aadr_1240k_ground_truth.tsv"
    tree_path = project_root / "data" / "yfull_tree.json"
    snp_db_path = project_root / "data" / "validation" / "ybrowse_snps_hg19.csv"
    cache_path = project_root / "data" / "ancient" / "aadr_variant_density.json"

    for path in [vcf_path, ground_truth_path, tree_path, snp_db_path]:
        if not path.exists():
            print(f"ERROR: {path} not found")
            return 1

    # Variant density bins (percentage of called variants)
    if args.range:
        # Custom range with 1% bins, e.g., "1-10" -> 1-2%, 2-3%, ..., 9-10%
        try:
            start, end = map(int, args.range.split("-"))
        except ValueError:
            print(f"ERROR: Invalid range format '{args.range}'. Use 'START-END' like '1-10'")
            return 1
        bins = [
            (f"{i}-{i+1}%", float(i), float(i + 1) + (0.1 if i == end - 1 else 0.0))
            for i in range(start, end)
        ]
    elif args.deciles:
        # 10 equal-width bins: 0-10%, 10-20%, ..., 90-100%
        bins = [
            (f"{i*10}-{(i+1)*10}%", i * 10.0, (i + 1) * 10.0 + (0.1 if i == 9 else 0.0))
            for i in range(10)
        ]
    else:
        bins = [
            ("<1%", 0.0, 1.0),
            ("1-10%", 1.0, 10.0),
            ("10-50%", 10.0, 50.0),
            ("≥50%", 50.0, 100.1),  # 100.1 to include 100%
        ]

    # Modes to test
    modes = {
        "heuristic_tv": "transversions",
        "bayesian_ancient": "bayesian",
        "bayesian_tv": "bayesian_tv",
    }

    print("=" * 70)
    print("AADR Stratified Validation (by Variant Density)")
    print("=" * 70)

    # Load samples
    print("\nLoading samples...")
    effective_cache = None if args.no_cache else cache_path
    samples = load_aadr_samples(vcf_path, ground_truth_path, effective_cache)
    print(f"  Found {len(samples)} samples with variant density and ground truth")

    # Stratify
    print(f"\nStratifying into {len(bins)} bins ({args.samples_per_bin} samples each)...")
    stratified = stratify_samples(samples, bins, args.samples_per_bin)
    for name, samples_list in stratified.items():
        if samples_list:
            mean_density = sum(s.variant_density for s in samples_list) / len(samples_list)
            print(f"  {name}: {len(samples_list)} samples (mean density: {mean_density:.1f}%)")

    # Flatten for classification
    all_samples = []
    for samples_list in stratified.values():
        all_samples.extend(samples_list)
    print(f"\nTotal samples to classify: {len(all_samples)}")

    # Run each mode
    mode_results: dict[str, dict[str, BinResult]] = {}

    for mode_name, mode_type in modes.items():
        print(f"\n[{mode_name}] Running classification...")
        start = time.time()

        results = run_yallhap_on_samples(
            all_samples,
            vcf_path,
            tree_path,
            snp_db_path,
            mode_type,
            threads=args.threads,
        )

        elapsed = time.time() - start
        print(f"  Completed in {elapsed:.1f}s")

        # Group by bin
        result_map = {r.sample: r for r in results}
        bin_results: dict[str, BinResult] = {}

        for bin_name, low, high in bins:
            bin_samples = stratified[bin_name]
            if not bin_samples:
                continue

            bin_classifications = [result_map[s.genetic_id] for s in bin_samples]
            correct = sum(1 for r in bin_classifications if r.correct_major)
            confs = [r.confidence for r in bin_classifications if r.confidence]
            mean_conf = sum(confs) / len(confs) if confs else None
            mean_density = sum(s.variant_density for s in bin_samples) / len(bin_samples)

            bin_results[bin_name] = BinResult(
                bin_name=bin_name,
                density_range=(low, high),
                samples=len(bin_samples),
                correct=correct,
                mean_confidence=mean_conf,
                mean_density=mean_density,
                results=bin_classifications,
            )

            print(f"    {bin_name}: {correct}/{len(bin_samples)} ({100*correct/len(bin_samples):.1f}%)")

        mode_results[mode_name] = bin_results

    # Generate report
    print(f"\nGenerating report: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    generate_report(mode_results, args.output, bins)

    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(main())

