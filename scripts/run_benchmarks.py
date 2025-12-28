#!/usr/bin/env python3
"""
Run comprehensive benchmarks comparing yallHap performance across datasets.

This script validates yallHap against three primary datasets:
1. 1000 Genomes Phase 3 (GRCh37) - 1,233 modern males
2. AADR v54 (GRCh37) - Ancient DNA samples with transversions-only mode
3. gnomAD HGDP/1000G (GRCh38) - High-coverage with AD fields for Bayesian mode

Output:
- Console summary table
- JSON results in results/benchmark_results.json

Usage:
    # Full benchmark (all samples)
    python scripts/run_benchmarks.py

    # Quick benchmark with subsampling and threading
    python scripts/run_benchmarks.py --subsample 50 --threads 16
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Progress bar support (optional)
try:
    from tqdm import tqdm

    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    tqdm = None  # type: ignore


# Variant density bin constants for stratified AADR analysis
# Based on actual called variants / total variants in Y-only VCF
# 5 bins with 1-4% / 4-10% split to highlight where Bayesian mode improves
DENSITY_BINS = ["<1%", "1-4%", "4-10%", "10-50%", ">=50%"]

# Minimum variant density (%) for reliable ISOGG comparison
# Below this threshold, samples are excluded from ISOGG metrics
ISOGG_MIN_DENSITY = 4.0


# =============================================================================
# ISOGG Comparison Functions
# =============================================================================


def normalize_isogg(hg: str) -> str:
    """
    Normalize ISOGG haplogroup for comparison.

    - Remove tildes (~) indicating uncertainty
    - Uppercase
    - Strip whitespace
    """
    if not hg:
        return ""
    return hg.replace("~", "").strip().upper()


def get_major_clade(hg: str) -> str:
    """Extract major clade (first letter) from haplogroup."""
    if not hg:
        return ""
    hg = normalize_isogg(hg)
    if not hg:
        return ""
    return hg[0] if hg[0].isalpha() else ""


def is_isogg_prefix_match(hg1: str, hg2: str) -> bool:
    """Check if one haplogroup is a prefix of the other."""
    h1 = normalize_isogg(hg1)
    h2 = normalize_isogg(hg2)
    if not h1 or not h2:
        return False
    return h1.startswith(h2) or h2.startswith(h1)


def classify_isogg_match(predicted: str, ground_truth: str) -> str:
    """
    Classify the match type between predicted and ground truth ISOGG.

    Returns one of: exact, prefix, major_clade, mismatch
    """
    pred_norm = normalize_isogg(predicted)
    gt_norm = normalize_isogg(ground_truth)

    if not pred_norm or not gt_norm:
        return "mismatch"

    # Exact match
    if pred_norm == gt_norm:
        return "exact"

    # Prefix match (one is ancestor of the other)
    if is_isogg_prefix_match(pred_norm, gt_norm):
        return "prefix"

    # Major clade match (same first letter)
    if get_major_clade(pred_norm) == get_major_clade(gt_norm):
        return "major_clade"

    return "mismatch"


def get_density_bin(density_pct: float) -> str:
    """Assign variant density percentage to bin."""
    if density_pct < 1.0:
        return "<1%"
    elif density_pct < 4.0:
        return "1-4%"
    elif density_pct < 10.0:
        return "4-10%"
    elif density_pct < 50.0:
        return "10-50%"
    else:
        return ">=50%"


def calculate_variant_density(
    vcf_path: Path,
    sample_ids: list[str],
) -> dict[str, float]:
    """
    Calculate variant density for each sample from the VCF.

    Variant density = (called variants / total variants) * 100

    Args:
        vcf_path: Path to VCF file
        sample_ids: List of sample IDs to calculate density for

    Returns:
        Dictionary mapping sample_id -> density percentage (0-100)
    """
    import pysam

    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)

    # Filter to samples in VCF
    target_samples = [s for s in sample_ids if s in vcf_samples]

    # Initialize counts
    called_counts: dict[str, int] = {s: 0 for s in target_samples}
    total_variants = 0

    # Count called variants per sample
    for rec in vcf.fetch():
        total_variants += 1
        for sample in target_samples:
            gt = rec.samples[sample]["GT"]
            if gt is not None and gt != (None,) and gt != (None, None):
                called_counts[sample] += 1

    vcf.close()

    # Calculate density percentages
    density_map: dict[str, float] = {}
    for sample, called in called_counts.items():
        density_map[sample] = (called / total_variants * 100) if total_variants > 0 else 0.0

    return density_map


def run_with_progress(
    cmd: list[str],
    desc: str,
    output_file: Path | None = None,
    source_file: Path | None = None,
) -> None:
    """
    Run a subprocess with a progress indicator.

    Shows elapsed time and optionally output file size growth.
    Uses tqdm if available, otherwise prints status updates.

    Args:
        cmd: Command to run
        desc: Description to show
        output_file: Optional output file to monitor for size progress
        source_file: Optional source file to estimate progress from
    """
    import subprocess
    import threading

    # Get source file size for progress estimation
    source_size = source_file.stat().st_size if source_file and source_file.exists() else None

    # Start the process
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if TQDM_AVAILABLE and output_file:
        # Show progress bar based on output file size
        # Estimate: output will be ~1-10% of input for filtered extraction
        estimated_size = source_size // 50 if source_size else None  # ~2% estimate

        with tqdm(
            total=estimated_size,
            desc=f"    {desc}",
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            leave=False,
            bar_format="{desc}: {n_fmt}/{total_fmt} [{elapsed}<{remaining}]" if estimated_size else "{desc}: {n_fmt} [{elapsed}]",
        ) as pbar:
            last_size = 0
            while process.poll() is None:
                if output_file.exists():
                    current_size = output_file.stat().st_size
                    pbar.update(current_size - last_size)
                    last_size = current_size
                time.sleep(0.5)

            # Final update
            if output_file.exists():
                final_size = output_file.stat().st_size
                pbar.update(final_size - last_size)
                # Update total to actual if we underestimated
                if estimated_size and final_size > estimated_size:
                    pbar.total = final_size
                    pbar.refresh()
    else:
        # Fallback: print status with elapsed time
        start = time.time()
        while process.poll() is None:
            elapsed = time.time() - start
            size_str = ""
            if output_file and output_file.exists():
                size_mb = output_file.stat().st_size / (1024 * 1024)
                size_str = f" ({size_mb:.1f} MB)"
            print(f"\r    {desc}: {elapsed:.0f}s{size_str}...", end="", flush=True)
            time.sleep(1)
        print()  # Newline after progress

    # Check for errors
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, cmd, stdout, stderr)


@dataclass
class StratifiedResult:
    """Results for a single variant density bin."""

    bin_name: str
    accuracy: float
    correct: int
    total: int
    samples_tested: int


@dataclass
class ISOGGMetrics:
    """ISOGG nomenclature comparison metrics."""

    total: int = 0
    exact: int = 0
    prefix: int = 0  # One is ancestor of other (compatible)
    major_clade: int = 0  # Same major clade but different subclade
    mismatch: int = 0  # Different major clades

    @property
    def compatible(self) -> int:
        """Exact + prefix matches (phylogenetically compatible assignments)."""
        return self.exact + self.prefix

    @property
    def compatible_rate(self) -> float:
        """Percentage of compatible assignments."""
        return 100 * self.compatible / self.total if self.total > 0 else 0.0

    @property
    def exact_rate(self) -> float:
        """Percentage of exact matches."""
        return 100 * self.exact / self.total if self.total > 0 else 0.0


@dataclass
class BenchmarkResult:
    """Results for a single tool/dataset combination."""

    tool_name: str
    dataset: str
    total_samples: int
    same_major_lineage: int
    exact_match: int
    runtime_seconds: float
    mean_confidence: float | None = None
    mean_derived_snps: float | None = None
    notes: list[str] = field(default_factory=list)
    stratified: list[StratifiedResult] = field(default_factory=list)
    isogg_metrics: ISOGGMetrics | None = None

    @property
    def major_lineage_rate(self) -> float:
        """Percentage with same major haplogroup letter."""
        if self.total_samples == 0:
            return 0.0
        return 100 * self.same_major_lineage / self.total_samples

    @property
    def exact_match_rate(self) -> float:
        """Percentage with exact haplogroup match."""
        if self.total_samples == 0:
            return 0.0
        return 100 * self.exact_match / self.total_samples


def load_ground_truth(path: Path, haplogroup_col: str = "haplogroup") -> dict[str, str]:
    """Load ground truth haplogroups from TSV file."""
    ground_truth = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_id = row.get("sample_id", row.get("sample", ""))
            haplogroup = row.get(haplogroup_col, row.get("haplogroup", ""))
            if sample_id and haplogroup:
                ground_truth[sample_id] = haplogroup
    return ground_truth


def compare_haplogroups(
    ground_truth: dict[str, str],
    predictions: dict[str, str],
) -> tuple[int, int, int]:
    """
    Compare predictions to ground truth.

    Returns (same_major_lineage_count, exact_match_count, total_compared)
    """
    same_major = 0
    exact = 0
    total = 0

    for sample_id, gt_hg in ground_truth.items():
        if sample_id not in predictions:
            continue

        pred_hg = predictions[sample_id]
        if not pred_hg or pred_hg == "NA":
            continue

        total += 1

        # Extract first letter for major lineage comparison
        gt_major = gt_hg[0].upper() if gt_hg else ""
        pred_major = pred_hg[0].upper() if pred_hg else ""

        if gt_major == pred_major:
            same_major += 1

        # Exact match (case-insensitive)
        if gt_hg.lower() == pred_hg.lower():
            exact += 1

    return same_major, exact, total


def run_yallhap_classification(
    vcf_path: Path,
    tree_path: Path,
    snp_db_path: Path,
    samples: list[str],
    reference: str = "grch37",
    bayesian: bool = False,
    ancient_mode: bool = False,
    transversions_only: bool = False,
    min_depth: int = 1,
    threads: int = 1,
    estimate_contamination: bool = False,
    max_tolerance: int = 3,
) -> tuple[dict[str, str], float, float, float, dict[str, float] | None]:
    """
    Run yallHap classification on samples.

    Args:
        vcf_path: Path to VCF file
        tree_path: Path to tree JSON
        snp_db_path: Path to SNP database CSV
        samples: List of sample IDs to classify
        reference: Reference genome (grch37, grch38, t2t)
        bayesian: Use Bayesian classifier
        ancient_mode: Enable ancient DNA mode
        transversions_only: Only use transversions
        min_depth: Minimum read depth
        threads: Number of parallel threads for batch processing
        estimate_contamination: Estimate Y-chromosome contamination
        max_tolerance: Maximum ancestral calls for path traversal

    Returns (predictions, runtime, mean_confidence, mean_derived_snps, contamination_rates)
    """
    from yallhap.classifier import HaplogroupClassifier
    from yallhap.snps import SNPDatabase
    from yallhap.tree import Tree

    # Load resources
    tree = Tree.from_json(tree_path)
    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)

    classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference=reference,
        bayesian=bayesian,
        ancient_mode=ancient_mode,
        transversions_only=transversions_only,
        min_depth=min_depth,
    )

    # Run classification with optional parallelism and progress bar
    start_time = time.time()

    mode_str = "Bayesian" if bayesian else "Heuristic"
    desc = f"    {mode_str} ({threads} threads)"

    if TQDM_AVAILABLE:
        pbar = tqdm(total=len(samples), desc=desc, unit="samples", leave=False)

        def progress_callback() -> None:
            pbar.update(1)

        results = classifier.classify_batch(
            vcf_path, samples, threads=threads, progress_callback=progress_callback
        )
        pbar.close()
    else:
        print(f"{desc}: classifying {len(samples)} samples...")
        results = classifier.classify_batch(vcf_path, samples, threads=threads)

    runtime = time.time() - start_time

    # Extract predictions and statistics
    predictions = {}
    confidences = []
    derived_counts = []

    for result in results:
        predictions[result.sample] = result.haplogroup
        if result.confidence is not None:
            confidences.append(result.confidence)
        if result.snp_stats and result.snp_stats.derived > 0:
            derived_counts.append(result.snp_stats.derived)

    mean_conf = sum(confidences) / len(confidences) if confidences else 0.0
    mean_derived = sum(derived_counts) / len(derived_counts) if derived_counts else 0.0

    # Collect contamination rates if requested
    contamination_rates: dict[str, float] | None = None
    if estimate_contamination:
        contamination_rates = {}
        for result in results:
            if hasattr(result, "contamination") and result.contamination is not None:
                contamination_rates[result.sample] = result.contamination.rate

    return predictions, runtime, mean_conf, mean_derived, contamination_rates


def benchmark_1kg(
    base_dir: Path,
    subsample: int | None = None,
    threads: int = 1,
    shared_samples: list[str] | None = None,
) -> tuple[BenchmarkResult | None, BenchmarkResult | None]:
    """Run 1000 Genomes Phase 3 benchmark with heuristic and Bayesian modes."""
    validation_dir = base_dir / "data" / "validation"
    data_dir = base_dir / "data"

    ground_truth_path = validation_dir / "poznik2016_haplogroups.tsv"
    vcf_path = validation_dir / "1kg_chrY_phase3.vcf.gz"
    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = validation_dir / "ybrowse_snps_hg19.csv"

    # Check required files
    for path in [ground_truth_path, vcf_path, tree_path, snp_db_path]:
        if not path.exists():
            print(f"  SKIP: {path.name} not found")
            return None, None

    # Load ground truth
    ground_truth = load_ground_truth(ground_truth_path)

    # Use shared_samples if provided (for cross-dataset comparison)
    if shared_samples:
        samples = [s for s in shared_samples if s in ground_truth]
    else:
        samples = list(ground_truth.keys())
        # Subsample if requested
        if subsample and len(samples) > subsample:
            random.seed(42)
            samples = random.sample(samples, subsample)

    print(f"  Loaded {len(samples)} samples")

    # Run yallHap heuristic mode
    heur_preds, heur_runtime, heur_conf, heur_derived, _ = run_yallhap_classification(
        vcf_path=vcf_path,
        tree_path=tree_path,
        snp_db_path=snp_db_path,
        samples=samples,
        reference="grch37",
        bayesian=False,
        threads=threads,
    )
    heur_major, heur_exact, heur_total = compare_haplogroups(ground_truth, heur_preds)

    heuristic_result = BenchmarkResult(
        tool_name="yallHap",
        dataset="1000G Phase 3",
        total_samples=heur_total,
        same_major_lineage=heur_major,
        exact_match=heur_exact,
        runtime_seconds=heur_runtime,
        mean_confidence=heur_conf,
        mean_derived_snps=heur_derived,
        notes=["GRCh37 reference", "Heuristic mode"],
    )

    # Run yallHap Bayesian mode
    bayes_preds, bayes_runtime, bayes_conf, bayes_derived, _ = run_yallhap_classification(
        vcf_path=vcf_path,
        tree_path=tree_path,
        snp_db_path=snp_db_path,
        samples=samples,
        reference="grch37",
        bayesian=True,
        threads=threads,
    )
    bayes_major, bayes_exact, bayes_total = compare_haplogroups(ground_truth, bayes_preds)

    bayesian_result = BenchmarkResult(
        tool_name="yallHap-Bayes",
        dataset="1000G Phase 3",
        total_samples=bayes_total,
        same_major_lineage=bayes_major,
        exact_match=bayes_exact,
        runtime_seconds=bayes_runtime,
        mean_confidence=bayes_conf,
        mean_derived_snps=bayes_derived,
        notes=["GRCh37 reference", "Bayesian mode"],
    )

    return heuristic_result, bayesian_result


def benchmark_aadr(
    base_dir: Path,
    samples_per_bin: int = 300,
    threads: int = 1,
    estimate_contamination: bool = False,
    max_tolerance: int = 3,
    isogg: bool = False,
    isogg_db_path: Path | None = None,
) -> tuple[BenchmarkResult | None, BenchmarkResult | None]:
    """
    Run AADR ancient DNA benchmark with stratified variant density analysis.

    Samples are stratified by variant density into bins (<1%, 1-4%, 4-10%, 10-50%, >=50%)
    based on the percentage of called variants in the Y-only VCF.

    Args:
        base_dir: Project base directory
        samples_per_bin: Maximum samples per density bin (default: 300)
        threads: Parallel threads for classification
        estimate_contamination: Estimate Y-chromosome contamination
        max_tolerance: Maximum ancestral calls for path traversal
        isogg: Output ISOGG haplogroup nomenclature
        isogg_db_path: Path to ISOGG SNP database

    Returns:
        Tuple of (heuristic_result, bayesian_result) with stratified results
    """
    ancient_dir = base_dir / "data" / "ancient"
    data_dir = base_dir / "data"
    validation_dir = base_dir / "data" / "validation"

    ground_truth_path = ancient_dir / "aadr_1240k_ground_truth.tsv"
    vcf_path = ancient_dir / "aadr_chrY_v2.vcf.gz"
    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = validation_dir / "ybrowse_snps_hg19.csv"

    # Check required files (no longer need anno file for coverage)
    for path in [ground_truth_path, vcf_path, tree_path, snp_db_path]:
        if not path.exists():
            print(f"  SKIP: {path.name} not found")
            return None, None

    # Load ground truth - filter to proper haplogroup format
    raw_ground_truth = load_ground_truth(ground_truth_path, "haplogroup_terminal")
    ground_truth = {}
    for sample_id, hg in raw_ground_truth.items():
        # Only keep proper haplogroup names (X-YYYY format)
        if hg and "-" in hg and hg[0].isalpha():
            ground_truth[sample_id] = hg

    # Load ISOGG ground truth if available
    isogg_ground_truth = load_ground_truth(ground_truth_path, "haplogroup_isogg")

    # Get samples in VCF
    import pysam

    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)
    vcf.close()

    # Filter ground truth to samples in VCF
    valid_samples = [s for s in ground_truth if s in vcf_samples]
    print(f"  {len(valid_samples)} samples with ground truth in VCF")

    # Calculate variant density for all samples
    print("  Calculating variant density from VCF (this may take a minute)...")
    density_map = calculate_variant_density(vcf_path, valid_samples)
    print(f"    Calculated density for {len(density_map)} samples")

    # Group samples by density bin
    samples_by_bin: dict[str, list[str]] = defaultdict(list)
    for sample_id in valid_samples:
        density = density_map.get(sample_id, 0.0)
        bin_name = get_density_bin(density)
        samples_by_bin[bin_name].append(sample_id)

    print("  Samples by variant density bin:")
    for bin_name in DENSITY_BINS:
        count = len(samples_by_bin.get(bin_name, []))
        print(f"    {bin_name}: {count}")

    # Load classifier resources
    from yallhap.classifier import HaplogroupClassifier
    from yallhap.snps import SNPDatabase
    from yallhap.tree import Tree

    tree = Tree.from_json(tree_path)
    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)

    # Load ISOGG database if requested
    isogg_mapper = None
    if isogg:
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper

        isogg_path = isogg_db_path or (base_dir / "data" / "isogg_snps_grch38.txt")
        if isogg_path.exists():
            isogg_db = ISOGGDatabase.from_file(isogg_path)
            isogg_mapper = ISOGGMapper(tree, isogg_db)
            print(f"  Loaded {len(isogg_db)} ISOGG SNPs")
        else:
            print(f"  Warning: ISOGG database not found at {isogg_path}")

    # Create classifiers with enhanced ancient mode options
    # Note: estimate_contamination and max_tolerance are post-processing options,
    # not classifier constructor arguments
    heur_classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference="grch37",
        transversions_only=True,
    )

    bayes_classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference="grch37",
        transversions_only=True,
        bayesian=True,
        ancient_mode=True,
    )

    # Store config for potential post-processing
    _ = estimate_contamination  # Will be used in future enhancement
    _ = max_tolerance  # Will be used in future enhancement

    # Run stratified validation
    heur_stratified: list[StratifiedResult] = []
    bayes_stratified: list[StratifiedResult] = []
    heur_total_correct = 0
    heur_total_tested = 0
    bayes_total_correct = 0
    bayes_total_tested = 0
    total_runtime = 0.0
    all_confidences: list[float] = []
    all_derived: list[int] = []

    # ISOGG metrics if mapper is available
    heur_isogg = ISOGGMetrics() if isogg_mapper else None
    bayes_isogg = ISOGGMetrics() if isogg_mapper else None

    for bin_name in DENSITY_BINS:
        samples = samples_by_bin.get(bin_name, [])
        if not samples:
            continue

        # Limit samples per bin
        if len(samples) > samples_per_bin:
            random.seed(42)
            samples = random.sample(samples, samples_per_bin)

        print(f"  {bin_name}: classifying {len(samples)} samples...")
        start = time.time()

        # Heuristic mode
        heur_results = heur_classifier.classify_batch(vcf_path, samples, threads=threads)

        heur_correct = 0
        heur_bin_total = 0
        for result in heur_results:
            expected = ground_truth.get(result.sample, "")
            called = result.haplogroup

            expected_major = expected[0].upper() if expected else ""
            called_major = called[0].upper() if called and called != "NA" else ""

            if called and called != "NA":
                heur_bin_total += 1
                if expected_major == called_major:
                    heur_correct += 1
                if result.confidence is not None:
                    all_confidences.append(result.confidence)
                if result.snp_stats and result.snp_stats.derived > 0:
                    all_derived.append(result.snp_stats.derived)

                # ISOGG comparison (only for samples with sufficient coverage)
                if heur_isogg and isogg_mapper:
                    sample_density = density_map.get(result.sample, 0)
                    if sample_density >= ISOGG_MIN_DENSITY:
                        isogg_gt = isogg_ground_truth.get(result.sample, "")
                        if isogg_gt and isogg_gt not in (".", "..", "NA", "n/a"):
                            isogg_pred = isogg_mapper.to_isogg(result.haplogroup)
                            match = classify_isogg_match(isogg_pred, isogg_gt)
                            heur_isogg.total += 1
                            if match == "exact":
                                heur_isogg.exact += 1
                            elif match == "prefix":
                                heur_isogg.prefix += 1
                            elif match == "major_clade":
                                heur_isogg.major_clade += 1
                            else:
                                heur_isogg.mismatch += 1

        heur_accuracy = heur_correct / heur_bin_total if heur_bin_total > 0 else 0
        heur_stratified.append(StratifiedResult(
            bin_name=bin_name,
            accuracy=heur_accuracy,
            correct=heur_correct,
            total=heur_bin_total,
            samples_tested=len(samples),
        ))
        heur_total_correct += heur_correct
        heur_total_tested += heur_bin_total

        # Bayesian mode
        bayes_results = bayes_classifier.classify_batch(vcf_path, samples, threads=threads)

        bayes_correct = 0
        bayes_bin_total = 0
        for result in bayes_results:
            expected = ground_truth.get(result.sample, "")
            called = result.haplogroup

            expected_major = expected[0].upper() if expected else ""
            called_major = called[0].upper() if called and called != "NA" else ""

            if called and called != "NA":
                bayes_bin_total += 1
                if expected_major == called_major:
                    bayes_correct += 1

                # ISOGG comparison (only for samples with sufficient coverage)
                if bayes_isogg and isogg_mapper:
                    sample_density = density_map.get(result.sample, 0)
                    if sample_density >= ISOGG_MIN_DENSITY:
                        isogg_gt = isogg_ground_truth.get(result.sample, "")
                        if isogg_gt and isogg_gt not in (".", "..", "NA", "n/a"):
                            isogg_pred = isogg_mapper.to_isogg(result.haplogroup)
                            match = classify_isogg_match(isogg_pred, isogg_gt)
                            bayes_isogg.total += 1
                            if match == "exact":
                                bayes_isogg.exact += 1
                            elif match == "prefix":
                                bayes_isogg.prefix += 1
                            elif match == "major_clade":
                                bayes_isogg.major_clade += 1
                            else:
                                bayes_isogg.mismatch += 1

        bayes_accuracy = bayes_correct / bayes_bin_total if bayes_bin_total > 0 else 0
        bayes_stratified.append(StratifiedResult(
            bin_name=bin_name,
            accuracy=bayes_accuracy,
            correct=bayes_correct,
            total=bayes_bin_total,
            samples_tested=len(samples),
        ))
        bayes_total_correct += bayes_correct
        bayes_total_tested += bayes_bin_total

        elapsed = time.time() - start
        total_runtime += elapsed
        print(f"    Accuracy: {heur_accuracy:.1%} ({heur_correct}/{heur_bin_total}) in {elapsed:.1f}s")

    # Compute overall weighted accuracy
    heur_overall = heur_total_correct / heur_total_tested if heur_total_tested > 0 else 0
    bayes_overall = bayes_total_correct / bayes_total_tested if bayes_total_tested > 0 else 0
    mean_conf = sum(all_confidences) / len(all_confidences) if all_confidences else 0.0
    mean_derived = sum(all_derived) / len(all_derived) if all_derived else 0.0

    print(f"  Overall: {heur_overall:.1%} ({heur_total_correct}/{heur_total_tested})")

    # Print ISOGG metrics if available
    if heur_isogg and heur_isogg.total > 0:
        print(f"  ISOGG compatible: {heur_isogg.compatible_rate:.1f}% "
              f"({heur_isogg.compatible}/{heur_isogg.total})")

    heuristic_result = BenchmarkResult(
        tool_name="yallHap",
        dataset="AADR v54 (stratified)",
        total_samples=heur_total_tested,
        same_major_lineage=heur_total_correct,
        exact_match=0,  # Not computed for stratified
        runtime_seconds=total_runtime,
        mean_confidence=mean_conf,
        mean_derived_snps=mean_derived,
        notes=["GRCh37 reference", "Transversions-only", "Stratified by variant density"],
        stratified=heur_stratified,
        isogg_metrics=heur_isogg,
    )

    bayesian_result = BenchmarkResult(
        tool_name="yallHap-Bayes",
        dataset="AADR v54 (stratified)",
        total_samples=bayes_total_tested,
        same_major_lineage=bayes_total_correct,
        exact_match=0,
        runtime_seconds=total_runtime,
        mean_confidence=mean_conf,
        mean_derived_snps=mean_derived,
        notes=["GRCh37 reference", "Bayesian + transversions-only", "Stratified by variant density"],
        stratified=bayes_stratified,
        isogg_metrics=bayes_isogg,
    )

    return heuristic_result, bayesian_result


def benchmark_gnomad_highcov(
    base_dir: Path,
    subsample: int | None = None,
    threads: int = 1,
    shared_samples: list[str] | None = None,
) -> tuple[BenchmarkResult | None, BenchmarkResult | None]:
    """
    Run gnomAD high-coverage benchmark with heuristic and Bayesian modes.

    Args:
        base_dir: Project base directory
        subsample: Max samples to use
        threads: Parallel threads
        shared_samples: If provided, use these specific samples (for cross-dataset comparison)

    Returns (heuristic_result, bayesian_result)
    """
    highcov_dir = base_dir / "data" / "validation_highcov"
    data_dir = base_dir / "data"
    validation_dir = base_dir / "data" / "validation"

    # Use 1KG Phase 3 ground truth for cross-dataset comparison
    ground_truth_path = validation_dir / "poznik2016_haplogroups.tsv"
    full_vcf_path = highcov_dir / "vcf" / "gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz"

    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = data_dir / "ybrowse_snps.csv"

    # Check required files
    for path in [tree_path, snp_db_path]:
        if not path.exists():
            print(f"  SKIP: {path.name} not found")
            return None, None, None

    if not full_vcf_path.exists():
        print("  SKIP: gnomAD VCF not found")
        return None, None, None

    # Load ground truth from 1KG Phase 3 (same as benchmark_1kg)
    if not ground_truth_path.exists():
        print("  SKIP: Ground truth not found")
        return None, None, None

    ground_truth = load_ground_truth(ground_truth_path)

    # Get available samples from gnomAD VCF
    import pysam

    vcf = pysam.VariantFile(str(full_vcf_path))
    vcf_samples = set(vcf.header.samples)
    vcf.close()

    # Use shared_samples if provided (for cross-dataset comparison)
    if shared_samples:
        samples = [s for s in shared_samples if s in vcf_samples]
    else:
        # Find overlap between ground truth and gnomAD samples
        samples = [s for s in ground_truth.keys() if s in vcf_samples]

        # Subsample if requested
        if subsample and len(samples) > subsample:
            random.seed(42)
            samples = random.sample(samples, subsample)

    # Create or reuse subset VCF with exactly the samples we need
    import subprocess

    # Prefer the comprehensive diagnostic VCF if available (all shared samples, filtered)
    comprehensive_vcf = highcov_dir / "vcf" / "gnomad_1kg_shared_diagnostic.vcf.gz"
    subset_vcf_path = highcov_dir / "vcf" / f"gnomad_subset_{len(samples)}_filtered.vcf.gz"

    if comprehensive_vcf.exists():
        vcf_path = comprehensive_vcf
        print(f"  Using pre-extracted diagnostic VCF: {comprehensive_vcf.name}")
    elif subset_vcf_path.exists():
        vcf_path = subset_vcf_path
        print(f"  Using cached subset: {subset_vcf_path.name}")
    else:
        subset_start = time.time()
        sample_list = ",".join(samples)

        # Create targets file with diagnostic SNP positions (if not exists)
        # Format: chrY<TAB>position (1-based)
        from yallhap.snps import SNPDatabase

        targets_path = highcov_dir / "vcf" / "diagnostic_positions.tsv"
        if not targets_path.exists():
            print(f"  Creating diagnostic positions file...")
            snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)
            positions = set()
            for snp in snp_db:
                pos = snp.get_position("grch38")
                if pos and 2700000 <= pos <= 60000000:
                    positions.add(pos)
            # Write sorted positions
            with open(targets_path, "w") as f:
                for pos in sorted(positions):
                    f.write(f"chrY\t{pos}\n")
            print(f"    Created {len(positions)} target positions")

        desc = f"Extracting {len(samples)} samples + filtering"
        # -T streams through file (fast), -R uses index seeking (slow for many regions)
        bcftools_cmd = [
            "bcftools", "view",
            "-s", sample_list,
            "-T", str(targets_path),
            "-Oz", "-o", str(subset_vcf_path),
            str(full_vcf_path),
        ]

        try:
            run_with_progress(
                bcftools_cmd,
                desc=desc,
                output_file=subset_vcf_path,
                source_file=full_vcf_path,
            )
            print(f"    Indexing subset VCF...")
            subprocess.run(
                ["tabix", "-f", str(subset_vcf_path)],
                check=True,
                capture_output=True,
            )
            vcf_path = subset_vcf_path
            subset_time = time.time() - subset_start
            print(f"  Created {subset_vcf_path.name} ({subset_time:.1f}s)")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"  Warning: Could not create subset VCF: {e}")
            print("  Using full VCF (slow)...")
            vcf_path = full_vcf_path

    print(f"  Loaded {len(samples)} samples with ground truth")

    if not samples:
        print("  SKIP: No samples found in VCF")
        return None, None, None

    # Run heuristic mode
    heur_predictions, heur_runtime, heur_conf, heur_derived, _ = run_yallhap_classification(
        vcf_path=vcf_path,
        tree_path=tree_path,
        snp_db_path=snp_db_path,
        samples=samples,
        reference="grch38",
        bayesian=False,
        threads=threads,
    )

    heur_same, heur_exact, heur_total = compare_haplogroups(ground_truth, heur_predictions)

    heuristic_result = BenchmarkResult(
        tool_name="yallHap",
        dataset="gnomAD HGDP/1KG",
        total_samples=heur_total,
        same_major_lineage=heur_same,
        exact_match=heur_exact,
        runtime_seconds=heur_runtime,
        mean_confidence=heur_conf,
        mean_derived_snps=heur_derived,
        notes=["GRCh38 reference", "Heuristic mode", "High-coverage WGS"],
    )

    # Run Bayesian mode
    bayes_predictions, bayes_runtime, bayes_conf, bayes_derived, _ = run_yallhap_classification(
        vcf_path=vcf_path,
        tree_path=tree_path,
        snp_db_path=snp_db_path,
        samples=samples,
        reference="grch38",
        bayesian=True,
        threads=threads,
    )

    bayes_same, bayes_exact, bayes_total = compare_haplogroups(ground_truth, bayes_predictions)

    bayesian_result = BenchmarkResult(
        tool_name="yallHap-Bayesian",
        dataset="gnomAD HGDP/1KG",
        total_samples=bayes_total,
        same_major_lineage=bayes_same,
        exact_match=bayes_exact,
        runtime_seconds=bayes_runtime,
        mean_confidence=bayes_conf,
        mean_derived_snps=bayes_derived,
        notes=["GRCh38 reference", "Bayesian mode with AD", "High-coverage WGS"],
    )

    # Note: Yleaf comparison runs on 1KG Phase 3 (hg19) separately
    # gnomAD VCF format causes issues with Yleaf's parsing

    return heuristic_result, bayesian_result


def run_yleaf(
    vcf_path: Path,
    samples: list[str],
    reference: str = "grch38",
) -> dict[str, str] | None:
    """
    Run Yleaf on samples and return predictions.

    Returns None if Yleaf is not installed.
    """
    import shutil
    import subprocess
    import tempfile

    # Check if Yleaf is installed
    yleaf_path = shutil.which("Yleaf")
    if not yleaf_path:
        return None

    predictions = {}

    # Yleaf processes one sample at a time
    for sample in samples:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Extract single sample VCF
            sample_vcf = Path(tmpdir) / f"{sample}.vcf.gz"
            try:
                subprocess.run(
                    ["bcftools", "view", "-s", sample, "-Oz", "-o", str(sample_vcf), str(vcf_path)],
                    check=True,
                    capture_output=True,
                )
                subprocess.run(
                    ["tabix", "-f", str(sample_vcf)],
                    check=True,
                    capture_output=True,
                )

                # Run Yleaf
                output_dir = Path(tmpdir) / "yleaf_out"
                result = subprocess.run(
                    [
                        "Yleaf",
                        "-vcf", str(sample_vcf),
                        "-o", str(output_dir),
                        "-rg", "hg38" if reference == "grch38" else "hg19",
                    ],
                    capture_output=True,
                    timeout=30,  # 30s should be plenty for a single chrY VCF
                )

                # Parse Yleaf output - hg_prediction.hg
                # Format: Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score...
                hg_file = output_dir / "hg_prediction.hg"
                if hg_file.exists():
                    with open(hg_file) as f:
                        for line in f:
                            if line.startswith("Sample_name"):
                                continue  # Skip header
                            parts = line.strip().split("\t")
                            if len(parts) >= 2:
                                predictions[sample] = parts[1]
                                break
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
                continue

    return predictions if predictions else None


def benchmark_yleaf(
    base_dir: Path,
    vcf_path: Path,
    ground_truth: dict[str, str],
    samples: list[str],
    reference: str = "grch38",
    dataset_name: str = "gnomAD HGDP/1KG",
) -> BenchmarkResult | None:
    """Run Yleaf benchmark."""
    start_time = time.time()
    predictions = run_yleaf(vcf_path, samples, reference=reference)
    runtime = time.time() - start_time

    if predictions is None:
        return None

    same_major, exact, total = compare_haplogroups(ground_truth, predictions)

    return BenchmarkResult(
        tool_name="Yleaf",
        dataset=dataset_name,
        total_samples=total,
        same_major_lineage=same_major,
        exact_match=exact,
        runtime_seconds=runtime,
        notes=[f"{reference.upper()} reference", "External tool comparison"],
    )


def benchmark_yleaf_1kg(
    base_dir: Path,
    subsample: int | None = None,
    shared_samples: list[str] | None = None,
) -> BenchmarkResult | None:
    """Run Yleaf benchmark on 1000 Genomes Phase 3 (hg19)."""
    validation_dir = base_dir / "data" / "validation"

    ground_truth_path = validation_dir / "poznik2016_haplogroups.tsv"
    vcf_path = validation_dir / "1kg_chrY_phase3.vcf.gz"

    if not ground_truth_path.exists() or not vcf_path.exists():
        return None

    ground_truth = load_ground_truth(ground_truth_path)

    if shared_samples:
        samples = [s for s in shared_samples if s in ground_truth]
    else:
        samples = list(ground_truth.keys())
        if subsample and len(samples) > subsample:
            random.seed(42)
            samples = random.sample(samples, subsample)

    return benchmark_yleaf(
        base_dir=base_dir,
        vcf_path=vcf_path,
        ground_truth={s: ground_truth[s] for s in samples if s in ground_truth},
        samples=samples,
        reference="grch37",
        dataset_name="1KG Phase 3",
    )


def print_results_table(results: list[BenchmarkResult]) -> None:
    """Print formatted results table."""
    # Check if any results have ISOGG metrics
    has_isogg = any(r.isogg_metrics and r.isogg_metrics.total > 0 for r in results)

    print("\n" + "=" * 110)
    print("BENCHMARK RESULTS")
    print("=" * 110)
    print()

    if has_isogg:
        print(f"{'Dataset':<20} {'Tool':<18} {'Samples':>8} {'Major':>8} {'Exact':>8} {'ISOGG':>10} {'Conf':>6}")
        print("-" * 110)
    else:
    print(f"{'Dataset':<20} {'Tool':<18} {'Samples':>8} {'Major Match':>12} {'Exact':>8} {'Conf':>6} {'SNPs':>6}")
        print("-" * 110)

    for r in results:
        conf_str = f"{r.mean_confidence:.2f}" if r.mean_confidence else "N/A"

        if has_isogg:
            if r.isogg_metrics and r.isogg_metrics.total > 0:
                isogg_str = f"{r.isogg_metrics.compatible_rate:.1f}%"
            else:
                isogg_str = "N/A"
            print(
                f"{r.dataset:<20} {r.tool_name:<18} {r.total_samples:>8} "
                f"{r.major_lineage_rate:>6.1f}% {r.exact_match_rate:>6.1f}% "
                f"{isogg_str:>10} {conf_str:>6}"
            )
        else:
        snps_str = f"{r.mean_derived_snps:.1f}" if r.mean_derived_snps else "N/A"
        print(
            f"{r.dataset:<20} {r.tool_name:<18} {r.total_samples:>8} "
            f"{r.major_lineage_rate:>10.2f}% {r.exact_match_rate:>6.2f}% "
            f"{conf_str:>6} {snps_str:>6}"
        )

    print("-" * 110)

    # Print ISOGG legend if applicable
    if has_isogg:
        print("\nISOGG column = Compatible rate (exact + prefix matches)")
        print("  - Exact: Identical ISOGG haplogroups")
        print("  - Prefix: One is ancestor of the other (e.g., R1b vs R1b1a1b)")


def save_results(base_dir: Path, results: list[BenchmarkResult]) -> None:
    """Save results to JSON file."""
    results_dir = base_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # JSON output
    json_path = results_dir / "benchmark_results.json"
    json_data: dict = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "results": [],
    }

    for r in results:
        result_dict: dict = {
            "dataset": r.dataset,
            "tool": r.tool_name,
            "total_samples": r.total_samples,
            "same_major_lineage": r.same_major_lineage,
            "same_major_lineage_pct": round(r.major_lineage_rate, 2),
            "exact_match": r.exact_match,
            "exact_match_pct": round(r.exact_match_rate, 2),
            "runtime_seconds": round(r.runtime_seconds, 1),
            "mean_confidence": round(r.mean_confidence, 3) if r.mean_confidence else None,
            "mean_derived_snps": round(r.mean_derived_snps, 1) if r.mean_derived_snps else None,
            "notes": r.notes,
        }

        # Add stratified results if present
        if r.stratified:
            result_dict["stratified"] = [
                {
                    "bin": s.bin_name,
                    "accuracy_pct": round(s.accuracy * 100, 1),
                    "correct": s.correct,
                    "total": s.total,
                    "samples_tested": s.samples_tested,
                }
                for s in r.stratified
            ]

        # Add ISOGG metrics if present
        if r.isogg_metrics and r.isogg_metrics.total > 0:
            result_dict["isogg"] = {
                "total": r.isogg_metrics.total,
                "exact": r.isogg_metrics.exact,
                "exact_pct": round(r.isogg_metrics.exact_rate, 2),
                "prefix": r.isogg_metrics.prefix,
                "compatible": r.isogg_metrics.compatible,
                "compatible_pct": round(r.isogg_metrics.compatible_rate, 2),
                "major_clade": r.isogg_metrics.major_clade,
                "mismatch": r.isogg_metrics.mismatch,
            }

        json_data["results"].append(result_dict)

    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"\nResults saved to: {json_path}")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run yallHap benchmark suite across validation datasets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Full benchmark (all samples)
    python scripts/run_benchmarks.py

    # Quick benchmark with 50 samples per dataset
    python scripts/run_benchmarks.py --subsample 50

    # Fast benchmark with 50 samples and 16 threads
    python scripts/run_benchmarks.py --subsample 50 --threads 16
""",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=None,
        metavar="N",
        help="Limit each dataset to N samples (default: use all)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        metavar="N",
        help="Number of parallel threads for classification (default: 1)",
    )
    parser.add_argument(
        "--gnomad-samples",
        type=int,
        default=10,
        metavar="N",
        help="Number of gnomAD high-coverage samples, balanced by superpopulation (default: 10)",
    )
    parser.add_argument(
        "--aadr-samples-per-bin",
        type=int,
        default=300,
        metavar="N",
        help="Max samples per variant density bin for AADR stratified analysis (default: 300)",
    )
    parser.add_argument(
        "--skip-1kg",
        action="store_true",
        help="Skip 1000 Genomes Phase 3 benchmark",
    )
    parser.add_argument(
        "--skip-aadr",
        action="store_true",
        help="Skip AADR ancient DNA benchmark",
    )
    parser.add_argument(
        "--skip-gnomad",
        action="store_true",
        help="Skip gnomAD high-coverage benchmark",
    )
    # Enhanced ancient mode options
    parser.add_argument(
        "--estimate-contamination",
        action="store_true",
        help="Estimate Y-chromosome contamination for ancient samples",
    )
    parser.add_argument(
        "--max-tolerance",
        type=int,
        default=3,
        metavar="N",
        help="Max ancestral calls for path traversal (default: 3)",
    )
    parser.add_argument(
        "--isogg",
        action="store_true",
        help="Output ISOGG haplogroup nomenclature alongside YFull",
    )
    parser.add_argument(
        "--isogg-db",
        type=str,
        default=None,
        metavar="PATH",
        help="Path to ISOGG SNP database (default: data/isogg_snps_grch38.txt)",
    )
    parser.add_argument(
        "--create-shared-diagnostic",
        action="store_true",
        help="Create the gnomAD shared diagnostic VCF and exit (extracts 1KG-shared samples at diagnostic positions)",
    )
    return parser.parse_args()


def create_shared_diagnostic_vcf(base_dir: Path) -> int:
    """
    Create the gnomAD shared diagnostic VCF containing all 1KG-shared samples
    at diagnostic SNP positions.

    This creates a pre-filtered VCF that can be reused for benchmarks without
    re-extracting from the full gnomAD VCF each time.

    Args:
        base_dir: Project base directory

    Returns:
        Exit code (0 = success, 1 = error)
    """
    import subprocess

    import pysam

    from yallhap.snps import SNPDatabase

    validation_dir = base_dir / "data" / "validation"
    highcov_dir = base_dir / "data" / "validation_highcov"
    data_dir = base_dir / "data"

    ground_truth_path = validation_dir / "poznik2016_haplogroups.tsv"
    full_vcf_path = highcov_dir / "vcf" / "gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz"
    snp_db_path = data_dir / "ybrowse_snps.csv"
    output_vcf_path = highcov_dir / "vcf" / "gnomad_1kg_shared_diagnostic.vcf.gz"

    print("=" * 70)
    print("Creating gnomAD Shared Diagnostic VCF")
    print("=" * 70)
    print()

    # Check required files
    print("Checking required files...")
    missing = []
    if not ground_truth_path.exists():
        missing.append(f"  - Ground truth: {ground_truth_path}")
    if not full_vcf_path.exists():
        missing.append(f"  - gnomAD VCF: {full_vcf_path}")
    if not snp_db_path.exists():
        missing.append(f"  - SNP database: {snp_db_path}")

    if missing:
        print("ERROR: Missing required files:")
        for m in missing:
            print(m)
        print("\nPlease download the required data first:")
        print("  python scripts/download_1kg_validation.py")
        print("  python scripts/download_gnomad_highcov.py")
        return 1

    print(f"  ✓ Ground truth: {ground_truth_path.name}")
    print(f"  ✓ gnomAD VCF: {full_vcf_path.name}")
    print(f"  ✓ SNP database: {snp_db_path.name}")
    print()

    # Load ground truth samples
    print("Loading ground truth samples...")
    ground_truth = load_ground_truth(ground_truth_path)
    print(f"  Found {len(ground_truth)} samples with haplogroup assignments")

    # Find samples in gnomAD VCF
    print("Finding samples in gnomAD VCF...")
    vcf = pysam.VariantFile(str(full_vcf_path))
    gnomad_samples = set(vcf.header.samples)
    vcf.close()
    print(f"  gnomAD VCF contains {len(gnomad_samples)} samples")

    # Find overlap
    shared_samples = [s for s in ground_truth.keys() if s in gnomad_samples]
    print(f"  Overlap (1KG males in gnomAD): {len(shared_samples)} samples")
    print()

    if not shared_samples:
        print("ERROR: No shared samples found between ground truth and gnomAD VCF")
        return 1

    # Create targets file with diagnostic SNP positions
    targets_path = highcov_dir / "vcf" / "diagnostic_positions.tsv"
    print("Creating diagnostic positions file...")
    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)
    positions = set()
    for snp in snp_db:
        pos = snp.get_position("grch38")
        if pos and 2700000 <= pos <= 60000000:  # Y chromosome bounds
            positions.add(pos)

    with open(targets_path, "w") as f:
        for pos in sorted(positions):
            f.write(f"chrY\t{pos}\n")
    print(f"  Created {len(positions)} target positions at: {targets_path.name}")
    print()

    # Extract samples and filter to diagnostic positions
    print(f"Extracting {len(shared_samples)} samples at {len(positions)} diagnostic positions...")
    print("  This may take several minutes for the full gnomAD VCF...")
    print()

    sample_list = ",".join(shared_samples)
    bcftools_cmd = [
        "bcftools", "view",
        "-s", sample_list,
        "-T", str(targets_path),
        "-Oz", "-o", str(output_vcf_path),
        str(full_vcf_path),
    ]

    start_time = time.time()
    try:
        run_with_progress(
            bcftools_cmd,
            desc="Extracting samples + filtering",
            output_file=output_vcf_path,
            source_file=full_vcf_path,
        )
    except subprocess.CalledProcessError as e:
        print(f"ERROR: bcftools failed: {e}")
        return 1
    except FileNotFoundError:
        print("ERROR: bcftools not found. Please install bcftools.")
        return 1

    extraction_time = time.time() - start_time

    # Index the output VCF
    print("Indexing output VCF...")
    try:
        subprocess.run(
            ["tabix", "-f", str(output_vcf_path)],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"ERROR: tabix failed: {e}")
        return 1
    except FileNotFoundError:
        print("ERROR: tabix not found. Please install htslib/tabix.")
        return 1

    total_time = time.time() - start_time

    # Report results
    output_size_mb = output_vcf_path.stat().st_size / (1024 * 1024)
    index_path = Path(str(output_vcf_path) + ".tbi")
    index_size_kb = index_path.stat().st_size / 1024 if index_path.exists() else 0

    print()
    print("=" * 70)
    print("SUCCESS: Shared diagnostic VCF created")
    print("=" * 70)
    print()
    print("Output files:")
    print(f"  VCF:   {output_vcf_path}")
    print(f"         ({output_size_mb:.1f} MB)")
    print(f"  Index: {index_path}")
    print(f"         ({index_size_kb:.1f} KB)")
    print()
    print("Contents:")
    print(f"  Samples:   {len(shared_samples)} (1KG Phase 3 males in gnomAD)")
    print(f"  Positions: {len(positions)} diagnostic SNPs (GRCh38)")
    print()
    print(f"Time: {total_time:.1f}s ({total_time/60:.1f}m)")
    print()
    print("This VCF will be automatically used by future benchmark runs.")
    print("To re-run benchmarks with gnomAD high-coverage data:")
    print("  python scripts/run_benchmarks.py --gnomad-samples 100 --threads 8")

    return 0


def find_shared_samples(
    base_dir: Path,
    subsample: int | None = None,
    gnomad_samples: int = 10,
) -> tuple[list[str], list[str]]:
    """
    Find samples for 1KG and gnomAD benchmarks.

    For gnomAD, selects samples balanced by superpopulation from those
    that exist in both 1KG Phase 3 ground truth AND gnomAD VCF.

    Args:
        base_dir: Project base directory
        subsample: Max samples for 1KG benchmark (None = all)
        gnomad_samples: Number of gnomAD samples, balanced by superpopulation

    Returns:
        (1kg_samples, gnomad_shared_samples)
    """
    import pysam

    validation_dir = base_dir / "data" / "validation"
    highcov_dir = base_dir / "data" / "validation_highcov"

    ground_truth_path = validation_dir / "poznik2016_haplogroups.tsv"
    panel_path = validation_dir / "1kg_sample_panel.txt"
    gnomad_vcf_path = highcov_dir / "vcf" / "gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz"

    if not ground_truth_path.exists():
        return [], []

    # Load ground truth samples
    ground_truth = load_ground_truth(ground_truth_path)

    # Load population info
    sample_to_superpop: dict[str, str] = {}
    if panel_path.exists():
        with open(panel_path) as f:
            for line in f:
                if line.startswith("sample"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4 and parts[3] == "male":
                    sample_to_superpop[parts[0]] = parts[2]  # super_pop

    # 1KG samples (all males with ground truth, optionally subsampled)
    kg_samples = list(ground_truth.keys())
    if subsample and len(kg_samples) > subsample:
        random.seed(42)
        kg_samples = random.sample(kg_samples, subsample)

    # gnomAD shared samples - balanced by superpopulation
    gnomad_shared: list[str] = []
    if gnomad_vcf_path.exists():
        vcf = pysam.VariantFile(str(gnomad_vcf_path))
        gnomad_vcf_samples = set(vcf.header.samples)
        vcf.close()

        # Find overlap with ground truth
        shared = [s for s in ground_truth.keys() if s in gnomad_vcf_samples]

        # Group by superpopulation
        superpops = ["AFR", "AMR", "EAS", "EUR", "SAS"]
        by_superpop: dict[str, list[str]] = {sp: [] for sp in superpops}
        for s in shared:
            sp = sample_to_superpop.get(s, "")
            if sp in by_superpop:
                by_superpop[sp].append(s)

        # Select ~equal from each superpopulation
        per_pop = max(1, gnomad_samples // len(superpops))
        remainder = gnomad_samples % len(superpops)

        random.seed(42)
        for i, sp in enumerate(superpops):
            n = per_pop + (1 if i < remainder else 0)
            available = by_superpop[sp]
            if available:
                gnomad_shared.extend(random.sample(available, min(n, len(available))))

    return kg_samples, gnomad_shared


def main() -> int:
    """Run all benchmarks."""
    args = parse_args()
    base_dir = Path(__file__).parent.parent

    # Handle --create-shared-diagnostic early exit
    if args.create_shared_diagnostic:
        return create_shared_diagnostic_vcf(base_dir)

    total_start = time.time()

    print("=" * 70)
    print("yallHap Comprehensive Benchmark Suite")
    print("=" * 70)
    if args.subsample:
        print(f"  Subsampling: {args.subsample} samples per dataset")
    if args.threads > 1:
        print(f"  Threads: {args.threads}")

    results: list[BenchmarkResult] = []

    # Find samples for benchmarks
    print("\nFinding samples for benchmarks...")
    kg_samples, gnomad_samples = find_shared_samples(
        base_dir,
        subsample=args.subsample,
        gnomad_samples=args.gnomad_samples,
    )
    if kg_samples:
        print(f"  1KG Phase 3: {len(kg_samples)} samples")
    if gnomad_samples:
        print(f"  gnomAD (shared with 1KG): {len(gnomad_samples)} samples, balanced by superpopulation")

    # 1. 1000 Genomes Phase 3 (GRCh37)
    if args.skip_1kg:
        print("\n[1/3] 1000 Genomes Phase 3 (GRCh37)... SKIPPED")
    else:
        print("\n[1/3] 1000 Genomes Phase 3 (GRCh37)...")
        step_start = time.time()
        heur_result, bayes_result = benchmark_1kg(
            base_dir,
            subsample=args.subsample,
            threads=args.threads,
            shared_samples=kg_samples if kg_samples else None,
        )
        step_time = time.time() - step_start
        if heur_result:
            results.append(heur_result)
            print(f"  + Heuristic: {heur_result.major_lineage_rate:.2f}% major lineage accuracy")
        if bayes_result:
            results.append(bayes_result)
            print(f"  + Bayesian: {bayes_result.major_lineage_rate:.2f}% major lineage accuracy")
        print(f"  Time: {step_time:.1f}s")

    # 2. AADR Ancient DNA (stratified by coverage)
    if args.skip_aadr:
        print("\n[2/3] AADR Ancient DNA (stratified by coverage)... SKIPPED")
    else:
        print("\n[2/3] AADR Ancient DNA (stratified by coverage)...")
        if args.estimate_contamination:
            print("  Enhanced ancient mode: contamination estimation enabled")
        if args.isogg:
            print("  Enhanced ancient mode: ISOGG output enabled")
        step_start = time.time()
        heur_result, bayes_result = benchmark_aadr(
            base_dir,
            samples_per_bin=args.aadr_samples_per_bin,
            threads=args.threads,
            estimate_contamination=args.estimate_contamination,
            max_tolerance=args.max_tolerance,
            isogg=args.isogg,
            isogg_db_path=Path(args.isogg_db) if args.isogg_db else None,
        )
        step_time = time.time() - step_start
        if heur_result:
            results.append(heur_result)
            print(f"  + Heuristic: {heur_result.major_lineage_rate:.2f}% major lineage accuracy (weighted)")
            if heur_result.stratified:
                for s in heur_result.stratified:
                    print(f"      {s.bin_name}: {s.accuracy*100:.1f}% ({s.correct}/{s.total})")
        if bayes_result:
            results.append(bayes_result)
            print(f"  + Bayesian: {bayes_result.major_lineage_rate:.2f}% major lineage accuracy (weighted)")
        print(f"  Time: {step_time:.1f}s")

    # 3. gnomAD High-Coverage (GRCh38) - use balanced superpopulation samples
    if args.skip_gnomad:
        print("\n[3/3] gnomAD HGDP/1000G High-Coverage (GRCh38)... SKIPPED")
    else:
        print("\n[3/3] gnomAD HGDP/1000G High-Coverage (GRCh38)...")
        step_start = time.time()
        heur_result, bayes_result = benchmark_gnomad_highcov(
            base_dir,
            subsample=None,  # Don't subsample further, use gnomad_samples list
            threads=args.threads,
            shared_samples=gnomad_samples if gnomad_samples else None,
        )
        step_time = time.time() - step_start
        if heur_result:
            results.append(heur_result)
            print(f"  + Heuristic: {heur_result.major_lineage_rate:.2f}% major lineage accuracy")
        if bayes_result:
            results.append(bayes_result)
            print(f"  + Bayesian: {bayes_result.major_lineage_rate:.2f}% major lineage accuracy")
        print(f"  Time: {step_time:.1f}s")


    if not results:
        print("\nNo benchmarks completed. Please download validation data first:")
        print("  python scripts/download_1kg_validation.py")
        print("  python scripts/download_ancient_test_data.py")
        print("  python scripts/download_gnomad_highcov.py")
        return 1

    # Print and save results
    print_results_table(results)
    save_results(base_dir, results)

    total_time = time.time() - total_start
    print(f"\n⏱ Total benchmark time: {total_time:.1f}s ({total_time/60:.1f}m)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
