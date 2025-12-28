#!/usr/bin/env python3
"""
Comprehensive ISOGG nomenclature validation.

Compares yallHap's ISOGG output against ground truth ISOGG annotations,
using fuzzy matching to assess quality:
- Exact match: Identical strings
- Prefix match: One is a prefix of the other (e.g., R1b vs R1b1a1b)
- Major clade match: Same first letter (e.g., R1b vs R1a)
- Mismatch: Different major clades (e.g., I2 vs Q)

Usage:
    python scripts/validate_isogg.py [--samples N] [--threads N]

Output:
    results/isogg_validation_full.json - Detailed results
    results/isogg_validation_full.md - Human-readable report
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from yallhap.classifier import HaplogroupClassifier
from yallhap.isogg import ISOGGDatabase, ISOGGMapper
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


@dataclass
class ISOGGComparison:
    """Result of comparing yallHap ISOGG to ground truth ISOGG."""

    sample: str
    yfull_predicted: str
    isogg_predicted: str
    isogg_ground_truth: str
    match_type: str  # exact, prefix, major_clade, mismatch
    confidence: float


@dataclass
class ISOGGValidationResult:
    """Summary of ISOGG validation."""

    total_samples: int
    exact_matches: int
    prefix_matches: int
    major_clade_matches: int
    mismatches: int
    no_ground_truth: int
    runtime_seconds: float
    comparisons: list[ISOGGComparison] = field(default_factory=list)

    @property
    def exact_rate(self) -> float:
        return self.exact_matches / self.total_samples * 100 if self.total_samples > 0 else 0

    @property
    def prefix_rate(self) -> float:
        return self.prefix_matches / self.total_samples * 100 if self.total_samples > 0 else 0

    @property
    def major_clade_rate(self) -> float:
        return self.major_clade_matches / self.total_samples * 100 if self.total_samples > 0 else 0

    @property
    def mismatch_rate(self) -> float:
        return self.mismatches / self.total_samples * 100 if self.total_samples > 0 else 0

    @property
    def compatible_rate(self) -> float:
        """Exact + prefix matches (compatible assignments)."""
        return (
            (self.exact_matches + self.prefix_matches) / self.total_samples * 100
            if self.total_samples > 0
            else 0
        )


def load_aadr_ground_truth(tsv_path: Path) -> dict[str, tuple[str, str]]:
    """
    Load AADR ground truth with both YFull and ISOGG haplogroups.

    Returns:
        Dict of sample_id -> (yfull_haplogroup, isogg_haplogroup)
    """
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
            if len(parts) >= 3:
                sample_id = parts[0]
                yfull_hg = parts[1]
                isogg_hg = parts[2]
                # Only include samples with valid ISOGG annotations
                if isogg_hg and isogg_hg not in (".", "..", "NA", "n/a", ""):
                    ground_truth[sample_id] = (yfull_hg, isogg_hg)
    return ground_truth


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
    # First letter
    return hg[0] if hg[0].isalpha() else ""


def is_prefix_match(hg1: str, hg2: str) -> bool:
    """Check if one haplogroup is a prefix of the other."""
    h1 = normalize_isogg(hg1)
    h2 = normalize_isogg(hg2)
    if not h1 or not h2:
        return False
    return h1.startswith(h2) or h2.startswith(h1)


def classify_match(predicted: str, ground_truth: str) -> str:
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
    if is_prefix_match(pred_norm, gt_norm):
        return "prefix"

    # Major clade match (same first letter)
    if get_major_clade(pred_norm) == get_major_clade(gt_norm):
        return "major_clade"

    return "mismatch"


def run_validation(
    vcf_path: Path,
    tree: Tree,
    snp_db: SNPDatabase,
    isogg_mapper: ISOGGMapper,
    ground_truth: dict[str, tuple[str, str]],
    sample_ids: list[str],
    reference: str = "grch37",
    threads: int = 16,
) -> ISOGGValidationResult:
    """
    Run ISOGG validation on samples.
    """
    start_time = time.time()

    # Create classifier
    classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference=reference,  # type: ignore
    )

    # Classify all samples
    print(f"  Classifying {len(sample_ids)} samples...")
    results = classifier.classify_batch(vcf_path, sample_ids, threads=threads)
    result_map = {r.sample: r for r in results}

    # Compare ISOGG outputs
    comparisons: list[ISOGGComparison] = []
    exact = prefix = major = mismatch = no_gt = 0

    for sample in sample_ids:
        if sample not in result_map:
            no_gt += 1
            continue

        result = result_map[sample]
        gt_yfull, gt_isogg = ground_truth.get(sample, ("", ""))

        if not gt_isogg:
            no_gt += 1
            continue

        # Get ISOGG mapping from yallHap
        isogg_predicted = isogg_mapper.to_isogg(result.haplogroup)

        # Classify match
        match_type = classify_match(isogg_predicted, gt_isogg)

        comparisons.append(
            ISOGGComparison(
                sample=sample,
                yfull_predicted=result.haplogroup,
                isogg_predicted=isogg_predicted,
                isogg_ground_truth=gt_isogg,
                match_type=match_type,
                confidence=result.confidence,
            )
        )

        if match_type == "exact":
            exact += 1
        elif match_type == "prefix":
            prefix += 1
        elif match_type == "major_clade":
            major += 1
        else:
            mismatch += 1

    runtime = time.time() - start_time

    return ISOGGValidationResult(
        total_samples=len(comparisons),
        exact_matches=exact,
        prefix_matches=prefix,
        major_clade_matches=major,
        mismatches=mismatch,
        no_ground_truth=no_gt,
        runtime_seconds=runtime,
        comparisons=comparisons,
    )


def generate_report(result: ISOGGValidationResult, output_path: Path) -> None:
    """Generate human-readable markdown report."""
    lines = [
        "# ISOGG Validation Report",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Summary",
        "",
        f"| Metric | Count | Percentage |",
        f"|--------|-------|------------|",
        f"| Total samples | {result.total_samples} | 100.0% |",
        f"| Exact match | {result.exact_matches} | {result.exact_rate:.1f}% |",
        f"| Prefix match | {result.prefix_matches} | {result.prefix_rate:.1f}% |",
        f"| Major clade match | {result.major_clade_matches} | {result.major_clade_rate:.1f}% |",
        f"| Mismatch | {result.mismatches} | {result.mismatch_rate:.1f}% |",
        f"| No ground truth | {result.no_ground_truth} | - |",
        "",
        f"**Compatible (exact + prefix):** {result.exact_matches + result.prefix_matches} ({result.compatible_rate:.1f}%)",
        "",
        f"Runtime: {result.runtime_seconds:.1f} seconds",
        "",
        "## Match Type Definitions",
        "",
        "- **Exact match**: Normalized ISOGG strings are identical",
        "- **Prefix match**: One is a prefix of the other (e.g., R1b vs R1b1a1b1a1a)",
        "- **Major clade match**: Same first letter but different subclades (e.g., R1a vs R1b)",
        "- **Mismatch**: Different major clades (e.g., I2 vs R1b)",
        "",
        "## Sample Results by Match Type",
        "",
    ]

    # Group by match type
    for match_type in ["exact", "prefix", "major_clade", "mismatch"]:
        samples = [c for c in result.comparisons if c.match_type == match_type]
        lines.append(f"### {match_type.replace('_', ' ').title()} ({len(samples)} samples)")
        lines.append("")

        if samples:
            lines.append("| Sample | YFull | ISOGG (yallHap) | ISOGG (Ground Truth) | Conf |")
            lines.append("|--------|-------|-----------------|----------------------|------|")

            # Show up to 20 examples per category
            for c in samples[:20]:
                lines.append(
                    f"| {c.sample} | {c.yfull_predicted} | {c.isogg_predicted} | "
                    f"{c.isogg_ground_truth} | {c.confidence:.3f} |"
                )

            if len(samples) > 20:
                lines.append(f"| ... | ({len(samples) - 20} more) | ... | ... | ... |")

        lines.append("")

    # Analysis of mismatches
    if result.mismatches > 0:
        lines.append("## Mismatch Analysis")
        lines.append("")
        lines.append("Major clade transitions in mismatches:")
        lines.append("")

        # Count transitions
        transitions: dict[str, int] = {}
        for c in result.comparisons:
            if c.match_type == "mismatch":
                pred_major = get_major_clade(c.isogg_predicted)
                gt_major = get_major_clade(c.isogg_ground_truth)
                key = f"{gt_major} → {pred_major}"
                transitions[key] = transitions.get(key, 0) + 1

        lines.append("| Ground Truth Major | Predicted Major | Count |")
        lines.append("|-------------------|-----------------|-------|")
        for trans, count in sorted(transitions.items(), key=lambda x: -x[1])[:20]:
            lines.append(f"| {trans.split(' → ')[0]} | {trans.split(' → ')[1]} | {count} |")
        lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def calculate_variant_density(
    vcf_path: Path,
    sample_ids: list[str],
    cache_path: Path | None = None,
) -> dict[str, float]:
    """
    Calculate variant density for each sample from the VCF.

    Variant density = (called variants / total variants) * 100

    Returns:
        Dictionary mapping sample_id -> density percentage (0-100)
    """
    import pysam

    # Check cache
    if cache_path and cache_path.exists():
        vcf_mtime = vcf_path.stat().st_mtime
        cache_mtime = cache_path.stat().st_mtime
        if cache_mtime > vcf_mtime:
            print("    Loading cached variant density...")
            with open(cache_path) as f:
                cached = json.load(f)
            return {s: cached[s] for s in sample_ids if s in cached}

    print("    Calculating variant density (this may take a few minutes)...")

    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)
    target_samples = [s for s in sample_ids if s in vcf_samples]

    called_counts: dict[str, int] = {s: 0 for s in target_samples}
    total_variants = 0

    for rec in vcf.fetch():
        total_variants += 1
        for sample in target_samples:
            gt = rec.samples[sample]["GT"]
            if gt is not None and gt != (None,) and gt != (None, None):
                called_counts[sample] += 1

    vcf.close()

    density_map: dict[str, float] = {}
    for sample, called in called_counts.items():
        density_map[sample] = (called / total_variants * 100) if total_variants > 0 else 0.0

    # Cache results
    if cache_path:
        all_densities = {}
        if cache_path.exists():
            with open(cache_path) as f:
                all_densities = json.load(f)
        all_densities.update(density_map)
        with open(cache_path, "w") as f:
            json.dump(all_densities, f)

    return density_map


def main() -> None:
    """Run ISOGG validation."""
    parser = argparse.ArgumentParser(description="Validate ISOGG nomenclature output")
    parser.add_argument(
        "--samples",
        type=int,
        default=None,
        help="Maximum samples to test (default: all)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads (default: 16)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for subsampling (default: 42)",
    )
    parser.add_argument(
        "--min-density",
        type=float,
        default=4.0,
        help="Minimum variant density %% to include sample (default: 4.0)",
    )
    parser.add_argument(
        "--no-density-filter",
        action="store_true",
        help="Disable variant density filtering (include all samples)",
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
    isogg_path = data_dir / "isogg_snps_grch38.txt"
    ground_truth_path = ancient_dir / "aadr_1240k_ground_truth.tsv"

    output_json = results_dir / "isogg_validation_full.json"
    output_md = results_dir / "isogg_validation_full.md"

    print("=" * 70)
    print("ISOGG Nomenclature Validation")
    print("=" * 70)

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

    ground_truth = load_aadr_ground_truth(ground_truth_path)
    print(f"  Ground truth with ISOGG: {len(ground_truth)} samples")

    # Get sample list
    import pysam

    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)
    vcf.close()

    # Filter to samples in VCF with ISOGG ground truth
    valid_samples = [s for s in ground_truth.keys() if s in vcf_samples]
    print(f"  Valid samples in VCF: {len(valid_samples)}")

    # Calculate variant density and filter
    if not args.no_density_filter:
        cache_path = ancient_dir / "aadr_variant_density.json"
        density_map = calculate_variant_density(vcf_path, valid_samples, cache_path)

        # Filter by minimum density
        before_count = len(valid_samples)
        valid_samples = [s for s in valid_samples if density_map.get(s, 0) >= args.min_density]
        print(f"  After density filter (>={args.min_density}%): {len(valid_samples)} "
              f"(removed {before_count - len(valid_samples)} low-coverage samples)")

    # Subsample if requested
    if args.samples and len(valid_samples) > args.samples:
        import random

        random.seed(args.seed)
        valid_samples = random.sample(valid_samples, args.samples)
        print(f"  Subsampled to: {len(valid_samples)}")

    # Run validation
    print("\nRunning validation...")
    result = run_validation(
        vcf_path=vcf_path,
        tree=tree,
        snp_db=snp_db,
        isogg_mapper=isogg_mapper,
        ground_truth=ground_truth,
        sample_ids=valid_samples,
        reference="grch37",
        threads=args.threads,
    )

    # Print summary
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"\nTotal samples: {result.total_samples}")
    print(f"\nMatch breakdown:")
    print(f"  Exact match:       {result.exact_matches:>5} ({result.exact_rate:>5.1f}%)")
    print(f"  Prefix match:      {result.prefix_matches:>5} ({result.prefix_rate:>5.1f}%)")
    print(f"  Major clade match: {result.major_clade_matches:>5} ({result.major_clade_rate:>5.1f}%)")
    print(f"  Mismatch:          {result.mismatches:>5} ({result.mismatch_rate:>5.1f}%)")
    print(f"\nCompatible (exact + prefix): {result.compatible_rate:.1f}%")
    print(f"Runtime: {result.runtime_seconds:.1f} seconds")

    # Save JSON results
    results_dir.mkdir(exist_ok=True)

    json_data = {
        "timestamp": datetime.now().isoformat(),
        "min_density_filter": None if args.no_density_filter else args.min_density,
        "total_samples": result.total_samples,
        "exact_matches": result.exact_matches,
        "exact_rate": round(result.exact_rate, 2),
        "prefix_matches": result.prefix_matches,
        "prefix_rate": round(result.prefix_rate, 2),
        "major_clade_matches": result.major_clade_matches,
        "major_clade_rate": round(result.major_clade_rate, 2),
        "mismatches": result.mismatches,
        "mismatch_rate": round(result.mismatch_rate, 2),
        "compatible_rate": round(result.compatible_rate, 2),
        "runtime_seconds": round(result.runtime_seconds, 2),
        "comparisons": [
            {
                "sample": c.sample,
                "yfull_predicted": c.yfull_predicted,
                "isogg_predicted": c.isogg_predicted,
                "isogg_ground_truth": c.isogg_ground_truth,
                "match_type": c.match_type,
                "confidence": round(c.confidence, 4),
            }
            for c in result.comparisons
        ],
    }

    with open(output_json, "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"\nJSON results: {output_json}")

    # Generate markdown report
    generate_report(result, output_md)
    print(f"Markdown report: {output_md}")


if __name__ == "__main__":
    main()

