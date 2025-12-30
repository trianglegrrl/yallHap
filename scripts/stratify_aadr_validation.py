#!/usr/bin/env python3
"""
Stratify AADR validation by coverage and publication.

Reads the AADR annotation file to extract coverage and publication info,
then runs yallHap validation stratified by coverage bins.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd


def load_anno_file(anno_path: Path) -> pd.DataFrame:
    """Load AADR annotation file."""
    # Read with tab separator
    df = pd.read_csv(anno_path, sep='\t', low_memory=False)

    # Rename columns for easier access
    df = df.rename(columns={
        'Genetic ID': 'sample_id',
        'Publication': 'publication',
        '1240k coverage (taken from original pulldown where possible)': 'coverage',
        'Y haplogroup (manual curation in terminal mutation format)': 'y_hg_terminal',
        'Y haplogroup (manual curation in ISOGG format)': 'y_hg_isogg',
    })

    return df


def load_ground_truth(gt_path: Path) -> dict[str, dict]:
    """Load ground truth file."""
    ground_truth = {}
    with open(gt_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_id = row['sample_id']
            ground_truth[sample_id] = {
                'haplogroup_terminal': row.get('haplogroup_terminal', ''),
                'haplogroup_isogg': row.get('haplogroup_isogg', ''),
            }
    return ground_truth


def load_validation_results(results_path: Path) -> dict[str, dict]:
    """Load existing validation results."""
    with open(results_path) as f:
        data = json.load(f)

    results = {}
    for r in data.get('results', []):
        results[r['sample_id']] = r
    return results


def parse_coverage(cov_str: str) -> float | None:
    """Parse coverage string to float."""
    if pd.isna(cov_str) or cov_str == '..' or cov_str == '' or cov_str == 'n/a':
        return None
    try:
        return float(cov_str)
    except ValueError:
        return None


def get_coverage_bin(coverage: float | None) -> str:
    """Assign coverage to bin."""
    if coverage is None:
        return 'unknown'
    elif coverage < 0.1:
        return '<0.1x'
    elif coverage < 0.5:
        return '0.1-0.5x'
    elif coverage < 1.0:
        return '0.5-1x'
    else:
        return '>=1x'


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Stratify AADR validation by coverage"
    )
    parser.add_argument(
        "--anno",
        type=Path,
        default=Path("data/ancient/v54.1.p1_HO_public.anno"),
        help="Path to AADR annotation file",
    )
    parser.add_argument(
        "--ground-truth",
        type=Path,
        default=Path("data/ancient/aadr_1240k_ground_truth.tsv"),
        help="Path to ground truth file",
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=None,
        help="Path to existing validation results JSON (optional)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON file",
    )

    args = parser.parse_args()

    # Load annotation file
    print("Loading annotation file...", file=sys.stderr)
    anno_df = load_anno_file(args.anno)
    print(f"  Loaded {len(anno_df)} samples", file=sys.stderr)

    # Load ground truth
    print("Loading ground truth...", file=sys.stderr)
    ground_truth = load_ground_truth(args.ground_truth)
    print(f"  Loaded {len(ground_truth)} samples", file=sys.stderr)

    # Create sample_id to coverage mapping
    # Sample IDs in ground truth may have different suffixes than anno file
    coverage_map = {}
    publication_map = {}

    for _, row in anno_df.iterrows():
        sample_id = str(row['sample_id'])
        coverage = parse_coverage(str(row.get('coverage', '')))
        publication = str(row.get('publication', ''))

        coverage_map[sample_id] = coverage
        publication_map[sample_id] = publication

        # Also map without suffix
        base_id = sample_id.rsplit('.', 1)[0] if '.' in sample_id else sample_id
        if base_id not in coverage_map:
            coverage_map[base_id] = coverage
            publication_map[base_id] = publication

    # Match ground truth samples to coverage
    matched = 0
    coverage_bins: dict[str, list] = defaultdict(list)
    publication_counts: dict[str, int] = defaultdict(int)

    for sample_id in ground_truth:
        coverage = None
        publication = ''

        # Try exact match first
        if sample_id in coverage_map:
            coverage = coverage_map[sample_id]
            publication = publication_map.get(sample_id, '')
        else:
            # Try base ID
            base_id = sample_id.rsplit('.', 1)[0] if '.' in sample_id else sample_id
            if base_id in coverage_map:
                coverage = coverage_map[base_id]
                publication = publication_map.get(base_id, '')

        if coverage is not None:
            matched += 1

        bin_name = get_coverage_bin(coverage)
        coverage_bins[bin_name].append({
            'sample_id': sample_id,
            'coverage': coverage,
            'publication': publication,
            'expected_hg': ground_truth[sample_id]['haplogroup_terminal'],
        })

        if publication:
            publication_counts[publication] += 1

    print(f"\nMatched {matched}/{len(ground_truth)} samples to coverage", file=sys.stderr)

    print("\nCoverage distribution:", file=sys.stderr)
    for bin_name in ['<0.1x', '0.1-0.5x', '0.5-1x', '>=1x', 'unknown']:
        count = len(coverage_bins.get(bin_name, []))
        print(f"  {bin_name}: {count} samples", file=sys.stderr)

    print(f"\nTop 10 publications:", file=sys.stderr)
    sorted_pubs = sorted(publication_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    for pub, count in sorted_pubs:
        print(f"  {pub}: {count}", file=sys.stderr)

    # If results file provided, compute stratified accuracy
    if args.results and args.results.exists():
        print("\nLoading validation results...", file=sys.stderr)
        results = load_validation_results(args.results)
        print(f"  Loaded {len(results)} results", file=sys.stderr)

        print("\nStratified accuracy by coverage:", file=sys.stderr)
        stratified_results = {}

        for bin_name in ['<0.1x', '0.1-0.5x', '0.5-1x', '>=1x', 'unknown']:
            samples = coverage_bins.get(bin_name, [])
            if not samples:
                continue

            correct = 0
            total = 0

            for s in samples:
                sample_id = s['sample_id']
                if sample_id in results:
                    r = results[sample_id]
                    if r.get('major_match', False):
                        correct += 1
                    total += 1

            if total > 0:
                accuracy = correct / total
                print(f"  {bin_name}: {accuracy:.1%} ({correct}/{total})", file=sys.stderr)
                stratified_results[bin_name] = {
                    'accuracy': accuracy,
                    'correct': correct,
                    'total': total,
                }

        # Output
        output_data = {
            'coverage_bins': {
                bin_name: len(samples)
                for bin_name, samples in coverage_bins.items()
            },
            'stratified_accuracy': stratified_results,
            'matched_samples': matched,
            'total_samples': len(ground_truth),
        }

        if args.output:
            with open(args.output, 'w') as f:
                json.dump(output_data, f, indent=2)
            print(f"\nResults written to {args.output}", file=sys.stderr)
        else:
            print(json.dumps(output_data, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())





