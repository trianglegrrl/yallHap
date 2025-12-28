#!/usr/bin/env python3
"""
Calculate two-proportion z-test p-values for Supplementary Tables S3 and S4.

Replaces the conservative non-overlapping CI test with proper statistical inference.

Usage:
    python scripts/calculate_significance.py
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from scipy import stats
import numpy as np


@dataclass
class BinComparison:
    """Comparison between Heuristic and Bayesian modes for a density bin."""

    bin_name: str
    n: int
    heuristic_pct: float
    bayesian_pct: float
    delta_pp: float

    # Computed fields
    heuristic_correct: int = 0
    bayesian_correct: int = 0
    z_stat: float = 0.0
    p_value: float = 1.0
    significant_005: bool = False
    significant_001: bool = False

    def compute_test(self) -> None:
        """Compute two-proportion z-test."""
        # Convert percentages to counts
        self.heuristic_correct = round(self.n * self.heuristic_pct / 100)
        self.bayesian_correct = round(self.n * self.bayesian_pct / 100)

        # Two-proportion z-test
        # H0: p1 = p2 (no difference between modes)
        # Ha: p1 != p2 (two-tailed)

        p1 = self.heuristic_correct / self.n
        p2 = self.bayesian_correct / self.n

        # Pooled proportion under H0
        p_pooled = (self.heuristic_correct + self.bayesian_correct) / (2 * self.n)

        # Standard error
        se = np.sqrt(p_pooled * (1 - p_pooled) * (2 / self.n))

        if se > 0:
            self.z_stat = (p2 - p1) / se
            # Two-tailed p-value
            self.p_value = 2 * (1 - stats.norm.cdf(abs(self.z_stat)))
        else:
            self.z_stat = 0.0
            self.p_value = 1.0

        self.significant_005 = self.p_value < 0.05
        self.significant_001 = self.p_value < 0.01


def main() -> None:
    """Calculate significance tests for supplementary tables."""

    # Table 3: AADR Ancient DNA Validation by Variant Density (main paper)
    table_3_data = [
        ("<1%", 101, 28.7, 33.7, +5.0),
        ("1–4%", 478, 37.2, 37.9, +0.7),
        ("4–10%", 736, 52.4, 71.7, +19.3),
        ("10–50%", 4101, 97.0, 97.8, +0.8),
        ("≥50%", 1927, 99.1, 99.0, -0.1),
        ("Overall", 7333, 88.3, 90.7, +2.4),
    ]

    # Supplementary Table S3: 1-10% Variant Density Breakdown
    table_s3_data = [
        ("1–2%", 167, 37.7, 36.5, -1.2),
        ("2–3%", 156, 39.1, 34.0, -5.1),
        ("3–4%", 155, 34.8, 43.2, +8.4),
        ("4–5%", 151, 35.8, 59.6, +23.8),
        ("5–6%", 172, 44.2, 68.0, +23.8),
        ("6–7%", 122, 47.5, 68.9, +21.4),
        ("7–8%", 98, 56.1, 76.5, +20.4),
        ("8–9%", 101, 69.3, 82.2, +12.9),
        ("9–10%", 92, 79.3, 85.9, +6.6),
        ("Overall", 1214, 46.5, 58.4, +11.9),
    ]

    # Supplementary Table S4: Decile Analysis
    table_s4_data = [
        ("0–10%", 400, 43.2, 57.2, +14.0),
        ("10–20%", 400, 84.0, 92.0, +8.0),
        ("20–30%", 400, 98.2, 98.5, +0.3),
        ("30–40%", 400, 99.5, 99.5, 0.0),
        ("40–50%", 400, 98.8, 97.5, -1.3),
        ("50–60%", 400, 99.0, 99.0, 0.0),
        ("60–70%", 400, 98.8, 98.8, 0.0),
        ("70–80%", 192, 97.4, 96.9, -0.5),
        ("80–90%", 77, 100.0, 100.0, 0.0),
        ("90–100%", 148, 100.0, 99.3, -0.7),
        ("Overall", 3217, 90.1, 92.6, +2.5),
    ]

    # Process Table 3 (main paper)
    print("=" * 70)
    print("Table 3: Two-Proportion Z-Test Results")
    print("=" * 70)
    print()

    t3_results = []
    for bin_name, n, heur, bayes, delta in table_3_data:
        comparison = BinComparison(bin_name, n, heur, bayes, delta)
        comparison.compute_test()
        t3_results.append(comparison)

        sig_marker = "***" if comparison.significant_001 else ("*" if comparison.significant_005 else "")
        print(f"{bin_name:>10}: Δ = {comparison.delta_pp:+5.1f} pp, "
              f"z = {comparison.z_stat:+5.2f}, p = {comparison.p_value:.4f} {sig_marker}")

    print()
    print("* p < 0.05, *** p < 0.01")
    print()

    # Process Table S3
    print("=" * 70)
    print("Supplementary Table S3: Two-Proportion Z-Test Results")
    print("=" * 70)
    print()

    s3_results = []
    for bin_name, n, heur, bayes, delta in table_s3_data:
        comparison = BinComparison(bin_name, n, heur, bayes, delta)
        comparison.compute_test()
        s3_results.append(comparison)

        sig_marker = "***" if comparison.significant_001 else ("*" if comparison.significant_005 else "")
        print(f"{bin_name:>10}: Δ = {comparison.delta_pp:+5.1f} pp, "
              f"z = {comparison.z_stat:+5.2f}, p = {comparison.p_value:.4f} {sig_marker}")

    print()
    print("* p < 0.05, *** p < 0.01")
    print()

    # Process Table S4
    print("=" * 70)
    print("Supplementary Table S4: Two-Proportion Z-Test Results")
    print("=" * 70)
    print()

    s4_results = []
    for bin_name, n, heur, bayes, delta in table_s4_data:
        comparison = BinComparison(bin_name, n, heur, bayes, delta)
        comparison.compute_test()
        s4_results.append(comparison)

        sig_marker = "***" if comparison.significant_001 else ("*" if comparison.significant_005 else "")
        print(f"{bin_name:>10}: Δ = {comparison.delta_pp:+5.1f} pp, "
              f"z = {comparison.z_stat:+5.2f}, p = {comparison.p_value:.4f} {sig_marker}")

    print()
    print("* p < 0.05, *** p < 0.01")

    # Save results to JSON
    output = {
        "table_3": [
            {
                "bin": r.bin_name,
                "n": r.n,
                "heuristic_pct": r.heuristic_pct,
                "bayesian_pct": r.bayesian_pct,
                "delta_pp": r.delta_pp,
                "z_statistic": round(r.z_stat, 3),
                "p_value": round(r.p_value, 6),
                "significant_005": bool(r.significant_005),
                "significant_001": bool(r.significant_001),
            }
            for r in t3_results
        ],
        "table_s3": [
            {
                "bin": r.bin_name,
                "n": r.n,
                "heuristic_pct": r.heuristic_pct,
                "bayesian_pct": r.bayesian_pct,
                "delta_pp": r.delta_pp,
                "z_statistic": round(r.z_stat, 3),
                "p_value": round(r.p_value, 6),
                "significant_005": bool(r.significant_005),
                "significant_001": bool(r.significant_001),
            }
            for r in s3_results
        ],
        "table_s4": [
            {
                "bin": r.bin_name,
                "n": r.n,
                "heuristic_pct": r.heuristic_pct,
                "bayesian_pct": r.bayesian_pct,
                "delta_pp": r.delta_pp,
                "z_statistic": round(r.z_stat, 3),
                "p_value": round(r.p_value, 6),
                "significant_005": bool(r.significant_005),
                "significant_001": bool(r.significant_001),
            }
            for r in s4_results
        ],
    }

    output_path = Path(__file__).parent.parent / "results" / "significance_tests.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)

    print()
    print(f"Results saved to: {output_path}")

    # Generate markdown table for S3
    print()
    print("=" * 70)
    print("Updated Supplementary Table S3 (with p-values)")
    print("=" * 70)
    print()
    print("| Density Bin | n | Heuristic TV | Bayesian Ancient | Δ | p-value |")
    print("|-------------|---|--------------|------------------|---|---------|")
    for r in s3_results:
        p_str = f"{r.p_value:.4f}" if r.p_value >= 0.0001 else "<0.0001"
        sig = "**" if r.significant_001 else ("*" if r.significant_005 else "")
        delta_str = f"**{r.delta_pp:+.1f} pp**" if r.significant_005 else f"{r.delta_pp:+.1f} pp"
        print(f"| {r.bin_name} | {r.n} | {r.heuristic_pct:.1f}% | {r.bayesian_pct:.1f}% | {delta_str} | {p_str}{sig} |")


if __name__ == "__main__":
    main()

