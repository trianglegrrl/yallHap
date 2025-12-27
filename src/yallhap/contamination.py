"""
Y-chromosome contamination estimation.

Estimates contamination levels from allele depths at informative SNP sites.
Based on the approach used by pathPhynder's contamin.py.
"""

from __future__ import annotations

from dataclasses import dataclass

from yallhap.ancient import is_transversion
from yallhap.snps import ReferenceGenome, SNPDatabase
from yallhap.vcf import Variant


@dataclass
class ContaminationResult:
    """
    Result of contamination estimation.

    Attributes:
        rate: Estimated contamination rate (0.0 - 1.0)
        n_sites: Number of informative sites used in estimation
    """

    rate: float
    n_sites: int


def estimate_contamination(
    variants: dict[int, Variant],
    min_depth: int = 3,
    transversions_only: bool = False,
) -> tuple[float, int]:
    """
    Estimate Y-chromosome contamination from allele depths.

    Contamination is estimated by looking at the minor allele frequency
    at each site. For a haploid Y chromosome, we expect 100% of reads
    to support the true allele. Minor allele reads indicate contamination.

    Args:
        variants: Dictionary of position -> Variant with allele depth data
        min_depth: Minimum total depth to include a site
        transversions_only: If True, exclude C>T and G>A transitions
                           (which may be ancient DNA damage)

    Returns:
        Tuple of (contamination_rate, number_of_sites_used)
    """
    contamination_estimates: list[float] = []

    for _pos, variant in variants.items():
        # Skip variants without allele depth data
        if variant.allele_depth is None:
            continue

        # Check minimum depth
        total_depth = sum(variant.allele_depth)
        if total_depth < min_depth:
            continue

        # Skip if transversions_only and this is a transition
        if transversions_only and variant.alt and len(variant.alt) > 0:
            # Check if it's a C>T, T>C, G>A, or A>G transition
            ref = variant.ref.upper()
            alt = variant.alt[0].upper()
            if not is_transversion(ref, alt):
                continue

        # Calculate minor allele frequency
        # Minor allele = smaller of the two allele counts
        if len(variant.allele_depth) >= 2:
            ref_count = variant.allele_depth[0]
            alt_count = variant.allele_depth[1]

            if ref_count + alt_count > 0:
                minor_count = min(ref_count, alt_count)
                major_count = max(ref_count, alt_count)
                contamination_at_site = minor_count / (minor_count + major_count)
                contamination_estimates.append(contamination_at_site)

    if not contamination_estimates:
        return 0.0, 0

    # Return mean contamination estimate
    mean_rate = sum(contamination_estimates) / len(contamination_estimates)
    return mean_rate, len(contamination_estimates)


def estimate_contamination_with_snpdb(
    variants: dict[int, Variant],
    snp_db: SNPDatabase,
    reference: ReferenceGenome = "grch38",
    min_depth: int = 3,
    transversions_only: bool = False,
) -> tuple[float, int]:
    """
    Estimate contamination using only sites from SNP database.

    This restricts estimation to known informative sites, reducing
    noise from non-variant positions.

    Args:
        variants: Dictionary of position -> Variant
        snp_db: SNP database with informative positions
        reference: Reference genome for position lookup
        min_depth: Minimum total depth to include a site
        transversions_only: If True, exclude C>T and G>A transitions

    Returns:
        Tuple of (contamination_rate, number_of_sites_used)
    """
    # Filter to only positions in SNP database
    db_positions = snp_db.positions[reference]
    filtered_variants = {pos: var for pos, var in variants.items() if pos in db_positions}

    return estimate_contamination(
        filtered_variants,
        min_depth=min_depth,
        transversions_only=transversions_only,
    )


def is_contaminated(
    result: ContaminationResult,
    threshold: float = 0.10,
) -> bool:
    """
    Check if sample exceeds contamination threshold.

    Args:
        result: ContaminationResult from estimation
        threshold: Maximum acceptable contamination rate (default 10%)

    Returns:
        True if contamination rate exceeds threshold
    """
    return result.rate > threshold
