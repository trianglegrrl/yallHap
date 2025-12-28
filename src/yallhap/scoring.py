"""
Variant scoring with allelic depth and damage awareness.

Provides probabilistic scoring of variants incorporating:
- Allelic depth (AD) for read-level support
- Genotype quality (GQ) for call confidence
- Ancient DNA damage modeling for transitions
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from yallhap.snps import SNP
from yallhap.vcf import Variant

# Damage-typical transitions (C>T on forward strand, G>A on reverse)
DAMAGE_TRANSITIONS = frozenset([("C", "T"), ("G", "A")])

# All transitions (including reverse complement damage patterns)
ALL_TRANSITIONS = frozenset([("C", "T"), ("G", "A"), ("T", "C"), ("A", "G")])


@dataclass
class VariantScore:
    """
    Probabilistic score for a variant call.

    Attributes:
        log_likelihood: Log-likelihood contribution to path
        weight: Combined confidence weight (0.0-1.0)
        is_derived: Whether the derived allele was called
        support_ratio: Fraction of reads supporting called allele (from AD)
        quality_factor: Quality-based confidence factor (from GQ)
        damage_probability: Probability this call is due to DNA damage
    """

    log_likelihood: float
    weight: float
    is_derived: bool
    support_ratio: float | None = None
    quality_factor: float = 1.0
    damage_probability: float = 0.0


def compute_damage_probability(
    ref: str,
    alt: str,
    damage_rate: float = 0.1,
) -> float:
    """
    Compute probability that a mutation is due to ancient DNA damage.

    Ancient DNA damage primarily causes C>T and G>A transitions due to
    cytosine deamination. This function returns a probability that accounts
    for the damage_rate parameter.

    Args:
        ref: Reference allele
        alt: Alternative allele
        damage_rate: Expected rate of damage-induced mutations (0.0-1.0)

    Returns:
        Probability (0.0-1.0) that this mutation is damage-induced
    """
    if (ref, alt) in DAMAGE_TRANSITIONS:
        # Primary damage pattern: C>T or G>A
        return damage_rate
    elif (ref, alt) in ALL_TRANSITIONS:
        # Reverse complement damage (T>C, A>G) - less common
        return damage_rate * 0.3
    else:
        # Transversions are not caused by deamination damage
        return 0.0


def compute_variant_weight(
    quality: int | None,
    support_ratio: float | None,
    error_rate: float = 0.001,
) -> float:
    """
    Compute combined weight for a variant call.

    Combines genotype quality and read support into a single weight.

    Args:
        quality: Genotype quality (GQ) score, typically 0-99
        support_ratio: Fraction of reads supporting called allele
        error_rate: Expected sequencing error rate

    Returns:
        Weight between 0.0 and 1.0
    """
    # Quality-based factor
    if quality is not None:
        # Convert GQ to probability of correct call
        # GQ = -10 * log10(P(error))
        # P(correct) = 1 - 10^(-GQ/10)
        p_correct = 1.0 - math.pow(10, -quality / 10.0)
        quality_factor = max(0.0, min(1.0, p_correct))
    else:
        quality_factor = 0.5  # Unknown quality = 50% confidence

    # Support-based factor
    if support_ratio is not None:
        # High support = high confidence
        # Penalize heavily if support < 50%
        if support_ratio < 0.5:
            support_factor = support_ratio * 0.5  # Heavy penalty
        elif support_ratio < 0.7:
            support_factor = 0.5 + (support_ratio - 0.5) * 1.0  # Gradual recovery
        else:
            support_factor = 0.7 + (support_ratio - 0.7) * 1.0  # Near-linear
    else:
        support_factor = 1.0  # No AD data = assume good support

    # Error rate adjustment
    error_factor = 1.0 - error_rate

    # Combine factors geometrically
    if support_ratio is not None:
        weight = (quality_factor * support_factor * error_factor) ** 0.5
    else:
        weight = quality_factor * error_factor

    return max(0.0, min(1.0, weight))


def score_variant_with_ad(
    variant: Variant,
    snp: SNP,
    error_rate: float = 0.001,
    damage_rate: float = 0.1,
    min_support: float = 0.7,
    is_ancient: bool = False,
) -> VariantScore:
    """
    Score a variant call incorporating allelic depth and damage probability.

    This is the main scoring function for Bayesian classification. It produces
    a VariantScore that encapsulates all confidence information about a call.

    Args:
        variant: The variant call from VCF
        snp: The SNP definition from database
        error_rate: Expected sequencing error rate
        damage_rate: Expected ancient DNA damage rate
        min_support: Minimum read support for confident calls
        is_ancient: Whether to apply ancient DNA damage modeling

    Returns:
        VariantScore with all scoring components
    """
    # Determine if derived allele was called
    called_allele = variant.called_allele
    is_derived = called_allele == snp.derived

    # Get support ratio from variant
    support_ratio = variant.support_ratio

    # Compute quality factor from GQ
    quality = variant.quality
    if quality is not None:
        p_correct = 1.0 - math.pow(10, -quality / 10.0)
        quality_factor = max(0.0, min(1.0, p_correct))
    else:
        quality_factor = 0.5

    # Compute damage probability if ancient mode
    if is_ancient and is_derived:
        damage_prob = compute_damage_probability(snp.ancestral, snp.derived, damage_rate)
    else:
        damage_prob = 0.0

    # Compute base weight
    weight = compute_variant_weight(quality, support_ratio, error_rate)

    # Apply damage penalty if applicable
    if damage_prob > 0:
        # Reduce weight by the probability this is damage
        weight = weight * (1.0 - damage_prob)

    # Apply minimum support penalty
    if support_ratio is not None and support_ratio < min_support:
        # Heavy penalty for low support
        penalty = (support_ratio / min_support) ** 2
        weight = weight * penalty

    # Compute log-likelihood contribution
    # For derived calls: log(weight)
    # For ancestral calls: log(weight) (also positive - we're scoring the call quality)
    if weight > 0:
        log_likelihood = math.log(weight)
    else:
        log_likelihood = -100.0  # Very negative for zero weight

    return VariantScore(
        log_likelihood=log_likelihood,
        weight=weight,
        is_derived=is_derived,
        support_ratio=support_ratio,
        quality_factor=quality_factor,
        damage_probability=damage_prob,
    )


