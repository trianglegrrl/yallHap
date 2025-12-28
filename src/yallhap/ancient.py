"""
Ancient DNA damage filtering for Y-chromosome haplogroup inference.

Implements filters for post-mortem DNA damage patterns commonly seen in
ancient DNA samples, including cytosine deamination (C>T and G>A transitions).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from yallhap.vcf import Variant


class MutationType(Enum):
    """Classification of mutation types."""

    TRANSITION = "transition"  # Purine<->Purine or Pyrimidine<->Pyrimidine
    TRANSVERSION = "transversion"  # Purine<->Pyrimidine


# Transition pairs (within-class mutations)
TRANSITIONS = frozenset([
    ("A", "G"), ("G", "A"),  # Purine <-> Purine
    ("C", "T"), ("T", "C"),  # Pyrimidine <-> Pyrimidine
])

# Transversion pairs (between-class mutations)
TRANSVERSIONS = frozenset([
    ("A", "C"), ("C", "A"),
    ("A", "T"), ("T", "A"),
    ("G", "C"), ("C", "G"),
    ("G", "T"), ("T", "G"),
])


def get_mutation_type(ref: str, alt: str) -> MutationType | None:
    """
    Determine if a mutation is a transition or transversion.

    Args:
        ref: Reference allele (single nucleotide)
        alt: Alternate allele (single nucleotide)

    Returns:
        MutationType.TRANSITION, MutationType.TRANSVERSION, or None if invalid
    """
    ref = ref.upper()
    alt = alt.upper()

    if (ref, alt) in TRANSITIONS:
        return MutationType.TRANSITION
    elif (ref, alt) in TRANSVERSIONS:
        return MutationType.TRANSVERSION
    return None


def is_transversion(ref: str, alt: str) -> bool:
    """Check if mutation is a transversion (not affected by deamination damage)."""
    return get_mutation_type(ref, alt) == MutationType.TRANSVERSION


def is_transition(ref: str, alt: str) -> bool:
    """Check if mutation is a transition (potentially affected by deamination)."""
    return get_mutation_type(ref, alt) == MutationType.TRANSITION


@dataclass
class DamagePattern:
    """
    A specific DNA damage pattern to filter.

    Attributes:
        ancestral: The ancestral (undamaged) allele
        derived: The derived (damaged) allele that mimics a real mutation
        description: Human-readable description of the damage
    """

    ancestral: str
    derived: str
    description: str

    def matches(self, ref: str, alt: str) -> bool:
        """Check if a mutation matches this damage pattern."""
        return ref.upper() == self.ancestral and alt.upper() == self.derived


# Standard ancient DNA damage patterns
# C>T: Cytosine deamination on + strand
# G>A: Cytosine deamination on - strand (reverse complement)
DEAMINATION_DAMAGE = [
    DamagePattern("C", "T", "C>T deamination (5' end bias)"),
    DamagePattern("G", "A", "G>A deamination (3' end bias)"),
]


class DamageFilter:
    """
    Filter for ancient DNA damage patterns.

    Filters out variants that match known post-mortem damage signatures,
    particularly cytosine deamination which causes C>T and G>A transitions.
    """

    def __init__(
        self,
        damage_patterns: list[DamagePattern] | None = None,
        strict: bool = False,
    ):
        """
        Initialize damage filter.

        Args:
            damage_patterns: List of damage patterns to filter.
                            Defaults to standard deamination patterns.
            strict: If True, filter all transitions (not just deamination).
        """
        self.damage_patterns = damage_patterns or DEAMINATION_DAMAGE
        self.strict = strict

    def is_damage_like(
        self,
        ref: str,
        alt: str,
        ancestral: str | None = None,
        derived: str | None = None,
    ) -> bool:
        """
        Check if a variant looks like ancient DNA damage.

        For haplogroup classification, we care about whether the called
        derived allele could be damage rather than a real mutation.

        Args:
            ref: Reference allele
            alt: Alternate allele (the called variant)
            ancestral: Expected ancestral allele for the SNP (optional)
            derived: Expected derived allele for the SNP (optional)

        Returns:
            True if the variant could be damage and should be filtered
        """
        # In strict mode, filter all transitions
        if self.strict:
            return is_transition(ref, alt)

        # Check against known damage patterns
        for pattern in self.damage_patterns:
            # Check if the mutation matches a damage pattern
            if pattern.matches(ref, alt):
                return True

            # If we know the SNP's expected alleles, check if damage
            # could explain an apparent derived call
            if ancestral and derived:
                # Damage is when: ancestral allele gets damaged to look like derived
                if pattern.ancestral == ancestral.upper() and pattern.derived == derived.upper():
                    return True

        return False

    def filter_variants(
        self,
        variants: dict[int, Variant],  # Forward reference to avoid circular import
    ) -> tuple[dict[int, Variant], int]:
        """
        Filter variants that look like damage.

        Args:
            variants: Dictionary of position -> Variant

        Returns:
            Tuple of (filtered variants dict, number filtered)
        """
        filtered: dict[int, Variant] = {}
        count_filtered = 0

        for pos, variant in variants.items():
            if self.is_damage_like(variant.ref, variant.alt[0] if variant.alt else ""):
                count_filtered += 1
            else:
                filtered[pos] = variant

        return filtered, count_filtered


class TransversionOnlyFilter:
    """
    Filter that only allows transversion mutations.

    Transversions are not affected by cytosine deamination and are thus
    more reliable for ancient DNA haplogroup calling, though at the cost
    of reduced resolution (fewer informative SNPs).
    """

    def is_allowed(self, ref: str, alt: str) -> bool:
        """Check if mutation is a transversion (allowed)."""
        return is_transversion(ref, alt)

    def filter_variants(
        self,
        variants: dict[int, Variant],
    ) -> tuple[dict[int, Variant], int]:
        """
        Filter variants to only keep transversions.

        Args:
            variants: Dictionary of position -> Variant

        Returns:
            Tuple of (filtered variants dict, number filtered)
        """
        filtered: dict[int, Variant] = {}
        count_filtered = 0

        for pos, variant in variants.items():
            alt = variant.alt[0] if variant.alt else ""
            if self.is_allowed(variant.ref, alt):
                filtered[pos] = variant
            else:
                count_filtered += 1

        return filtered, count_filtered


# Type alias for damage rescaling mode
DamageRescaleMode = Literal["none", "moderate", "aggressive"]


def apply_damage_rescale(
    quality: float | int,
    ref: str,
    alt: str,
    mode: DamageRescaleMode = "moderate",
) -> float:
    """
    Rescale quality scores for potentially damaged variants.

    Applies a penalty to quality scores for transitions that could be damage.
    Transversions are not penalized.

    Args:
        quality: Original quality score
        ref: Reference allele
        alt: Alternate allele
        mode: Rescaling aggressiveness
            - "none": No rescaling
            - "moderate": 50% penalty for C>T/G>A, 25% for other transitions
            - "aggressive": 75% penalty for C>T/G>A, 50% for other transitions

    Returns:
        Rescaled quality score
    """
    if mode == "none":
        return float(quality)

    # Transversions are not affected by damage
    if is_transversion(ref, alt):
        return float(quality)

    # Define penalties based on mode
    if mode == "moderate":
        ct_ga_penalty = 0.5  # C>T and G>A get 50% penalty
        other_transition_penalty = 0.25
    else:  # aggressive
        ct_ga_penalty = 0.75
        other_transition_penalty = 0.5

    # Apply appropriate penalty
    ref = ref.upper()
    alt = alt.upper()

    if (ref == "C" and alt == "T") or (ref == "G" and alt == "A"):
        return float(quality) * (1 - ct_ga_penalty)
    elif is_transition(ref, alt):
        return float(quality) * (1 - other_transition_penalty)

    return float(quality)

