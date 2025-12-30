"""
Tests for variant scoring with allelic depth support.

TDD tests for:
- AD-aware variant scoring
- Quality-dependent error modeling
- Damage probability adjustments for ancient DNA
"""

from __future__ import annotations

import pytest

from yallhap.scoring import (
    VariantScore,
    compute_damage_probability,
    compute_variant_weight,
    score_variant_with_ad,
)
from yallhap.snps import SNP
from yallhap.vcf import Variant


class TestVariantScore:
    """Tests for VariantScore dataclass."""

    def test_variant_score_creation(self) -> None:
        """VariantScore can be created with all fields."""
        score = VariantScore(
            log_likelihood=-0.5,
            weight=0.95,
            is_derived=True,
            support_ratio=0.9,
            quality_factor=0.99,
            damage_probability=0.0,
        )

        assert score.log_likelihood == -0.5
        assert score.weight == 0.95
        assert score.is_derived is True
        assert score.support_ratio == 0.9

    def test_variant_score_defaults(self) -> None:
        """VariantScore has sensible defaults."""
        score = VariantScore(
            log_likelihood=0.0,
            weight=1.0,
            is_derived=False,
        )

        assert score.support_ratio is None
        assert score.quality_factor == 1.0
        assert score.damage_probability == 0.0


class TestScoreVariantWithAD:
    """Tests for score_variant_with_ad function."""

    @pytest.fixture
    def sample_snp(self) -> SNP:
        """A sample SNP for testing."""
        return SNP(
            name="M269",
            position_grch38=22749853,
            ancestral="C",
            derived="T",
            haplogroup="R1b",
        )

    def test_derived_high_support(self, sample_snp: SNP) -> None:
        """Derived call with high support gets high weight."""
        variant = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=1,  # Derived (alt) called
            quality=99,
            allele_depth=(0, 30),  # 100% support
        )

        score = score_variant_with_ad(variant, sample_snp)

        assert score.is_derived is True
        assert score.weight > 0.9
        assert score.support_ratio == 1.0

    def test_ancestral_high_support(self, sample_snp: SNP) -> None:
        """Ancestral call with high support gets high weight."""
        variant = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=0,  # Ancestral (ref) called
            quality=99,
            allele_depth=(30, 0),  # 100% support
        )

        score = score_variant_with_ad(variant, sample_snp)

        assert score.is_derived is False
        assert score.weight > 0.9
        assert score.support_ratio == 1.0

    def test_low_support_reduces_weight(self, sample_snp: SNP) -> None:
        """Low read support reduces variant weight."""
        # High support variant
        high_support = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=1,
            quality=99,
            allele_depth=(0, 30),
        )

        # Low support variant (50/50 reads)
        low_support = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=1,
            quality=99,
            allele_depth=(15, 15),
        )

        high_score = score_variant_with_ad(high_support, sample_snp)
        low_score = score_variant_with_ad(low_support, sample_snp)

        assert high_score.weight > low_score.weight

    def test_missing_ad_uses_quality(self, sample_snp: SNP) -> None:
        """Missing AD falls back to quality-based scoring."""
        variant = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=1,
            quality=99,
            allele_depth=None,
        )

        score = score_variant_with_ad(variant, sample_snp)

        assert score.is_derived is True
        assert score.support_ratio is None
        assert score.weight > 0  # Still has some weight from quality

    def test_low_quality_reduces_weight(self, sample_snp: SNP) -> None:
        """Low genotype quality reduces weight."""
        high_quality = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=1,
            quality=99,
            allele_depth=(0, 30),
        )

        low_quality = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=1,
            quality=10,
            allele_depth=(0, 30),
        )

        high_score = score_variant_with_ad(high_quality, sample_snp)
        low_score = score_variant_with_ad(low_quality, sample_snp)

        assert high_score.weight > low_score.weight
        assert low_score.quality_factor < high_score.quality_factor


class TestComputeVariantWeight:
    """Tests for compute_variant_weight function."""

    def test_perfect_quality_and_support(self) -> None:
        """Perfect quality and support gives weight near 1.0."""
        weight = compute_variant_weight(
            quality=99,
            support_ratio=1.0,
            error_rate=0.001,
        )

        assert weight > 0.95

    def test_no_support_data(self) -> None:
        """Missing support data uses quality only."""
        weight = compute_variant_weight(
            quality=99,
            support_ratio=None,
            error_rate=0.001,
        )

        assert 0.5 < weight < 1.0

    def test_low_support_ratio(self) -> None:
        """Low support ratio significantly reduces weight."""
        weight = compute_variant_weight(
            quality=99,
            support_ratio=0.5,  # 50% support
            error_rate=0.001,
        )

        # 50% support should give reduced weight compared to 100%
        assert weight < 0.75

    def test_missing_quality(self) -> None:
        """Missing quality uses support ratio only."""
        weight = compute_variant_weight(
            quality=None,
            support_ratio=0.9,
            error_rate=0.001,
        )

        assert 0.5 < weight < 1.0

    def test_error_rate_effect(self) -> None:
        """Higher error rate reduces weight."""
        low_error = compute_variant_weight(
            quality=99,
            support_ratio=0.95,
            error_rate=0.001,
        )

        high_error = compute_variant_weight(
            quality=99,
            support_ratio=0.95,
            error_rate=0.05,
        )

        assert low_error > high_error


class TestDamageProbability:
    """Tests for ancient DNA damage probability."""

    def test_ct_transition_high_damage(self) -> None:
        """C>T transition has high damage probability."""
        prob = compute_damage_probability(
            ref="C",
            alt="T",
            damage_rate=0.1,
        )

        assert prob > 0.05

    def test_ga_transition_high_damage(self) -> None:
        """G>A transition has high damage probability."""
        prob = compute_damage_probability(
            ref="G",
            alt="A",
            damage_rate=0.1,
        )

        assert prob > 0.05

    def test_transversion_no_damage(self) -> None:
        """Transversions have no damage probability."""
        prob = compute_damage_probability(
            ref="C",
            alt="A",
            damage_rate=0.1,
        )

        assert prob == 0.0

    def test_tc_transition_lower_damage(self) -> None:
        """T>C transition (not typical damage) has lower probability."""
        ct_prob = compute_damage_probability(ref="C", alt="T", damage_rate=0.1)
        tc_prob = compute_damage_probability(ref="T", alt="C", damage_rate=0.1)

        # C>T is typical damage, T>C is reverse complement damage (less common)
        assert ct_prob >= tc_prob

    def test_damage_rate_scaling(self) -> None:
        """Damage probability scales with damage rate."""
        low_rate = compute_damage_probability(ref="C", alt="T", damage_rate=0.05)
        high_rate = compute_damage_probability(ref="C", alt="T", damage_rate=0.2)

        assert high_rate > low_rate


class TestScoreWithDamageAdjustment:
    """Tests for scoring with ancient DNA damage adjustments."""

    @pytest.fixture
    def damage_snp(self) -> SNP:
        """SNP with C>T mutation (typical damage pattern)."""
        return SNP(
            name="L2",
            position_grch38=2850442,
            ancestral="C",
            derived="T",  # C>T is typical damage
            haplogroup="A",
        )

    @pytest.fixture
    def non_damage_snp(self) -> SNP:
        """SNP with C>A mutation (transversion, not damage)."""
        return SNP(
            name="M168",
            position_grch38=12866459,
            ancestral="C",
            derived="A",  # Transversion
            haplogroup="CT",
        )

    def test_damage_transition_penalized_in_ancient_mode(self, damage_snp: SNP) -> None:
        """C>T derived call is penalized in ancient mode."""
        variant = Variant(
            chrom="Y",
            position=2850442,
            ref="C",
            alt=("T",),
            genotype=1,
            quality=99,
            allele_depth=(0, 30),
        )

        modern_score = score_variant_with_ad(variant, damage_snp, is_ancient=False, damage_rate=0.0)
        ancient_score = score_variant_with_ad(variant, damage_snp, is_ancient=True, damage_rate=0.1)

        # Ancient mode should reduce confidence in C>T calls
        assert modern_score.weight > ancient_score.weight
        assert ancient_score.damage_probability > 0

    def test_transversion_not_penalized(self, non_damage_snp: SNP) -> None:
        """Transversion calls are not penalized in ancient mode."""
        variant = Variant(
            chrom="Y",
            position=12866459,
            ref="C",
            alt=("A",),
            genotype=1,
            quality=99,
            allele_depth=(0, 30),
        )

        modern_score = score_variant_with_ad(variant, non_damage_snp, is_ancient=False)
        ancient_score = score_variant_with_ad(
            variant, non_damage_snp, is_ancient=True, damage_rate=0.1
        )

        # Weights should be equal (no damage penalty for transversions)
        assert abs(modern_score.weight - ancient_score.weight) < 0.01
        assert ancient_score.damage_probability == 0.0


class TestMinSupportThreshold:
    """Tests for minimum support threshold handling."""

    @pytest.fixture
    def snp(self) -> SNP:
        """A sample SNP."""
        return SNP(
            name="M269",
            position_grch38=22749853,
            ancestral="C",
            derived="T",
            haplogroup="R1b",
        )

    def test_below_min_support_very_low_weight(self, snp: SNP) -> None:
        """Support below min_support threshold gets very low weight."""
        variant = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=1,
            quality=99,
            allele_depth=(25, 5),  # Only ~17% support for called allele
        )

        score = score_variant_with_ad(variant, snp, min_support=0.7)

        assert score.weight < 0.3
        assert score.support_ratio < 0.7

    def test_above_min_support_normal_weight(self, snp: SNP) -> None:
        """Support above min_support gets normal weight."""
        variant = Variant(
            chrom="Y",
            position=22749853,
            ref="C",
            alt=("T",),
            genotype=1,
            quality=99,
            allele_depth=(5, 25),  # ~83% support
        )

        score = score_variant_with_ad(variant, snp, min_support=0.7)

        assert score.weight > 0.7
        assert score.support_ratio > 0.7
