"""Tests for ancient DNA damage filtering."""

from __future__ import annotations

import pytest

from yallhap.ancient import (
    DamageFilter,
    DamagePattern,
    MutationType,
    TransversionOnlyFilter,
    apply_damage_rescale,
    get_mutation_type,
    is_transition,
    is_transversion,
)


class TestMutationType:
    """Tests for mutation type classification."""

    def test_transitions(self) -> None:
        """Test transition identification."""
        # Purine <-> Purine
        assert get_mutation_type("A", "G") == MutationType.TRANSITION
        assert get_mutation_type("G", "A") == MutationType.TRANSITION
        # Pyrimidine <-> Pyrimidine
        assert get_mutation_type("C", "T") == MutationType.TRANSITION
        assert get_mutation_type("T", "C") == MutationType.TRANSITION

    def test_transversions(self) -> None:
        """Test transversion identification."""
        # All purine <-> pyrimidine combinations
        assert get_mutation_type("A", "C") == MutationType.TRANSVERSION
        assert get_mutation_type("A", "T") == MutationType.TRANSVERSION
        assert get_mutation_type("G", "C") == MutationType.TRANSVERSION
        assert get_mutation_type("G", "T") == MutationType.TRANSVERSION
        assert get_mutation_type("C", "A") == MutationType.TRANSVERSION
        assert get_mutation_type("T", "A") == MutationType.TRANSVERSION
        assert get_mutation_type("C", "G") == MutationType.TRANSVERSION
        assert get_mutation_type("T", "G") == MutationType.TRANSVERSION

    def test_case_insensitive(self) -> None:
        """Test case insensitivity."""
        assert get_mutation_type("a", "g") == MutationType.TRANSITION
        assert get_mutation_type("A", "t") == MutationType.TRANSVERSION

    def test_invalid_alleles(self) -> None:
        """Test invalid allele handling."""
        assert get_mutation_type("A", "A") is None  # Same allele
        assert get_mutation_type("X", "Y") is None  # Invalid bases

    def test_is_transition_helper(self) -> None:
        """Test is_transition helper function."""
        assert is_transition("C", "T") is True
        assert is_transition("A", "C") is False

    def test_is_transversion_helper(self) -> None:
        """Test is_transversion helper function."""
        assert is_transversion("A", "T") is True
        assert is_transversion("C", "T") is False


class TestDamagePattern:
    """Tests for DamagePattern class."""

    def test_matches_exact(self) -> None:
        """Test exact pattern matching."""
        pattern = DamagePattern("C", "T", "C>T deamination")

        assert pattern.matches("C", "T") is True
        assert pattern.matches("G", "A") is False

    def test_matches_case_insensitive(self) -> None:
        """Test case insensitive matching."""
        pattern = DamagePattern("C", "T", "C>T deamination")

        assert pattern.matches("c", "t") is True
        assert pattern.matches("C", "t") is True


class TestDamageFilter:
    """Tests for DamageFilter class."""

    def test_default_patterns_c_to_t(self) -> None:
        """Test C>T damage detection with default patterns."""
        filt = DamageFilter()

        # C>T is damage-like
        assert filt.is_damage_like("C", "T") is True

    def test_default_patterns_g_to_a(self) -> None:
        """Test G>A damage detection with default patterns."""
        filt = DamageFilter()

        # G>A is damage-like (reverse complement)
        assert filt.is_damage_like("G", "A") is True

    def test_transversions_not_damage(self) -> None:
        """Test that transversions are not flagged as damage."""
        filt = DamageFilter()

        assert filt.is_damage_like("A", "T") is False
        assert filt.is_damage_like("G", "C") is False

    def test_other_transitions_not_damage(self) -> None:
        """Test that non-deamination transitions are not flagged by default."""
        filt = DamageFilter()

        # T>C and A>G are transitions but not deamination
        assert filt.is_damage_like("T", "C") is False
        assert filt.is_damage_like("A", "G") is False

    def test_strict_mode_all_transitions(self) -> None:
        """Test strict mode flags all transitions."""
        filt = DamageFilter(strict=True)

        # All transitions should be flagged
        assert filt.is_damage_like("C", "T") is True
        assert filt.is_damage_like("T", "C") is True
        assert filt.is_damage_like("A", "G") is True
        assert filt.is_damage_like("G", "A") is True

        # Transversions still not flagged
        assert filt.is_damage_like("A", "T") is False

    def test_with_snp_context(self) -> None:
        """Test damage detection with SNP ancestral/derived context."""
        filt = DamageFilter()

        # If ancestral=C and derived=T, and we see a variant, it could be damage
        assert filt.is_damage_like("C", "T", ancestral="C", derived="T") is True

        # If ancestral=A and derived=G, C>T damage won't explain it
        assert filt.is_damage_like("A", "G", ancestral="A", derived="G") is False


class TestTransversionOnlyFilter:
    """Tests for TransversionOnlyFilter class."""

    def test_allows_transversions(self) -> None:
        """Test that transversions are allowed."""
        filt = TransversionOnlyFilter()

        assert filt.is_allowed("A", "T") is True
        assert filt.is_allowed("G", "C") is True
        assert filt.is_allowed("C", "A") is True

    def test_blocks_transitions(self) -> None:
        """Test that transitions are blocked."""
        filt = TransversionOnlyFilter()

        assert filt.is_allowed("C", "T") is False
        assert filt.is_allowed("A", "G") is False


class TestDamageRescale:
    """Tests for damage rescaling function."""

    def test_no_rescale_mode(self) -> None:
        """Test no rescaling mode."""
        assert apply_damage_rescale(30, "C", "T", mode="none") == 30.0
        assert apply_damage_rescale(30, "A", "G", mode="none") == 30.0

    def test_transversions_not_penalized(self) -> None:
        """Test that transversions are never penalized."""
        assert apply_damage_rescale(30, "A", "T", mode="moderate") == 30.0
        assert apply_damage_rescale(30, "A", "T", mode="aggressive") == 30.0

    def test_moderate_ct_ga_penalty(self) -> None:
        """Test moderate penalty for C>T and G>A."""
        # 50% penalty in moderate mode
        assert apply_damage_rescale(30, "C", "T", mode="moderate") == 15.0
        assert apply_damage_rescale(30, "G", "A", mode="moderate") == 15.0

    def test_aggressive_ct_ga_penalty(self) -> None:
        """Test aggressive penalty for C>T and G>A."""
        # 75% penalty in aggressive mode
        assert apply_damage_rescale(40, "C", "T", mode="aggressive") == 10.0
        assert apply_damage_rescale(40, "G", "A", mode="aggressive") == 10.0

    def test_other_transitions_lesser_penalty(self) -> None:
        """Test that other transitions get lesser penalty."""
        # T>C and A>G get smaller penalty
        # Moderate: 25% penalty
        assert apply_damage_rescale(40, "T", "C", mode="moderate") == 30.0
        # Aggressive: 50% penalty
        assert apply_damage_rescale(40, "T", "C", mode="aggressive") == 20.0

