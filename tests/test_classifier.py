"""
Unit tests for yclade.classifier module.
"""

import pytest

from yclade.classifier import (
    HaplogroupCall,
    HaplogroupClassifier,
    QCScores,
    SNPStats,
)
from yclade.snps import SNPDatabase
from yclade.tree import Tree


class TestQCScores:
    """Tests for QCScores dataclass."""

    def test_qc_scores_combined(self) -> None:
        """Test combined score calculation."""
        scores = QCScores(
            qc1_backbone=1.0,
            qc2_terminal=1.0,
            qc3_path=1.0,
            qc4_posterior=1.0,
        )
        assert scores.combined == 1.0

    def test_qc_scores_combined_partial(self) -> None:
        """Test combined score with partial values."""
        scores = QCScores(
            qc1_backbone=0.8,
            qc2_terminal=0.9,
            qc3_path=0.85,
            qc4_posterior=0.95,
        )
        # Geometric mean
        expected = (0.8 * 0.9 * 0.85 * 0.95) ** 0.25
        assert abs(scores.combined - expected) < 0.001

    def test_qc_scores_combined_zero(self) -> None:
        """Test combined score with zero value."""
        scores = QCScores(
            qc1_backbone=0.0,
            qc2_terminal=1.0,
            qc3_path=1.0,
            qc4_posterior=1.0,
        )
        assert scores.combined == 0.0


class TestSNPStats:
    """Tests for SNPStats dataclass."""

    def test_snp_stats_total_called(self) -> None:
        """Test total_called property."""
        stats = SNPStats(
            informative_tested=100,
            derived=50,
            ancestral=40,
            missing=10,
        )
        assert stats.total_called == 90


class TestHaplogroupCall:
    """Tests for HaplogroupCall dataclass."""

    def test_to_dict(self) -> None:
        """Test JSON serialization."""
        call = HaplogroupCall(
            sample="TEST",
            haplogroup="R-L21",
            confidence=0.95,
            qc_scores=QCScores(0.9, 1.0, 0.95, 0.95),
            path=["ROOT", "R", "R-L21"],
            defining_snps=["L21"],
            alternatives=[("R1b", 0.8)],
            snp_stats=SNPStats(100, 50, 40, 10),
            reference="grch38",
            tree_version="YFull v13",
        )

        d = call.to_dict()
        assert d["sample"] == "TEST"
        assert d["haplogroup"] == "R-L21"
        assert d["confidence"] == 0.95
        assert "quality_scores" in d
        assert d["quality_scores"]["qc1_backbone"] == 0.9


class TestHaplogroupClassifier:
    """Tests for HaplogroupClassifier class."""

    @pytest.fixture
    def classifier(self, sample_tree_dict: dict, sample_snps_csv) -> HaplogroupClassifier:
        """Create a classifier with test data."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        return HaplogroupClassifier(
            tree=tree,
            snp_db=snp_db,
            reference="grch38",
        )

    def test_classifier_creation(self, classifier: HaplogroupClassifier) -> None:
        """Test classifier initialization."""
        assert classifier.tree is not None
        assert classifier.snp_db is not None
        assert classifier.reference == "grch38"

    def test_classifier_ancient_mode(self, sample_tree_dict: dict, sample_snps_csv) -> None:
        """Test ancient DNA mode configuration."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)

        classifier = HaplogroupClassifier(
            tree=tree,
            snp_db=snp_db,
            ancient_mode=True,
            min_depth=1,
        )

        assert classifier.ancient_mode is True
        assert classifier.min_depth == 1


class TestClassifierDamageFiltering:
    """Tests for ancient DNA damage filtering."""

    @pytest.fixture
    def classifier_ancient(self, sample_tree_dict: dict, sample_snps_csv) -> HaplogroupClassifier:
        """Create classifier in ancient mode."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        return HaplogroupClassifier(
            tree=tree,
            snp_db=snp_db,
            ancient_mode=True,
        )

    def test_damage_like_c_to_t(self, classifier_ancient: HaplogroupClassifier) -> None:
        """Test C>T is flagged as damage-like."""
        from yclade.snps import SNP

        snp = SNP(name="TEST", ancestral="C", derived="T")
        assert classifier_ancient._is_damage_like(snp, "T")

    def test_damage_like_g_to_a(self, classifier_ancient: HaplogroupClassifier) -> None:
        """Test G>A is flagged as damage-like."""
        from yclade.snps import SNP

        snp = SNP(name="TEST", ancestral="G", derived="A")
        assert classifier_ancient._is_damage_like(snp, "A")

    def test_not_damage_like_transversion(self, classifier_ancient: HaplogroupClassifier) -> None:
        """Test transversions are not flagged as damage."""
        from yclade.snps import SNP

        snp = SNP(name="TEST", ancestral="C", derived="A")
        assert not classifier_ancient._is_damage_like(snp, "A")


# Integration tests would go here but require real VCF data
class TestClassifierIntegration:
    """Integration tests for classifier (marked slow)."""

    @pytest.mark.skip(reason="Requires real test data")
    def test_classify_real_sample(self) -> None:
        """Test classification on real 1000 Genomes sample."""
        pass
