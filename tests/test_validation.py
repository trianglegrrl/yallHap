"""
Unit tests for yallhap.validation module.

Tests for loading ground truth data, comparing haplogroups,
and running validation against 1000 Genomes samples.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


class TestGroundTruth:
    """Tests for loading Poznik 2016 ground truth data."""

    def test_load_ground_truth_tsv(self, tmp_path: Path) -> None:
        """Test loading ground truth haplogroup assignments."""
        tsv_content = """sample_id\thaplogroup\tpopulation
HG00096\tR1b1a2a1a2c1\tGBR
HG00097\tJ2a1\tGBR
NA19239\tE1b1a1a1g1a\tYRI"""

        tsv_path = tmp_path / "ground_truth.tsv"
        tsv_path.write_text(tsv_content)

        from yallhap.validation import load_ground_truth

        ground_truth = load_ground_truth(tsv_path)

        assert len(ground_truth) == 3
        assert ground_truth["HG00096"] == "R1b1a2a1a2c1"
        assert ground_truth["NA19239"] == "E1b1a1a1g1a"

    def test_load_ground_truth_missing_file_raises(self, tmp_path: Path) -> None:
        """Loading non-existent file raises FileNotFoundError."""
        from yallhap.validation import load_ground_truth

        with pytest.raises(FileNotFoundError):
            load_ground_truth(tmp_path / "nonexistent.tsv")

    def test_load_ground_truth_empty_file_raises(self, tmp_path: Path) -> None:
        """Loading empty file raises ValueError."""
        tsv_path = tmp_path / "empty.tsv"
        tsv_path.write_text("sample_id\thaplogroup\tpopulation\n")

        from yallhap.validation import load_ground_truth

        with pytest.raises(ValueError, match="No samples found"):
            load_ground_truth(tsv_path)


class TestHaplogroupComparison:
    """Tests for haplogroup comparison functions."""

    def test_compare_haplogroups_exact_match(self) -> None:
        """Test exact haplogroup match detection."""
        from yallhap.validation import compare_haplogroups

        result = compare_haplogroups("R1b1a2a1a2c1", "R1b1a2a1a2c1")

        assert result.exact_match is True
        assert result.major_match is True
        assert result.depth_difference == 0

    def test_compare_haplogroups_major_match(self) -> None:
        """Test major haplogroup match when terminal differs."""
        from yallhap.validation import compare_haplogroups

        result = compare_haplogroups("R1b1a2a1a2c1", "R1b1a2")

        assert result.exact_match is False
        assert result.major_match is True  # Both are R1b

    def test_compare_haplogroups_mismatch(self) -> None:
        """Test detection of haplogroup mismatch."""
        from yallhap.validation import compare_haplogroups

        result = compare_haplogroups("R1b1a2a1a2c1", "J2a1")

        assert result.exact_match is False
        assert result.major_match is False

    def test_compare_haplogroups_with_snp_suffix(self) -> None:
        """Test comparison handles SNP suffixes like R-L21."""
        from yallhap.validation import compare_haplogroups

        result = compare_haplogroups("R-L21", "R1b1a2a1a2c")

        assert result.exact_match is False
        assert result.major_match is True  # Both are R

    def test_compare_haplogroups_na_handling(self) -> None:
        """Test that NA haplogroup is handled as no match."""
        from yallhap.validation import compare_haplogroups

        result = compare_haplogroups("R1b", "NA")

        assert result.exact_match is False
        assert result.major_match is False


class TestValidationResult:
    """Tests for ValidationResult dataclass."""

    def test_validation_result_creation(self) -> None:
        """Test creating a ValidationResult."""
        from yallhap.validation import ComparisonResult, ValidationResult

        result = ValidationResult(
            sample_id="HG00096",
            expected_haplogroup="R1b1a2a1a2c1",
            called_haplogroup="R-L21",
            confidence=0.95,
            comparison=ComparisonResult(
                exact_match=False,
                major_match=True,
                depth_difference=2,
            ),
        )

        assert result.sample_id == "HG00096"
        assert result.expected_haplogroup == "R1b1a2a1a2c1"
        assert result.called_haplogroup == "R-L21"
        assert result.confidence == 0.95
        assert result.comparison.major_match is True


class TestValidationRunner:
    """Tests for ValidationRunner class."""

    def test_validation_runner_initialization(
        self, sample_tree_dict: dict, sample_snps_csv: Path
    ) -> None:
        """Test ValidationRunner can be initialized."""
        from yallhap.validation import ValidationRunner

        ground_truth = {"SAMPLE1": "R-L21", "SAMPLE2": "J2a"}

        runner = ValidationRunner(
            tree=Tree.from_dict(sample_tree_dict),
            snp_db=SNPDatabase.from_csv(sample_snps_csv),
            ground_truth=ground_truth,
        )

        assert runner is not None
        assert len(runner.ground_truth) == 2

    def test_validation_runner_get_expected_haplogroup(
        self, sample_tree_dict: dict, sample_snps_csv: Path
    ) -> None:
        """Test getting expected haplogroup for a sample."""
        from yallhap.validation import ValidationRunner

        ground_truth = {"SAMPLE1": "R-L21", "SAMPLE2": "J2a"}

        runner = ValidationRunner(
            tree=Tree.from_dict(sample_tree_dict),
            snp_db=SNPDatabase.from_csv(sample_snps_csv),
            ground_truth=ground_truth,
        )

        assert runner.get_expected_haplogroup("SAMPLE1") == "R-L21"
        assert runner.get_expected_haplogroup("SAMPLE2") == "J2a"

    def test_validation_runner_missing_sample_raises(
        self, sample_tree_dict: dict, sample_snps_csv: Path
    ) -> None:
        """Test that missing sample raises KeyError."""
        from yallhap.validation import ValidationRunner

        ground_truth = {"SAMPLE1": "R-L21"}

        runner = ValidationRunner(
            tree=Tree.from_dict(sample_tree_dict),
            snp_db=SNPDatabase.from_csv(sample_snps_csv),
            ground_truth=ground_truth,
        )

        with pytest.raises(KeyError, match="UNKNOWN"):
            runner.get_expected_haplogroup("UNKNOWN")


class TestValidationMetrics:
    """Tests for validation accuracy metrics."""

    def test_compute_metrics_perfect_accuracy(self) -> None:
        """Test metrics with 100% accuracy."""
        from yallhap.validation import ComparisonResult, ValidationResult, compute_metrics

        results = [
            ValidationResult(
                sample_id="S1",
                expected_haplogroup="R1b",
                called_haplogroup="R1b",
                confidence=0.99,
                comparison=ComparisonResult(
                    exact_match=True, major_match=True, depth_difference=0
                ),
            ),
            ValidationResult(
                sample_id="S2",
                expected_haplogroup="J2a",
                called_haplogroup="J2a",
                confidence=0.98,
                comparison=ComparisonResult(
                    exact_match=True, major_match=True, depth_difference=0
                ),
            ),
        ]

        metrics = compute_metrics(results)

        assert metrics.exact_match_rate == 1.0
        assert metrics.major_match_rate == 1.0
        assert metrics.total_samples == 2

    def test_compute_metrics_partial_accuracy(self) -> None:
        """Test metrics with partial accuracy."""
        from yallhap.validation import ComparisonResult, ValidationResult, compute_metrics

        results = [
            ValidationResult(
                sample_id="S1",
                expected_haplogroup="R1b1a2a1a2c1",
                called_haplogroup="R1b",
                confidence=0.95,
                comparison=ComparisonResult(
                    exact_match=False, major_match=True, depth_difference=5
                ),
            ),
            ValidationResult(
                sample_id="S2",
                expected_haplogroup="J2a",
                called_haplogroup="E1b",
                confidence=0.80,
                comparison=ComparisonResult(
                    exact_match=False, major_match=False, depth_difference=10
                ),
            ),
        ]

        metrics = compute_metrics(results)

        assert metrics.exact_match_rate == 0.0
        assert metrics.major_match_rate == 0.5
        assert metrics.total_samples == 2


class TestIntegrationValidation:
    """Integration tests for validation (require downloaded data)."""

    @pytest.mark.integration
    def test_downloaded_vcf_has_y_chromosome(self, fixtures_dir: Path) -> None:
        """Verify downloaded VCF contains Y chromosome data."""
        from yallhap.vcf import VCFReader

        vcf_path = fixtures_dir / "1kg_chrY_subset.vcf.gz"

        if not vcf_path.exists():
            pytest.skip("Integration test data not downloaded")

        with VCFReader(vcf_path) as reader:
            samples = reader.samples
            assert len(samples) > 0

            # Should have Y chromosome variants
            variants = list(reader.iter_variants())
            assert len(variants) > 1000  # Phase 3 has ~65k variants
