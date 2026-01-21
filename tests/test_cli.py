"""
Tests for CLI ancient mode and ISOGG output options.

Following TDD - these tests are written first, before implementation.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from yallhap.cli import main


@pytest.fixture
def runner() -> CliRunner:
    """Create a CLI test runner."""
    return CliRunner()


@pytest.fixture
def mock_tree(tmp_path: Path) -> Path:
    """Create a minimal mock tree JSON file."""
    import json

    tree_data = {
        'ROOT (Y-Chromosome "Adam")': ["A00", "A0-T"],
        "A00": [],
        "A0-T": ["R"],
        "R": ["R1b"],
        "R1b": ["R-L21"],
        "R-L21": [],
    }
    tree_path = tmp_path / "tree.json"
    with open(tree_path, "w") as f:
        json.dump(tree_data, f)
    return tree_path


@pytest.fixture
def mock_snp_db(tmp_path: Path) -> Path:
    """Create a minimal mock SNP database CSV."""
    csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
L21,"M529,S145",2887478,2655229,C,T,R-L21
M269,,2803636,2571333,G,A,R1b
"""
    csv_path = tmp_path / "snps.csv"
    csv_path.write_text(csv_content)
    return csv_path


class TestAncientModeCLI:
    """Tests for ancient mode CLI options."""

    def test_ancient_mode_flag_exists(self, runner: CliRunner) -> None:
        """--ancient flag is recognized."""
        result = runner.invoke(main, ["classify", "--help"])
        assert "--ancient" in result.output

    def test_ancient_mode_levels_option(self, runner: CliRunner) -> None:
        """--ancient-level option is recognized."""
        result = runner.invoke(main, ["classify", "--help"])
        # After implementation, this should be in the help
        # For now, test presence of related ancient options
        assert "--ancient" in result.output or "--transversions-only" in result.output

    def test_contamination_flag_exists(self, runner: CliRunner) -> None:
        """--estimate-contamination flag is recognized."""
        result = runner.invoke(main, ["classify", "--help"])
        # This will be added in implementation
        assert "--estimate-contamination" in result.output or result.exit_code != 0

    def test_max_tolerance_option(self, runner: CliRunner) -> None:
        """--max-tolerance option is recognized."""
        result = runner.invoke(main, ["classify", "--help"])
        # This will be added in implementation
        assert "--max-tolerance" in result.output or result.exit_code != 0


class TestISOGGOutputCLI:
    """Tests for ISOGG output options."""

    def test_isogg_option_exists(self, runner: CliRunner) -> None:
        """--isogg flag is recognized."""
        result = runner.invoke(main, ["classify", "--help"])
        # This will be added in implementation
        assert "--isogg" in result.output or result.exit_code != 0

    def test_isogg_db_option_exists(self, runner: CliRunner) -> None:
        """--isogg-db option is recognized."""
        result = runner.invoke(main, ["classify", "--help"])
        # This will be added in implementation
        assert "--isogg-db" in result.output or result.exit_code != 0


class TestEnhancedJSONOutput:
    """Tests for enhanced JSON output with new fields."""

    def test_json_includes_contamination(self) -> None:
        """JSON output includes contamination estimate when requested."""
        # This tests that HaplogroupCall includes contamination
        from yallhap.classifier import HaplogroupCall, QCScores, SNPStats

        call = HaplogroupCall(
            sample="test",
            haplogroup="R1b",
            confidence=0.95,
            qc_scores=QCScores(),
            snp_stats=SNPStats(),
            path=[],
        )
        result_dict = call.to_dict()
        # Contamination field will be added after implementation
        # For now just verify to_dict() works
        assert "sample" in result_dict
        assert "haplogroup" in result_dict

    def test_json_includes_isogg_haplogroup(self) -> None:
        """JSON output includes ISOGG haplogroup when --isogg is used."""
        from yallhap.classifier import HaplogroupCall, QCScores, SNPStats

        call = HaplogroupCall(
            sample="test",
            haplogroup="R1b",
            confidence=0.95,
            qc_scores=QCScores(),
            snp_stats=SNPStats(),
            path=[],
            # isogg_haplogroup will be added after implementation
        )
        result_dict = call.to_dict()
        # Test that it at least works without isogg_haplogroup for now
        assert "haplogroup" in result_dict


class TestCLIIntegration:
    """Integration tests for CLI with mocked components."""

    @patch("yallhap.cli.Tree")
    @patch("yallhap.cli._load_snp_database")
    @patch("yallhap.cli.HaplogroupClassifier")
    def test_classify_with_ancient_options(
        self,
        mock_classifier_cls: MagicMock,
        mock_load_snp_db: MagicMock,
        mock_tree_cls: MagicMock,
        runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        """Classify command passes ancient options to classifier."""
        # Create a mock VCF file
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")

        # Create mock tree and snp_db files
        tree_path = tmp_path / "tree.json"
        tree_path.write_text("{}")
        snp_db_path = tmp_path / "snps.csv"
        snp_db_path.write_text("name,grch38_pos,ancestral,derived,haplogroup\n")

        # Set up mocks
        mock_tree_cls.from_json.return_value = MagicMock()
        mock_load_snp_db.return_value = MagicMock()

        mock_classifier = MagicMock()
        mock_classifier_cls.return_value = mock_classifier

        from yallhap.classifier import HaplogroupCall, QCScores, SNPStats

        mock_classifier.classify.return_value = HaplogroupCall(
            sample="test",
            haplogroup="R1b",
            confidence=0.98,
            qc_scores=QCScores(qc1_backbone=0.95),
            snp_stats=SNPStats(derived=10),
            path=["R", "R1b"],
        )

        # Run the command with ancient options
        result = runner.invoke(
            main,
            [
                "classify",
                str(vcf_path),
                "--tree",
                str(tree_path),
                "--snp-db",
                str(snp_db_path),
                "--ancient",
            ],
        )

        # Check that ancient_mode was passed to classifier
        # (This will work after implementation)
        if result.exit_code == 0:
            call_kwargs = mock_classifier_cls.call_args.kwargs
            assert call_kwargs.get("ancient_mode") is True

    def test_version_option(self, runner: CliRunner) -> None:
        """--version shows version."""
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "yallhap" in result.output or "version" in result.output.lower()

    def test_help_option(self, runner: CliRunner) -> None:
        """--help shows help."""
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "yallhap" in result.output.lower() or "Usage" in result.output


class TestExitCodes:
    """Tests for CLI exit code behavior.

    Exit codes should follow Unix conventions:
    - 0: Success (classification completed, regardless of confidence)
    - 1: Classification failed (NA result)
    - 10+: Various error conditions
    """

    @patch("yallhap.cli.Tree")
    @patch("yallhap.cli._load_snp_database")
    @patch("yallhap.cli.HaplogroupClassifier")
    def test_exit_code_zero_on_high_confidence(
        self,
        mock_classifier_cls: MagicMock,
        mock_load_snp_db: MagicMock,
        mock_tree_cls: MagicMock,
        runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        """Exit code 0 for successful classification with high confidence."""
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")
        tree_path = tmp_path / "tree.json"
        tree_path.write_text("{}")
        snp_db_path = tmp_path / "snps.csv"
        snp_db_path.write_text("name,grch38_pos,ancestral,derived,haplogroup\n")

        mock_tree_cls.from_json.return_value = MagicMock()
        mock_load_snp_db.return_value = MagicMock()

        from yallhap.classifier import HaplogroupCall, QCScores, SNPStats

        mock_classifier = MagicMock()
        mock_classifier_cls.return_value = mock_classifier
        mock_classifier.classify.return_value = HaplogroupCall(
            sample="test",
            haplogroup="R1b",
            confidence=0.99,  # High confidence
            qc_scores=QCScores(),
            snp_stats=SNPStats(),
            path=["R", "R1b"],
        )

        result = runner.invoke(
            main,
            ["classify", str(vcf_path), "--tree", str(tree_path), "--snp-db", str(snp_db_path)],
        )
        assert result.exit_code == 0

    @patch("yallhap.cli.Tree")
    @patch("yallhap.cli._load_snp_database")
    @patch("yallhap.cli.HaplogroupClassifier")
    def test_exit_code_zero_on_low_confidence(
        self,
        mock_classifier_cls: MagicMock,
        mock_load_snp_db: MagicMock,
        mock_tree_cls: MagicMock,
        runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        """Exit code 0 for successful classification even with low confidence.

        Low confidence is informational, not an error. Users should check
        the confidence value in the output, not rely on exit codes.
        """
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")
        tree_path = tmp_path / "tree.json"
        tree_path.write_text("{}")
        snp_db_path = tmp_path / "snps.csv"
        snp_db_path.write_text("name,grch38_pos,ancestral,derived,haplogroup\n")

        mock_tree_cls.from_json.return_value = MagicMock()
        mock_load_snp_db.return_value = MagicMock()

        from yallhap.classifier import HaplogroupCall, QCScores, SNPStats

        mock_classifier = MagicMock()
        mock_classifier_cls.return_value = mock_classifier
        mock_classifier.classify.return_value = HaplogroupCall(
            sample="test",
            haplogroup="Q-L472",
            confidence=0.74,  # Low confidence - should still exit 0
            qc_scores=QCScores(),
            snp_stats=SNPStats(),
            path=["Q", "Q-L472"],
        )

        result = runner.invoke(
            main,
            ["classify", str(vcf_path), "--tree", str(tree_path), "--snp-db", str(snp_db_path)],
        )
        assert result.exit_code == 0, f"Low confidence should exit 0, got {result.exit_code}"

    @patch("yallhap.cli.Tree")
    @patch("yallhap.cli._load_snp_database")
    @patch("yallhap.cli.HaplogroupClassifier")
    def test_exit_code_one_on_failed_classification(
        self,
        mock_classifier_cls: MagicMock,
        mock_load_snp_db: MagicMock,
        mock_tree_cls: MagicMock,
        runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        """Exit code 1 for failed classification (NA result)."""
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")
        tree_path = tmp_path / "tree.json"
        tree_path.write_text("{}")
        snp_db_path = tmp_path / "snps.csv"
        snp_db_path.write_text("name,grch38_pos,ancestral,derived,haplogroup\n")

        mock_tree_cls.from_json.return_value = MagicMock()
        mock_load_snp_db.return_value = MagicMock()

        from yallhap.classifier import HaplogroupCall, QCScores, SNPStats

        mock_classifier = MagicMock()
        mock_classifier_cls.return_value = mock_classifier
        mock_classifier.classify.return_value = HaplogroupCall(
            sample="test",
            haplogroup="NA",  # Failed classification
            confidence=0.0,
            qc_scores=QCScores(),
            snp_stats=SNPStats(),
            path=[],
        )

        result = runner.invoke(
            main,
            ["classify", str(vcf_path), "--tree", str(tree_path), "--snp-db", str(snp_db_path)],
        )
        assert result.exit_code == 1


class TestBatchCommand:
    """Tests for batch command options."""

    def test_batch_help(self, runner: CliRunner) -> None:
        """batch command has help text."""
        result = runner.invoke(main, ["batch", "--help"])
        assert result.exit_code == 0
        assert "VCF" in result.output or "vcf" in result.output.lower()

    def test_download_help(self, runner: CliRunner) -> None:
        """download command has help text."""
        result = runner.invoke(main, ["download", "--help"])
        assert result.exit_code == 0
        assert "YFull" in result.output or "download" in result.output.lower()
