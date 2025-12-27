"""
Integration tests for enhanced ancient mode.

These tests verify that all the new features work together correctly.
"""

from __future__ import annotations

from pathlib import Path

import pytest


class TestAncientModeIntegration:
    """Integration tests for ancient DNA mode."""

    def test_contamination_module_imports(self) -> None:
        """Contamination module imports successfully."""
        from yallhap.contamination import (
            ContaminationResult,
            estimate_contamination,
            estimate_contamination_with_snpdb,
            is_contaminated,
        )

        assert ContaminationResult is not None
        assert estimate_contamination is not None
        assert estimate_contamination_with_snpdb is not None
        assert is_contaminated is not None

    def test_path_traversal_module_imports(self) -> None:
        """Path traversal function imports successfully."""
        from yallhap.classifier import traverse_with_tolerance

        assert traverse_with_tolerance is not None

    def test_coalescent_prior_imports(self) -> None:
        """Coalescent prior function imports successfully."""
        from yallhap.bayesian import compute_coalescent_prior

        assert compute_coalescent_prior is not None

    def test_isogg_module_imports(self) -> None:
        """ISOGG module imports successfully."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper, ISOGGSNP

        assert ISOGGDatabase is not None
        assert ISOGGMapper is not None
        assert ISOGGSNP is not None

    def test_isogg_database_loads(self) -> None:
        """ISOGG database loads from data file."""
        from yallhap.isogg import ISOGGDatabase

        isogg_path = Path(__file__).parent.parent / "data" / "isogg_snps_grch38.txt"
        if not isogg_path.exists():
            pytest.skip("ISOGG database file not found")

        db = ISOGGDatabase.from_file(isogg_path)
        assert len(db) > 50000

    def test_full_workflow_contamination(self) -> None:
        """Test contamination estimation on simulated data."""
        from yallhap.contamination import ContaminationResult, estimate_contamination
        from yallhap.vcf import Variant

        # Simulate a clean sample (no contamination)
        clean_variants = {
            1000: Variant(
                chrom="Y",
                position=1000,
                ref="A",
                alt=("G",),
                genotype=0,
                depth=10,
                quality=30,
                allele_depth=(10, 0),
            ),
            2000: Variant(
                chrom="Y",
                position=2000,
                ref="C",
                alt=("T",),
                genotype=1,
                depth=10,
                quality=30,
                allele_depth=(0, 10),
            ),
        }

        rate, n_sites = estimate_contamination(clean_variants, min_depth=1)
        assert rate < 0.05
        assert n_sites == 2

        # Simulate a contaminated sample
        contaminated_variants = {
            1000: Variant(
                chrom="Y",
                position=1000,
                ref="A",
                alt=("G",),
                genotype=0,
                depth=10,
                quality=30,
                allele_depth=(7, 3),  # 30% minor allele
            ),
        }

        rate, n_sites = estimate_contamination(contaminated_variants, min_depth=1)
        assert rate > 0.20
        assert n_sites == 1

    def test_full_workflow_path_traversal(self) -> None:
        """Test path traversal with tolerance on tree."""
        from yallhap.classifier import traverse_with_tolerance
        from yallhap.tree import Tree

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
            "R1b": ["R-L21"],
        }
        tree = Tree.from_dict(tree_dict)

        # Good path - should go all the way
        good_scores = {
            "R": {"derived": 5, "ancestral": 0, "missing": 0},
            "R1": {"derived": 3, "ancestral": 0, "missing": 0},
            "R1b": {"derived": 2, "ancestral": 0, "missing": 0},
            "R-L21": {"derived": 1, "ancestral": 0, "missing": 0},
        }
        hg, path, stats = traverse_with_tolerance(tree, good_scores, max_tolerance=3)
        assert hg == "R-L21"
        assert len(path) == 4

        # Blocked path - should stop early
        blocked_scores = {
            "R": {"derived": 5, "ancestral": 0, "missing": 0},
            "R1": {"derived": 3, "ancestral": 5, "missing": 0},  # Too many ancestral
        }
        hg, path, stats = traverse_with_tolerance(tree, blocked_scores, max_tolerance=3)
        assert hg == "R"  # Stops before R1

    def test_full_workflow_isogg_mapping(self) -> None:
        """Test ISOGG mapping from YFull to ISOGG names."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        # Create tree
        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
            "R1b": ["R-L21"],
        }
        tree = Tree.from_dict(tree_dict)

        # Create minimal ISOGG database
        isogg_path = Path(__file__).parent / "fixtures" / "sample_isogg.txt"
        if not isogg_path.exists():
            pytest.skip("Sample ISOGG file not found")

        db = ISOGGDatabase.from_file(isogg_path)
        mapper = ISOGGMapper(tree, db)

        # Unmapped haplogroup returns original
        result = mapper.to_isogg("R-XYZ")
        assert result == "R-XYZ"


class TestCLIAncientModeIntegration:
    """Integration tests for CLI ancient mode options."""

    def test_cli_has_all_options(self) -> None:
        """CLI classify command has all enhanced ancient mode options."""
        from click.testing import CliRunner

        from yallhap.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["classify", "--help"])

        assert result.exit_code == 0
        assert "--ancient" in result.output
        assert "--transversions-only" in result.output
        assert "--estimate-contamination" in result.output
        assert "--max-tolerance" in result.output
        assert "--isogg" in result.output
        assert "--isogg-db" in result.output
        assert "--bayesian" in result.output


@pytest.mark.slow
class TestAADRValidation:
    """Validation tests for AADR ancient DNA samples."""

    def test_aadr_validation_files_exist(self) -> None:
        """Required AADR validation files exist."""
        base_path = Path(__file__).parent.parent / "data" / "ancient"

        files_needed = [
            "aadr_chrY_v2.vcf.gz",
            "aadr_1240k_ground_truth.tsv",
        ]

        for filename in files_needed:
            filepath = base_path / filename
            if not filepath.exists():
                pytest.skip(f"AADR file not found: {filepath}")

    def test_yallhap_tree_and_snpdb_exist(self) -> None:
        """Required yallhap data files exist."""
        base_path = Path(__file__).parent.parent / "data"

        if not (base_path / "yfull_tree.json").exists():
            pytest.skip("YFull tree not found")
        if not (base_path / "ybrowse_snps.csv").exists():
            pytest.skip("YBrowse SNPs not found")
