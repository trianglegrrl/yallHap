"""
Unit tests for yclade.snps module.
"""

from pathlib import Path

import pytest

from yclade.snps import SNP, SNPDatabase


class TestSNP:
    """Tests for SNP dataclass."""

    def test_snp_creation(self) -> None:
        """Test basic SNP creation."""
        snp = SNP(
            name="L21",
            aliases=["M529", "S145"],
            position_grch38=2655229,
            ancestral="C",
            derived="T",
            haplogroup="R-L21",
        )

        assert snp.name == "L21"
        assert "M529" in snp.aliases
        assert snp.position_grch38 == 2655229
        assert snp.ancestral == "C"
        assert snp.derived == "T"
        assert snp.haplogroup == "R-L21"

    def test_snp_get_position(self) -> None:
        """Test position lookup by reference."""
        snp = SNP(
            name="L21",
            position_grch37=2887478,
            position_grch38=2655229,
            position_t2t=2700000,
        )

        assert snp.get_position("grch37") == 2887478
        assert snp.get_position("grch38") == 2655229
        assert snp.get_position("t2t") == 2700000

    def test_snp_get_position_unknown_reference(self) -> None:
        """Unknown reference raises ValueError."""
        snp = SNP(name="L21")

        with pytest.raises(ValueError, match="Unknown reference"):
            snp.get_position("unknown")  # type: ignore

    def test_snp_all_names(self) -> None:
        """Test all_names property."""
        snp = SNP(name="L21", aliases=["M529", "S145"])

        all_names = snp.all_names
        assert "L21" in all_names
        assert "M529" in all_names
        assert "S145" in all_names
        assert len(all_names) == 3


class TestSNPDatabase:
    """Tests for SNPDatabase class."""

    def test_database_from_csv(self, sample_snps_csv: Path) -> None:
        """Test loading database from CSV."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        assert len(db) > 0
        assert "L21" in db

    def test_database_get_by_name(self, sample_snps_csv: Path) -> None:
        """Test getting SNP by name."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        snp = db.get_by_name("L21")
        assert snp.name == "L21"
        assert snp.haplogroup == "R-L21"

    def test_database_get_by_alias(self, sample_snps_csv: Path) -> None:
        """Test getting SNP by alias."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        snp = db.get_by_name("M529")  # Alias for L21
        assert snp.name == "L21"

    def test_database_get_missing_raises(self, sample_snps_csv: Path) -> None:
        """Getting missing SNP raises KeyError."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        with pytest.raises(KeyError, match="SNP not found"):
            db.get_by_name("NOT_EXISTS")

    def test_database_get_by_position(self, sample_snps_csv: Path) -> None:
        """Test getting SNPs by position."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        # L21 is at 2655229 in GRCh38
        snps = db.get_by_position(2655229, "grch38")
        assert len(snps) > 0
        assert any(s.name == "L21" for s in snps)

    def test_database_get_by_position_empty(self, sample_snps_csv: Path) -> None:
        """No SNP at position returns empty list."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        snps = db.get_by_position(99999999, "grch38")
        assert snps == []

    def test_database_contains(self, sample_snps_csv: Path) -> None:
        """Test __contains__ operator."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        assert "L21" in db
        assert "M529" in db  # Alias
        assert "NOT_EXISTS" not in db

    def test_database_iter(self, sample_snps_csv: Path) -> None:
        """Test iteration over SNPs."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        snps = list(db)
        assert len(snps) == len(db)
        assert all(isinstance(s, SNP) for s in snps)

    def test_database_positions(self, sample_snps_csv: Path) -> None:
        """Test positions property."""
        db = SNPDatabase.from_csv(sample_snps_csv)

        positions = db.positions
        assert "grch37" in positions
        assert "grch38" in positions
        assert "t2t" in positions
        assert len(positions["grch38"]) > 0


class TestSNPDatabaseEdgeCases:
    """Edge case tests for SNPDatabase."""

    def test_empty_database(self, tmp_path: Path) -> None:
        """Test handling of empty CSV."""
        csv_path = tmp_path / "empty.csv"
        csv_path.write_text("name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup\n")

        db = SNPDatabase.from_csv(csv_path)
        assert len(db) == 0

    def test_missing_positions(self, tmp_path: Path) -> None:
        """Test handling of missing position values."""
        csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
TEST1,,,2655229,C,T,R-L21
"""
        csv_path = tmp_path / "partial.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_csv(csv_path)
        snp = db.get_by_name("TEST1")
        assert snp.position_grch37 is None
        assert snp.position_grch38 == 2655229

    def test_multiple_snps_same_position(self, tmp_path: Path) -> None:
        """Test multiple SNPs at same position."""
        csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
SNP1,,1000,1000,C,T,HG1
SNP2,,1000,1000,C,T,HG2
"""
        csv_path = tmp_path / "multi.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_csv(csv_path)
        snps = db.get_by_position(1000, "grch38")
        assert len(snps) == 2
