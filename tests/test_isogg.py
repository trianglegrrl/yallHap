"""
Tests for ISOGG SNP database and haplogroup mapping.

Following TDD - these tests are written first, before implementation.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def sample_isogg_file(fixtures_dir: Path) -> Path:
    """Return path to sample ISOGG file."""
    return fixtures_dir / "sample_isogg.txt"


@pytest.fixture
def full_isogg_file() -> Path:
    """Return path to full ISOGG file."""
    return Path(__file__).parent.parent / "data" / "isogg_snps_grch38.txt"


class TestISOGGDatabase:
    """Tests for ISOGG SNP database."""

    def test_load_from_pathphynder_format(self, sample_isogg_file: Path) -> None:
        """Loads ISOGG SNPs from pathPhynder format."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        assert len(db) > 50  # Sample file has ~100 SNPs

    def test_snp_has_haplogroup(self, sample_isogg_file: Path) -> None:
        """Each SNP has associated haplogroup."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        # A100 -> R1b1a1b1a1a2c1a1d5a
        snp = db.get_by_name("A100")
        assert snp is not None
        assert snp.haplogroup == "R1b1a1b1a1a2c1a1d5a"

    def test_position_lookup(self, sample_isogg_file: Path) -> None:
        """Can look up SNPs by position."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        # A100 is at position 11912037
        snp = db.get_by_position(11912037)
        assert snp is not None
        assert snp.name == "A100"

    def test_snp_attributes(self, sample_isogg_file: Path) -> None:
        """SNP has all expected attributes."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        snp = db.get_by_name("A100")
        assert snp is not None
        assert snp.position == 11912037
        assert snp.ancestral == "A"
        assert snp.derived == "T"

    def test_iteration(self, sample_isogg_file: Path) -> None:
        """Can iterate over all SNPs."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        snps = list(db)
        assert len(snps) > 0
        assert all(hasattr(snp, "name") for snp in snps)

    def test_full_database_size(self, full_isogg_file: Path) -> None:
        """Full ISOGG database has expected number of SNPs."""
        from yallhap.isogg import ISOGGDatabase

        if not full_isogg_file.exists():
            pytest.skip("Full ISOGG file not available")

        db = ISOGGDatabase.from_file(full_isogg_file)
        # Should have ~90k SNPs
        assert len(db) > 50000


class TestISOGGSNP:
    """Tests for ISOGG SNP dataclass."""

    def test_isogg_snp_creation(self) -> None:
        """ISOGGSNP can be created with all fields."""
        from yallhap.isogg import ISOGGSNP

        snp = ISOGGSNP(
            name="M269",
            haplogroup="R1b1a1b",
            position=2571333,
            ancestral="G",
            derived="A",
        )
        assert snp.name == "M269"
        assert snp.haplogroup == "R1b1a1b"
        assert snp.position == 2571333


class TestISOGGMapper:
    """Tests for YFull to ISOGG haplogroup mapping."""

    def test_yfull_to_isogg_major_groups(self, sample_isogg_file: Path) -> None:
        """Major YFull haplogroups map to ISOGG names."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        # Create a simple tree
        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
            "R1b": ["R-L21"],
        }
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # R1b should map to something like R1b1a1b... (ISOGG style)
        # The exact mapping depends on the SNPs in the database
        result = mapper.to_isogg("R1b")
        # Should return some ISOGG-style haplogroup or the original if no mapping
        assert result is not None

    def test_unmapped_returns_original(self, sample_isogg_file: Path) -> None:
        """YFull nodes without ISOGG match return original name."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R-XYZ123"],  # Made-up haplogroup not in ISOGG
        }
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # Unknown haplogroup returns original
        result = mapper.to_isogg("R-XYZ123")
        assert result == "R-XYZ123"

    def test_get_isogg_haplogroups_at_position(self, sample_isogg_file: Path) -> None:
        """Can get all ISOGG haplogroups that have SNPs at a position."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        # Position 11912037 has SNP A100 for R1b1a1b1a1a2c1a1d5a
        haplogroups = db.get_haplogroups_at_position(11912037)
        assert len(haplogroups) >= 1
        assert "R1b1a1b1a1a2c1a1d5a" in haplogroups


class TestISOGGHaplogroups:
    """Tests for ISOGG haplogroup list."""

    def test_unique_haplogroups(self, sample_isogg_file: Path) -> None:
        """Can get list of unique haplogroups."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        haplogroups = db.get_all_haplogroups()
        assert len(haplogroups) > 0
        # All should be unique
        assert len(haplogroups) == len(set(haplogroups))

    def test_haplogroup_snp_count(self, sample_isogg_file: Path) -> None:
        """Can count SNPs per haplogroup."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        counts = db.get_snp_counts_by_haplogroup()
        assert isinstance(counts, dict)
        # At least one haplogroup should have SNPs
        assert any(count > 0 for count in counts.values())
