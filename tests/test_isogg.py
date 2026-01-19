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

    def test_unmapped_returns_original_with_tilde(self, sample_isogg_file: Path) -> None:
        """YFull nodes without ISOGG match return original name with tilde."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R-XYZ123"],  # Made-up haplogroup not in ISOGG
        }
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # Unknown haplogroup returns original with tilde
        result = mapper.to_isogg("R-XYZ123")
        assert result == "R-XYZ123~"

    def test_multi_haplogroup_snp_prefers_matching_clade(self, sample_isogg_file: Path) -> None:
        """When SNP maps to multiple haplogroups, prefers one matching major clade."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        # L596 maps to both I2a2 and N1a1a1a1a1a6~ in sample file
        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["I"],
            "I": ["I2"],
            "I2": ["I-L596"],
        }
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # I-L596 should map to I2a2 (matches I major clade), not N1a1a1a1a1a6~
        result = mapper.to_isogg("I-L596")
        assert result == "I2a2"

    def test_multi_haplogroup_snp_n_clade(self, sample_isogg_file: Path) -> None:
        """Multi-haplogroup SNP chooses N when that's the major clade."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        # Same L596 SNP but from N context
        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["N"],
            "N": ["N-L596"],
        }
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # N-L596 should map to N1a1a1a1a1a6~ (matches N major clade)
        result = mapper.to_isogg("N-L596")
        assert result == "N1a1a1a1a1a6~"

    def test_get_isogg_haplogroups_at_position(self, sample_isogg_file: Path) -> None:
        """Can get all ISOGG haplogroups that have SNPs at a position."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        # Position 11912037 has SNP A100 for R1b1a1b1a1a2c1a1d5a
        haplogroups = db.get_haplogroups_at_position(11912037)
        assert len(haplogroups) >= 1
        assert "R1b1a1b1a1a2c1a1d5a" in haplogroups

    def test_direct_snp_lookup_priority(self, sample_isogg_file: Path) -> None:
        """Direct SNP name lookup takes priority over tree traversal."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        # A100 maps to R1b1a1b1a1a2c1a1d5a in sample ISOGG file
        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R-A100"],  # SNP A100 in the haplogroup name
        }
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # R-A100 should map to R1b1a1b1a1a2c1a1d5a via direct SNP lookup
        result = mapper.to_isogg("R-A100")
        assert result == "R1b1a1b1a1a2c1a1d5a"

    def test_to_isogg_from_snps_finds_most_specific(self, sample_isogg_file: Path) -> None:
        """to_isogg_from_snps finds the most specific haplogroup from SNP list."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        tree_dict: dict[str, list[str]] = {'ROOT (Y-Chromosome "Adam")': []}
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # A100 -> R1b1a1b1a1a2c1a1d5a (very specific)
        # This should prefer the more specific haplogroup
        result = mapper.to_isogg_from_snps(["A100"])
        assert result is not None
        assert result == "R1b1a1b1a1a2c1a1d5a"

    def test_to_isogg_from_snps_empty_list(self, sample_isogg_file: Path) -> None:
        """to_isogg_from_snps returns None for empty SNP list."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        tree_dict: dict[str, list[str]] = {'ROOT (Y-Chromosome "Adam")': []}
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        result = mapper.to_isogg_from_snps([])
        assert result is None

    def test_to_isogg_from_snps_unknown_snps(self, sample_isogg_file: Path) -> None:
        """to_isogg_from_snps returns None for unknown SNPs."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        tree_dict: dict[str, list[str]] = {'ROOT (Y-Chromosome "Adam")': []}
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        result = mapper.to_isogg_from_snps(["FAKE_SNP_123", "UNKNOWN_456"])
        assert result is None

    def test_base_haplogroup_fallback(self, sample_isogg_file: Path) -> None:
        """Falls back to base haplogroup letter if in ISOGG database."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["I2"],
            "I2": [],
        }
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # I2 is in the ISOGG database as a haplogroup name
        result = mapper.to_isogg("I2")
        # Should return I2 or a more specific I2 haplogroup
        assert result is not None
        assert "I2" in result or result.startswith("I2")

    def test_direct_snp_lookup_uses_major_clade_hint(self, sample_isogg_file: Path) -> None:
        """Direct SNP lookup in to_isogg uses major clade from YFull name."""
        from yallhap.isogg import ISOGGDatabase, ISOGGMapper
        from yallhap.tree import Tree

        tree_dict: dict[str, list[str]] = {'ROOT (Y-Chromosome "Adam")': []}
        tree = Tree.from_dict(tree_dict)

        db = ISOGGDatabase.from_file(sample_isogg_file)
        mapper = ISOGGMapper(tree, db)

        # M253 is I1-defining SNP
        result = mapper.to_isogg("I-M253")
        assert result == "I1"

        # M269 is R1b1a1b-defining SNP
        result = mapper.to_isogg("R-M269")
        assert result == "R1b1a1b"


class TestISOGGDatabaseMultiHaplogroup:
    """Tests for ISOGG database handling of multi-haplogroup SNPs."""

    def test_get_all_by_name_returns_multiple(self, sample_isogg_file: Path) -> None:
        """get_all_by_name returns all haplogroups for recurrent SNPs."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        # L596 is in both I2a2 and N1a1a1a1a1a6~
        snps = db.get_all_by_name("L596")
        assert len(snps) == 2
        haplogroups = {s.haplogroup for s in snps}
        assert "I2a2" in haplogroups
        assert "N1a1a1a1a1a6~" in haplogroups

    def test_get_by_name_returns_first(self, sample_isogg_file: Path) -> None:
        """get_by_name returns first SNP for multi-haplogroup entries."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        snp = db.get_by_name("L596")
        assert snp is not None
        # Returns first one encountered
        assert snp.haplogroup in ("I2a2", "N1a1a1a1a1a6~")

    def test_len_counts_all_entries(self, sample_isogg_file: Path) -> None:
        """len() counts all SNP entries including duplicates."""
        from yallhap.isogg import ISOGGDatabase

        db = ISOGGDatabase.from_file(sample_isogg_file)
        # Should count L596 twice (once for each haplogroup)
        total = len(db)
        assert total > 100  # Sample file has ~100+ entries


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
