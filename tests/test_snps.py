"""
Unit tests for yallhap.snps module.
"""

from pathlib import Path

import pytest

from yallhap.snps import SNP, SNPDatabase


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


class TestYBrowseGFFCSVParsing:
    """Tests for parsing actual YBrowse GFF-style CSV format."""

    def test_from_ybrowse_gff_csv_real_format(self, tmp_path: Path) -> None:
        """Test parsing actual YBrowse CSV format with GFF-style columns."""
        csv_content = '''"seqid","source","type","start","end","score","strand","phase","Name","ID","allele_anc","allele_der","YCC_haplogroup","ISOGG_haplogroup","mutation","count_tested","count_derived","ref","comment"
"chrY","point","snp","2655229","2655229",".","+",".","L21","L21","C","T","R1b","R1b1a2a1a2c","C to T","100","50","Poznik","."'''

        csv_path = tmp_path / "ybrowse.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_ybrowse_gff_csv(csv_path)

        assert len(db) == 1
        snp = db.get_by_name("L21")
        assert snp.position_grch38 == 2655229
        assert snp.ancestral == "C"
        assert snp.derived == "T"
        assert snp.haplogroup == "R1b"

    def test_from_ybrowse_gff_csv_multiple_snps(self, tmp_path: Path) -> None:
        """Test parsing multiple SNPs from YBrowse format."""
        csv_content = '''"seqid","source","type","start","end","score","strand","phase","Name","ID","allele_anc","allele_der","YCC_haplogroup","ISOGG_haplogroup","mutation","count_tested","count_derived","ref","comment"
"chrY","point","snp","2655229","2655229",".","+",".","L21","L21","C","T","R1b","R1b1a2a1a2c","C to T","100","50","Poznik","."
"chrY","point","snp","2571333","2571333",".","+",".","M269","M269","G","A","R1b","R1b1a2","G to A","200","100","Poznik","."
"chrY","point","snp","14696931","14696931",".","+",".","M168","M168","C","T","CT","CT","C to T","50","25","Poznik","."'''

        csv_path = tmp_path / "ybrowse_multi.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_ybrowse_gff_csv(csv_path)

        assert len(db) == 3
        assert "L21" in db
        assert "M269" in db
        assert "M168" in db

        m269 = db.get_by_name("M269")
        assert m269.position_grch38 == 2571333
        assert m269.haplogroup == "R1b"

    def test_from_ybrowse_gff_csv_skips_indels(self, tmp_path: Path) -> None:
        """Test that indels are skipped (type != 'snp' or non-single-base alleles)."""
        csv_content = '''"seqid","source","type","start","end","score","strand","phase","Name","ID","allele_anc","allele_der","YCC_haplogroup","ISOGG_haplogroup","mutation","count_tested","count_derived","ref","comment"
"chrY","point","snp","2655229","2655229",".","+",".","L21","L21","C","T","R1b","R1b1a2a1a2c","C to T","100","50","Poznik","."
"chrY","indel","snp","1012648","1012650",".","+",".","FGC57219","FGC57219","ins","del","unknown","unknown","3T to 2T","0","0","Full genomes","."'''

        csv_path = tmp_path / "ybrowse_indel.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_ybrowse_gff_csv(csv_path)

        # Should only have L21, not the indel
        assert len(db) == 1
        assert "L21" in db
        assert "FGC57219" not in db


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


class TestT2TLiftover:
    """Tests for T2T liftover functionality."""

    @pytest.fixture
    def liftover_chain_path(self) -> Path:
        """Return path to GRCh38-to-T2T chain file."""
        chain_path = Path("data/liftover/grch38-chm13v2.chain")
        if not chain_path.exists():
            pytest.skip("Liftover chain files not downloaded")
        return chain_path

    @pytest.fixture
    def grch37_chain_path(self) -> Path:
        """Return path to GRCh37-to-T2T chain file."""
        chain_path = Path("data/liftover/hg19-chm13v2.chain")
        if not chain_path.exists():
            pytest.skip("Liftover chain files not downloaded")
        return chain_path

    def test_lift_to_t2t_basic(self, tmp_path: Path, liftover_chain_path: Path) -> None:
        """Test basic liftover from GRCh38 to T2T."""
        # Create database with known SNP
        csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
L21,,2887478,2655229,C,T,R-L21
"""
        csv_path = tmp_path / "snps.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_csv(csv_path)

        # Initially no T2T position
        snp = db.get_by_name("L21")
        assert snp.position_t2t is None

        # Lift to T2T
        lifted_count = db.lift_to_t2t(liftover_chain_path, "grch38")

        # Should have lifted the SNP
        assert lifted_count == 1
        snp = db.get_by_name("L21")
        assert snp.position_t2t is not None
        assert snp.position_t2t > 0

    def test_lift_to_t2t_from_grch37(self, tmp_path: Path, grch37_chain_path: Path) -> None:
        """Test liftover from GRCh37 to T2T."""
        csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
L21,,2887478,,C,T,R-L21
"""
        csv_path = tmp_path / "snps.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_csv(csv_path)
        lifted_count = db.lift_to_t2t(grch37_chain_path, "grch37")

        assert lifted_count == 1
        snp = db.get_by_name("L21")
        assert snp.position_t2t is not None

    def test_lift_to_t2t_chain_not_found(self, tmp_path: Path) -> None:
        """Liftover with missing chain file raises FileNotFoundError."""
        csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
L21,,2887478,2655229,C,T,R-L21
"""
        csv_path = tmp_path / "snps.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_csv(csv_path)

        with pytest.raises(FileNotFoundError, match="Chain file not found"):
            db.lift_to_t2t("/nonexistent/chain.chain", "grch38")

    def test_lift_to_t2t_from_t2t_raises(self, tmp_path: Path) -> None:
        """Liftover from T2T to T2T raises ValueError."""
        csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
L21,,2887478,2655229,C,T,R-L21
"""
        csv_path = tmp_path / "snps.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_csv(csv_path)

        with pytest.raises(ValueError, match="Cannot liftover from T2T to T2T"):
            db.lift_to_t2t("dummy.chain", "t2t")

    def test_position_index_updated_after_liftover(
        self, tmp_path: Path, liftover_chain_path: Path
    ) -> None:
        """T2T position index is updated after liftover."""
        csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
L21,,2887478,2655229,C,T,R-L21
"""
        csv_path = tmp_path / "snps.csv"
        csv_path.write_text(csv_content)

        db = SNPDatabase.from_csv(csv_path)
        db.lift_to_t2t(liftover_chain_path, "grch38")

        snp = db.get_by_name("L21")
        t2t_pos = snp.position_t2t
        assert t2t_pos is not None

        # Should be able to find SNP by T2T position
        snps_at_pos = db.get_by_position(t2t_pos, "t2t")
        assert len(snps_at_pos) >= 1
        assert any(s.name == "L21" for s in snps_at_pos)
