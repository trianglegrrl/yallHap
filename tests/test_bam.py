"""
Tests for BAM file reading for Y-chromosome haplogroup classification.

TDD tests for BAMReader class that reads directly from BAM files,
enabling Bayesian classification without requiring pre-called VCF files.
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from yallhap.vcf import Variant


class TestBAMReaderBasics:
    """Tests for basic BAMReader functionality."""

    @pytest.fixture
    def simple_bam(self, tmp_path: Path) -> Path:
        """Create a simple BAM file with Y chromosome reads for testing."""
        bam_path = tmp_path / "test.bam"

        # Create BAM header with Y chromosome
        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
                "RG": [{"ID": "rg1", "SM": "TestSample"}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            # Create a read at position 1000 with base 'T'
            read1 = pysam.AlignedSegment()
            read1.query_name = "read1"
            read1.query_sequence = "ACGTACGTACGTACGTACGT"
            read1.flag = 0
            read1.reference_id = 0  # Y chromosome
            read1.reference_start = 990  # 0-based, so covers positions 991-1010
            read1.mapping_quality = 60
            read1.cigar = [(0, 20)]  # 20M
            read1.query_qualities = pysam.qualitystring_to_array("I" * 20)
            read1.set_tag("RG", "rg1")
            bam.write(read1)

            # Create another read at same position
            read2 = pysam.AlignedSegment()
            read2.query_name = "read2"
            read2.query_sequence = "ACGTACGTACGTACGTACGT"
            read2.flag = 0
            read2.reference_id = 0
            read2.reference_start = 990
            read2.mapping_quality = 60
            read2.cigar = [(0, 20)]
            read2.query_qualities = pysam.qualitystring_to_array("I" * 20)
            read2.set_tag("RG", "rg1")
            bam.write(read2)

        # Index the BAM
        pysam.index(str(bam_path))

        return bam_path

    def test_bam_reader_opens_file(self, simple_bam: Path) -> None:
        """BAMReader successfully opens a BAM file."""
        from yallhap.bam import BAMReader

        with BAMReader(simple_bam) as reader:
            assert reader._bam is not None

    def test_bam_reader_detects_y_chromosome(self, simple_bam: Path) -> None:
        """BAMReader detects Y chromosome naming convention."""
        from yallhap.bam import BAMReader

        with BAMReader(simple_bam) as reader:
            assert reader.y_chrom == "Y"

    def test_bam_reader_detects_sample_name(self, simple_bam: Path) -> None:
        """BAMReader extracts sample name from read groups."""
        from yallhap.bam import BAMReader

        with BAMReader(simple_bam) as reader:
            assert reader.sample == "TestSample"

    def test_bam_reader_falls_back_to_filename(self, tmp_path: Path) -> None:
        """BAMReader uses filename as sample when no read groups."""
        bam_path = tmp_path / "my_sample.bam"

        # Create BAM without read groups
        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header):
            pass  # Empty BAM

        pysam.index(str(bam_path))

        from yallhap.bam import BAMReader

        with BAMReader(bam_path) as reader:
            assert reader.sample == "my_sample"


class TestBAMReaderYChromDetection:
    """Tests for Y chromosome detection with various naming conventions."""

    @pytest.fixture
    def make_bam_with_chrom(self, tmp_path: Path):
        """Factory fixture to create BAM with specific chromosome name."""

        def _make_bam(chrom_name: str) -> Path:
            bam_path = tmp_path / f"test_{chrom_name}.bam"

            header = pysam.AlignmentHeader.from_dict(
                {
                    "HD": {"VN": "1.6", "SO": "coordinate"},
                    "SQ": [{"SN": chrom_name, "LN": 59373566}],
                }
            )

            with pysam.AlignmentFile(str(bam_path), "wb", header=header):
                pass

            pysam.index(str(bam_path))
            return bam_path

        return _make_bam

    def test_detects_Y(self, make_bam_with_chrom) -> None:
        """Detects 'Y' chromosome name."""
        from yallhap.bam import BAMReader

        bam_path = make_bam_with_chrom("Y")
        with BAMReader(bam_path) as reader:
            assert reader.y_chrom == "Y"

    def test_detects_chrY(self, make_bam_with_chrom) -> None:
        """Detects 'chrY' chromosome name."""
        from yallhap.bam import BAMReader

        bam_path = make_bam_with_chrom("chrY")
        with BAMReader(bam_path) as reader:
            assert reader.y_chrom == "chrY"

    def test_raises_on_missing_y(self, make_bam_with_chrom) -> None:
        """Raises error when no Y chromosome found."""
        from yallhap.bam import BAMReader

        bam_path = make_bam_with_chrom("chr1")
        with pytest.raises(ValueError, match="No Y chromosome found"), BAMReader(bam_path):
            pass


class TestPileupAtPosition:
    """Tests for pileup_at_position method."""

    @pytest.fixture
    def bam_with_coverage(self, tmp_path: Path) -> Path:
        """Create BAM with known coverage at specific positions."""
        bam_path = tmp_path / "coverage.bam"

        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
                "RG": [{"ID": "rg1", "SM": "TestSample"}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            # Position 1000 (1-based): 3 reads with 'T', 1 read with 'C'
            # Read sequences designed so position 1000 has specific bases
            # Position 1000 is at index 9 in a read starting at 991

            # 3 reads with T at position 1000
            for i in range(3):
                read = pysam.AlignedSegment()
                read.query_name = f"read_T_{i}"
                # Position 1000 (1-based) = index 9 in read (991 + 9 = 1000)
                read.query_sequence = "AAAAAAAAAT" + "A" * 10  # T at position 10 (0-based index 9)
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 990  # 0-based, position 991 in 1-based
                read.mapping_quality = 60
                read.cigar = [(0, 20)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                read.set_tag("RG", "rg1")
                bam.write(read)

            # 1 read with C at position 1000
            read = pysam.AlignedSegment()
            read.query_name = "read_C_0"
            read.query_sequence = (
                "AAAAAAAAAC" + "A" * 10
            )  # C at index 9 (0-based) for position 1000
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 990
            read.mapping_quality = 60
            read.cigar = [(0, 20)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 20)
            read.set_tag("RG", "rg1")
            bam.write(read)

            # Position 2000: No coverage (no reads)

            # Position 3000: 5 reads all with 'G'
            for i in range(5):
                read = pysam.AlignedSegment()
                read.query_name = f"read_G_{i}"
                read.query_sequence = "AAAAAAAAAG" + "A" * 10  # G at position 10 (0-based index 9)
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 2990  # Position 2991-3010
                read.mapping_quality = 60
                read.cigar = [(0, 20)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                read.set_tag("RG", "rg1")
                bam.write(read)

        pysam.index(str(bam_path))
        return bam_path

    def test_pileup_returns_allele_counts(self, bam_with_coverage: Path) -> None:
        """pileup_at_position returns correct allele counts."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_with_coverage) as reader:
            result = reader.pileup_at_position(1000)

        assert result is not None
        assert result.position == 1000
        assert result.allele_counts.get("T", 0) == 3
        assert result.allele_counts.get("C", 0) == 1
        assert result.total_depth == 4

    def test_pileup_returns_none_for_no_coverage(self, bam_with_coverage: Path) -> None:
        """pileup_at_position returns None when no coverage."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_with_coverage) as reader:
            result = reader.pileup_at_position(2000)

        assert result is None

    def test_pileup_at_different_position(self, bam_with_coverage: Path) -> None:
        """pileup_at_position works at different positions."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_with_coverage) as reader:
            result = reader.pileup_at_position(3000)

        assert result is not None
        assert result.allele_counts.get("G", 0) == 5
        assert result.total_depth == 5


class TestPileupAtPositions:
    """Tests for pileup_at_positions method (batch query)."""

    @pytest.fixture
    def bam_multi_positions(self, tmp_path: Path) -> Path:
        """Create BAM with coverage at multiple positions."""
        bam_path = tmp_path / "multi.bam"

        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
                "RG": [{"ID": "rg1", "SM": "TestSample"}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            # Create reads spanning positions 1000, 1005, 1010
            for i in range(10):
                read = pysam.AlignedSegment()
                read.query_name = f"read_{i}"
                # T at pos 1000 (index 4), C at pos 1005 (index 9), G at pos 1010 (index 14)
                read.query_sequence = "AAAATAAAACAAAAGA" + "A" * 4
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 995  # 0-based, positions 996-1015
                read.mapping_quality = 60
                read.cigar = [(0, 20)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                read.set_tag("RG", "rg1")
                bam.write(read)

        pysam.index(str(bam_path))
        return bam_path

    def test_pileup_at_multiple_positions(self, bam_multi_positions: Path) -> None:
        """pileup_at_positions returns results for all covered positions."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_multi_positions) as reader:
            positions = {1000, 1005, 1010}
            results = list(reader.pileup_at_positions(positions))

        assert len(results) == 3
        positions_found = {r.position for r in results}
        assert positions_found == {1000, 1005, 1010}

    def test_pileup_skips_uncovered_positions(self, bam_multi_positions: Path) -> None:
        """pileup_at_positions skips positions with no coverage."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_multi_positions) as reader:
            positions = {1000, 2000, 3000}  # 2000 and 3000 have no coverage
            results = list(reader.pileup_at_positions(positions))

        assert len(results) == 1
        assert results[0].position == 1000

    def test_pileup_empty_positions(self, bam_multi_positions: Path) -> None:
        """pileup_at_positions handles empty position set."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_multi_positions) as reader:
            results = list(reader.pileup_at_positions(set()))

        assert len(results) == 0


class TestGetVariantAtPosition:
    """Tests for converting pileup to Variant objects."""

    @pytest.fixture
    def bam_for_variants(self, tmp_path: Path) -> Path:
        """Create BAM for testing variant conversion."""
        bam_path = tmp_path / "variants.bam"

        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
                "RG": [{"ID": "rg1", "SM": "TestSample"}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            # Position 1000: 8 T reads, 2 C reads (T is derived)
            for i in range(8):
                read = pysam.AlignedSegment()
                read.query_name = f"read_T_{i}"
                read.query_sequence = "AAAAAAAAATAAAAAAAAAA"  # T at index 9
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 990
                read.mapping_quality = 60
                read.cigar = [(0, 20)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                read.set_tag("RG", "rg1")
                bam.write(read)

            for i in range(2):
                read = pysam.AlignedSegment()
                read.query_name = f"read_C_{i}"
                read.query_sequence = "AAAAAAAAACAAAAAAAAAA"  # C at index 9
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 990
                read.mapping_quality = 60
                read.cigar = [(0, 20)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                read.set_tag("RG", "rg1")
                bam.write(read)

            # Position 2000: 10 G reads (all reference)
            for i in range(10):
                read = pysam.AlignedSegment()
                read.query_name = f"read_G_{i}"
                read.query_sequence = "AAAAAAAAAGAAAAAAAAAA"  # G at index 9
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 1990
                read.mapping_quality = 60
                read.cigar = [(0, 20)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                read.set_tag("RG", "rg1")
                bam.write(read)

        pysam.index(str(bam_path))
        return bam_path

    def test_get_variant_returns_variant_object(self, bam_for_variants: Path) -> None:
        """get_variant_at_position returns a Variant object."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_for_variants) as reader:
            variant = reader.get_variant_at_position(1000, ref="C", alt="T")

        assert variant is not None
        assert isinstance(variant, Variant)
        assert variant.position == 1000
        assert variant.ref == "C"
        assert variant.alt == ("T",)

    def test_get_variant_allele_depth(self, bam_for_variants: Path) -> None:
        """get_variant_at_position correctly sets allele_depth."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_for_variants) as reader:
            variant = reader.get_variant_at_position(1000, ref="C", alt="T")

        assert variant is not None
        assert variant.allele_depth == (2, 8)  # 2 C (ref), 8 T (alt)
        assert variant.depth == 10

    def test_get_variant_genotype_alt_majority(self, bam_for_variants: Path) -> None:
        """Genotype is alt (1) when alt reads are majority."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_for_variants) as reader:
            variant = reader.get_variant_at_position(1000, ref="C", alt="T")

        assert variant is not None
        assert variant.genotype == 1  # Alt called (8 T > 2 C)

    def test_get_variant_genotype_ref_majority(self, bam_for_variants: Path) -> None:
        """Genotype is ref (0) when ref reads are majority."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_for_variants) as reader:
            variant = reader.get_variant_at_position(2000, ref="G", alt="A")

        assert variant is not None
        assert variant.genotype == 0  # Ref called (10 G, 0 A)
        assert variant.allele_depth == (10, 0)

    def test_get_variant_none_for_no_coverage(self, bam_for_variants: Path) -> None:
        """get_variant_at_position returns None when no coverage."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_for_variants) as reader:
            variant = reader.get_variant_at_position(5000, ref="A", alt="G")

        assert variant is None


class TestGetVariantsAtSNPPositions:
    """Tests for batch variant retrieval at SNP positions."""

    @pytest.fixture
    def bam_for_snps(self, tmp_path: Path) -> Path:
        """Create BAM for testing SNP position queries."""
        bam_path = tmp_path / "snps.bam"

        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
                "RG": [{"ID": "rg1", "SM": "TestSample"}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            # Create long reads covering multiple SNP positions
            # Positions: 1000 (T), 1010 (G), 1020 (C)
            for i in range(5):
                read = pysam.AlignedSegment()
                read.query_name = f"read_{i}"
                # Build sequence with specific bases at positions
                # Position 1000 at index 9, 1010 at index 19, 1020 at index 29
                seq = "A" * 9 + "T" + "A" * 9 + "G" + "A" * 9 + "C" + "A" * 22  # 52 chars
                read.query_sequence = seq
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 990  # 0-based
                read.mapping_quality = 60
                read.cigar = [(0, len(seq))]
                read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                read.set_tag("RG", "rg1")
                bam.write(read)

        pysam.index(str(bam_path))
        return bam_path

    def test_get_variants_at_snp_positions(self, bam_for_snps: Path) -> None:
        """get_variants_at_snp_positions returns dict of variants."""
        from yallhap.bam import BAMReader

        snp_positions = {
            1000: ("C", "T"),  # ref C, alt T
            1010: ("A", "G"),  # ref A, alt G
            1020: ("G", "C"),  # ref G, alt C
        }

        with BAMReader(bam_for_snps) as reader:
            variants = reader.get_variants_at_snp_positions(snp_positions)

        assert len(variants) == 3
        assert 1000 in variants
        assert 1010 in variants
        assert 1020 in variants

    def test_get_variants_correct_allele_depths(self, bam_for_snps: Path) -> None:
        """Variants have correct allele depths from BAM."""
        from yallhap.bam import BAMReader

        snp_positions = {
            1000: ("C", "T"),  # All reads have T
        }

        with BAMReader(bam_for_snps) as reader:
            variants = reader.get_variants_at_snp_positions(snp_positions)

        variant = variants[1000]
        assert variant.allele_depth[1] == 5  # 5 alt (T) reads
        assert variant.allele_depth[0] == 0  # 0 ref (C) reads

    def test_get_variants_skips_no_coverage(self, bam_for_snps: Path) -> None:
        """Positions without coverage are not in result dict."""
        from yallhap.bam import BAMReader

        snp_positions = {
            1000: ("C", "T"),  # Has coverage
            5000: ("A", "G"),  # No coverage
        }

        with BAMReader(bam_for_snps) as reader:
            variants = reader.get_variants_at_snp_positions(snp_positions)

        assert 1000 in variants
        assert 5000 not in variants


class TestReadBamVariantsConvenience:
    """Tests for the read_bam_variants convenience function."""

    @pytest.fixture
    def simple_bam_for_convenience(self, tmp_path: Path) -> Path:
        """Create simple BAM for convenience function test."""
        bam_path = tmp_path / "simple.bam"

        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
                "RG": [{"ID": "rg1", "SM": "ConvenienceSample"}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "AAAAAAAAATAAAAAAAAAA"  # T at position 1000
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 990
            read.mapping_quality = 60
            read.cigar = [(0, 20)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 20)
            read.set_tag("RG", "rg1")
            bam.write(read)

        pysam.index(str(bam_path))
        return bam_path

    def test_read_bam_variants_returns_tuple(self, simple_bam_for_convenience: Path) -> None:
        """read_bam_variants returns (sample_name, variants_dict)."""
        from yallhap.bam import read_bam_variants

        snp_positions = {1000: ("C", "T")}
        sample, variants = read_bam_variants(simple_bam_for_convenience, snp_positions)

        assert sample == "ConvenienceSample"
        assert isinstance(variants, dict)
        assert 1000 in variants


class TestBAMReaderQualityFiltering:
    """Tests for quality filtering in BAM reading."""

    @pytest.fixture
    def bam_with_quality_variation(self, tmp_path: Path) -> Path:
        """Create BAM with reads of varying quality."""
        bam_path = tmp_path / "quality.bam"

        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
                "RG": [{"ID": "rg1", "SM": "TestSample"}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            # High quality read with T
            read1 = pysam.AlignedSegment()
            read1.query_name = "high_qual"
            read1.query_sequence = "AAAAAAAAATAAAAAAAAAA"
            read1.flag = 0
            read1.reference_id = 0
            read1.reference_start = 990
            read1.mapping_quality = 60
            read1.cigar = [(0, 20)]
            read1.query_qualities = pysam.qualitystring_to_array("I" * 20)  # High quality
            read1.set_tag("RG", "rg1")
            bam.write(read1)

            # Low quality read with C
            read2 = pysam.AlignedSegment()
            read2.query_name = "low_qual"
            read2.query_sequence = "AAAAAAAAACAAAAAAAAAA"
            read2.flag = 0
            read2.reference_id = 0
            read2.reference_start = 990
            read2.mapping_quality = 60
            read2.cigar = [(0, 20)]
            read2.query_qualities = pysam.qualitystring_to_array("!" * 20)  # Low quality (0)
            read2.set_tag("RG", "rg1")
            bam.write(read2)

        pysam.index(str(bam_path))
        return bam_path

    def test_filters_low_base_quality(self, bam_with_quality_variation: Path) -> None:
        """Low base quality reads are filtered out."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_with_quality_variation, min_base_quality=20) as reader:
            result = reader.pileup_at_position(1000)

        assert result is not None
        # Only high quality T read should be counted
        assert result.allele_counts.get("T", 0) == 1
        assert result.allele_counts.get("C", 0) == 0
        assert result.total_depth == 1

    def test_includes_all_with_low_threshold(self, bam_with_quality_variation: Path) -> None:
        """All reads included when quality threshold is 0."""
        from yallhap.bam import BAMReader

        with BAMReader(bam_with_quality_variation, min_base_quality=0) as reader:
            result = reader.pileup_at_position(1000)

        assert result is not None
        assert result.total_depth == 2  # Both reads included


class TestPileupResultHelpers:
    """Tests for PileupResult helper methods."""

    def test_get_allele_depth(self) -> None:
        """get_allele_depth returns correct depths for ref and alt."""
        from yallhap.bam import PileupResult

        result = PileupResult(
            position=1000,
            ref_allele="C",
            allele_counts={"C": 5, "T": 15, "A": 2},
            total_depth=22,
        )

        ref_depth, alt_depth = result.get_allele_depth("C", "T")
        assert ref_depth == 5
        assert alt_depth == 15

    def test_get_allele_depth_missing_allele(self) -> None:
        """get_allele_depth returns 0 for missing alleles."""
        from yallhap.bam import PileupResult

        result = PileupResult(
            position=1000,
            ref_allele="C",
            allele_counts={"C": 10},
            total_depth=10,
        )

        ref_depth, alt_depth = result.get_allele_depth("C", "T")
        assert ref_depth == 10
        assert alt_depth == 0

    def test_get_allele_depth_case_insensitive(self) -> None:
        """get_allele_depth is case insensitive."""
        from yallhap.bam import PileupResult

        result = PileupResult(
            position=1000,
            ref_allele="C",
            allele_counts={"C": 5, "T": 15},
            total_depth=20,
        )

        ref_depth, alt_depth = result.get_allele_depth("c", "t")
        assert ref_depth == 5
        assert alt_depth == 15


class TestClassifierBAMIntegration:
    """Tests for HaplogroupClassifier BAM support."""

    def test_classifier_has_classify_from_bam_method(self) -> None:
        """HaplogroupClassifier has classify_from_bam method."""
        from yallhap.classifier import HaplogroupClassifier

        assert hasattr(HaplogroupClassifier, "classify_from_bam")

    def test_classifier_accepts_bam_path(self, tmp_path: Path) -> None:
        """Classifier can classify directly from BAM file."""
        import json

        import pysam

        from yallhap.classifier import HaplogroupClassifier
        from yallhap.snps import SNPDatabase
        from yallhap.tree import Tree

        # Create minimal tree
        tree_data = {
            'ROOT (Y-Chromosome "Adam")': ["A0-T"],
            "A0-T": ["BT"],
            "BT": ["CT"],
            "CT": [],
        }
        tree_path = tmp_path / "tree.json"
        with open(tree_path, "w") as f:
            json.dump(tree_data, f)

        # Create minimal SNP database with SNP at position 1000
        csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
TEST1,,1000,1000,C,T,CT
"""
        snp_path = tmp_path / "snps.csv"
        snp_path.write_text(csv_content)

        # Create BAM with derived allele at position 1000
        bam_path = tmp_path / "test.bam"
        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "Y", "LN": 59373566}],
                "RG": [{"ID": "rg1", "SM": "TestSample"}],
            }
        )

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            for i in range(10):
                read = pysam.AlignedSegment()
                read.query_name = f"read_{i}"
                read.query_sequence = "AAAAAAAAATAAAAAAAAAA"  # T at position 1000
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 990
                read.mapping_quality = 60
                read.cigar = [(0, 20)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                read.set_tag("RG", "rg1")
                bam.write(read)

        pysam.index(str(bam_path))

        # Create classifier
        tree = Tree.from_json(tree_path)
        snp_db = SNPDatabase.from_csv(snp_path)
        classifier = HaplogroupClassifier(
            tree=tree,
            snp_db=snp_db,
            reference="grch37",
            min_depth=1,
            min_quality=0,
        )

        # Classify from BAM
        result = classifier.classify_from_bam(bam_path)

        assert result is not None
        assert result.sample == "TestSample"


class TestCLIBAMSupport:
    """Tests for CLI BAM file support."""

    def test_cli_detects_bam_extension(self) -> None:
        """CLI recognizes .bam extension."""
        from yallhap.cli import _is_bam_file

        assert _is_bam_file(Path("sample.bam")) is True
        assert _is_bam_file(Path("sample.BAM")) is True
        assert _is_bam_file(Path("sample.vcf.gz")) is False
        assert _is_bam_file(Path("sample.vcf")) is False

    def test_cli_help_mentions_bam(self) -> None:
        """CLI help text mentions BAM support."""
        from click.testing import CliRunner

        from yallhap.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["classify", "--help"])

        # After implementation, help should mention BAM
        assert "BAM" in result.output or "bam" in result.output
