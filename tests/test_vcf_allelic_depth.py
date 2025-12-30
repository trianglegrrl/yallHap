"""
Tests for allelic depth (AD) support in VCF parsing.

TDD tests for:
- Parsing AD field from VCF
- Calculating support_ratio for variant calls
- Handling missing/ambiguous AD data
"""

from __future__ import annotations

from pathlib import Path

import pytest

from yallhap.vcf import Variant, VCFReader


class TestVariantAllelicDepth:
    """Tests for Variant allelic depth handling."""

    def test_variant_with_ad_field(self) -> None:
        """Variant correctly stores AD field."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            depth=30,
            quality=99,
            allele_depth=(5, 25),  # 5 ref reads, 25 alt reads
        )

        assert variant.allele_depth == (5, 25)
        assert variant.allele_depth[0] == 5  # ref reads
        assert variant.allele_depth[1] == 25  # alt reads

    def test_variant_without_ad_field(self) -> None:
        """Variant without AD has None allele_depth."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            depth=30,
            quality=99,
        )

        assert variant.allele_depth is None


class TestVariantSupportRatio:
    """Tests for support_ratio property."""

    def test_support_ratio_homozygous_alt(self) -> None:
        """Support ratio is 1.0 for clean hom-alt call with all alt reads."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,  # Alt called
            allele_depth=(0, 30),  # 0 ref, 30 alt
        )

        assert variant.support_ratio == 1.0

    def test_support_ratio_homozygous_ref(self) -> None:
        """Support ratio is 1.0 for clean hom-ref call with all ref reads."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=0,  # Ref called
            allele_depth=(30, 0),  # 30 ref, 0 alt
        )

        assert variant.support_ratio == 1.0

    def test_support_ratio_mixed_reads_alt_called(self) -> None:
        """Support ratio reflects actual read distribution when alt called."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,  # Alt called
            allele_depth=(10, 20),  # 10 ref, 20 alt
        )

        # Alt is called, so support ratio = alt_reads / total = 20/30
        expected = 20 / 30
        assert abs(variant.support_ratio - expected) < 0.001

    def test_support_ratio_mixed_reads_ref_called(self) -> None:
        """Support ratio reflects actual read distribution when ref called."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=0,  # Ref called
            allele_depth=(20, 10),  # 20 ref, 10 alt
        )

        # Ref is called, so support ratio = ref_reads / total = 20/30
        expected = 20 / 30
        assert abs(variant.support_ratio - expected) < 0.001

    def test_support_ratio_missing_ad(self) -> None:
        """Missing AD returns None for support_ratio."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            allele_depth=None,
        )

        assert variant.support_ratio is None

    def test_support_ratio_zero_depth(self) -> None:
        """Zero total depth returns None for support_ratio."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            allele_depth=(0, 0),
        )

        assert variant.support_ratio is None

    def test_support_ratio_missing_genotype(self) -> None:
        """Missing genotype returns None for support_ratio."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=None,
            allele_depth=(10, 20),
        )

        assert variant.support_ratio is None

    def test_support_ratio_multiallelic(self) -> None:
        """Support ratio works with multiallelic sites."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T", "G"),  # Two alts
            genotype=2,  # Second alt (G) called
            allele_depth=(5, 10, 15),  # 5 ref, 10 T, 15 G
        )

        # G is called (genotype=2), so support ratio = 15 / 30
        expected = 15 / 30
        assert abs(variant.support_ratio - expected) < 0.001


class TestVariantIsAmbiguous:
    """Tests for is_ambiguous property based on support ratio."""

    def test_is_ambiguous_high_support(self) -> None:
        """High support ratio (>=0.7) is not ambiguous."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            allele_depth=(3, 27),  # 90% support
        )

        assert not variant.is_ambiguous

    def test_is_ambiguous_low_support(self) -> None:
        """Low support ratio (<0.7) is ambiguous."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            allele_depth=(15, 15),  # 50% support
        )

        assert variant.is_ambiguous

    def test_is_ambiguous_borderline(self) -> None:
        """Exactly 70% support is not ambiguous."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            allele_depth=(30, 70),  # Exactly 70% support
        )

        assert not variant.is_ambiguous

    def test_is_ambiguous_missing_ad(self) -> None:
        """Missing AD means not ambiguous (unknown, not low support)."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            allele_depth=None,
        )

        # When we don't have AD data, we can't determine ambiguity
        assert not variant.is_ambiguous


class TestVCFReaderAD:
    """Tests for VCFReader AD field parsing."""

    @pytest.fixture
    def vcf_with_ad(self, tmp_path: Path) -> Path:
        """Create a VCF file with AD field for testing."""
        import pysam

        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=Y,length=59373566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
Y\t1000\t.\tC\tT\t100\tPASS\t.\tGT:DP:GQ:AD\t1:30:99:5,25
Y\t2000\t.\tG\tA\t100\tPASS\t.\tGT:DP:GQ:AD\t1:20:80:0,20
Y\t3000\t.\tA\tC\t100\tPASS\t.\tGT:DP:GQ:AD\t0:25:99:25,0
"""
        vcf_path = tmp_path / "test_ad.vcf"
        with open(vcf_path, "w") as f:
            f.write(vcf_content)

        # Compress and index
        vcf_gz = tmp_path / "test_ad.vcf.gz"
        pysam.tabix_compress(str(vcf_path), str(vcf_gz), force=True)
        pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)

        return vcf_gz

    def test_reader_parses_ad_field(self, vcf_with_ad: Path) -> None:
        """VCFReader correctly parses AD field from VCF."""
        with VCFReader(vcf_with_ad) as reader:
            variants = list(reader.iter_variants())

        assert len(variants) == 3

        # First variant: 5 ref, 25 alt
        assert variants[0].allele_depth == (5, 25)
        assert variants[0].position == 1000

        # Second variant: 0 ref, 20 alt
        assert variants[1].allele_depth == (0, 20)

        # Third variant: 25 ref, 0 alt
        assert variants[2].allele_depth == (25, 0)

    def test_reader_handles_missing_ad(self, tmp_path: Path) -> None:
        """VCFReader handles VCF without AD field."""
        import pysam

        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=Y,length=59373566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
Y\t1000\t.\tC\tT\t100\tPASS\t.\tGT\t1
"""
        vcf_path = tmp_path / "no_ad.vcf"
        with open(vcf_path, "w") as f:
            f.write(vcf_content)

        vcf_gz = tmp_path / "no_ad.vcf.gz"
        pysam.tabix_compress(str(vcf_path), str(vcf_gz), force=True)
        pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)

        with VCFReader(vcf_gz) as reader:
            variants = list(reader.iter_variants())

        assert len(variants) == 1
        assert variants[0].allele_depth is None


class TestTotalAllelicDepth:
    """Tests for total_allelic_depth property."""

    def test_total_allelic_depth_present(self) -> None:
        """Total allelic depth is sum of all allele depths."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            allele_depth=(10, 20),
        )

        assert variant.total_allelic_depth == 30

    def test_total_allelic_depth_multiallelic(self) -> None:
        """Total allelic depth works for multiallelic sites."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T", "G"),
            genotype=1,
            allele_depth=(5, 10, 15),
        )

        assert variant.total_allelic_depth == 30

    def test_total_allelic_depth_missing(self) -> None:
        """Total allelic depth is None when AD missing."""
        variant = Variant(
            chrom="Y",
            position=1000,
            ref="C",
            alt=("T",),
            genotype=1,
            allele_depth=None,
        )

        assert variant.total_allelic_depth is None
