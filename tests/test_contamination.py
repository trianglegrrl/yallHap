"""
Tests for Y-chromosome contamination estimation.

Following TDD - these tests are written first, before implementation.
"""

from __future__ import annotations

from yallhap.vcf import Variant


class TestContaminationEstimation:
    """Tests for Y-chromosome contamination estimation."""

    def test_contamination_from_mixed_reads(self) -> None:
        """Mixed allele depths indicate contamination."""
        from yallhap.contamination import estimate_contamination

        # Variant with 8 ref reads, 2 alt reads at informative site
        # Called genotype is ref (0), but we see ~20% alt reads = contamination
        variants = {
            1000: Variant(
                chrom="Y",
                position=1000,
                ref="A",
                alt=("G",),
                genotype=0,
                depth=10,
                quality=30,
                allele_depth=(8, 2),
            )
        }
        rate, n_sites = estimate_contamination(variants, min_depth=1)
        assert 0.15 <= rate <= 0.25
        assert n_sites >= 1

    def test_no_contamination_clean_sample(self) -> None:
        """Clean sample with only major allele shows 0% contamination."""
        from yallhap.contamination import estimate_contamination

        variants = {
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
                depth=15,
                quality=30,
                allele_depth=(0, 15),
            ),
        }
        rate, n_sites = estimate_contamination(variants, min_depth=1)
        assert rate < 0.05
        assert n_sites >= 2

    def test_transversions_only_mode(self) -> None:
        """Transversions-only excludes C>T/G>A sites."""
        from yallhap.contamination import estimate_contamination

        variants = {
            # C>T transition - should be excluded in transversions_only mode
            1000: Variant(
                chrom="Y",
                position=1000,
                ref="C",
                alt=("T",),
                genotype=1,
                depth=10,
                quality=30,
                allele_depth=(5, 5),  # High contamination
            ),
            # A>C transversion - should be included
            2000: Variant(
                chrom="Y",
                position=2000,
                ref="A",
                alt=("C",),
                genotype=0,
                depth=10,
                quality=30,
                allele_depth=(10, 0),  # No contamination
            ),
        }
        rate, n_sites = estimate_contamination(variants, min_depth=1, transversions_only=True)
        # Only transversion site counted, which shows no contamination
        assert rate < 0.05
        assert n_sites == 1

    def test_minimum_depth_filter(self) -> None:
        """Sites below min_depth are excluded."""
        from yallhap.contamination import estimate_contamination

        variants = {
            # Low depth - should be excluded
            1000: Variant(
                chrom="Y",
                position=1000,
                ref="A",
                alt=("G",),
                genotype=0,
                depth=2,
                quality=30,
                allele_depth=(1, 1),  # 50% contamination
            ),
            # High depth - should be included
            2000: Variant(
                chrom="Y",
                position=2000,
                ref="A",
                alt=("G",),
                genotype=0,
                depth=10,
                quality=30,
                allele_depth=(10, 0),  # No contamination
            ),
        }
        rate, n_sites = estimate_contamination(variants, min_depth=5)
        # Only high-depth site counted
        assert rate < 0.05
        assert n_sites == 1

    def test_no_allele_depth_returns_none(self) -> None:
        """Sites without AD field are excluded."""
        from yallhap.contamination import estimate_contamination

        variants = {
            1000: Variant(
                chrom="Y",
                position=1000,
                ref="A",
                alt=("G",),
                genotype=0,
                depth=10,
                quality=30,
                allele_depth=None,  # No AD field
            ),
        }
        rate, n_sites = estimate_contamination(variants, min_depth=1)
        assert n_sites == 0
        # Rate should be 0.0 when no sites available
        assert rate == 0.0

    def test_high_contamination_detected(self) -> None:
        """High contamination (>10%) is correctly estimated."""
        from yallhap.contamination import estimate_contamination

        variants = {
            1000: Variant(
                chrom="Y",
                position=1000,
                ref="A",
                alt=("G",),
                genotype=0,
                depth=20,
                quality=30,
                allele_depth=(14, 6),  # 30% minor allele
            ),
            2000: Variant(
                chrom="Y",
                position=2000,
                ref="T",
                alt=("C",),
                genotype=1,
                depth=20,
                quality=30,
                allele_depth=(5, 15),  # 25% minor allele
            ),
        }
        rate, n_sites = estimate_contamination(variants, min_depth=1)
        # Average should be around 27.5%
        assert 0.20 <= rate <= 0.35
        assert n_sites == 2

    def test_empty_variants_returns_zero(self) -> None:
        """Empty variants dict returns 0% contamination."""
        from yallhap.contamination import estimate_contamination

        rate, n_sites = estimate_contamination({}, min_depth=1)
        assert rate == 0.0
        assert n_sites == 0


class TestContaminationFilter:
    """Tests for contamination filtering threshold."""

    def test_filter_threshold_flags_contaminated(self) -> None:
        """Samples exceeding threshold are flagged."""
        from yallhap.contamination import ContaminationResult, is_contaminated

        result = ContaminationResult(rate=0.15, n_sites=50)
        assert is_contaminated(result, threshold=0.10)
        assert not is_contaminated(result, threshold=0.20)

    def test_contamination_result_dataclass(self) -> None:
        """ContaminationResult stores rate and site count."""
        from yallhap.contamination import ContaminationResult

        result = ContaminationResult(rate=0.05, n_sites=100)
        assert result.rate == 0.05
        assert result.n_sites == 100


class TestContaminationWithSNPDB:
    """Tests for contamination using SNP database filtering."""

    def test_only_informative_sites_used(self) -> None:
        """Only sites from SNP database are used for estimation."""
        from yallhap.contamination import estimate_contamination_with_snpdb
        from yallhap.snps import SNP, SNPDatabase

        # Create minimal SNP database
        db = SNPDatabase()
        # Add one SNP at position 1000
        snp = SNP(
            name="TEST1",
            position_grch38=1000,
            ancestral="A",
            derived="G",
            haplogroup="R",
        )
        db._add_snp(snp)

        variants = {
            # In database - should be used
            1000: Variant(
                chrom="Y",
                position=1000,
                ref="A",
                alt=("G",),
                genotype=0,
                depth=10,
                quality=30,
                allele_depth=(8, 2),
            ),
            # Not in database - should be excluded
            2000: Variant(
                chrom="Y",
                position=2000,
                ref="C",
                alt=("T",),
                genotype=1,
                depth=10,
                quality=30,
                allele_depth=(5, 5),
            ),
        }
        rate, n_sites = estimate_contamination_with_snpdb(
            variants, db, reference="grch38", min_depth=1
        )
        # Only position 1000 used
        assert n_sites == 1
        assert 0.15 <= rate <= 0.25
