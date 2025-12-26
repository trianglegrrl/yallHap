"""
Unit tests for yallhap.classifier module.
"""

import pytest

from yallhap.classifier import (
    HaplogroupCall,
    HaplogroupClassifier,
    QCScores,
    SNPStats,
)
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


class TestQCScores:
    """Tests for QCScores dataclass."""

    def test_qc_scores_combined(self) -> None:
        """Test combined score calculation."""
        scores = QCScores(
            qc1_backbone=1.0,
            qc2_terminal=1.0,
            qc3_path=1.0,
            qc4_posterior=1.0,
        )
        assert scores.combined == 1.0

    def test_qc_scores_combined_partial(self) -> None:
        """Test combined score with partial values."""
        scores = QCScores(
            qc1_backbone=0.8,
            qc2_terminal=0.9,
            qc3_path=0.85,
            qc4_posterior=0.95,
        )
        # Geometric mean
        expected = (0.8 * 0.9 * 0.85 * 0.95) ** 0.25
        assert abs(scores.combined - expected) < 0.001

    def test_qc_scores_combined_zero(self) -> None:
        """Test combined score with zero value."""
        scores = QCScores(
            qc1_backbone=0.0,
            qc2_terminal=1.0,
            qc3_path=1.0,
            qc4_posterior=1.0,
        )
        assert scores.combined == 0.0


class TestSNPStats:
    """Tests for SNPStats dataclass."""

    def test_snp_stats_total_called(self) -> None:
        """Test total_called property."""
        stats = SNPStats(
            informative_tested=100,
            derived=50,
            ancestral=40,
            missing=10,
        )
        assert stats.total_called == 90


class TestHaplogroupCall:
    """Tests for HaplogroupCall dataclass."""

    def test_to_dict(self) -> None:
        """Test JSON serialization."""
        call = HaplogroupCall(
            sample="TEST",
            haplogroup="R-L21",
            confidence=0.95,
            qc_scores=QCScores(0.9, 1.0, 0.95, 0.95),
            path=["ROOT", "R", "R-L21"],
            defining_snps=["L21"],
            alternatives=[("R1b", 0.8)],
            snp_stats=SNPStats(100, 50, 40, 10),
            reference="grch38",
            tree_version="YFull v13",
        )

        d = call.to_dict()
        assert d["sample"] == "TEST"
        assert d["haplogroup"] == "R-L21"
        assert d["confidence"] == 0.95
        assert "quality_scores" in d
        assert d["quality_scores"]["qc1_backbone"] == 0.9


class TestHaplogroupClassifier:
    """Tests for HaplogroupClassifier class."""

    @pytest.fixture
    def classifier(self, sample_tree_dict: dict, sample_snps_csv) -> HaplogroupClassifier:
        """Create a classifier with test data."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        return HaplogroupClassifier(
            tree=tree,
            snp_db=snp_db,
            reference="grch38",
        )

    def test_classifier_creation(self, classifier: HaplogroupClassifier) -> None:
        """Test classifier initialization."""
        assert classifier.tree is not None
        assert classifier.snp_db is not None
        assert classifier.reference == "grch38"

    def test_classifier_ancient_mode(self, sample_tree_dict: dict, sample_snps_csv) -> None:
        """Test ancient DNA mode configuration."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)

        classifier = HaplogroupClassifier(
            tree=tree,
            snp_db=snp_db,
            ancient_mode=True,
            min_depth=1,
        )

        assert classifier.ancient_mode is True
        assert classifier.min_depth == 1


class TestClassifierDamageFiltering:
    """Tests for ancient DNA damage filtering."""

    @pytest.fixture
    def classifier_ancient(self, sample_tree_dict: dict, sample_snps_csv) -> HaplogroupClassifier:
        """Create classifier in ancient mode."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        return HaplogroupClassifier(
            tree=tree,
            snp_db=snp_db,
            ancient_mode=True,
        )

    def test_damage_like_c_to_t(self, classifier_ancient: HaplogroupClassifier) -> None:
        """Test C>T is flagged as damage-like."""
        from yallhap.snps import SNP

        snp = SNP(name="TEST", ancestral="C", derived="T")
        assert classifier_ancient._is_damage_like(snp, "T")

    def test_damage_like_g_to_a(self, classifier_ancient: HaplogroupClassifier) -> None:
        """Test G>A is flagged as damage-like."""
        from yallhap.snps import SNP

        snp = SNP(name="TEST", ancestral="G", derived="A")
        assert classifier_ancient._is_damage_like(snp, "A")

    def test_not_damage_like_transversion(self, classifier_ancient: HaplogroupClassifier) -> None:
        """Test transversions are not flagged as damage."""
        from yallhap.snps import SNP

        snp = SNP(name="TEST", ancestral="C", derived="A")
        assert not classifier_ancient._is_damage_like(snp, "A")


class TestClassifyBatch:
    """Tests for batch classification functionality."""

    @pytest.fixture
    def multi_sample_vcf(self, tmp_path) -> str:
        """Create a minimal multi-sample VCF for testing."""
        import pysam

        # Write uncompressed VCF first
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=Y,length=59373566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
Y	2887824	M343	C	A	.	PASS	.	GT	1/1	0/0	1/1
Y	15654428	L21	C	G	.	PASS	.	GT	1/1	0/0	0/0
Y	22739367	M269	T	C	.	PASS	.	GT	1/1	0/0	./.
"""
        vcf_path = tmp_path / "test.vcf"
        with open(vcf_path, "w") as f:
            f.write(vcf_content)

        # Compress with bgzip and index
        vcf_gz = tmp_path / "test.vcf.gz"
        pysam.tabix_compress(str(vcf_path), str(vcf_gz), force=True)
        pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)

        return str(vcf_gz)

    @pytest.fixture
    def batch_classifier(self, sample_tree_dict: dict, sample_snps_csv) -> HaplogroupClassifier:
        """Create classifier for batch tests."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        return HaplogroupClassifier(tree=tree, snp_db=snp_db, reference="grch38")

    def test_classify_batch_returns_list(
        self, batch_classifier: HaplogroupClassifier, multi_sample_vcf: str
    ) -> None:
        """Test batch classification returns list of results."""
        samples = ["SAMPLE1", "SAMPLE2"]
        results = batch_classifier.classify_batch(multi_sample_vcf, samples)

        assert isinstance(results, list)
        assert len(results) == 2

    def test_classify_batch_result_samples_match(
        self, batch_classifier: HaplogroupClassifier, multi_sample_vcf: str
    ) -> None:
        """Test batch results have correct sample names."""
        samples = ["SAMPLE1", "SAMPLE2", "SAMPLE3"]
        results = batch_classifier.classify_batch(multi_sample_vcf, samples)

        assert results[0].sample == "SAMPLE1"
        assert results[1].sample == "SAMPLE2"
        assert results[2].sample == "SAMPLE3"

    def test_classify_batch_invalid_sample_raises(
        self, batch_classifier: HaplogroupClassifier, multi_sample_vcf: str
    ) -> None:
        """Test batch with invalid sample name raises ValueError."""
        with pytest.raises(ValueError, match="Sample INVALID not found"):
            batch_classifier.classify_batch(multi_sample_vcf, ["INVALID"])

    def test_classify_batch_empty_samples(
        self, batch_classifier: HaplogroupClassifier, multi_sample_vcf: str
    ) -> None:
        """Test batch with empty sample list returns empty list."""
        results = batch_classifier.classify_batch(multi_sample_vcf, [])
        assert results == []

    def test_classify_batch_with_ancient_mode(
        self, sample_tree_dict: dict, sample_snps_csv, multi_sample_vcf: str
    ) -> None:
        """Test batch classification with transversions-only mode."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        classifier = HaplogroupClassifier(
            tree=tree, snp_db=snp_db, reference="grch38", transversions_only=True
        )

        samples = ["SAMPLE1", "SAMPLE2"]
        results = classifier.classify_batch(multi_sample_vcf, samples)

        assert len(results) == 2
        # Results should have fewer derived SNPs due to transition filtering
        for r in results:
            assert isinstance(r, HaplogroupCall)


# Integration tests would go here but require real VCF data
class TestClassifierIntegration:
    """Integration tests for classifier (marked slow)."""

    @pytest.mark.skip(reason="Requires real test data")
    def test_classify_real_sample(self) -> None:
        """Test classification on real 1000 Genomes sample."""
        pass
