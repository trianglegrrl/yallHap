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
            tree_version="YFull (185780 SNPs, hash: a1b2c3d4)",
        )

        d = call.to_dict()
        assert d["sample"] == "TEST"
        assert d["haplogroup"] == "R-L21"
        assert d["confidence"] == 0.95
        assert "quality_scores" in d
        assert d["quality_scores"]["qc1_backbone"] == 0.9
        assert "YFull" in d["tree_version"]

    def test_tree_version_contains_metadata(self) -> None:
        """Test that tree_version contains meaningful metadata."""
        call = HaplogroupCall(
            sample="TEST",
            haplogroup="R-L21",
            confidence=0.95,
            qc_scores=QCScores(0.9, 1.0, 0.95, 0.95),
            path=["ROOT", "R", "R-L21"],
            defining_snps=["L21"],
            alternatives=[],
            snp_stats=SNPStats(100, 50, 40, 10),
            reference="grch38",
            tree_version="YFull (185780 SNPs, hash: a1b2c3d4)",
        )

        d = call.to_dict()
        # tree_version should contain source, SNP count, and hash
        assert "YFull" in d["tree_version"]
        assert "SNPs" in d["tree_version"]
        assert "hash:" in d["tree_version"]


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


class TestParallelBatchProcessing:
    """Tests for parallel batch processing in CLI."""

    @pytest.fixture
    def vcf_files(self, tmp_path) -> list:
        """Create multiple VCF files for parallel processing tests."""
        import pysam

        vcf_paths = []

        # VCF 1 - Sample with R-L21 markers (positions must be sorted!)
        vcf_content_1 = """##fileformat=VCFv4.2
##contig=<ID=Y,length=59373566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE_R1
Y	2571333	M269	G	A	.	PASS	.	GT	1/1
Y	2655229	L21	C	T	.	PASS	.	GT	1/1
Y	2656183	M343	G	A	.	PASS	.	GT	1/1
"""
        vcf1_path = tmp_path / "sample1.vcf"
        with open(vcf1_path, "w") as f:
            f.write(vcf_content_1)
        vcf1_gz = tmp_path / "sample1.vcf.gz"
        pysam.tabix_compress(str(vcf1_path), str(vcf1_gz), force=True)
        pysam.tabix_index(str(vcf1_gz), preset="vcf", force=True)
        vcf_paths.append(vcf1_gz)

        # VCF 2 - Different sample
        vcf_content_2 = """##fileformat=VCFv4.2
##contig=<ID=Y,length=59373566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE_R2
Y	2571333	M269	G	A	.	PASS	.	GT	1/1
Y	2655229	L21	C	T	.	PASS	.	GT	0/0
Y	2656183	M343	G	A	.	PASS	.	GT	1/1
"""
        vcf2_path = tmp_path / "sample2.vcf"
        with open(vcf2_path, "w") as f:
            f.write(vcf_content_2)
        vcf2_gz = tmp_path / "sample2.vcf.gz"
        pysam.tabix_compress(str(vcf2_path), str(vcf2_gz), force=True)
        pysam.tabix_index(str(vcf2_gz), preset="vcf", force=True)
        vcf_paths.append(vcf2_gz)

        # VCF 3 - Third sample
        vcf_content_3 = """##fileformat=VCFv4.2
##contig=<ID=Y,length=59373566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE_R3
Y	2571333	M269	G	A	.	PASS	.	GT	0/0
Y	2655229	L21	C	T	.	PASS	.	GT	1/1
Y	2656183	M343	G	A	.	PASS	.	GT	1/1
"""
        vcf3_path = tmp_path / "sample3.vcf"
        with open(vcf3_path, "w") as f:
            f.write(vcf_content_3)
        vcf3_gz = tmp_path / "sample3.vcf.gz"
        pysam.tabix_compress(str(vcf3_path), str(vcf3_gz), force=True)
        pysam.tabix_index(str(vcf3_gz), preset="vcf", force=True)
        vcf_paths.append(vcf3_gz)

        return vcf_paths

    @pytest.fixture
    def tree_json(self, tmp_path, sample_tree_dict: dict) -> str:
        """Create tree JSON file for CLI tests."""
        import json

        tree_path = tmp_path / "tree.json"
        with open(tree_path, "w") as f:
            json.dump(sample_tree_dict, f)
        return str(tree_path)

    def test_worker_init_and_classify(
        self, vcf_files: list, tree_json: str, sample_snps_csv
    ) -> None:
        """Test worker initialization and classification functions."""
        from pathlib import Path

        from yallhap.cli import _classify_file, _init_worker

        # Initialize worker
        _init_worker(
            Path(tree_json),
            sample_snps_csv,
            {"reference": "grch38"},
        )

        # Classify a file
        result = _classify_file(vcf_files[0])

        assert result is not None
        assert result.sample == "SAMPLE_R1"
        assert isinstance(result, HaplogroupCall)

    def test_parallel_results_match_sequential(
        self, vcf_files: list, sample_tree_dict: dict, sample_snps_csv
    ) -> None:
        """Test that parallel processing returns same results as sequential."""
        from pathlib import Path

        # Run sequential classification
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        classifier = HaplogroupClassifier(tree=tree, snp_db=snp_db, reference="grch38")

        sequential_results = []
        for vcf in vcf_files:
            result = classifier.classify(vcf)
            sequential_results.append(result)

        # Run parallel classification using ProcessPoolExecutor
        from concurrent.futures import ProcessPoolExecutor

        from yallhap.cli import _classify_file, _init_worker

        with ProcessPoolExecutor(
            max_workers=2,
            initializer=_init_worker,
            initargs=(
                Path(sample_snps_csv).parent / "tree.json",
                sample_snps_csv,
                {"reference": "grch38"},
            ),
        ) as _executor:
            # Need to create tree.json first
            import json

            tree_path = Path(sample_snps_csv).parent / "tree.json"
            with open(tree_path, "w") as f:
                json.dump(sample_tree_dict, f)

            # Now run with actual init
            with ProcessPoolExecutor(
                max_workers=2,
                initializer=_init_worker,
                initargs=(tree_path, sample_snps_csv, {"reference": "grch38"}),
            ) as executor2:
                parallel_results = list(executor2.map(_classify_file, vcf_files))

        # Compare results - haplogroups should match
        sequential_hgs = sorted([r.haplogroup for r in sequential_results])
        parallel_hgs = sorted([r.haplogroup for r in parallel_results])
        assert sequential_hgs == parallel_hgs

        # Sample names should match
        sequential_samples = sorted([r.sample for r in sequential_results])
        parallel_samples = sorted([r.sample for r in parallel_results])
        assert sequential_samples == parallel_samples

    def test_parallel_single_file_with_multiple_threads(
        self, vcf_files: list, tree_json: str, sample_snps_csv
    ) -> None:
        """Test parallel processing with more threads than files."""
        from concurrent.futures import ProcessPoolExecutor
        from pathlib import Path

        from yallhap.cli import _classify_file, _init_worker

        # Use only 1 file but request 4 threads
        single_file = [vcf_files[0]]

        with ProcessPoolExecutor(
            max_workers=4,  # More workers than files
            initializer=_init_worker,
            initargs=(Path(tree_json), sample_snps_csv, {"reference": "grch38"}),
        ) as executor:
            results = list(executor.map(_classify_file, single_file))

        assert len(results) == 1
        assert results[0].sample == "SAMPLE_R1"

    @pytest.fixture
    def multi_sample_vcf(self, tmp_path) -> str:
        """Create a multi-sample VCF for batch parallel tests."""
        import pysam

        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=Y,length=59373566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE_A	SAMPLE_B	SAMPLE_C	SAMPLE_D
Y	2571333	M269	G	A	.	PASS	.	GT	1/1	1/1	0/0	1/1
Y	2655229	L21	C	T	.	PASS	.	GT	1/1	0/0	1/1	0/0
Y	2656183	M343	G	A	.	PASS	.	GT	1/1	1/1	1/1	0/0
"""
        vcf_path = tmp_path / "multi_sample.vcf"
        with open(vcf_path, "w") as f:
            f.write(vcf_content)
        vcf_gz = tmp_path / "multi_sample.vcf.gz"
        pysam.tabix_compress(str(vcf_path), str(vcf_gz), force=True)
        pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)
        return str(vcf_gz)

    def test_classify_batch_parallel_matches_sequential(
        self, multi_sample_vcf: str, sample_tree_dict: dict, sample_snps_csv
    ) -> None:
        """Test that parallel batch classification matches sequential results."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        classifier = HaplogroupClassifier(tree=tree, snp_db=snp_db, reference="grch38")

        samples = ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C", "SAMPLE_D"]

        # Run sequential (threads=1)
        sequential_results = classifier.classify_batch(multi_sample_vcf, samples, threads=1)

        # Run parallel (threads=4)
        parallel_results = classifier.classify_batch(multi_sample_vcf, samples, threads=4)

        # Results should match exactly
        assert len(sequential_results) == len(parallel_results) == 4

        for seq, par in zip(sequential_results, parallel_results, strict=True):
            assert seq.sample == par.sample
            assert seq.haplogroup == par.haplogroup
            assert seq.confidence == par.confidence
            # SNP stats should also match
            if seq.snp_stats and par.snp_stats:
                assert seq.snp_stats.derived == par.snp_stats.derived
                assert seq.snp_stats.ancestral == par.snp_stats.ancestral

    def test_classify_batch_parallel_preserves_order(
        self, multi_sample_vcf: str, sample_tree_dict: dict, sample_snps_csv
    ) -> None:
        """Test that parallel batch classification preserves sample order."""
        tree = Tree.from_dict(sample_tree_dict)
        snp_db = SNPDatabase.from_csv(sample_snps_csv)
        classifier = HaplogroupClassifier(tree=tree, snp_db=snp_db, reference="grch38")

        samples = ["SAMPLE_D", "SAMPLE_A", "SAMPLE_C", "SAMPLE_B"]  # Non-alphabetical order

        results = classifier.classify_batch(multi_sample_vcf, samples, threads=4)

        # Results should be in the same order as input samples
        assert [r.sample for r in results] == samples


# Integration tests would go here but require real VCF data
class TestClassifierIntegration:
    """Integration tests for classifier (marked slow)."""

    @pytest.mark.skip(reason="Requires real test data")
    def test_classify_real_sample(self) -> None:
        """Test classification on real 1000 Genomes sample."""
        pass
