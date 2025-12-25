"""
Pytest configuration and fixtures for yclade tests.
"""

from pathlib import Path

import pytest

# Path to test fixtures
FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    """Return path to fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture
def sample_tree_dict() -> dict:
    """
    Return a minimal tree structure for testing.

    Structure:
        ROOT
        ├── A00
        │   ├── A00a
        │   └── A00b
        └── A0-T
            └── A1
                └── BT
                    └── CT
                        └── R
                            └── R1
                                └── R1b
                                    └── R-L21
    """
    return {
        'ROOT (Y-Chromosome "Adam")': ["A00", "A0-T"],
        "A00": ["A00a", "A00b"],
        "A0-T": ["A1"],
        "A1": ["BT"],
        "BT": ["CT"],
        "CT": ["R"],
        "R": ["R1"],
        "R1": ["R1b"],
        "R1b": ["R-L21"],
    }


@pytest.fixture
def sample_snps_csv(tmp_path: Path) -> Path:
    """Create a sample SNP database CSV for testing."""
    csv_content = """name,aliases,grch37_pos,grch38_pos,ancestral,derived,haplogroup
L21,M529;S145,2887478,2655229,C,T,R-L21
M269,,2803636,2571333,G,A,R1b
M343,,2888436,2656183,G,A,R1
M207,,14524515,14279666,A,G,R
M168,,14941780,14696931,C,T,CT
"""
    csv_path = tmp_path / "snps.csv"
    csv_path.write_text(csv_content)
    return csv_path


@pytest.fixture
def sample_vcf(tmp_path: Path) -> Path:
    """Create a sample VCF for testing."""
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chrY,length=59373566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chrY	2655229	.	C	T	100	PASS	.	GT:DP:GQ	1:30:99
chrY	2571333	.	G	A	100	PASS	.	GT:DP:GQ	1:25:99
chrY	2656183	.	G	A	100	PASS	.	GT:DP:GQ	1:28:99
"""
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text(vcf_content)
    return vcf_path
