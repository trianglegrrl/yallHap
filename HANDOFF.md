# YClade: Modern Y-Chromosome Haplogroup Caller

## Handoff Document for Claude Code Agent

**Project Name:** yclade  
**Purpose:** A modern, pipeline-friendly Y-chromosome haplogroup inference tool  
**Target Users:** Population geneticists, forensic geneticists, paleogenomicists  
**License:** MIT (or GPL-3.0 for compatibility with dependent tools)

---

## 1. Project Vision

Build a state-of-the-art Y-chromosome haplogroup caller that:

1. **Supports modern and ancient DNA** with built-in damage modeling
2. **Uses the YFull tree** (185,780+ SNPs, monthly updates) as primary reference
3. **Accepts VCF input** (avoids variant calling complexity)
4. **Provides probabilistic confidence scores** (not just SNP counting)
5. **Supports T2T-CHM13v2.0** reference (62Mb vs 28Mb in GRCh38)
6. **Is pipeline-friendly** (Nextflow/Snakemake compatible, Bioconda/Docker)

---

## 2. Why Not Existing Tools?

| Tool | Limitation |
|------|------------|
| **Yleaf** | Not pipeline-friendly, interactive prompts, no probabilistic scoring |
| **yhaplo** (23andMe) | ISOGG tree (deprecated since 2020), genotyping array focused |
| **pathPhynder** | R-based, complex setup, designed for placement not classification |
| **Y-mer** | K-mer only, no tree-based resolution |

**Our gap:** No existing tool combines YFull tree + probabilistic scoring + ancient DNA support + T2T reference + pipeline compatibility.

---

## 3. Technical Architecture

### Language Choice: Python with Rust Extensions

**Primary:** Python 3.10+  
**Hot paths:** Rust via PyO3 (optional optimization later)  
**Rationale:** 
- Python for rapid development and bioinformatics ecosystem (pysam, pandas)
- Rust extensions can be added later for performance-critical tree traversal
- Bioconda/Docker compatibility is straightforward

### Core Algorithm

Implement a **likelihood-based approach** inspired by pathPhynder:

```
For each tree branch b:
    L(b) = P(observed variants | sample descended from b)
    
    For each informative SNP at position p:
        If derived allele observed:
            L(b) *= P(derived | descended from b)
        If ancestral allele observed:
            L(b) *= P(ancestral | descended from b)
        If missing:
            L(b) *= 1  (no information)
```

**Quality Scores (from Yleaf, improved):**

- **QC1:** Backbone consistency (intermediate markers match expected states)
- **QC2:** Terminal marker consistency (defining markers for final haplogroup)
- **QC3:** Within-haplogroup consistency (path from major haplogroup to terminal)
- **QC4 (new):** Posterior probability from likelihood calculation

### Data Flow

```
VCF Input
    ↓
Parse Y-chromosome variants (filter to chrY/Y)
    ↓
Map positions to reference (GRCh37/38/T2T via liftover)
    ↓
Load YFull tree (JSON)
    ↓
For each variant:
    - Look up in SNP database
    - Determine ancestral/derived state
    - Apply damage filter (if ancient DNA mode)
    ↓
Tree traversal with likelihood calculation
    ↓
Output: Haplogroup + confidence + path metrics
```

---

## 4. Input/Output Specification

### Input

**Required:**
- VCF file (`.vcf` or `.vcf.gz`, indexed)
- Reference genome version (`--reference grch37|grch38|t2t`)

**Optional:**
- `--ancient` flag for damage-aware mode
- `--min-depth N` minimum read depth (default: 1 for ancient, 10 for modern)
- `--min-qual N` minimum base quality (default: 20)
- `--threads N` parallel processing

### Output (JSON)

```json
{
  "sample": "SAMPLE_ID",
  "haplogroup": "R-L21",
  "confidence": 0.97,
  "method": "likelihood",
  "reference": "GRCh38",
  "tree_version": "YFull v13.06",
  "snp_stats": {
    "informative_tested": 1247,
    "derived": 145,
    "ancestral": 1089,
    "missing": 13,
    "filtered_damage": 0
  },
  "quality_scores": {
    "qc1_backbone": 0.98,
    "qc2_terminal": 1.0,
    "qc3_path": 0.95,
    "qc4_posterior": 0.97
  },
  "path": ["Y-Adam", "A0-T", "A1", "A1b", "BT", "CT", "CF", "F", "GHIJK", "HIJK", "IJK", "K", "K2", "K2b", "P", "R", "R1", "R1b", "R1b1a1a2", "R1b1a1a2a", "R1b1a1a2a1", "R1b1a1a2a1a", "R1b1a1a2a1a2", "R1b1a1a2a1a2c", "R-L21"],
  "defining_snps": ["L21/M529/S145"],
  "alternative_calls": [
    {"haplogroup": "R-DF13", "posterior": 0.02},
    {"haplogroup": "R-L21*", "posterior": 0.01}
  ]
}
```

### TSV Output (for batch processing)

```
sample	haplogroup	confidence	qc1	qc2	qc3	qc4	derived	ancestral	missing
SAMPLE1	R-L21	0.97	0.98	1.0	0.95	0.97	145	1089	13
```

---

## 5. Data Sources

### YFull Tree (PRIMARY)

**Source:** https://github.com/YFullTeam/YTree  
**File:** `current_tree.json`  
**Update frequency:** Monthly  
**Format:** JSON with nested structure

```json
{
  "A00": {
    "snps": ["AF6/V148", "V149", ...],
    "formed": 254300,
    "tmrca": 226400,
    "children": ["A00a", "A00b", "A00c"]
  },
  ...
}
```

### YBrowse SNP Database

**Source:** http://ybrowse.org/gbrowse2/gff/  
**Files:**
- `snps_hg38.vcf.gz` - VCF format positions
- `snps_hg38.csv` - CSV with all metadata
- Cross-reference between naming systems (M-series, L-series, Y-series, CTS, rs#)

### Reference Genomes

| Reference | Source | Y-chr size |
|-----------|--------|------------|
| GRCh37 | UCSC/Ensembl | ~28 Mb |
| GRCh38 | UCSC/Ensembl | ~28 Mb |
| T2T-CHM13v2.0 | `s3://human-pangenomics/T2T/CHM13/assemblies/` | 62.4 Mb |

### Validation Data

**1000 Genomes Phase 3:** 1,233 males with published haplogroups (Poznik et al. 2016)  
**Download:** ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

---

## 6. Implementation Roadmap

### Phase 1: Core Infrastructure (Week 1-2)

- [ ] Project scaffolding (pyproject.toml, CI/CD)
- [ ] YFull tree parser and data structures
- [ ] SNP database loader (YBrowse CSV)
- [ ] VCF parser (pysam-based)
- [ ] Unit tests for all parsers

### Phase 2: Classification Algorithm (Week 3-4)

- [ ] Simple tree traversal (Yleaf-style) as baseline
- [ ] QC1/QC2/QC3 score calculation
- [ ] Likelihood-based scoring (pathPhynder-style)
- [ ] Posterior probability calculation
- [ ] Integration tests with 1000 Genomes samples

### Phase 3: Ancient DNA Support (Week 5-6)

- [ ] Damage pattern detection
- [ ] Transversion-only mode
- [ ] Damage-adjusted quality scores
- [ ] Validation with published ancient samples

### Phase 4: Pipeline Integration (Week 7-8)

- [ ] CLI with proper exit codes
- [ ] Bioconda recipe
- [ ] Dockerfile
- [ ] Nextflow module example
- [ ] Snakemake rule example

### Phase 5: T2T Support (Week 9-10)

- [ ] T2T coordinate positions database
- [ ] Liftover chain integration
- [ ] Validation on T2T-aligned samples

---

## 7. Development Guidelines

### TDD Workflow

1. Write test first (in `tests/`)
2. Run test, confirm it fails
3. Implement minimum code to pass
4. Refactor if needed
5. Commit with descriptive message

### Code Style

- **Formatter:** Black (line length 100)
- **Linter:** Ruff
- **Type hints:** Required on all public functions
- **Docstrings:** Google style

### Error Handling

```python
# CORRECT: Let exceptions bubble up
def load_tree(path: Path) -> Tree:
    with open(path) as f:
        data = json.load(f)
    return Tree.from_dict(data)

# WRONG: Don't swallow errors
def load_tree(path: Path) -> Tree | None:
    try:
        with open(path) as f:
            return Tree.from_dict(json.load(f))
    except Exception:
        return None  # NO! This hides the problem
```

### Commit Messages

```
feat(tree): Add YFull tree parser with node lookup

- Implement Tree and Node dataclasses
- Add from_json() class method
- Support depth-first traversal
- 100% test coverage

Closes #12
```

---

## 8. Key Data Structures

### Tree Node

```python
@dataclass
class Node:
    name: str                    # e.g., "R-L21"
    snps: list[str]              # Defining SNPs
    parent: str | None           # Parent node name
    children: list[str]          # Child node names
    depth: int                   # Distance from root
    formed: int | None           # Years before present (estimated)
    tmrca: int | None            # Time to most recent common ancestor
```

### SNP Record

```python
@dataclass
class SNP:
    name: str                    # Primary name (e.g., "L21")
    aliases: list[str]           # Other names (M529, S145)
    position_grch37: int | None
    position_grch38: int | None
    position_t2t: int | None
    ancestral: str               # Ancestral allele
    derived: str                 # Derived allele
    haplogroup: str              # Associated haplogroup
```

### Call Result

```python
@dataclass
class HaplogroupCall:
    sample: str
    haplogroup: str
    confidence: float
    qc_scores: QCScores
    path: list[str]
    defining_snps: list[str]
    alternatives: list[tuple[str, float]]
    snp_stats: SNPStats
```

---

## 9. Testing Strategy

### Unit Tests

```python
# tests/test_tree.py
def test_tree_loads_yfull_json():
    tree = Tree.from_json(FIXTURES / "yfull_sample.json")
    assert tree.root.name == "Y-Adam"
    assert len(tree.nodes) > 100

def test_node_depth_calculation():
    tree = Tree.from_json(FIXTURES / "yfull_sample.json")
    r_l21 = tree.get("R-L21")
    assert r_l21.depth > 0
    assert r_l21.depth > tree.get("R").depth

def test_path_to_root():
    tree = Tree.from_json(FIXTURES / "yfull_sample.json")
    path = tree.path_to_root("R-L21")
    assert path[0] == "R-L21"
    assert path[-1] == "Y-Adam"
```

### Integration Tests

```python
# tests/test_integration.py
def test_classify_1000g_sample():
    """Test against known 1000 Genomes haplogroup."""
    vcf_path = FIXTURES / "NA12878_chrY.vcf.gz"
    result = classify(vcf_path, reference="grch38")
    # NA12878 is female, but for a male sample:
    # assert result.haplogroup == "R-L21"  # Known haplogroup
    assert result.confidence > 0.9
```

### Validation Tests

```python
# tests/test_validation.py
@pytest.mark.slow
def test_1000g_accuracy():
    """Validate against all 1000 Genomes males."""
    results = []
    for sample in THOUSAND_GENOMES_MALES:
        call = classify(sample.vcf, reference="grch38")
        results.append(call.haplogroup == sample.known_haplogroup)
    accuracy = sum(results) / len(results)
    assert accuracy > 0.94  # Yleaf benchmark
```

---

## 10. Dependencies

### Runtime

```toml
[project]
dependencies = [
    "pysam>=0.22.0",
    "pandas>=2.0.0",
    "numpy>=1.24.0",
    "click>=8.0.0",
]
```

### Development

```toml
[project.optional-dependencies]
dev = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "black>=23.0.0",
    "ruff>=0.1.0",
    "mypy>=1.0.0",
]
```

---

## 11. External Resources

### Documentation
- YFull API: https://www.yfull.com/api/
- YBrowse: http://ybrowse.org/
- T2T Consortium: https://github.com/marbl/CHM13

### Papers
- Yleaf 3.0: https://doi.org/10.1016/j.fsigen.2023.102870
- pathPhynder: https://doi.org/10.1093/bioinformatics/btab359
- Y-mer: (2025, check bioRxiv)
- T2T-Y: https://doi.org/10.1038/s41586-023-06457-y

### Test Data
- 1000 Genomes: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/
- HGDP: https://www.internationalgenome.org/data-portal/data-collection/hgdp

---

## 12. Quick Start for Agent

```bash
# Clone and setup
cd yclade
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Download YFull tree (do this first!)
python scripts/download_yfull_tree.py

# Download test data
python scripts/download_test_data.py

# Run classification
yclade classify sample.vcf.gz --reference grch38 --output result.json
```

---

## 13. File Structure

```
yclade/
├── HANDOFF.md              # This document
├── README.md               # User-facing documentation
├── RESEARCH.md             # Background research (from Overview)
├── pyproject.toml          # Project configuration
├── src/
│   └── yclade/
│       ├── __init__.py
│       ├── cli.py          # Click-based CLI
│       ├── tree.py         # YFull tree parser
│       ├── snps.py         # SNP database
│       ├── vcf.py          # VCF parsing
│       ├── classifier.py   # Main classification logic
│       ├── likelihood.py   # Probabilistic scoring
│       ├── ancient.py      # Ancient DNA handling
│       └── output.py       # Result formatting
├── tests/
│   ├── conftest.py         # Pytest fixtures
│   ├── fixtures/           # Test data
│   ├── test_tree.py
│   ├── test_snps.py
│   ├── test_vcf.py
│   ├── test_classifier.py
│   └── test_integration.py
├── scripts/
│   ├── download_yfull_tree.py
│   └── download_test_data.py
├── data/
│   └── .gitkeep            # Data downloaded at runtime
├── docker/
│   └── Dockerfile
└── conda/
    └── meta.yaml
```

---

## 14. First Tasks for Agent

1. **Implement `tree.py`** - Parse YFull JSON into Tree/Node structure
2. **Write tests for tree** - Ensure correct depth calculation, path finding
3. **Implement `snps.py`** - Load YBrowse CSV, map positions
4. **Write tests for snps** - Position lookup, alias resolution
5. **Implement basic `classifier.py`** - Simple derived/ancestral counting
6. **Integration test** - Classify one known sample correctly

Start with TDD: write the test, then implement the code.

---

## 15. Questions to Resolve

1. **Naming convention priority:** When SNP has multiple names (L21/M529/S145), which to display?
2. **Haplogroup notation:** YFull style (R-L21) vs ISOGG style (R1b1a1a2a1a2c)?
3. **Minimum confidence threshold:** Default 0.95 (like Yleaf) or configurable?
4. **Ancient DNA default settings:** Transversion-only or damage-filtered transitions?
5. **Multi-sample VCF:** Process all samples or require single-sample?

---

**End of Handoff Document**

*Generated: December 2024*  
*For questions: Review RESEARCH.md for detailed background*
