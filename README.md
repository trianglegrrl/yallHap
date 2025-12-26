# yallHap

Modern, pipeline-friendly Y-chromosome haplogroup inference.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Features

- **YFull tree**: Uses the most comprehensive Y-chromosome phylogeny (185,780+ SNPs)
- **Probabilistic scoring**: Likelihood-based confidence scores, not just SNP counting
- **Ancient DNA support**: Built-in damage filtering, transversions-only mode, quality rescaling
- **Multiple references**: Supports GRCh37, GRCh38, and T2T-CHM13v2.0 with automatic liftover
- **Multi-threaded**: Parallel sample processing with `--threads N` for population-scale studies
- **Batch processing**: Classify thousands of samples efficiently with `classify_batch()`
- **Pipeline-friendly**: Proper exit codes, JSON/TSV output, Nextflow/Snakemake examples
- **Bioconda/Docker**: Easy installation and containerized execution

## Accuracy

Validated against established datasets:

| Dataset | Samples | Same Major Lineage | Reference | Notes |
|---------|---------|-------------------|-----------|-------|
| 1000 Genomes Phase 3 | 1,233 | **99.76%** | GRCh37 | Modern WGS |
| AADR Ancient DNA | 1,991 | **88.65%** | GRCh37 | Transversions-only mode |
| gnomAD HGDP/1KG | 10 | **100.00%** | GRCh38 | High-coverage WGS |

**1000 Genomes details:**
- Only 3 misclassified samples (2 rare A0 haplogroups, 1 NO/K confusion)
- Mean confidence: 0.994
- Mean derived SNPs: 15.4

**AADR Ancient DNA details:**
- Validated with properly formatted ground truth (X-YYYY haplogroup names)
- Transversions-only mode for maximum damage resistance
- Mean derived SNPs: 11.8 (low coverage typical for aDNA)

**gnomAD High-Coverage details:**
- 30× high-coverage whole-genome sequencing
- Mean derived SNPs: 79.4 (5× more than Phase 3 due to improved coverage)

See [VALIDATION_TESTING.md](VALIDATION_TESTING.md) for reproducible validation protocols.

## Installation

### pip (recommended)

```bash
pip install yallhap
```

### Conda

```bash
conda install -c bioconda yallhap
```

### Docker

```bash
docker pull trianglegrrl/yallhap
```

### From source

```bash
git clone https://github.com/trianglegrrl/yallhap.git
cd yallhap
pip install -e ".[dev]"
```

## Quick Start

### 1. Download reference data

```bash
yallhap download --output-dir data/
```

This downloads:
- YFull tree JSON (~14 MB)
- YBrowse SNP database (~400 MB)

### 2. Classify a sample

```bash
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --reference grch38 \
    --output result.json
```

### 3. View results

```bash
cat result.json | jq '.haplogroup, .confidence'
# "R-L21"
# 0.97
```

## Usage

### Single Sample Classification

```bash
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --reference grch38 \
    --output result.json
```

### Multi-Sample VCF

For VCFs containing multiple samples, specify which sample to classify:

```bash
yallhap classify multi_sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --sample NA12878 \
    --output result.json
```

### Batch Processing

Process multiple VCF files into a single TSV:

```bash
yallhap batch sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --output results.tsv
```

### Parallel Processing

Use multiple threads for faster batch processing:

```bash
yallhap batch samples/*.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --threads 16 \
    --output results.tsv
```

With 16 threads, processing 1,000+ samples takes approximately 10 minutes.

### TSV Output Format

Use `--format tsv` for tab-separated output (useful for pipelines):

```bash
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --format tsv \
    --output result.tsv
```

### Reference Genomes

yallHap supports three reference genomes:

```bash
# GRCh37 (hg19) - default for 1000 Genomes
yallhap classify sample.vcf.gz --reference grch37 ...

# GRCh38 (hg38) - current standard
yallhap classify sample.vcf.gz --reference grch38 ...

# T2T-CHM13v2.0 - complete Y chromosome (62 Mb)
yallhap classify sample.vcf.gz --reference t2t ...
```

**T2T Note**: T2T coordinates are computed automatically via liftover from GRCh37/38 positions. Ensure liftover chain files are available (run `python scripts/download_liftover_chains.py`).

## Ancient DNA Mode

yallHap includes specialized handling for ancient DNA samples with post-mortem damage.

### Basic Ancient Mode

Filters C>T and G>A transitions at read termini:

```bash
yallhap classify ancient.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --ancient \
    --min-depth 1 \
    --output result.json
```

### Transversions-Only Mode

Strictest mode for heavily damaged samples (ignores all transitions):

```bash
yallhap classify ancient.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --transversions-only \
    --output result.json
```

### Damage Rescaling

Downweight potentially damaged variants without excluding them:

```bash
yallhap classify ancient.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --ancient \
    --damage-rescale moderate \
    --output result.json
```

Options for `--damage-rescale`:
- `none` (default): No rescaling
- `moderate`: 50% weight reduction for damage-like transitions
- `aggressive`: 80% weight reduction

## Python API

### Single Sample

```python
from yallhap.tree import Tree
from yallhap.snps import SNPDatabase
from yallhap.classifier import HaplogroupClassifier

# Load resources
tree = Tree.from_json("data/yfull_tree.json")
snp_db = SNPDatabase.from_csv("data/ybrowse_snps.csv")

# Create classifier
classifier = HaplogroupClassifier(
    tree=tree,
    snp_db=snp_db,
    reference="grch38",
)

# Classify
result = classifier.classify("sample.vcf.gz")
print(f"{result.sample}: {result.haplogroup} (confidence: {result.confidence:.2f})")
```

### Batch Classification (Multi-Sample VCF)

For multi-sample VCFs, `classify_batch()` is 10x faster than calling `classify()` repeatedly:

```python
# Get list of sample names to classify
samples = ["NA12878", "NA12891", "NA12892"]

# Classify all samples in one pass
results = classifier.classify_batch("multi_sample.vcf.gz", samples)

for result in results:
    print(f"{result.sample}: {result.haplogroup}")
```

### Ancient DNA Mode

```python
classifier = HaplogroupClassifier(
    tree=tree,
    snp_db=snp_db,
    reference="grch37",
    ancient_mode=True,
    transversions_only=True,  # Strictest mode
    damage_rescale="moderate",
)
```

## Output Format

### JSON (default)

```json
{
  "sample": "SAMPLE1",
  "haplogroup": "R-L21",
  "confidence": 0.97,
  "reference": "grch38",
  "tree_version": "YFull",
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
  "path": ["ROOT", "A0-T", "A1", "...", "R-L21"],
  "defining_snps": ["L21"]
}
```

### TSV (for batch processing)

```
sample	haplogroup	confidence	qc1	qc2	qc3	qc4	derived	ancestral	missing
SAMPLE1	R-L21	0.9700	0.9800	1.0000	0.9500	0.9700	145	1089	13
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success (high confidence, ≥0.95) |
| 1 | Classification failed (no haplogroup) |
| 2 | Low confidence (<0.95) |
| 10 | File not found |
| 11 | Invalid input |
| 99 | Unexpected error |

## Quality Scores

| Score | Name | Description |
|-------|------|-------------|
| QC1 | Backbone | Intermediate markers on path to haplogroup match expected states |
| QC2 | Terminal | Defining markers for called haplogroup are present |
| QC3 | Path | Consistency within the called haplogroup branch |
| QC4 | Posterior | Overall posterior probability from likelihood calculation |

## CLI Reference

### `yallhap classify`

Classify a single VCF file.

```
Usage: yallhap classify [OPTIONS] VCF

Options:
  -t, --tree PATH          Path to YFull tree JSON [required]
  -s, --snp-db PATH        Path to SNP database CSV [required]
  -r, --reference TEXT     Reference genome: grch37, grch38, t2t [default: grch38]
  --sample TEXT            Sample name (for multi-sample VCFs)
  --ancient                Enable ancient DNA mode
  --transversions-only     Only use transversions (strictest aDNA mode)
  --damage-rescale TEXT    Rescale quality: none, moderate, aggressive
  --min-depth INTEGER      Minimum read depth [default: 1]
  --min-quality INTEGER    Minimum base quality [default: 20]
  -o, --output PATH        Output file (stdout if omitted)
  --format TEXT            Output format: json, tsv [default: json]
```

### `yallhap batch`

Batch process multiple VCF files.

```
Usage: yallhap batch [OPTIONS] VCF_FILES...

Options:
  -t, --tree PATH          Path to YFull tree JSON [required]
  -s, --snp-db PATH        Path to SNP database CSV [required]
  -r, --reference TEXT     Reference genome: grch37, grch38, t2t [default: grch38]
  --ancient                Enable ancient DNA mode
  --transversions-only     Only use transversions
  --damage-rescale TEXT    Rescale quality: none, moderate, aggressive
  -o, --output PATH        Output TSV file [required]
  --threads INTEGER        Parallel threads [default: 1]
```

### `yallhap download`

Download reference data.

```
Usage: yallhap download [OPTIONS]

Options:
  -o, --output-dir PATH    Output directory [default: data/]
  -f, --force              Overwrite existing files
```

## Pipeline Integration

### Nextflow

See [pipelines/nextflow/](pipelines/nextflow/) for a complete example.

```nextflow
process YALLHAP {
    input:
    path vcf

    output:
    path "*.json"

    script:
    """
    yallhap classify ${vcf} \
        --tree ${params.tree} \
        --snp-db ${params.snp_db} \
        --reference ${params.reference} \
        --output ${vcf.baseName}.json
    """
}
```

### Snakemake

See [pipelines/snakemake/](pipelines/snakemake/) for a complete example.

```python
rule yallhap:
    input:
        vcf="{sample}.vcf.gz"
    output:
        json="{sample}.haplogroup.json"
    params:
        tree=config["yallhap_tree"],
        snp_db=config["yallhap_snps"]
    shell:
        """
        yallhap classify {input.vcf} \
            --tree {params.tree} \
            --snp-db {params.snp_db} \
            --output {output.json}
        """
```

## Experimental Features

### Bayesian Mode

An experimental Bayesian classification mode is available that computes true posterior probabilities over tree paths:

```bash
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --bayesian \
    --output result.json
```

This mode incorporates allelic depth (AD) information when available. However, validation testing showed **no accuracy improvement** over the default heuristic approach on well-characterized samples (identical results across 3,200+ samples from 1000 Genomes, AADR, and gnomAD). Bayesian mode is disabled by default and intended for research purposes.

## Development

```bash
# Clone repository
git clone https://github.com/trianglegrrl/yallhap.git
cd yallhap

# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Run linters
black src/ tests/
ruff check src/ tests/
mypy src/
```

## Citation

If you use yallHap in your research, please cite:

```bibtex
@software{yallhap,
  title = {yallHap: Modern Y-chromosome haplogroup inference},
  year = {2025},
  url = {https://github.com/trianglegrrl/yallhap}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [YFull](https://www.yfull.com/) for maintaining the comprehensive Y-chromosome phylogeny
- [YBrowse](http://ybrowse.org/) for the SNP database
- Yleaf and pathPhynder for algorithmic inspiration
- 1000 Genomes Project and AADR for validation data
