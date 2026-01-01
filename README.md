# yallHap

Modern, pipeline-friendly Y-chromosome haplogroup inference.

[![License: PolyForm Noncommercial](https://img.shields.io/badge/License-PolyForm%20Noncommercial-blue.svg)](https://polyformproject.org/licenses/noncommercial/1.0.0)
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
| 1000 Genomes Phase 3 | 1,233 | **99.8%** (95% CI: 99.3-100%) | GRCh37 | Modern WGS, heuristic mode |
| AADR Ancient DNA | 7,333 | **90.7%** Bayesian / **88.3%** Heuristic | GRCh37 | Full dataset, stratified by variant density |
| gnomAD HGDP/1KG | 1,231 | **99.9%** (95% CI: 99.5-100%) | GRCh38 | High-coverage WGS |

**1000 Genomes details:**
- Only 3 misclassified samples (2 rare A0 haplogroups, 1 NO/K confusion)
- Mean confidence: 0.994
- Mean derived SNPs: 15.4

**AADR Ancient DNA details (7,333 samples):**
- Overall: 90.7% accuracy with Bayesian ancient mode vs 88.3% with heuristic transversions-only
- Stratified by variant density: <1% (33.7%), 1-4% (37.9%), 4-10% (71.7%), 10-50% (97.8%), ≥50% (99.0%)
- At ≥10% variant density, both modes achieve 97-99% accuracy, comparable to modern WGS
- **Bayesian mode recommended for 4-10% variant density** (+12-24 pp improvement)
- Variant density = (called variants / total variants in chrY VCF) × 100%

**gnomAD High-Coverage details:**
- 200 samples randomly selected from 1,231 overlapping with 1000 Genomes
- 30× high-coverage whole-genome sequencing
- Mean derived SNPs: 26.7
- 95% confidence interval: 98.17-100%

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
- YBrowse SNP database for GRCh38 (~430 MB)
- YBrowse SNP database for GRCh37 (~50 MB)

### 2. Classify a sample

Use the SNP database matching your VCF's reference genome:

```bash
# For GRCh38/hg38 VCFs
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --reference grch38 \
    --output result.json

# For GRCh37/hg19 VCFs
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch37.csv \
    --reference grch37 \
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
    --snp-db data/ybrowse_snps_grch38.csv \
    --reference grch38 \
    --output result.json
```

### Multi-Sample VCF

For VCFs containing multiple samples, specify which sample to classify:

```bash
yallhap classify multi_sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --sample NA12878 \
    --output result.json
```

### Batch Processing

Process multiple VCF files into a single TSV:

```bash
yallhap batch sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --output results.tsv
```

### Parallel Processing

Use multiple threads for faster batch processing:

```bash
yallhap batch samples/*.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --threads 16 \
    --output results.tsv
```

With 16 threads, processing 1,000+ samples takes approximately 10 minutes.

### TSV Output Format

Use `--format tsv` for tab-separated output (useful for pipelines):

```bash
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --format tsv \
    --output result.tsv
```

### Reference Genomes

yallHap supports three reference genomes. **Use the SNP database matching your VCF's reference**:

| VCF Reference | SNP Database | `-r` flag |
|---------------|--------------|-----------|
| GRCh37/hg19 | `ybrowse_snps_grch37.csv` | `grch37` |
| GRCh38/hg38 | `ybrowse_snps_grch38.csv` | `grch38` |
| T2T-CHM13v2.0 | `ybrowse_snps_grch38.csv` | `t2t` |

```bash
# GRCh37 (hg19) - 1000 Genomes Phase 3, many ancient DNA datasets
yallhap classify sample.vcf.gz \
    -s data/ybrowse_snps_grch37.csv -r grch37 ...

# GRCh38 (hg38) - current standard, gnomAD, most modern studies
yallhap classify sample.vcf.gz \
    -s data/ybrowse_snps_grch38.csv -r grch38 ...

# T2T-CHM13v2.0 - complete Y chromosome (62 Mb)
yallhap classify sample.vcf.gz \
    -s data/ybrowse_snps_grch38.csv -r t2t ...
```

**T2T Note**: T2T coordinates are computed automatically via liftover from GRCh38 positions. Ensure liftover chain files are available (run `python scripts/download_liftover_chains.py`).

## Ancient DNA Mode

yallHap includes specialized handling for ancient DNA samples with post-mortem damage.

### Recommended: Bayesian Ancient Mode

For ancient DNA samples with moderate variant density (4–10%), **Bayesian ancient mode is recommended**, achieving +12–24 percentage point improvement over heuristic mode in this range:

```bash
yallhap classify ancient.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --ancient \
    --bayesian \
    --output result.json
```

**Variant density** is calculated as `(called variants / total variants in chrY VCF) × 100%`. You can estimate this from your VCF or calculate it directly. At ≥10% variant density, both modes achieve comparable accuracy (97–99%); below 4%, classification is unreliable regardless of mode.

### Basic Ancient Mode

Filters C>T and G>A transitions at read termini:

```bash
yallhap classify ancient.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --ancient \
    --min-depth 1 \
    --output result.json
```

### Transversions-Only Mode

Strictest mode for heavily damaged samples (ignores all transitions):

```bash
yallhap classify ancient.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --transversions-only \
    --output result.json
```

### Damage Rescaling

Downweight potentially damaged variants without excluding them:

```bash
yallhap classify ancient.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
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
snp_db = SNPDatabase.from_csv("data/ybrowse_snps_grch38.csv")

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
# Recommended: Bayesian ancient mode for moderate variant density (4-10%)
classifier = HaplogroupClassifier(
    tree=tree,
    snp_db=snp_db,
    reference="grch37",
    ancient_mode=True,
    bayesian=True,  # Recommended for 4-10% variant density
)

# Alternative: Transversions-only mode (strictest filtering)
classifier = HaplogroupClassifier(
    tree=tree,
    snp_db=snp_db,
    reference="grch37",
    ancient_mode=True,
    transversions_only=True,
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
  "tree_version": "YFull (185780 SNPs, hash: a1b2c3d4)",
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

### Reproducibility

The `tree_version` field includes a hash of the tree file content, enabling exact reproducibility. When citing yallHap results, include the `tree_version` value to document the exact phylogeny version used. The format is:

```
YFull (<snp_count> SNPs, hash: <8-char SHA256>)
```

Example: `"YFull (185780 SNPs, hash: a1b2c3d4)"`

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

Download reference data (YFull tree + SNP databases for all reference genomes).

```
Usage: yallhap download [OPTIONS]

Options:
  -o, --output-dir PATH    Output directory [default: data/]
  -f, --force              Overwrite existing files
```

Downloads:
- `yfull_tree.json` - YFull phylogenetic tree (~14 MB)
- `ybrowse_snps_grch38.csv` - SNP positions for GRCh38/hg38 (~430 MB)
- `ybrowse_snps_grch37.csv` - SNP positions for GRCh37/hg19 (~50 MB)

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

A Bayesian classification mode is available that computes posterior probabilities over tree paths using log-likelihood ratios:

```bash
# For modern samples
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --bayesian \
    --output result.json

# For ancient DNA (recommended for 4-10% variant density)
yallhap classify ancient.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps_grch38.csv \
    --ancient \
    --bayesian \
    --output result.json
```

**Performance:** On modern high-coverage samples (1000 Genomes, gnomAD), Bayesian mode produces identical results to heuristic mode—no accuracy improvement. However, **for ancient DNA with moderate variant density (4–10%)**, Bayesian ancient mode achieves +12–24 percentage point improvement over heuristic mode (71.7% vs 52.4% accuracy). On the full AADR ancient DNA dataset (7,333 samples), Bayesian ancient mode achieves 90.7% accuracy vs 88.3% for heuristic transversions-only mode.

This mode incorporates allelic depth (AD) information when available and uses adjusted error rates for ancient DNA damage modeling. For modern samples, heuristic mode is recommended for speed; for ancient DNA at 4–10% variant density, Bayesian mode is recommended for improved accuracy.

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

PolyForm Noncommercial License 1.0.0 - see [LICENSE](LICENSE) for details.

This license allows use for noncommercial purposes, including research, education, and personal projects. For commercial use, please contact the maintainers.

## Acknowledgments

- [YFull](https://www.yfull.com/) for maintaining the comprehensive Y-chromosome phylogeny
- [YBrowse](http://ybrowse.org/) for the SNP database
- Yleaf and pathPhynder for algorithmic inspiration
- 1000 Genomes Project and AADR for validation data
