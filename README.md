# YClade

Modern, pipeline-friendly Y-chromosome haplogroup inference.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Features

- **YFull tree**: Uses the most comprehensive Y-chromosome phylogeny (185,780+ SNPs)
- **Probabilistic scoring**: Likelihood-based confidence scores, not just SNP counting
- **Ancient DNA support**: Built-in damage filtering for aDNA samples
- **Multiple references**: Supports GRCh37, GRCh38, and T2T-CHM13v2.0
- **Pipeline-friendly**: Proper exit codes, JSON/TSV output, batch processing
- **Bioconda/Docker**: Easy installation and containerized execution

## Installation

### pip (recommended)

```bash
pip install yclade
```

### Conda

```bash
conda install -c bioconda yclade
```

### Docker

```bash
docker pull yourusername/yclade
```

### From source

```bash
git clone https://github.com/yourusername/yclade.git
cd yclade
pip install -e ".[dev]"
```

## Quick Start

```bash
# Download reference data
yclade download --output-dir data/

# Classify a sample
yclade classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --reference grch38 \
    --output result.json

# Ancient DNA mode
yclade classify ancient_sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --ancient \
    --min-depth 1 \
    --output result.json
```

## Output Format

### JSON (default)

```json
{
  "sample": "SAMPLE1",
  "haplogroup": "R-L21",
  "confidence": 0.97,
  "reference": "grch38",
  "tree_version": "YFull v13.06",
  "quality_scores": {
    "qc1_backbone": 0.98,
    "qc2_terminal": 1.0,
    "qc3_path": 0.95,
    "qc4_posterior": 0.97
  },
  "path": ["Y-Adam", "A0-T", "...", "R-L21"],
  "defining_snps": ["L21/M529/S145"]
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
| 0 | Success (high confidence) |
| 1 | Classification failed (no haplogroup) |
| 2 | Low confidence (<0.95) |
| 10 | File not found |
| 11 | Invalid input |
| 99 | Unexpected error |

## Quality Scores

- **QC1 (backbone)**: Intermediate markers match expected states
- **QC2 (terminal)**: Defining markers for called haplogroup
- **QC3 (path)**: Within-haplogroup consistency
- **QC4 (posterior)**: Overall posterior probability

## Pipeline Integration

### Nextflow

```nextflow
process YCLADE {
    input:
    path vcf

    output:
    path "*.json"

    script:
    """
    yclade classify ${vcf} \
        --tree ${params.tree} \
        --snp-db ${params.snp_db} \
        --output ${vcf.baseName}.json
    """
}
```

### Snakemake

```python
rule yclade:
    input:
        vcf="{sample}.vcf.gz"
    output:
        json="{sample}.haplogroup.json"
    params:
        tree=config["yclade_tree"],
        snp_db=config["yclade_snps"]
    shell:
        """
        yclade classify {input.vcf} \
            --tree {params.tree} \
            --snp-db {params.snp_db} \
            --output {output.json}
        """
```

## Development

```bash
# Clone repository
git clone https://github.com/yourusername/yclade.git
cd yclade

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

If you use YClade in your research, please cite:

```bibtex
@software{yclade,
  title = {YClade: Modern Y-chromosome haplogroup inference},
  year = {2024},
  url = {https://github.com/yourusername/yclade}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- YFull for maintaining the comprehensive Y-chromosome phylogeny
- YBrowse for the SNP database
- Yleaf and pathPhynder for algorithmic inspiration
