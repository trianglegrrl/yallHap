# yallHap

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
pip install yallhap
```

### Conda

```bash
conda install -c bioconda yallhap
```

### Docker

```bash
docker pull yourusername/yallhap
```

### From source

```bash
git clone https://github.com/yourusername/yallhap.git
cd yallhap
pip install -e ".[dev]"
```

## Quick Start

```bash
# Download reference data
yallhap download --output-dir data/

# Classify a sample
yallhap classify sample.vcf.gz \
    --tree data/yfull_tree.json \
    --snp-db data/ybrowse_snps.csv \
    --reference grch38 \
    --output result.json

# Ancient DNA mode
yallhap classify ancient_sample.vcf.gz \
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
        --output ${vcf.baseName}.json
    """
}
```

### Snakemake

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

## Development

```bash
# Clone repository
git clone https://github.com/yourusername/yallhap.git
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

- YFull for maintaining the comprehensive Y-chromosome phylogeny
- YBrowse for the SNP database
- Yleaf and pathPhynder for algorithmic inspiration
