# yallHap Nextflow Pipeline

Nextflow pipeline for Y-chromosome haplogroup inference using yallHap.

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04.0
- One of: Docker, Singularity, or Conda

## Quick Start

```bash
# Clone the repository
git clone https://github.com/trianglegrrl/yallhap.git
cd yallhap/pipelines/nextflow

# Download reference data
yallhap download --output-dir data/

# Run the pipeline
nextflow run main.nf \
    --input samples.csv \
    --tree data/yfull_tree.json \
    --snp_db data/ybrowse_snps.csv \
    --outdir results/ \
    -profile docker
```

## Input Files

### Sample Sheet (CSV)

Create a CSV file with sample information:

```csv
sample_id,vcf_path
SAMPLE1,/path/to/sample1.vcf.gz
SAMPLE2,/path/to/sample2.vcf.gz
```

### Reference Data

- `yfull_tree.json` - YFull phylogenetic tree
- `ybrowse_snps.csv` - YBrowse SNP database

Download using: `yallhap download`

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Sample sheet CSV (required) | - |
| `--tree` | YFull tree JSON file (required) | - |
| `--snp_db` | SNP database CSV file (required) | - |
| `--reference` | Reference genome: grch37, grch38, t2t | grch38 |
| `--ancient` | Enable ancient DNA mode | false |
| `--outdir` | Output directory | ./results |

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `conda` | Run with Conda environment |
| `slurm` | Submit jobs to SLURM cluster |
| `sge` | Submit jobs to SGE cluster |
| `test` | Run with minimal test resources |

## Output

```
results/
├── individual/           # Per-sample JSON results
│   ├── SAMPLE1.json
│   └── SAMPLE2.json
├── yallhap_results.tsv   # Combined TSV results
└── pipeline_info/        # Nextflow reports
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

## Usage Examples

### Basic Usage

```bash
nextflow run main.nf \
    --input samples.csv \
    --tree data/yfull_tree.json \
    --snp_db data/ybrowse_snps.csv
```

### Ancient DNA Mode

```bash
nextflow run main.nf \
    --input ancient_samples.csv \
    --tree data/yfull_tree.json \
    --snp_db data/ybrowse_snps.csv \
    --ancient true \
    -profile docker
```

### SLURM Cluster

```bash
nextflow run main.nf \
    --input samples.csv \
    --tree data/yfull_tree.json \
    --snp_db data/ybrowse_snps.csv \
    -profile slurm
```

### Resume Failed Run

```bash
nextflow run main.nf -resume
```

## Docker Container

Build the yallHap container:

```bash
cd ../../docker
docker build -t yallhap/yallhap:latest .
```

## License

MIT License

