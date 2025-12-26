# yallHap Snakemake Workflow

Snakemake workflow for Y-chromosome haplogroup inference using yallHap.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) >= 7.0
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)

## Quick Start

```bash
# Clone the repository
git clone https://github.com/yourusername/yallhap.git
cd yallhap/pipelines/snakemake

# Create sample sheet
cat > samples.csv << EOF
sample_id,vcf_path
SAMPLE1,/path/to/sample1.vcf.gz
SAMPLE2,/path/to/sample2.vcf.gz
EOF

# Download reference data
snakemake download_data --use-conda --cores 1

# Run the workflow
snakemake --use-conda --cores 4
```

## Configuration

Edit `config.yaml` to configure the workflow:

```yaml
# Sample sheet
samples: "samples.csv"

# Reference data
tree: "data/yfull_tree.json"
snp_db: "data/ybrowse_snps.csv"

# Reference genome
reference: "grch38"  # grch37, grch38, or t2t

# Output directory
outdir: "results"

# Ancient DNA mode
ancient: false
```

## Input Files

### Sample Sheet (CSV)

Create a CSV file with sample information:

```csv
sample_id,vcf_path
SAMPLE1,/path/to/sample1.vcf.gz
SAMPLE2,/path/to/sample2.vcf.gz
```

## Output

```
results/
├── individual/           # Per-sample results
│   ├── SAMPLE1.json
│   ├── SAMPLE1.tsv
│   ├── SAMPLE2.json
│   └── SAMPLE2.tsv
├── yallhap_results.tsv   # Combined results
└── logs/                 # Log files
    ├── SAMPLE1.log
    └── SAMPLE2.log
```

## Usage Examples

### Basic Usage

```bash
snakemake --use-conda --cores 4
```

### Dry Run (Preview)

```bash
snakemake --use-conda --cores 4 -n
```

### Ancient DNA Mode

Edit `config.yaml`:
```yaml
ancient: true
```

Then run:
```bash
snakemake --use-conda --cores 4
```

### SLURM Cluster

Create a cluster profile or use command line:

```bash
snakemake --use-conda \
    --cluster "sbatch -p standard -c {threads} --mem={resources.mem_mb}M" \
    --jobs 100
```

### With Mamba (faster)

```bash
snakemake --use-conda --conda-frontend mamba --cores 4
```

## Rules

| Rule | Description |
|------|-------------|
| `all` | Default target: run all samples |
| `classify_sample` | Classify a single sample |
| `merge_results` | Merge all TSV results |
| `batch_classify` | Batch classify (alternative) |
| `download_data` | Download reference data |

## DAG Visualization

Generate a workflow diagram:

```bash
snakemake --dag | dot -Tpng > workflow.png
```

## License

MIT License

