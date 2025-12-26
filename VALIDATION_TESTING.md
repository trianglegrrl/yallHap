# Validation Testing

This document describes how to reproduce yallHap's validation results and run validation tests on new datasets.

## Overview

yallHap is validated against three primary datasets:

| Dataset | Samples | Reference | Type | Same Major Lineage |
|---------|---------|-----------|------|-------------------|
| 1000 Genomes Phase 3 | 1,233 males | GRCh37 | Modern WGS | **99.76%** |
| AADR v54 | 1,991 | GRCh37 | Ancient DNA | **88.65%** |
| gnomAD HGDP/1000G | 10 | GRCh38 | High-coverage WGS | **100.00%** |

Notes:
- AADR accuracy is measured on samples with properly formatted ground truth (X-YYYY haplogroup names)
- gnomAD samples are selected from the overlap with 1000 Genomes Phase 3, balanced by superpopulation
- Both heuristic and Bayesian modes produce identical results on all datasets

## Benchmark Results Summary

The following results were generated on 2025-12-26 using yallHap with 16 threads:

| Dataset | Tool | Samples | Major Match | Exact | Confidence | SNPs |
|---------|------|---------|-------------|-------|------------|------|
| 1000G Phase 3 | yallHap | 1,233 | 99.76% | 8.84% | 0.994 | 15.4 |
| 1000G Phase 3 | yallHap-Bayes | 1,233 | 99.76% | 8.84% | 0.994 | 15.4 |
| AADR v54 | yallHap | 1,991 | 88.65% | 7.18% | 0.987 | 11.8 |
| AADR v54 | yallHap-Bayes | 1,991 | 88.65% | 7.18% | 0.987 | 11.8 |
| gnomAD HGDP/1KG | yallHap | 10 | 100.00% | 0.00% | 0.988 | 79.4 |
| gnomAD HGDP/1KG | yallHap-Bayesian | 10 | 100.00% | 0.00% | 0.988 | 79.4 |

**Key Finding**: Heuristic and Bayesian modes produce identical classifications on all validation datasets.

## Quick Start

### Option 1: Download from Zenodo (Recommended)

Download the pre-packaged validation data:

```bash
# Download validation bundle (~10 GB with gnomAD data)
wget https://zenodo.org/record/XXXXXX/files/yallhap-validation-v1.tar.gz

# Verify checksum
wget https://zenodo.org/record/XXXXXX/files/yallhap-validation-v1.md5
md5sum -c yallhap-validation-v1.md5

# Extract to data directory
tar -xzf yallhap-validation-v1.tar.gz -C data/
```

### Option 2: Download from Original Sources

```bash
# 1000 Genomes Phase 3 data (GRCh37)
python scripts/download_1kg_validation.py

# Ancient DNA data (requires AADR access)
python scripts/download_ancient_test_data.py

# gnomAD high-coverage data (GRCh38, ~9GB)
python scripts/download_gnomad_highcov.py

# T2T liftover chains (optional)
python scripts/download_liftover_chains.py
```

## Prerequisites

```bash
# Install yallhap with dev dependencies
pip install -e ".[dev]"

# Verify installation
yallhap --version

# Required external tools
bcftools --version  # For VCF manipulation
tabix --version     # For VCF indexing
```

## Validation Data Contents

After downloading, your `data/` directory should contain:

```
data/
├── yfull_tree.json              # YFull phylogenetic tree (~14 MB)
├── ybrowse_snps.csv             # GRCh38 SNP database (~400 MB)
├── validation/
│   ├── 1kg_chrY_phase3.vcf.gz   # 1000 Genomes Y-chr VCF (5.4 MB)
│   ├── 1kg_chrY_phase3.vcf.gz.tbi
│   ├── 1kg_sample_panel.txt     # Sample population mappings
│   ├── poznik2016_haplogroups.tsv  # Ground truth (24 KB)
│   └── ybrowse_snps_hg19.csv    # GRCh37 SNP positions (48 MB)
├── ancient/
│   ├── aadr_chrY_v2.vcf.gz      # Ancient DNA VCF (97 MB)
│   ├── aadr_chrY_v2.vcf.gz.tbi
│   └── aadr_1240k_ground_truth.tsv  # Ancient ground truth (252 KB)
├── validation_highcov/
│   └── vcf/
│       ├── gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz  # Full VCF (~9 GB)
│       ├── gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz.tbi
│       ├── gnomad_subset_10_filtered.vcf.gz  # Benchmark subset
│       └── diagnostic_positions.tsv  # SNP positions for filtering
└── liftover/
    ├── grch38-chm13v2.chain     # T2T liftover (6 MB)
    └── hg19-chm13v2.chain       # T2T liftover (6 MB)
```

## Running Validations

### All-in-One Benchmark (Reproducing Published Results)

Run the comprehensive benchmark suite exactly as used for the paper:

```bash
# Full benchmark with 16 threads (recommended)
python scripts/run_benchmarks.py --threads 16

# Quick development benchmark with subsampling
python scripts/run_benchmarks.py --subsample 50 --threads 16
```

This runs validation against all available datasets and generates:
- Console summary table
- `paper/benchmark_results.json` - Machine-readable results
- `paper/benchmark_summary.tsv` - Summary table
- `paper/benchmark_detailed.tsv` - Per-sample results

**Expected runtime:** ~30 minutes with 16 threads for full benchmark (both heuristic and Bayesian modes on all samples).

### 1000 Genomes Validation

This validates yallHap against 1,233 male samples from 1000 Genomes Phase 3.

```bash
python scripts/validate_1kg.py
```

**Expected output:**

```
Classifying 1233 samples...

Results:
  Same major lineage: 1230/1233 (99.76%)
  Wrong lineage: 3/1233 (0.24%)
  Exact match: 109/1233 (8.84%)

Mean confidence: 0.994
Mean derived SNPs: 15.4

Notes:
- "Exact match" requires identical haplogroup names
- Low exact match rate is expected due to different tree versions
  (yallHap uses YFull 2024, ground truth uses ISOGG 2016)
- "Same major lineage" means the first letter matches (R, J, I, etc.)
- 3 wrong samples: 2 rare A0 haplogroups, 1 NO/K confusion
```

#### Quick Validation (100 samples)

For faster testing during development:

```python
python << 'EOF'
from yallhap.classifier import HaplogroupClassifier
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree
import csv

tree = Tree.from_json('data/yfull_tree.json')
db = SNPDatabase.from_ybrowse_gff_csv('data/validation/ybrowse_snps_hg19.csv')
c = HaplogroupClassifier(tree=tree, snp_db=db, reference='grch37')

# Load ground truth
gt = {}
with open('data/validation/poznik2016_haplogroups.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gt[row['sample_id']] = row['haplogroup']

# Classify 100 samples
samples = list(gt.keys())[:100]
results = c.classify_batch('data/validation/1kg_chrY_phase3.vcf.gz', samples)

# Check accuracy
correct = sum(1 for i, s in enumerate(samples) if gt[s][0] == results[i].haplogroup[0])
print(f"Accuracy: {correct}/100 ({correct}%) same major lineage")
EOF
```

### Ancient DNA Validation

This validates yallHap against ancient samples from the Allen Ancient DNA Resource (AADR).

**Important:** AADR ground truth uses mixed nomenclature - some samples have proper haplogroup names (e.g., "R-L21") while others have SNP-only names (e.g., "M458"). The benchmark script automatically filters to samples with proper X-YYYY format (1,991 of ~2,000 samples).

```bash
python scripts/validate_ancient.py
```

**Expected results:**
- Same major lineage: 88.65% (1,991 samples)
- Mean derived SNPs: 11.8 (low coverage typical for aDNA)
- Mean confidence: 0.987

### gnomAD High-Coverage Validation

This validates yallHap against gnomAD HGDP/1000G high-coverage samples with allelic depth (AD) information.

#### Downloading gnomAD Data

```bash
# Full dataset (~9GB, required for comprehensive validation)
python scripts/download_gnomad_highcov.py

# Or extract just the validation subset (faster)
python scripts/download_gnomad_highcov.py --sample-only
```

**Data source:**
- URL: `https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz`
- Reference: GRCh38
- Samples: 4,151 (HGDP + 1000 Genomes high-coverage)
- Key feature: AD (allelic depth) fields for Bayesian classification

**Note:** The source index file may be outdated. The download script rebuilds the index locally using `tabix`.

#### How gnomAD Samples Are Selected

The benchmark script selects 10 samples (default, configurable via `--gnomad-samples N`) that:
1. Exist in both gnomAD VCF AND 1000 Genomes Phase 3 ground truth
2. Are balanced by superpopulation (AFR, AMR, EAS, EUR, SAS)

This allows cross-reference validation using the well-established Poznik 2016 haplogroup assignments.

#### Running gnomAD Validation

```bash
# Quick test with 5 samples
python << 'EOF'
from yallhap.classifier import HaplogroupClassifier
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree

tree = Tree.from_json('data/yfull_tree.json')
db = SNPDatabase.from_ybrowse_gff_csv('data/ybrowse_snps.csv')

# Heuristic mode (default)
c = HaplogroupClassifier(
    tree=tree,
    snp_db=db,
    reference='grch38',
    min_depth=1,
)

# Sample validation (first 5 from benchmark)
samples = ["HG00096", "HG00101", "HG00103", "HG00105", "HG00107"]
vcf = 'data/validation_highcov/vcf/gnomad_subset_10_filtered.vcf.gz'

print("Sample       | yallHap         | Confidence | SNPs")
print("-------------|-----------------|------------|-----")
for sample in samples:
    result = c.classify(vcf, sample)
    print(f"{sample:12} | {result.haplogroup:15} | {result.confidence:.3f}      | {result.snp_stats.derived}")
EOF
```

**Expected results:**
- Same major lineage: 100% (10 validation samples)
- All R haplogroups correctly identified as R-*
- All I haplogroups correctly identified as I-*
- High confidence (>0.98) with 60-100 derived SNPs per sample

## Understanding Results

### Major Lineage vs Exact Match

yallHap uses the YFull 2024 tree, while many ground truth datasets use older ISOGG nomenclature. This causes low "exact match" rates even when classifications are correct.

**Example:**
- Ground truth: `R1b1a1a2a1a2c` (ISOGG 2016)
- yallHap call: `R-L21` (YFull 2024)
- These represent the **same haplogroup** with different naming conventions

**What to check:**
1. **Same major lineage** (first letter): This should be very high (>99%)
2. **Phylogenetic consistency**: The called haplogroup should be on the correct branch
3. **Specificity**: yallHap may call a more or less specific haplogroup depending on coverage

### Interpreting Low Accuracy

If you see low accuracy:

1. **Check reference genome**: Ensure VCF reference matches `--reference` flag
2. **Check sample names**: Sample IDs must match between VCF and ground truth
3. **Check coverage**: Low-coverage samples may get less specific calls
4. **Check for damage**: Ancient DNA samples need `--ancient` or `--transversions-only`

### Heuristic vs Bayesian Mode

| Mode | Use Case | AD Required | Accuracy | Default |
|------|----------|-------------|----------|---------|
| Heuristic | General use | No | 99.76% (1KG) | Yes |
| Bayesian | Research | Optional | 99.76% (1KG) | No |

**Key finding from validation:** Both modes produce identical results on all tested datasets. Bayesian mode is experimental and disabled by default.

## Adding New Validation Sets

### Creating Ground Truth Files

Ground truth files should be TSV format with at least these columns:

```
sample_id	haplogroup
SAMPLE1	R-L21
SAMPLE2	J-M172
SAMPLE3	I-M253
```

### Running Custom Validation

```python
from yallhap.validation import ValidationRunner, load_ground_truth

# Load your ground truth
ground_truth = load_ground_truth("my_ground_truth.tsv")

# Run validation
runner = ValidationRunner(
    tree_path="data/yfull_tree.json",
    snp_db_path="data/ybrowse_snps.csv",
    vcf_path="my_samples.vcf.gz",
    reference="grch38",
)

results = runner.run(ground_truth)

# Print metrics
print(f"Exact match: {results.exact_match_rate:.1%}")
print(f"Same lineage: {results.same_lineage_rate:.1%}")
```

## Validation Scripts Reference

| Script | Purpose |
|--------|---------|
| `scripts/download_1kg_validation.py` | Downloads 1000 Genomes Phase 3 VCF and panel |
| `scripts/download_ancient_test_data.py` | Creates ancient DNA ground truth TSV |
| `scripts/download_gnomad_highcov.py` | Downloads gnomAD HGDP/1KG high-coverage VCF |
| `scripts/download_liftover_chains.py` | Downloads T2T liftover chain files |
| `scripts/run_benchmarks.py` | Runs comprehensive benchmark suite |
| `scripts/validate_1kg.py` | Runs validation against 1000 Genomes |
| `scripts/eigenstrat_to_vcf.py` | Converts EIGENSTRAT format to VCF |
| `scripts/prep_for_zenodo.py` | Packages validation data for Zenodo upload |

## Troubleshooting

### "No Y chromosome found in VCF"

The VCF must contain a Y chromosome contig. Check with:

```bash
bcftools view -h sample.vcf.gz | grep "^##contig.*Y"
```

### "Position index is empty"

Ensure the SNP database matches your reference:
- `ybrowse_snps.csv` for GRCh38
- `ybrowse_snps_hg19.csv` for GRCh37

### "Sample not found in VCF"

List samples in VCF with:

```bash
bcftools query -l sample.vcf.gz
```

### Low accuracy on ancient samples

Try different ancient DNA modes:

```python
# Start with transversions-only (strictest)
classifier = HaplogroupClassifier(..., transversions_only=True)

# If too few variants, try damage rescaling instead
classifier = HaplogroupClassifier(..., ancient_mode=True, damage_rescale="moderate")
```

### Slow classification on multi-sample VCF

For large multi-sample VCFs (like gnomAD with 4,151 samples), the benchmark script automatically:
1. Creates a diagnostic positions file (`diagnostic_positions.tsv`)
2. Uses bcftools to extract only relevant positions and samples
3. Creates a cached subset VCF for faster subsequent runs

If you need to manually extract:

```bash
# Extract single sample
bcftools view -s HG00096 input.vcf.gz -Oz -o HG00096.vcf.gz
tabix HG00096.vcf.gz

# Then classify
yallhap classify --vcf HG00096.vcf.gz
```

### gnomAD index file outdated

If you see warnings about outdated index, rebuild locally:

```bash
tabix -f data/validation_highcov/vcf/gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz
```

## Reproducing Paper Results

To exactly reproduce the benchmark results from the yallHap paper:

```bash
# 1. Ensure all data is downloaded
python scripts/download_1kg_validation.py
python scripts/download_ancient_test_data.py
python scripts/download_gnomad_highcov.py

# 2. Run full benchmarks with 16 threads
python scripts/run_benchmarks.py --threads 16

# 3. Check results match paper
cat paper/benchmark_results.json
```

Expected file outputs:
- `paper/benchmark_results.json` - Full machine-readable results
- `paper/benchmark_summary.tsv` - Summary table matching paper Table 2
