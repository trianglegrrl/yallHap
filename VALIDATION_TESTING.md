# Validation Testing

This document describes how to reproduce yallHap's validation results and run validation tests on new datasets.

## Overview

yallHap is validated against three primary datasets:

| Dataset | Samples | Reference | Type | Same Major Lineage |
|---------|---------|-----------|------|-------------------|
| 1000 Genomes Phase 3 | 1,233 males | GRCh37 | Modern WGS | **99.8%** (95% CI: 99.2-99.9%) |
| gnomAD HGDP/1000G | 1,231 | GRCh38 | High-coverage WGS | **99.9%** (95% CI: 99.5-100%) |
| AADR v54 (full) | 7,333 | GRCh37 | Ancient DNA | **90.7%** Bayesian / **88.3%** Heuristic |

Notes:
- AADR accuracy is stratified by **variant density** (% of called variants in chrY VCF):
  - <1%: 33.7% Bayesian / 28.7% Heuristic (n=101)
  - 1-10%: 58.4% Bayesian / 46.3% Heuristic (n=1,204)
  - 10-50%: 97.8% Bayesian / 97.0% Heuristic (n=4,101)
  - ≥50%: 99.0% Bayesian / 99.1% Heuristic (n=1,927)
- Full AADR dataset used (no subsampling)
- gnomAD samples are selected from the overlap with 1000 Genomes Phase 3 ground truth
- **For ancient DNA, Bayesian mode is recommended for 4-10% variant density** (+12-24 pp improvement in that range)

## Benchmark Results Summary

The following results were generated on 2025-12-27 using yallHap v0.2.0 on Apple M3 Max (16 cores, 64GB RAM):

| Dataset | Tool | Samples | Major Match | Confidence | Runtime |
|---------|------|---------|-------------|------------|---------|
| 1000G Phase 3 | yallHap | 1,233 | 99.8% | 0.994 | 314s |
| 1000G Phase 3 | yallHap-Bayes | 1,233 | 99.8% | 0.994 | 314s |
| gnomAD High-Cov | yallHap | 1,231 | 99.9% | 0.991 | 1,508s |
| gnomAD High-Cov | yallHap-Bayes | 1,231 | 99.9% | 0.991 | 1,508s |
| AADR v54 (full) | yallHap-Heuristic-TV | 7,333 | 88.3% | 0.983 | 1,367s |
| AADR v54 (full) | yallHap-Bayesian-Ancient | 7,333 | 90.7% | 0.962 | 1,367s |

**Key Finding**: For modern high-coverage samples, both modes produce identical results. For ancient DNA, **Bayesian ancient mode is recommended for samples with 4-10% variant density** (+12-24 pp improvement in that range). At very low (<4%) or high (>10%) densities, mode choice matters less.

## Quick Start

### Option 1: Download from Zenodo (Recommended)

Download the pre-packaged validation data:

```bash
# Download validation bundle (~700 MB, excludes full gnomAD VCF)
wget https://zenodo.org/record/XXXXXX/files/yallhap-validation-v3.tar.gz

# Verify checksum
wget https://zenodo.org/record/XXXXXX/files/yallhap-validation-v3.md5
md5sum -c yallhap-validation-v3.md5

# Extract to data directory
tar -xzf yallhap-validation-v3.tar.gz -C data/
```

The Zenodo bundle includes:
- Core reference data (YFull tree, YBrowse SNPs, ISOGG SNPs)
- 1000 Genomes Phase 3 validation data
- AADR v54 ancient DNA validation data
- gnomAD pre-extracted diagnostic subset (not full 9GB VCF)
- 9 ancient genome BAM/VCF files for tool comparison
- T2T liftover chain files
- Benchmark results and validation reports

### Option 2: Download from Original Sources

```bash
# 1000 Genomes Phase 3 data (GRCh37)
python scripts/download_1kg_validation.py

# Ancient DNA data (requires AADR access)
python scripts/download_ancient_test_data.py

# gnomAD high-coverage data (GRCh38, ~9GB full VCF)
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
├── isogg_snps_grch38.txt        # ISOGG SNP database (~3 MB)
├── validation/
│   ├── 1kg_chrY_phase3.vcf.gz   # 1000 Genomes Y-chr VCF (5.4 MB)
│   ├── 1kg_chrY_phase3.vcf.gz.tbi
│   ├── 1kg_sample_panel.txt     # Sample population mappings
│   ├── poznik2016_haplogroups.tsv  # Ground truth (24 KB)
│   └── ybrowse_snps_hg19.csv    # GRCh37 SNP positions (48 MB)
├── ancient/
│   ├── aadr_chrY_v2.vcf.gz      # Ancient DNA VCF (97 MB)
│   ├── aadr_chrY_v2.vcf.gz.tbi
│   ├── aadr_1240k_ground_truth.tsv  # Ancient ground truth
│   └── aadr_variant_density.json    # Cached variant density
├── ancient_genomes/             # Individual ancient samples for tool comparison
│   ├── I0231.390k.chrY.{bam,vcf.gz}      # Yamnaya
│   ├── I0443.390k.chrY.{bam,vcf.gz}      # Yamnaya
│   ├── Kennewick_*.chrY.{bam,vcf.gz}     # Paleoamerican
│   ├── SB524A_*.chrY.{bam,vcf.gz}        # Cheddar Man (lib 1)
│   ├── SB524A2_*.chrY.{bam,vcf.gz}       # Cheddar Man (lib 2)
│   ├── VK287.final.chrY.{bam,vcf.gz}     # Viking
│   ├── VK292.final.chrY.{bam,vcf.gz}     # Viking
│   ├── VK296.final.chrY.{bam,vcf.gz}     # Viking
│   └── VK582.final.chrY.{bam,vcf.gz}     # Viking
├── validation_highcov/
│   └── vcf/
│       ├── gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz  # Full VCF (~9 GB)
│       ├── gnomad_1kg_shared_diagnostic.vcf.gz  # Pre-extracted subset
│       └── diagnostic_positions.tsv  # SNP positions for filtering
└── liftover/
    ├── grch38-chm13v2.chain     # T2T liftover (6 MB)
    └── hg19-chm13v2.chain       # T2T liftover (6 MB)
```

## Running Validations

### All-in-One Benchmark (Reproducing Published Results)

Run the comprehensive benchmark suite exactly as used for the paper:

```bash
# Full benchmark with 12 threads (recommended)
python scripts/run_benchmarks.py --threads 12

# Quick development benchmark with subsampling
python scripts/run_benchmarks.py --subsample 50 --threads 12
```

This runs validation against all available datasets and generates:
- Console summary table
- `results/benchmark_results.json` - Machine-readable results

**Expected runtime:** ~53 minutes with 12 threads for full benchmark on Apple M3 Max.

### 1000 Genomes Validation

This validates yallHap against 1,233 male samples from 1000 Genomes Phase 3.

```bash
python scripts/validate_1kg.py
```

**Expected output:**

```
Classifying 1233 samples...

Results:
  Same major lineage: 1230/1233 (99.8%)
  Wrong lineage: 3/1233 (0.2%)
  Exact match: 109/1233 (8.84%)

Mean confidence: 0.994
Mean derived SNPs: 15.4

Notes:
- "Exact match" requires identical haplogroup names
- Low exact match rate is expected due to different tree versions
  (yallHap uses YFull 2025, ground truth uses ISOGG 2016)
- "Same major lineage" means the first letter matches (R, J, I, etc.)
- 3 wrong samples: 2 rare A0 haplogroups (YBrowse database issue), 1 NO/K confusion
```

### Ancient DNA Validation

This validates yallHap against ancient samples from the Allen Ancient DNA Resource (AADR).

**Important:** AADR samples are stratified by **variant density** (percentage of called variants in the chrY VCF), not by sequencing coverage metadata. This provides a more accurate measure of data quality directly from the VCF.

#### Stratified Validation

```bash
# Run stratified validation with up to 500 samples per density bin
python scripts/validate_aadr_stratified.py \
  --samples-per-bin 500 \
  --threads 12 \
  -o results/aadr_density_stratified.md
```

**Expected results (full dataset, by variant density):**

| Density Bin | Samples | Heuristic TV | Bayesian Ancient | Δ |
|-------------|---------|--------------|------------------|---|
| <1% | 101 | 28.7% | 33.7% | +5.0 pp |
| 1-10% | 1,204 | 46.3% | 58.4% | +12.1 pp |
| 10-50% | 4,101 | 97.0% | 97.8% | +0.8 pp |
| ≥50% | 1,927 | 99.1% | 99.0% | −0.1 pp |
| **Overall** | **7,333** | **88.3%** | **90.7%** | **+2.4 pp** |

Fine-grained analysis of the 1-10% range (see `results/aadr_1_10_pct_stratified_full.md`) shows Bayesian gains concentrated at 4-10% density (+12 to +24 pp).

#### Decile Analysis

For finer-grained analysis:

```bash
python scripts/validate_aadr_stratified.py \
  --deciles \
  --samples-per-bin 200 \
  --threads 12 \
  -o results/aadr_deciles.md
```

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

**Expected results:**
- Same major lineage: 99.9% (1,231 samples with ground truth overlap)
- 1 A0 misclassification (YBrowse database issue, same as 1KG)
- High confidence (>0.99) with 20-30 derived SNPs per sample on average

### Ancient Tool Comparison

Compare yallHap against yhaplo, Yleaf, and pathPhynder on the 9 ancient samples:

```bash
python scripts/gather_validation_and_comparative_data.py \
  -o results/validation_report_ancient_local.md \
  --threads 12
```

**Key findings from tool comparison:**
- yhaplo misclassifies Kennewick Man (calls CT instead of Q)
- Yleaf fails to produce calls for all 9 ancient samples
- pathPhynder achieves deepest resolution but requires BAM input and is slower
- yallHap Bayesian mode matches ground truth in 4/9 samples vs 1/9 for heuristic

## ISOGG Output

yallHap v0.2.0 supports ISOGG haplogroup nomenclature output:

```bash
# CLI usage
yallhap classify --vcf sample.vcf.gz --isogg

# Python API
result = classifier.classify(vcf_path, sample_id)
if result.isogg_haplogroup:
    print(f"ISOGG: {result.isogg_haplogroup}")
```

## Ancient DNA Modes

yallHap provides four approaches for ancient DNA:

1. **Standard mode**: No special handling; suitable for well-preserved samples
2. **Transversions-only mode** (`--transversions-only`): Excludes all transitions; strictest
3. **Damage rescaling** (`--ancient --damage-rescale moderate`): Downweights transitions
4. **Bayesian ancient mode** (`--ancient --bayesian`): Probabilistic scoring with adjusted error rates

**Recommended for ancient DNA:** Use `--ancient --bayesian` for best accuracy.

## Understanding Results

### Major Lineage vs Exact Match

yallHap uses the YFull 2025 tree, while many ground truth datasets use older ISOGG nomenclature. This causes low "exact match" rates even when classifications are correct.

**Example:**
- Ground truth: `R1b1a1a2a1a2c` (ISOGG 2016)
- yallHap call: `R-L21` (YFull 2025)
- These represent the **same haplogroup** with different naming conventions

### Heuristic vs Bayesian Mode

| Mode | Use Case | AD Required | Modern Accuracy | Ancient Accuracy |
|------|----------|-------------|-----------------|------------------|
| Heuristic | General use | No | 99.8% (1KG) | 88.3% (AADR) |
| Bayesian | Research, ancient DNA | Optional | 99.8% (1KG) | **90.7%** (AADR) |

**Key finding:** For modern samples, both modes produce identical results. For ancient DNA, Bayesian mode provides +2.4 pp improvement overall, with the largest gains at 4-10% variant density (+12-24 pp). At <4% or >10% density, mode choice matters less.

### Variant Density Metric

Variant density is calculated as:

```
variant_density = (called_variants / total_variants_in_VCF) × 100%
```

This is computed directly from the chrY VCF and cached in `data/ancient/aadr_variant_density.json` for efficiency. It provides a more reliable stratification metric than coverage metadata, which is only available for a subset of AADR samples.

## Validation Scripts Reference

| Script | Purpose |
|--------|---------|
| `scripts/download_1kg_validation.py` | Downloads 1000 Genomes Phase 3 VCF and panel |
| `scripts/download_ancient_test_data.py` | Creates ancient DNA ground truth TSV |
| `scripts/download_gnomad_highcov.py` | Downloads gnomAD HGDP/1KG high-coverage VCF |
| `scripts/download_liftover_chains.py` | Downloads T2T liftover chain files |
| `scripts/run_benchmarks.py` | Runs comprehensive benchmark suite |
| `scripts/validate_aadr_stratified.py` | AADR validation with variant density stratification |
| `scripts/gather_validation_and_comparative_data.py` | Multi-tool comparison on ancient samples |
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

Use Bayesian ancient mode for best results:

```python
classifier = HaplogroupClassifier(
    tree=tree,
    snp_db=snp_db,
    reference="grch37",
    bayesian=True,
    ancient_mode=True,
)
```

Or via CLI:

```bash
yallhap classify --vcf ancient.vcf.gz --ancient --bayesian
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

# 2. Run full benchmarks with 12 threads
python scripts/run_benchmarks.py --threads 12

# 3. Run AADR stratified validation
python scripts/validate_aadr_stratified.py \
  --samples-per-bin 500 \
  --threads 12 \
  -o results/aadr_density_stratified.md

# 4. Check results match paper
cat results/benchmark_results.json
```

Expected file outputs:
- `results/benchmark_results.json` - Full machine-readable results
- `results/aadr_density_stratified.md` - Stratified AADR results matching Table 3
