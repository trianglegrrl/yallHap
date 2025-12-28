# Reference Data for Y-chromosome Haplogroup Validation

This document describes the reference datasets used for yallHap validation and how to obtain them.

## Summary of Validation Datasets

| Dataset | Males | Coverage | Reference | Haplogroup Calls | yallHap Accuracy |
|---------|-------|----------|-----------|------------------|------------------|
| 1KG Phase 3 | 1,233 | ~4x | GRCh37 | Poznik 2016 | **99.76%** |
| AADR v54 | 1,991 | Variable | GRCh37 | AADR annotations | **88.65%** |
| gnomAD HGDP/1KG | 10* | 30x | GRCh38 | Poznik 2016 cross-ref | **100.00%** |

*gnomAD validation uses 10 samples balanced by superpopulation that overlap with 1KG Phase 3.

## 1000 Genomes Phase 3 (Primary Validation)

The 1000 Genomes Project provides the gold-standard validation dataset with published ground-truth haplogroup calls from Poznik et al. 2016.

### Phase 3 chrY VCF

**Primary download (EBI FTP):**
```bash
# Phase 3 chrY VCF and index (~13 MB compressed)
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi

# Sample panel file (maps sample IDs to populations)
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

**NCBI mirror (alternative):**
```bash
wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz
wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi
```

Key specifications: **1,244 males** from 26 populations, **~65,000 biallelic SNVs**, GRCh37 coordinates, Ensembl chromosome naming (uses "Y" not "chrY"), haploid genotype format. The .tbi tabix index enables random access.

### Poznik 2016 Ground-Truth Haplogroups

The seminal haplogroup reference for 1000 Genomes males comes from Poznik et al. 2016 (*Nature Genetics*), which assigned haplogroups to all **1,244 Phase 3 males** using ISOGG nomenclature.

**Paper and supplementary data:**
- **Full citation:** Poznik GD et al. "Punctuated bursts in human male demography inferred from 1,244 worldwide Y-chromosome sequences." *Nature Genetics* 48, 593–599 (2016)
- **DOI:** https://doi.org/10.1038/ng.3559
- **Supplementary Data:** Download from https://www.nature.com/articles/ng.3559#Sec14

**yhaplo software and reference data:**
```bash
# Clone the 23andMe yhaplo repository (contains ISOGG SNP data + example 1KG data)
git clone https://github.com/23andMe/yhaplo.git
cd yhaplo
pip install --editable .[vcf]

# Run on built-in 1000 Genomes example data
yhaplo --example_text
```

**Important note on naming convention:** Poznik 2016 uses the **ISOGG tree from January 2016**. The ISOGG tree has been updated substantially since then—yallHap uses the YFull 2024 tree with different nomenclature. This explains the low "exact match" rate (8.84%) despite high major lineage accuracy (99.76%).

### Automated Download

```bash
python scripts/download_1kg_validation.py
```

This creates:
- `data/validation/1kg_chrY_phase3.vcf.gz` - Filtered to males only
- `data/validation/poznik2016_haplogroups.tsv` - Ground truth mappings
- `data/validation/ybrowse_snps_hg19.csv` - GRCh37 SNP positions

## Allen Ancient DNA Resource (AADR)

The AADR provides ancient DNA samples with haplogroup annotations for archaeogenetic validation.

### Dataset Details

- **Version:** v54 (2024)
- **Total samples:** ~15,000+ (Y-chr annotations for subset)
- **Reference:** GRCh37
- **Format:** EIGENSTRAT (converted to VCF for yallHap)

### Ground Truth Quality

AADR ground truth uses mixed nomenclature:
- Some samples have proper haplogroup names (e.g., "R-L21") ✓
- Others have SNP-only names (e.g., "M458") - cannot compare directly

The benchmark script filters to samples with proper X-YYYY format, yielding **1,991 samples** for validation.

### Automated Download

```bash
python scripts/download_ancient_test_data.py
```

This creates:
- `data/ancient/aadr_chrY_v2.vcf.gz` - Y-chromosome variants
- `data/ancient/aadr_1240k_ground_truth.tsv` - Haplogroup annotations

### Ancient DNA Considerations

Ancient samples require special handling due to post-mortem damage:

1. **Transversions-only mode** (recommended): Ignores all C>T and G>A transitions
2. **Lower coverage**: Mean 11.8 derived SNPs vs 15.4 for modern samples
3. **Variable quality**: Some samples have very few callable positions

## gnomAD HGDP/1000G High-Coverage

The gnomAD v3.1 release includes high-coverage resequencing of HGDP and 1000 Genomes samples, enabling validation on GRCh38.

### Dataset Details

- **URL:** `https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz`
- **Reference:** GRCh38
- **Total samples:** 4,151 (HGDP + 1000 Genomes)
- **Coverage:** ~30x
- **Key feature:** AD (allelic depth) fields

### Sample Selection for Validation

The benchmark uses **10 samples** (configurable) that:
1. Exist in both gnomAD VCF AND 1000 Genomes Phase 3 ground truth
2. Are balanced by superpopulation (AFR, AMR, EAS, EUR, SAS - 2 each)

This allows cross-reference validation using Poznik 2016 haplogroup assignments.

### Automated Download

```bash
# Full VCF (~9 GB)
python scripts/download_gnomad_highcov.py

# Sample-only subset (faster)
python scripts/download_gnomad_highcov.py --sample-only
```

This creates:
- `data/validation_highcov/vcf/gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz`
- `data/validation_highcov/vcf/gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz.tbi`

The benchmark script automatically creates filtered subsets for efficiency:
- `data/validation_highcov/vcf/gnomad_subset_10_filtered.vcf.gz`
- `data/validation_highcov/vcf/diagnostic_positions.tsv`

### Benefits of gnomAD Validation

- **GRCh38 reference**: Validates coordinate mapping
- **High coverage**: More derived SNPs (mean 79.4 vs 15.4)
- **AD fields**: Enables Bayesian allelic depth integration
- **100% accuracy**: Perfect classification on all 10 samples

## Reference Phylogeny and SNP Database

### YFull Tree

The YFull consortium maintains the most comprehensive Y-chromosome phylogeny with monthly updates.

**Download:**
```bash
# Via yallHap CLI
yallhap download --output-dir data/

# Or manually
python scripts/download_yfull_tree.py
```

**Specifications:**
- **SNPs:** 185,780+
- **Subclades:** 15,000+
- **Format:** JSON
- **Size:** ~14 MB

### YBrowse SNP Database

Thomas Krahn (YSEQ) maintains the YBrowse database with SNP coordinates in multiple references.

**Files:**
- `data/ybrowse_snps.csv` - GRCh38 positions (~400 MB)
- `data/validation/ybrowse_snps_hg19.csv` - GRCh37 positions (~48 MB)

**Format:** GFF-style CSV with columns for position, ancestral/derived alleles, and aliases.

### T2T Liftover Chains

For T2T-CHM13v2.0 support, positions are computed via liftover:

```bash
python scripts/download_liftover_chains.py
```

This creates:
- `data/liftover/grch38-chm13v2.chain`
- `data/liftover/hg19-chm13v2.chain`

## Extracting chrY from Whole-Genome Files

If you need to extract Y-chromosome data from whole-genome VCFs:

```bash
# Extract chrY from Phase 3 (Ensembl naming - no "chr" prefix)
bcftools view -r Y input.vcf.gz -Oz -o chrY.vcf.gz
bcftools index -t chrY.vcf.gz

# Extract chrY from GRCh38 (UCSC naming)
bcftools view -r chrY input.vcf.gz -Oz -o chrY.vcf.gz

# Remote streaming without full download
bcftools view "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz" -Oz -o local_chrY.vcf.gz

# Extract specific samples (males only from Phase 3)
bcftools view -S male_samples.txt ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz -Oz -o males_only.vcf.gz
```

Required tools: **bcftools** (for VCF manipulation), **tabix** (for indexing), **htslib** (dependency).

## Complete Download Workflow

For comprehensive validation, download all datasets:

```bash
#!/bin/bash
# Complete download script for yallHap validation

# Create working directory
cd /path/to/yallhap

# 1. Reference data (tree + SNP database)
yallhap download --output-dir data/

# 2. 1000 Genomes Phase 3 (GRCh37, Poznik 2016 ground truth)
python scripts/download_1kg_validation.py

# 3. AADR Ancient DNA (GRCh37)
python scripts/download_ancient_test_data.py

# 4. gnomAD High-Coverage (GRCh38)
python scripts/download_gnomad_highcov.py

# 5. T2T liftover chains (optional)
python scripts/download_liftover_chains.py

echo "Download complete!"
```

**Expected file sizes:**
- YFull tree: ~14 MB
- YBrowse SNPs: ~400 MB (GRCh38) + ~48 MB (GRCh37)
- Phase 3 chrY VCF: ~5 MB
- AADR VCF: ~97 MB
- gnomAD chrY VCF: ~9 GB
- Total: **~10 GB**

## Key Considerations for Validation

### Coordinate Systems

| Dataset | Reference | Naming | SNP Database |
|---------|-----------|--------|--------------|
| 1KG Phase 3 | GRCh37 | "Y" | ybrowse_snps_hg19.csv |
| AADR | GRCh37 | "Y" | ybrowse_snps_hg19.csv |
| gnomAD | GRCh38 | "chrY" | ybrowse_snps.csv |

### Tree Version Compatibility

The Y-chromosome phylogeny evolves continuously:
- Poznik 2016 used ISOGG January 2016
- yallHap uses YFull 2024 (185K+ SNPs)
- This explains low "exact match" rates (8.84%) with high accuracy (99.76%)

### No Authentication Required

All datasets are fully open access:
- 1000 Genomes: Fort Lauderdale principles
- AADR: Reich Lab public data
- gnomAD: Broad Institute open access
- YFull/YBrowse: Public access

## Reproducing Benchmark Results

To exactly reproduce the published benchmark results:

```bash
# 1. Download all data
python scripts/download_1kg_validation.py
python scripts/download_ancient_test_data.py
python scripts/download_gnomad_highcov.py

# 2. Run benchmarks with 16 threads
python scripts/run_benchmarks.py --threads 16

# 3. Compare results
cat results/benchmark_results.json
```

Expected results:
| Dataset | Accuracy | Samples |
|---------|----------|---------|
| 1KG Phase 3 | 99.76% | 1,233 |
| AADR v54 | 88.65% | 1,991 |
| gnomAD | 100.00% | 10 |
