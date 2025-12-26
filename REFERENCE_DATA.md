# Obtaining 1000 Genomes Y-chromosome data for haplogroup caller validation

The 1000 Genomes Project provides **two distinct Y-chromosome datasets**: Phase 3 (low-coverage, GRCh37, ~65,000 variants from 1,244 males) and the NYGC high-coverage 30x resequencing (GRCh38, 3,202 samples). For haplogroup caller development, Phase 3 remains the gold standard with published ground-truth haplogroup calls from Poznik et al. 2016, while 30x data offers higher variant quality for novel discovery. Both are fully open access with no authentication required.

## Phase 3 chrY VCF: the primary validation dataset

The Phase 3 release includes a **dedicated chrY VCF file** separate from autosomes, specifically curated for haplogroup studies. This is the dataset underlying the Poznik 2016 haplogroup assignments.

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

## NYGC high-coverage 30x data for enhanced variant discovery

The 2020-2022 NYGC resequencing provides substantially higher quality variant calls. Y-chromosome VCFs are available per-chromosome.

**Raw genotype VCF (chrY, 30x):**
```bash
# Primary location - raw GT with annotations
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz.tbi
```

**AWS S3 access (faster for large downloads):**
```bash
# List available files
aws s3 ls s3://1000genomes/1000G_2504_high_coverage/ --no-sign-request

# Download chrY from AWS
aws s3 cp s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz . --no-sign-request
```

Key differences from Phase 3: **GRCh38 coordinates**, UCSC naming ("chrY"), 3,202 total samples including 698 related individuals completing 602 trios, **30x Illumina NovaSeq coverage**.

## Poznik 2016 ground-truth haplogroup assignments

The seminal haplogroup reference for 1000 Genomes males comes from Poznik et al. 2016 (*Nature Genetics*), which assigned haplogroups to all **1,244 Phase 3 males** using ISOGG nomenclature.

**Paper and supplementary data:**
- **Full citation:** Poznik GD et al. "Punctuated bursts in human male demography inferred from 1,244 worldwide Y-chromosome sequences." *Nature Genetics* 48, 593–599 (2016)
- **DOI:** https://doi.org/10.1038/ng.3559
- **Supplementary Data:** Download from https://www.nature.com/articles/ng.3559#Sec14
  - Supplementary Data ZIP (6.2 MB) contains sample-to-haplogroup mappings

**yhaplo software and reference data:**
```bash
# Clone the 23andMe yhaplo repository (contains ISOGG SNP data + example 1KG data)
git clone https://github.com/23andMe/yhaplo.git
cd yhaplo
pip install --editable .[vcf]

# Run on built-in 1000 Genomes example data
yhaplo --example_text
```

The yhaplo output format provides 4 columns: sample ID, haplogroup with observed SNP (e.g., "R1b-M343"), haplogroup with representative SNP, and long-form Y Chromosome Consortium nomenclature. The repository includes ISOGG SNP positions in `yhaplo/data/variants/` for GRCh37.

**Important note on naming convention:** Poznik 2016 uses the **ISOGG tree from January 2016**. The ISOGG tree has been updated substantially since then—current tools like Y-LineageTracker use ISOGG 2019, and Yleaf v3+ uses YFull v10.01.

## Extracting chrY from whole-genome files with bcftools

If you need to extract Y-chromosome data from whole-genome VCFs or process files remotely:

```bash
# Extract chrY from Phase 3 (Ensembl naming - no "chr" prefix)
bcftools view -r Y input.vcf.gz -Oz -o chrY.vcf.gz
bcftools index -t chrY.vcf.gz

# Extract chrY from 30x GRCh38 (UCSC naming)
bcftools view -r chrY input.vcf.gz -Oz -o chrY.vcf.gz

# Remote streaming without full download
bcftools view "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz" -Oz -o local_chrY.vcf.gz

# Extract specific samples (males only from Phase 3)
bcftools view -S male_samples.txt ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz -Oz -o males_only.vcf.gz
```

Required tools: **bcftools** (for VCF manipulation), **tabix** (for indexing), **htslib** (dependency).

## Alternative datasets for comprehensive validation

Beyond 1000 Genomes, three major datasets provide additional Y-haplogroup benchmarking opportunities:

| Dataset | Males | Coverage | Reference | Haplogroup calls available |
|---------|-------|----------|-----------|---------------------------|
| 1KG Phase 3 | 1,244 | ~4x | GRCh37 | Yes (Poznik 2016) |
| 1KG 30x NYGC | ~1,500 | 30x | GRCh38 | Must derive from VCF |
| HGDP Bergström | 492 | ~36x | GRCh38 | Community spreadsheets only |
| SGDP | ~160 | ~30x | GRCh37/38 | Must derive from VCF |

**HGDP high-coverage (Bergström 2020):**
```bash
# HGDP data from IGSR
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/

# gnomAD harmonized HGDP+1KG callset (includes Y-chr)
gsutil -m cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/* ./
```

**SGDP from Reich Lab:**
```bash
# SGDP data access
# Primary: https://reichdata.hms.harvard.edu/pub/datasets/sgdp/
# Google Cloud: gs://genomics-public-data/simons-genome-diversity-project
```

**Y-LineageTracker validation set:** The Y-LineageTracker paper (Chen et al. 2021, BMC Bioinformatics) validated against **1,233 1KG males** using the ISOGG 2019 tree, providing updated haplogroup calls in supplementary materials.

## Benchmark datasets from Y-haplogroup tools

Several published tools provide validation datasets:

- **yhaplo** (23andMe): Built-in 1KG example data at `yhaplo --example_text`
- **Y-LineageTracker**: Supplementary tables with 1,233 sample haplogroup calls (ISOGG 2019)
- **Yleaf**: Table 1 in Ralf et al. 2018 lists 51 validation samples with YFull assignments
- **Y-mer** (2025): V110 validation set with 110 individuals across 11 basal haplogroups

A 2023 benchmarking study (PMC10560978) compared 5 tools using 50 custom WGS/WES donors plus 54 1KG samples—this paper provides the most systematic cross-tool comparison.

## Complete download workflow for tool development

For a comprehensive validation dataset, download Phase 3 chrY (gold-standard calls) plus 30x data (higher quality variants):

```bash
#!/bin/bash
# Complete download script for Y-haplogroup tool validation

# Create working directory
mkdir -p 1kg_y_chromosome && cd 1kg_y_chromosome

# --- Phase 3 (GRCh37, Poznik 2016 ground truth) ---
echo "Downloading Phase 3 chrY data..."
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# --- 30x High-Coverage (GRCh38) ---
echo "Downloading 30x chrY data..."
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz.tbi

# --- Pedigree/sample information ---
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt

# --- yhaplo tool and ISOGG reference ---
git clone https://github.com/23andMe/yhaplo.git
cd yhaplo && pip install --editable .[vcf] && cd ..

echo "Download complete. Total size: ~200-300 MB for VCFs"
```

**Expected file sizes:**
- Phase 3 chrY VCF: ~13 MB compressed
- 30x chrY VCF: ~50-100 MB compressed
- Total download: **~200-300 MB** for core chrY data

## Key considerations for haplogroup caller development

**Coordinate systems matter:** Phase 3 uses GRCh37 with Ensembl naming ("Y"), while 30x uses GRCh38 with UCSC naming ("chrY"). Most ISOGG SNP databases provide both coordinate systems—yhaplo uses GRCh37 positions.

**Ground truth availability:** Only Phase 3 has published ground-truth haplogroup assignments (Poznik 2016). For 30x data, you must derive haplogroups using existing tools like yhaplo or Y-LineageTracker, then validate against those calls.

**ISOGG tree versioning:** The Y-chromosome phylogeny evolves continuously. Poznik 2016 used the January 2016 ISOGG tree; current tools use 2019+ versions with substantially refined branching. Consider which tree version your tool will target.

**No authentication required:** All 1000 Genomes data is fully open access under Fort Lauderdale principles. gnomAD notably does **not** include Y-chromosome calls in its releases. UK Biobank requires a formal application process for access to its ~152,000 male Y-genotypes.