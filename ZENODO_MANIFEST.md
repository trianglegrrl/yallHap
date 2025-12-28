# yallHap Validation Data Bundle

**Version:** v4
**DOI:** [To be assigned upon Zenodo upload]

This archive contains all data required to reproduce yallHap validation benchmarks and results published in the accompanying paper.

## Quick Start

```bash
# Extract archive
tar -xzf yallhap-validation-v4.tar.gz

# Copy data to yallHap data directory
cp -r yfull_tree.json ybrowse_snps.csv isogg_snps_grch38.txt /path/to/yallhap/data/
cp -r validation ancient ancient_genomes liftover /path/to/yallhap/data/

# Run benchmarks
python scripts/run_benchmarks.py --threads 16
```

---

## File Manifest

### Core Reference Data

| File | Size | Description |
|------|------|-------------|
| `yfull_tree.json` | 13.9 MB | YFull phylogenetic tree (Dec 2025, commit `7a6dc24`) |
| `ybrowse_snps.csv` | 410.8 MB | YBrowse SNP database (GRCh38 coordinates) |
| `isogg_snps_grch38.txt` | 2.7 MB | ISOGG SNP database (GRCh38 coordinates) |

### 1000 Genomes Phase 3 Validation

Ground truth from Poznik et al. 2016, 1,233 male samples with Y haplogroup calls.

| File | Size | Description |
|------|------|-------------|
| `validation/1kg_chrY_phase3.vcf.gz` | 5.4 MB | 1000 Genomes Phase 3 Y-chromosome VCF (GRCh37) |
| `validation/1kg_chrY_phase3.vcf.gz.tbi` | 0.0 MB | Tabix index |
| `validation/poznik2016_haplogroups.tsv` | 0.0 MB | Ground truth haplogroups from Poznik 2016 |
| `validation/1kg_sample_panel.txt` | 0.1 MB | Sample population/superpopulation metadata |
| `validation/ybrowse_snps_hg19.csv` | 48.0 MB | YBrowse SNP database (GRCh37 coordinates) |

### AADR Ancient DNA Validation

Allen Ancient DNA Resource v54, 7,333 ancient samples with published haplogroup assignments.

| File | Size | Description |
|------|------|-------------|
| `ancient/aadr_chrY_v2.vcf.gz` | 90.8 MB | AADR v54 Y-chromosome multi-sample VCF (GRCh37) |
| `ancient/aadr_chrY_v2.vcf.gz.tbi` | 0.0 MB | Tabix index |
| `ancient/aadr_1240k_ground_truth.tsv` | 0.2 MB | Ground truth haplogroups from AADR annotations |

### gnomAD High-Coverage Validation

Pre-extracted subset of gnomAD v3.1 HGDP/1KG high-coverage data for samples overlapping 1000 Genomes.

| File | Size | Description |
|------|------|-------------|
| `validation_highcov/vcf/gnomad_1kg_shared_diagnostic.vcf.gz` | 2.22 GB | Diagnostic SNP positions for shared samples (GRCh38) |
| `validation_highcov/vcf/gnomad_1kg_shared_diagnostic.vcf.gz.tbi` | 0.0 MB | Tabix index |
| `validation_highcov/vcf/diagnostic_positions.tsv` | 35.3 MB | List of diagnostic SNP positions used |

### T2T Liftover Chain Files

Chain files for lifting coordinates to T2T-CHM13v2 reference genome.

| File | Size | Description |
|------|------|-------------|
| `liftover/grch38-chm13v2.chain` | 6.0 MB | GRCh38 → T2T-CHM13v2 liftover chain |
| `liftover/hg19-chm13v2.chain` | 6.0 MB | hg19/GRCh37 → T2T-CHM13v2 liftover chain |

### Ancient Genome Test Cases

9 published ancient genomes for cross-tool comparison (BAM + VCF pairs).

| Sample | BAM Size | VCF Size | Description |
|--------|----------|----------|-------------|
| `I0231` | 2.6 MB | 6.2 MB | Bronze Age sample (390k capture) |
| `I0443` | 1.5 MB | 3.8 MB | Bronze Age sample (390k capture) |
| `Kennewick` | 6.3 MB | 17.1 MB | Kennewick Man (whole genome) |
| `SB524A` | 1.5 MB | 3.8 MB | Srubnaya sample |
| `SB524A2` | 2.2 MB | 3.2 MB | Srubnaya sample (replicate) |
| `VK287` | 4.1 MB | 11.4 MB | Viking Age sample |
| `VK292` | 0.8 MB | 2.4 MB | Viking Age sample |
| `VK296` | 3.7 MB | 10.0 MB | Viking Age sample |
| `VK582` | 1.0 MB | 2.8 MB | Viking Age sample |

Each sample includes:
- `ancient_genomes/{sample}.bam` - Aligned reads (Y chromosome only)
- `ancient_genomes/{sample}.bam.bai` - BAM index
- `ancient_genomes/{sample}.vcf.gz` - Variant calls
- `ancient_genomes/{sample}.vcf.gz.tbi` - VCF index

### Benchmark Results

Pre-computed results from running validation scripts.

| File | Description |
|------|-------------|
| `results/benchmark_results.json` | Machine-readable benchmark results |
| `results/coalescent_prior_comparison.json` | Uniform vs coalescent prior comparison (200 samples) |
| `results/isogg_validation.json` | ISOGG mapping validation summary |
| `results/isogg_validation_full.json` | Detailed ISOGG validation per-sample results |
| `results/isogg_mismatch_analysis.json` | Categorized ISOGG mismatch analysis |
| `results/isogg_mismatches_detailed.json` | Complete list of 148 ISOGG mismatches with sample details |
| `results/power_analysis_ancient.json` | Power analysis by variant density |
| `results/significance_tests.json` | Two-proportion z-test p-values for Tables S3/S4 |
| `results/validation_report_*.md` | Validation reports in Markdown format |
| `results/aadr_*.md` | AADR stratified analysis reports |
| `results/isogg_*.md` | ISOGG validation reports |

### Validation Scripts

Scripts for reproducing validation analyses from scratch.

| File | Description |
|------|-------------|
| `scripts/download_1kg_validation.py` | Download 1000 Genomes Phase 3 chrY data |
| `scripts/download_ancient_test_data.py` | Download AADR ancient DNA data |
| `scripts/download_gnomad_highcov.py` | Download gnomAD high-coverage data |
| `scripts/download_liftover_chains.py` | Download T2T liftover chain files |
| `scripts/download_yfull_tree.py` | Download YFull phylogenetic tree |
| `scripts/run_benchmarks.py` | Run comprehensive benchmark suite |
| `scripts/validate_1kg.py` | Validate against 1000 Genomes |
| `scripts/validate_aadr_stratified.py` | Stratified ancient DNA validation by variant density |
| `scripts/gather_validation_and_comparative_data.py` | Comprehensive validation and tool comparison |
| `scripts/validate_isogg.py` | ISOGG nomenclature mapping validation |
| `scripts/analyze_isogg_mismatches.py` | Categorize and analyze ISOGG mismatches |
| `scripts/generate_isogg_mismatch_table.py` | Generate Supplementary Table S5 (mismatch details) |
| `scripts/power_analysis_ancient.py` | Power analysis by variant density threshold |
| `scripts/calculate_significance.py` | Two-proportion z-tests for statistical significance |
| `scripts/test_coalescent_priors.py` | Compare uniform vs coalescent prior performance |
| `scripts/generate_assets.py` | Generate publication figures |
| `scripts/prep_for_zenodo.py` | Script used to create this bundle |

---

## Data Provenance

### YFull Tree
- **Source:** https://github.com/YFullTeam/YTree
- **Commit:** `7a6dc243a5535e239f2400d1d308760069f7d103`
- **Date:** December 25, 2025

### YBrowse SNP Database
- **Source:** https://ybrowse.org/
- **Format:** CSV export with GRCh37/GRCh38 coordinates

### 1000 Genomes Phase 3
- **Source:** https://www.internationalgenome.org/
- **Reference:** Poznik et al. 2016 (DOI: 10.1038/ng.3527)

### AADR (Allen Ancient DNA Resource)
- **Source:** https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
- **Version:** v54.1

### gnomAD High-Coverage
- **Source:** https://gnomad.broadinstitute.org/
- **Version:** v3.1 HGDP/1KG subset

### Ancient Genomes
- Various published studies; see individual sample publications

---

## License

| Component | License |
|-----------|---------|
| **yallHap software** | MIT License |
| **Derived data/results** (benchmark outputs, validation reports) | CC-BY-4.0 |
| **Original consortium data** | See respective provider terms |

### Original Data Terms

- **1000 Genomes:** Fort Lauderdale principles
- **AADR:** Reich Lab data access policies
- **gnomAD:** gnomAD terms of use
- **YFull/YBrowse:** See provider terms

---

## Citation

If you use this data bundle, please cite:

1. The yallHap paper (bioRxiv preprint / publication)
2. This Zenodo record (DOI to be assigned)
3. Original data sources as appropriate

