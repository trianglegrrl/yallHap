#!/bin/bash
# Comprehensive ancient DNA haplogroup comparison: yallHap vs Yleaf vs pathPhynder
# Usage: ./6_ancient_comparison.sh
#
# This script runs yallHap in multiple modes on ancient samples and compares
# results with Yleaf and pathPhynder outputs.

set -uo pipefail

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR="/20tb/yallHapTesting"
VCF_DIR="$BASE_DIR/vcf"
OUTPUT_DIR="$BASE_DIR/yallhap_comparison"
YLEAF_DIR="$BASE_DIR/yleaf_results"
PATHPHYNDER_DIR="$BASE_DIR/pathPhynder"

# yallHap paths (adjust if needed)
YALLHAP="/home/a/.local/bin/yallhap"
TREE="$BASE_DIR/data/yfull_tree.json"
SNP_DB="$BASE_DIR/data/ybrowse_snps_hg19.csv"
ISOGG_DB="$PATHPHYNDER_DIR/data/211108.snps_isogg_curated.b38.txt"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "Ancient DNA Haplogroup Comparison: yallHap vs Yleaf vs pathPhynder"
echo "============================================================"
echo ""

# ============================================================================
# Sample definitions: (sample_name, vcf_file, vcf_sample_id)
# ============================================================================
declare -A SAMPLES=(
    ["I0231"]="I0231.390k.chrY.vcf.gz|SM"
    ["I0443"]="I0443.390k.chrY.vcf.gz|SM"
    ["Kennewick"]="Kennewick_defaultMap1extr.realign.md.head.rmdup.chrY.vcf.gz|Ken19"
    ["SB524A"]="SB524A_lib.merged.markdup.chrY.vcf.gz|Cheddar"
    ["SB524A2"]="SB524A2_lib.merged.markdup.chrY.vcf.gz|Cheddar"
    ["VK287"]="VK287.final.chrY.vcf.gz|VK287"
    ["VK292"]="VK292.final.chrY.vcf.gz|VK292"
    ["VK296"]="VK296.final.chrY.vcf.gz|VK296"
    ["VK582"]="VK582.final.chrY.vcf.gz|IA_PD_22"
)

# ============================================================================
# Run yallHap in multiple modes
# ============================================================================

echo "Running yallHap classifications..."
echo "--------------------------------"

# Create output files
YALLHAP_REGULAR="$OUTPUT_DIR/yallhap_regular.tsv"
YALLHAP_ISOGG="$OUTPUT_DIR/yallhap_isogg.tsv"
YALLHAP_BAYESIAN="$OUTPUT_DIR/yallhap_bayesian.tsv"
YALLHAP_TRANSVERSIONS="$OUTPUT_DIR/yallhap_transversions.tsv"

# Header for TSV files
echo -e "sample\thaplogroup\tconfidence\tderived\tancestral\tnocall" > "$YALLHAP_REGULAR"
echo -e "sample\thaplogroup\tconfidence\tderived\tancestral\tnocall\tisogg_haplogroup" > "$YALLHAP_ISOGG"
echo -e "sample\thaplogroup\tconfidence\tderived\tancestral\tnocall" > "$YALLHAP_BAYESIAN"
echo -e "sample\thaplogroup\tconfidence\tderived\tancestral\tnocall" > "$YALLHAP_TRANSVERSIONS"

for sample in "${!SAMPLES[@]}"; do
    IFS='|' read -r vcf_file sample_id <<< "${SAMPLES[$sample]}"
    vcf_path="$VCF_DIR/$vcf_file"

    if [ ! -f "$vcf_path" ]; then
        echo "  WARNING: VCF not found: $vcf_path"
        continue
    fi

    echo "  Processing $sample..."

    # Mode 1: Regular ancient mode
    echo "    - Regular ancient mode"
    $YALLHAP classify "$vcf_path" \
        --tree "$TREE" \
        --snp-db "$SNP_DB" \
        --sample "$sample_id" \
        --ancient \
        --reference grch37 \
        --format tsv \
        --min-depth 1 \
        --min-quality 10 \
        2>/dev/null | tail -1 | awk -v s="$sample" '{print s"\t"$0}' >> "$YALLHAP_REGULAR"

    # Mode 2: Ancient + ISOGG
    echo "    - Ancient + ISOGG mode"
    if [ -f "$ISOGG_DB" ]; then
        $YALLHAP classify "$vcf_path" \
            --tree "$TREE" \
            --snp-db "$SNP_DB" \
            --sample "$sample_id" \
            --ancient \
            --isogg \
            --isogg-db "$ISOGG_DB" \
            --reference grch37 \
            --format tsv \
            --min-depth 1 \
            --min-quality 10 \
            2>/dev/null | tail -1 | awk -v s="$sample" '{print s"\t"$0}' >> "$YALLHAP_ISOGG"
    else
        echo "      ISOGG DB not found at $ISOGG_DB"
    fi

    # Mode 3: Bayesian ancient mode
    echo "    - Bayesian ancient mode"
    $YALLHAP classify "$vcf_path" \
        --tree "$TREE" \
        --snp-db "$SNP_DB" \
        --sample "$sample_id" \
        --ancient \
        --bayesian \
        --reference grch37 \
        --format tsv \
        --min-depth 1 \
        --min-quality 10 \
        2>/dev/null | tail -1 | awk -v s="$sample" '{print s"\t"$0}' >> "$YALLHAP_BAYESIAN"

    # Mode 4: Transversions-only (strictest ancient mode)
    echo "    - Transversions-only mode"
    $YALLHAP classify "$vcf_path" \
        --tree "$TREE" \
        --snp-db "$SNP_DB" \
        --sample "$sample_id" \
        --transversions-only \
        --reference grch37 \
        --format tsv \
        --min-depth 1 \
        --min-quality 10 \
        2>/dev/null | tail -1 | awk -v s="$sample" '{print s"\t"$0}' >> "$YALLHAP_TRANSVERSIONS"
done

echo ""

# ============================================================================
# Collect Yleaf results
# ============================================================================

echo "Collecting Yleaf results..."
echo "---------------------------"

YLEAF_RESULTS="$OUTPUT_DIR/yleaf_results.tsv"
echo -e "sample\thaplogroup\tqc_score" > "$YLEAF_RESULTS"

for sample in "${!SAMPLES[@]}"; do
    pred_file=$(find "$YLEAF_DIR/$sample" -name "*_hg_prediction.txt" 2>/dev/null | head -1)
    if [ -n "$pred_file" ] && [ -f "$pred_file" ]; then
        hg=$(tail -1 "$pred_file" | cut -f1)
        qc=$(tail -1 "$pred_file" | cut -f2)
        echo -e "$sample\t$hg\t$qc" >> "$YLEAF_RESULTS"
        echo "  $sample: $hg (QC: $qc)"
    else
        echo -e "$sample\tNA\tNA" >> "$YLEAF_RESULTS"
        echo "  $sample: No Yleaf result found"
    fi
done

echo ""

# ============================================================================
# Collect pathPhynder results (if available)
# ============================================================================

echo "Collecting pathPhynder results..."
echo "---------------------------------"

PATHPHYNDER_RESULTS="$OUTPUT_DIR/pathphynder_results.tsv"
echo -e "sample\thaplogroup\tisogg_haplogroup\tlikelihood" > "$PATHPHYNDER_RESULTS"

# Look for pathPhynder output files
for sample in "${!SAMPLES[@]}"; do
    # pathPhynder typically outputs to sample.results or similar
    result_file=$(find "$PATHPHYNDER_DIR" -name "*${sample}*" -type f 2>/dev/null | head -1)
    if [ -n "$result_file" ] && [ -f "$result_file" ]; then
        # Try to parse pathPhynder results (format may vary)
        echo "  $sample: Found at $result_file"
        # Add parsing logic based on pathPhynder output format
    else
        echo -e "$sample\tNA\tNA\tNA" >> "$PATHPHYNDER_RESULTS"
        echo "  $sample: No pathPhynder result found"
    fi
done

echo ""

# ============================================================================
# Generate summary table
# ============================================================================

echo "Generating summary..."
echo "--------------------"

SUMMARY_MD="$OUTPUT_DIR/99_output_summary.md"

cat > "$SUMMARY_MD" << 'EOF'
# Ancient DNA Haplogroup Classification Comparison

## Overview

This document summarizes Y-chromosome haplogroup classifications for ancient DNA samples
using multiple tools and modes:

- **yallHap Regular**: Standard ancient DNA mode (filters C>T, G>A transitions)
- **yallHap Bayesian**: Bayesian classification with ancient mode
- **yallHap Transversions**: Strictest mode using only transversions
- **yallHap ISOGG**: Ancient mode with ISOGG nomenclature output
- **Yleaf**: Reference tool with ancient DNA mode
- **pathPhynder**: Bayesian phylogenetic placement

## Sample Information

| Sample | Description | Period |
|--------|-------------|--------|
| I0231 | Yamnaya | Bronze Age |
| I0443 | Yamnaya | Bronze Age |
| Kennewick | Kennewick Man | Paleoamerican |
| SB524A | Cheddar Man (lib 1) | British Mesolithic |
| SB524A2 | Cheddar Man (lib 2) | British Mesolithic |
| VK287 | Viking | Viking Age |
| VK292 | Viking | Viking Age |
| VK296 | Viking | Viking Age |
| VK582 | Viking | Viking Age |

## Classification Results

EOF

# Build the comparison table
echo "" >> "$SUMMARY_MD"
echo "| Sample | yallHap Regular | yallHap Bayesian | yallHap Transversions | yallHap ISOGG | Yleaf | pathPhynder |" >> "$SUMMARY_MD"
echo "|--------|-----------------|------------------|----------------------|---------------|-------|-------------|" >> "$SUMMARY_MD"

for sample in I0231 I0443 Kennewick SB524A SB524A2 VK287 VK292 VK296 VK582; do
    # Get yallHap regular result
    yallhap_reg=$(grep "^$sample" "$YALLHAP_REGULAR" 2>/dev/null | cut -f2 || echo "NA")
    [ -z "$yallhap_reg" ] && yallhap_reg="NA"

    # Get yallHap Bayesian result
    yallhap_bayes=$(grep "^$sample" "$YALLHAP_BAYESIAN" 2>/dev/null | cut -f2 || echo "NA")
    [ -z "$yallhap_bayes" ] && yallhap_bayes="NA"

    # Get yallHap Transversions result
    yallhap_tv=$(grep "^$sample" "$YALLHAP_TRANSVERSIONS" 2>/dev/null | cut -f2 || echo "NA")
    [ -z "$yallhap_tv" ] && yallhap_tv="NA"

    # Get yallHap ISOGG result
    yallhap_isogg=$(grep "^$sample" "$YALLHAP_ISOGG" 2>/dev/null | awk '{print $2" ("$NF")"}' || echo "NA")
    [ -z "$yallhap_isogg" ] && yallhap_isogg="NA"

    # Get Yleaf result
    yleaf=$(grep "^$sample" "$YLEAF_RESULTS" 2>/dev/null | cut -f2 || echo "NA")
    [ -z "$yleaf" ] && yleaf="NA"

    # Get pathPhynder result
    pathphynder=$(grep "^$sample" "$PATHPHYNDER_RESULTS" 2>/dev/null | cut -f2 || echo "NA")
    [ -z "$pathphynder" ] && pathphynder="NA"

    echo "| $sample | $yallhap_reg | $yallhap_bayes | $yallhap_tv | $yallhap_isogg | $yleaf | $pathphynder |" >> "$SUMMARY_MD"
done

cat >> "$SUMMARY_MD" << 'EOF'

## Notes

- **Coverage**: Ancient samples typically have very low coverage (<1x)
- **Damage**: C>T and G>A transitions filtered due to deamination damage
- **Confidence**: Higher confidence values indicate more reliable calls
- **Derived SNPs**: Number of haplogroup-defining SNPs observed in derived state

## Method Details

### yallHap Modes

1. **Regular ancient mode** (`--ancient`): Filters damage-like transitions (C>T, G>A)
2. **Bayesian mode** (`--bayesian --ancient`): Uses allelic depth for probabilistic classification
3. **Transversions-only** (`--transversions-only`): Strictest mode, only uses transversions
4. **ISOGG mode** (`--isogg`): Adds ISOGG nomenclature to output

### Yleaf Settings

- `-aDNA`: Ancient DNA mode
- `-r 1`: Minimum 1 read
- `-q 10`: Minimum quality 10
- `-pq 0.5`: Prediction quality threshold 0.5

---

*Generated by 6_ancient_comparison.sh*
EOF

echo ""
echo "============================================================"
echo "Done! Results saved to:"
echo "  - $YALLHAP_REGULAR"
echo "  - $YALLHAP_BAYESIAN"
echo "  - $YALLHAP_TRANSVERSIONS"
echo "  - $YALLHAP_ISOGG"
echo "  - $YLEAF_RESULTS"
echo "  - $SUMMARY_MD"
echo "============================================================"
echo ""
echo "View summary: cat $SUMMARY_MD"

