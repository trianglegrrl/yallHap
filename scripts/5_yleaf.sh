#!/bin/bash
# Run Yleaf on ancient DNA samples for tool comparison
# Usage: ./5_yleaf.sh

set -uo pipefail

# Paths
YLEAF="/home/a/.local/bin/Yleaf"
VCF_DIR="/20tb/yallHapTesting/vcf"
OUTPUT_DIR="/20tb/yallHapTesting/yleaf_results"
FILTERED_DIR="/20tb/yallHapTesting/vcf_filtered"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FILTERED_DIR"

echo "Running Yleaf on ancient DNA samples..."
echo "================================================"

# Function to filter VCF to only variants (remove reference-only calls)
filter_vcf() {
    local input=$1
    local output=$2
    echo "  Filtering $input -> $output"
    bcftools view -v snps "$input" -Oz -o "$output"
    bcftools index -t "$output"
}

# Common options for ancient DNA:
# -rg hg19: reference genome (GRCh37/hg19)
# -aDNA: ancient DNA mode (ignores C>T and G>A transitions)
# -r 1: minimum reads threshold (relaxed for ancient DNA)
# -q 10: minimum quality (relaxed for ancient DNA)
# -pq 0.5: prediction quality threshold (relaxed for low coverage)
# -force: overwrite existing files without prompting

# I0231 (Yamnaya, Bronze Age)
echo "Processing I0231..."
filter_vcf "$VCF_DIR/I0231.390k.chrY.vcf.gz" "$FILTERED_DIR/I0231.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/I0231.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/I0231"

# I0443 (Yamnaya, Bronze Age)
echo "Processing I0443..."
filter_vcf "$VCF_DIR/I0443.390k.chrY.vcf.gz" "$FILTERED_DIR/I0443.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/I0443.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/I0443"

# Kennewick Man
echo "Processing Kennewick..."
filter_vcf "$VCF_DIR/Kennewick_defaultMap1extr.realign.md.head.rmdup.chrY.vcf.gz" "$FILTERED_DIR/Kennewick.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/Kennewick.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/Kennewick"

# SB524A (British Neolithic)
echo "Processing SB524A..."
filter_vcf "$VCF_DIR/SB524A_lib.merged.markdup.chrY.vcf.gz" "$FILTERED_DIR/SB524A.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/SB524A.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/SB524A"

# SB524A2 (British Neolithic)
echo "Processing SB524A2..."
filter_vcf "$VCF_DIR/SB524A2_lib.merged.markdup.chrY.vcf.gz" "$FILTERED_DIR/SB524A2.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/SB524A2.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/SB524A2"

# VK287 (Viking Age)
echo "Processing VK287..."
filter_vcf "$VCF_DIR/VK287.final.chrY.vcf.gz" "$FILTERED_DIR/VK287.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/VK287.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/VK287"

# VK292 (Viking Age)
echo "Processing VK292..."
filter_vcf "$VCF_DIR/VK292.final.chrY.vcf.gz" "$FILTERED_DIR/VK292.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/VK292.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/VK292"

# VK296 (Viking Age)
echo "Processing VK296..."
filter_vcf "$VCF_DIR/VK296.final.chrY.vcf.gz" "$FILTERED_DIR/VK296.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/VK296.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/VK296"

# VK582 (Viking Age)
echo "Processing VK582..."
filter_vcf "$VCF_DIR/VK582.final.chrY.vcf.gz" "$FILTERED_DIR/VK582.vcf.gz"
$YLEAF -vcf "$FILTERED_DIR/VK582.vcf.gz" \
    -rg hg19 -aDNA -r 1 -q 10 -pq 0.5 -force \
    -o "$OUTPUT_DIR/VK582"

echo "================================================"
echo "Done! Results in: $OUTPUT_DIR"
echo ""
echo "Summary:"
echo "--------"

# Print summary - Yleaf outputs results to *_hg_prediction.txt files
for dir in "$OUTPUT_DIR"/*/; do
    sample=$(basename "$dir")
    pred_file=$(find "$dir" -name "*_hg_prediction.txt" 2>/dev/null | head -1)
    if [ -n "$pred_file" ] && [ -f "$pred_file" ]; then
        # Get haplogroup from prediction file (usually second line, first column)
        hg=$(tail -1 "$pred_file" | cut -f1)
        qc=$(tail -1 "$pred_file" | cut -f2)
        echo "$sample: $hg (QC: $qc)"
    else
        echo "$sample: No prediction file found"
    fi
done

