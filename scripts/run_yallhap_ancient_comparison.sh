#!/bin/bash
# Run yallHap on ancient DNA samples for pathPhynder comparison
# Usage: ./run_yallhap_ancient_comparison.sh

# Note: Don't use -e because yallHap exits with code 2 for low-confidence calls
set -uo pipefail

# Paths
YALLHAP="/home/a/.local/bin/yallhap"
TREE="/references/human-y/yfull_tree.json"
SNPDB="/references/human-y/ybrowse_snps_hg19.csv"
VCF_DIR="/20tb/yallHapTesting/vcf"
OUTPUT_DIR="/20tb/yallHapTesting/yallhap_results"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Common options
# Ancient DNA needs relaxed depth/quality thresholds since coverage is typically <1x
OPTS="--tree $TREE --snp-db $SNPDB --reference grch37 --ancient --min-depth 1 --min-quality 0 --format json"

echo "Running yallHap on ancient DNA samples..."
echo "================================================"

# I0231 (Yamnaya, Bronze Age)
echo "Processing I0231..."
$YALLHAP classify "$VCF_DIR/I0231.390k.chrY.vcf.gz" \
    $OPTS --sample SM \
    -o "$OUTPUT_DIR/I0231.json"

# I0443 (Yamnaya, Bronze Age)
echo "Processing I0443..."
$YALLHAP classify "$VCF_DIR/I0443.390k.chrY.vcf.gz" \
    $OPTS --sample SM \
    -o "$OUTPUT_DIR/I0443.json"

# Kennewick Man
echo "Processing Kennewick..."
$YALLHAP classify "$VCF_DIR/Kennewick_defaultMap1extr.realign.md.head.rmdup.chrY.vcf.gz" \
    $OPTS --sample Ken19 \
    -o "$OUTPUT_DIR/Kennewick.json"

# SB524A (British Neolithic - Cheddar Man related)
echo "Processing SB524A..."
$YALLHAP classify "$VCF_DIR/SB524A_lib.merged.markdup.chrY.vcf.gz" \
    $OPTS --sample Cheddar \
    -o "$OUTPUT_DIR/SB524A.json"

# SB524A2 (British Neolithic - Cheddar Man related)
echo "Processing SB524A2..."
$YALLHAP classify "$VCF_DIR/SB524A2_lib.merged.markdup.chrY.vcf.gz" \
    $OPTS --sample Cheddar \
    -o "$OUTPUT_DIR/SB524A2.json"

# VK287 (Viking Age)
echo "Processing VK287..."
$YALLHAP classify "$VCF_DIR/VK287.final.chrY.vcf.gz" \
    $OPTS --sample VK287 \
    -o "$OUTPUT_DIR/VK287.json"

# VK292 (Viking Age)
echo "Processing VK292..."
$YALLHAP classify "$VCF_DIR/VK292.final.chrY.vcf.gz" \
    $OPTS --sample VK292 \
    -o "$OUTPUT_DIR/VK292.json"

# VK296 (Viking Age)
echo "Processing VK296..."
$YALLHAP classify "$VCF_DIR/VK296.final.chrY.vcf.gz" \
    $OPTS --sample VK296 \
    -o "$OUTPUT_DIR/VK296.json"

# VK582 (Viking Age)
echo "Processing VK582..."
$YALLHAP classify "$VCF_DIR/VK582.final.chrY.vcf.gz" \
    $OPTS --sample IA_PD_22 \
    -o "$OUTPUT_DIR/VK582.json"

echo "================================================"
echo "Done! Results in: $OUTPUT_DIR"
echo ""
echo "Summary:"
echo "--------"

# Print summary table
printf "%-12s %-25s %-10s %-10s\n" "Sample" "Haplogroup" "Confidence" "Derived"
printf "%-12s %-25s %-10s %-10s\n" "------" "----------" "----------" "-------"

for json in "$OUTPUT_DIR"/*.json; do
    sample=$(basename "$json" .json)
    if [ -f "$json" ]; then
        hg=$(jq -r '.haplogroup // "N/A"' "$json" 2>/dev/null || echo "N/A")
        conf=$(jq -r '.confidence // "N/A"' "$json" 2>/dev/null || echo "N/A")
        derived=$(jq -r '.derived_count // "N/A"' "$json" 2>/dev/null || echo "N/A")
        printf "%-12s %-25s %-10s %-10s\n" "$sample" "$hg" "$conf" "$derived"
    fi
done

