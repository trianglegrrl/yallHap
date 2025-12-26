# Validation Testing

This document describes how to reproduce yallHap's validation results and run validation tests on new datasets.

## Overview

yallHap is validated against two primary datasets:

| Dataset | Samples | Reference | Type | Accuracy |
|---------|---------|-----------|------|----------|
| 1000 Genomes Phase 3 | 1,244 males | GRCh37 | Modern WGS | 100% same major lineage |
| AADR v54 | 50 samples | GRCh37 | Ancient DNA | 96% |

## Quick Start

### Option 1: Download from Zenodo (Recommended)

Download the pre-packaged validation data:

```bash
# Download validation bundle (~165 MB)
wget https://zenodo.org/record/XXXXXX/files/yallhap-validation-v1.tar.gz

# Verify checksum
wget https://zenodo.org/record/XXXXXX/files/yallhap-validation-v1.md5
md5sum -c yallhap-validation-v1.md5

# Extract to data directory
tar -xzf yallhap-validation-v1.tar.gz -C data/
```

### Option 2: Download from Original Sources

```bash
# 1000 Genomes Phase 3 data
python scripts/download_1kg_validation.py

# Ancient DNA data (requires AADR access)
python scripts/download_ancient_test_data.py

# T2T liftover chains
python scripts/download_liftover_chains.py
```

## Prerequisites

```bash
# Install yallhap with dev dependencies
pip install -e ".[dev]"

# Verify installation
yallhap --version
```

## Validation Data Contents

After downloading, your `data/` directory should contain:

```
data/
├── yfull_tree.json              # YFull phylogenetic tree (~14 MB)
├── ybrowse_snps.csv             # Full SNP database (~400 MB)
├── validation/
│   ├── 1kg_chrY_phase3.vcf.gz   # 1000 Genomes Y-chr VCF (5.4 MB)
│   ├── 1kg_chrY_phase3.vcf.gz.tbi
│   ├── poznik2016_haplogroups.tsv  # Ground truth (24 KB)
│   └── ybrowse_snps_hg19.csv    # GRCh37 SNP positions (48 MB)
├── ancient/
│   ├── aadr_chrY_v2.vcf.gz      # Ancient DNA VCF (97 MB)
│   ├── aadr_chrY_v2.vcf.gz.tbi
│   └── aadr_1240k_ground_truth.tsv  # Ancient ground truth (252 KB)
└── liftover/
    ├── grch38-chm13v2.chain     # T2T liftover (6 MB)
    └── hg19-chm13v2.chain       # T2T liftover (6 MB)
```

## Running Validations

### 1000 Genomes Validation

This validates yallHap against 1,244 male samples from 1000 Genomes Phase 3.

```bash
python scripts/validate_1kg.py
```

**Expected output:**

```
Loading tree from data/yfull_tree.json...
Loading SNP database from data/validation/ybrowse_snps_hg19.csv...
Loading ground truth from data/validation/poznik2016_haplogroups.tsv...
Classifying 1244 samples from data/validation/1kg_chrY_phase3.vcf.gz...

Results:
  Same major lineage: 1244/1244 (100.0%)
  Exact match: 174/1244 (14.0%)
  More specific call: 892/1244 (71.7%)
  Less specific call: 178/1244 (14.3%)

Notes:
- "Exact match" requires identical haplogroup names
- Low exact match rate is expected due to different tree versions
  (yallHap uses YFull 2024, ground truth uses ISOGG 2016)
- "Same major lineage" means the first letter matches (R, J, I, etc.)
```

#### Quick Validation (100 samples)

For faster testing during development:

```bash
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

```bash
python scripts/validate_ancient.py
```

Or manually:

```bash
python << 'EOF'
from yallhap.classifier import HaplogroupClassifier
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree
import csv

tree = Tree.from_json('data/yfull_tree.json')
db = SNPDatabase.from_ybrowse_gff_csv('data/validation/ybrowse_snps_hg19.csv')

# Ancient DNA mode with transversions-only
c = HaplogroupClassifier(
    tree=tree,
    snp_db=db,
    reference='grch37',
    ancient_mode=True,
    transversions_only=True,
)

# Load ground truth
gt = {}
with open('data/ancient/aadr_1240k_ground_truth.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gt[row['sample_id']] = row['haplogroup']

# Classify 50 samples
samples = list(gt.keys())[:50]
results = c.classify_batch('data/ancient/aadr_chrY_v2.vcf.gz', samples)

# Check accuracy
correct = 0
for i, sample in enumerate(samples):
    expected = gt[sample]
    called = results[i].haplogroup
    # Check if same lineage (called is ancestor or descendant of expected)
    if expected[0] == called[0]:  # Same major lineage
        correct += 1

print(f"Accuracy: {correct}/50 ({correct*2}%) same major lineage")
EOF
```

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
| `scripts/validate_1kg.py` | Runs validation against 1000 Genomes |
| `scripts/download_ancient_test_data.py` | Creates ancient DNA ground truth TSV |
| `scripts/eigenstrat_to_vcf.py` | Converts EIGENSTRAT format to VCF |
| `scripts/download_liftover_chains.py` | Downloads T2T liftover chain files |
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

