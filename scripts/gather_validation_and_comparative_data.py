#!/usr/bin/env python3
"""
Comprehensive Y-chromosome haplogroup validation and tool comparison.

This script runs all validation benchmarks and tool comparisons:

1. **Modern DNA Validation**
   - 1000 Genomes Phase 3 (GRCh37) - yallHap, yhaplo, Yleaf
   - gnomAD High-Coverage (GRCh38) - yallHap (heuristic, Bayesian)

2. **Ancient DNA Validation**
   - AADR v54 (stratified by coverage) - yallHap (transversions-only)
   - Local ancient samples - yallHap, yhaplo, Yleaf, pathPhynder

3. **Tool Comparison Summary**
   - Comprehensive table comparing all tools across all samples
   - Concordance analysis
   - Runtime comparison

Output: A single markdown report with embedded TSV data.

Usage:
    # Full validation (all samples, slow)
    python scripts/gather_validation_and_comparative_data.py -o results/validation_report.md

    # Quick test run (5 samples per dataset)
    python scripts/gather_validation_and_comparative_data.py -o results/validation_report.md --subsample 5

    # Fast with threading
    python scripts/gather_validation_and_comparative_data.py -o results/validation_report.md --subsample 10 --threads 8

    # Skip slow tools
    python scripts/gather_validation_and_comparative_data.py -o results/validation_report.md --skip-pathphynder --skip-yleaf
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

# Add src to path for local imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


# Variant density bin constants for stratified AADR analysis
# Based on actual called variants / total variants in Y-only VCF
DENSITY_BINS = ["<1%", "1-10%", "10-50%", ">=50%"]


def get_density_bin(density_pct: float) -> str:
    """Assign variant density percentage to bin."""
    if density_pct < 1.0:
        return "<1%"
    elif density_pct < 10.0:
        return "1-10%"
    elif density_pct < 50.0:
        return "10-50%"
    else:
        return ">=50%"


def calculate_variant_density(
    vcf_path: Path,
    sample_ids: list[str],
    cache_path: Path | None = None,
) -> dict[str, float]:
    """
    Calculate variant density for each sample from the VCF.

    Variant density = (called variants / total variants) * 100

    Args:
        vcf_path: Path to VCF file
        sample_ids: List of sample IDs to calculate density for
        cache_path: Optional path to cache file

    Returns:
        Dictionary mapping sample_id -> density percentage (0-100)
    """
    import pysam

    # Check cache
    if cache_path and cache_path.exists():
        vcf_mtime = vcf_path.stat().st_mtime
        cache_mtime = cache_path.stat().st_mtime
        if cache_mtime > vcf_mtime:
            print("    Loading cached variant density...")
            with open(cache_path) as f:
                cached = json.load(f)
            # Filter to requested samples
            return {s: cached[s] for s in sample_ids if s in cached}

    print("    Calculating variant density from VCF (this may take a few minutes)...")

    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)

    # Filter to samples in VCF
    target_samples = [s for s in sample_ids if s in vcf_samples]

    # Initialize counts
    called_counts: dict[str, int] = {s: 0 for s in target_samples}
    total_variants = 0

    # Count called variants per sample
    for i, rec in enumerate(vcf.fetch()):
        total_variants += 1
        if i % 5000 == 0 and i > 0:
            print(f"      Processed {i} variants...")
        for sample in target_samples:
            gt = rec.samples[sample]["GT"]
            if gt is not None and gt != (None,) and gt != (None, None):
                called_counts[sample] += 1

    vcf.close()
    print(f"      Processed {total_variants} total variants")

    # Calculate density percentages
    density_map: dict[str, float] = {}
    for sample, called in called_counts.items():
        density_map[sample] = (called / total_variants * 100) if total_variants > 0 else 0.0

    # Cache results
    if cache_path:
        print(f"    Caching variant density to {cache_path.name}...")
        all_densities = {}
        if cache_path.exists():
            with open(cache_path) as f:
                all_densities = json.load(f)
        all_densities.update(density_map)
        with open(cache_path, "w") as f:
            json.dump(all_densities, f)

    return density_map


@dataclass
class SampleInfo:
    """Information about a sample for validation."""

    name: str
    vcf_path: Path
    vcf_sample_id: str
    bam_path: Path | None = None
    description: str = ""
    period: str = ""
    reference: str = "grch37"
    ground_truth_terminal: str = ""  # AADR terminal mutation format
    ground_truth_isogg: str = ""  # AADR ISOGG format


@dataclass
class HaplogroupCall:
    """A haplogroup call from any tool."""

    tool: str
    sample: str
    haplogroup: str
    confidence: float | None = None
    runtime_seconds: float = 0.0
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class ValidationResult:
    """Results from a validation run."""

    dataset: str
    tool: str
    total_samples: int
    successful_calls: int
    correct_major: int = 0  # Same major haplogroup as ground truth
    mean_confidence: float | None = None
    runtime_seconds: float = 0.0
    calls: list[HaplogroupCall] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def accuracy(self) -> float:
        """Major haplogroup accuracy percentage."""
        if self.total_samples == 0:
            return 0.0
        return 100 * self.correct_major / self.total_samples


def find_tool(name: str) -> str | None:
    """Find a tool in PATH or common locations."""
    path = shutil.which(name)
    if path:
        return path
    common_paths = [
        Path.home() / ".local" / "bin" / name,
        Path("/usr/local/bin") / name,
    ]
    for p in common_paths:
        if p.exists():
            return str(p)
    return None


def load_ground_truth(path: Path, haplogroup_col: str = "haplogroup") -> dict[str, str]:
    """Load ground truth haplogroups from TSV file."""
    ground_truth = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_id = row.get("sample_id", row.get("sample", ""))
            haplogroup = row.get(haplogroup_col, row.get("haplogroup", ""))
            if sample_id and haplogroup:
                ground_truth[sample_id] = haplogroup
    return ground_truth


def compare_major_haplogroup(expected: str, predicted: str) -> bool:
    """Compare major haplogroup letters (e.g., R, I, J)."""
    if not expected or not predicted:
        return False
    if predicted in ("NA", "ERROR", "TIMEOUT", "NO_BAM"):
        return False
    return expected[0].upper() == predicted[0].upper()


# =============================================================================
# yallHap runners
# =============================================================================


def run_yallhap_ancient(
    samples: list[SampleInfo],
    tree_path: Path,
    snp_db_path: Path,
    mode: str = "regular",
    threads: int = 1,
    isogg_db_path: Path | None = None,
) -> ValidationResult:
    """Run yallHap on ancient samples."""
    from yallhap.classifier import HaplogroupClassifier
    from yallhap.snps import SNPDatabase
    from yallhap.tree import Tree

    tree = Tree.from_json(tree_path)
    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)

    bayesian = mode == "bayesian"
    ancient_mode = mode in ("regular", "bayesian", "isogg")
    transversions_only = mode == "transversions"
    min_depth = 1 if (ancient_mode and bayesian) else 1
    min_quality = 0 if (ancient_mode and bayesian) else 10

    classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference="grch37",
        bayesian=bayesian,
        ancient_mode=ancient_mode,
        transversions_only=transversions_only,
        min_depth=min_depth,
        min_quality=min_quality,
    )

    # Load ISOGG mapper if in ISOGG mode
    isogg_mapper = None
    if mode == "isogg":
        if not isogg_db_path or not isogg_db_path.exists():
            raise ValueError(f"ISOGG mode requires ISOGG database, not found: {isogg_db_path}")

        from yallhap.isogg import ISOGGDatabase, ISOGGMapper

        isogg_db = ISOGGDatabase.from_file(isogg_db_path)
        isogg_mapper = ISOGGMapper(tree, isogg_db)

    calls: list[HaplogroupCall] = []
    start_time = time.time()

    for sample in samples:
        try:
            result = classifier.classify(sample.vcf_path, sample.vcf_sample_id)
            extra: dict[str, Any] = {}
            if result.snp_stats:
                extra["derived"] = result.snp_stats.derived
                extra["ancestral"] = result.snp_stats.ancestral

            # Add ISOGG haplogroup if mapper available
            if isogg_mapper and result.haplogroup not in ("NA", "ERROR"):
                try:
                    isogg_hg = isogg_mapper.to_isogg(result.haplogroup)
                    if isogg_hg:
                        extra["isogg_haplogroup"] = isogg_hg
                except Exception:
                    pass

            calls.append(
                HaplogroupCall(
                    tool=f"yallHap-{mode}",
                    sample=sample.name,
                    haplogroup=result.haplogroup,
                    confidence=result.confidence,
                    extra=extra,
                )
            )
        except Exception as e:
            calls.append(
                HaplogroupCall(
                    tool=f"yallHap-{mode}",
                    sample=sample.name,
                    haplogroup="ERROR",
                    extra={"error": str(e)},
                )
            )

    runtime = time.time() - start_time
    successful = sum(1 for c in calls if c.haplogroup not in ("NA", "ERROR"))
    mean_conf = (
        sum(c.confidence for c in calls if c.confidence is not None) / len(calls)
        if calls
        else None
    )

    return ValidationResult(
        dataset="ancient-local",
        tool=f"yallHap-{mode}",
        total_samples=len(samples),
        successful_calls=successful,
        mean_confidence=mean_conf,
        runtime_seconds=runtime,
        calls=calls,
    )


def run_yallhap_benchmark(
    vcf_path: Path,
    tree_path: Path,
    snp_db_path: Path,
    ground_truth: dict[str, str],
    samples: list[str],
    reference: str = "grch37",
    bayesian: bool = False,
    ancient_mode: bool = False,
    transversions_only: bool = False,
    threads: int = 1,
    dataset_name: str = "benchmark",
) -> ValidationResult:
    """Run yallHap benchmark on a VCF with multiple samples."""
    from yallhap.classifier import HaplogroupClassifier
    from yallhap.snps import SNPDatabase
    from yallhap.tree import Tree

    tree = Tree.from_json(tree_path)
    snp_db = SNPDatabase.from_ybrowse_gff_csv(snp_db_path)

    classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference=reference,
        bayesian=bayesian,
        ancient_mode=ancient_mode,
        transversions_only=transversions_only,
        min_depth=1,
    )

    start_time = time.time()
    results = classifier.classify_batch(vcf_path, samples, threads=threads)
    runtime = time.time() - start_time

    calls: list[HaplogroupCall] = []
    correct_major = 0
    confidences = []

    for result in results:
        calls.append(
            HaplogroupCall(
                tool="yallHap-Bayes" if bayesian else "yallHap",
                sample=result.sample,
                haplogroup=result.haplogroup,
                confidence=result.confidence,
            )
        )
        if result.confidence is not None:
            confidences.append(result.confidence)

        expected = ground_truth.get(result.sample, "")
        if compare_major_haplogroup(expected, result.haplogroup):
            correct_major += 1

    successful = sum(1 for c in calls if c.haplogroup not in ("NA", "ERROR"))
    mean_conf = sum(confidences) / len(confidences) if confidences else None

    tool_name = "yallHap-Bayes" if bayesian else "yallHap"
    notes = [f"{reference.upper()} reference"]
    if bayesian:
        notes.append("Bayesian mode")
    if ancient_mode:
        notes.append("Ancient mode")
    if transversions_only:
        notes.append("Transversions-only")

    return ValidationResult(
        dataset=dataset_name,
        tool=tool_name,
        total_samples=len(samples),
        successful_calls=successful,
        correct_major=correct_major,
        mean_confidence=mean_conf,
        runtime_seconds=runtime,
        calls=calls,
        notes=notes,
    )


# =============================================================================
# External tool runners
# =============================================================================


def run_yhaplo_ancient(samples: list[SampleInfo], output_dir: Path) -> ValidationResult:
    """Run yhaplo on ancient samples."""
    yhaplo_path = find_tool("yhaplo")
    if not yhaplo_path:
        return ValidationResult(
            dataset="ancient-local",
            tool="yhaplo",
            total_samples=len(samples),
            successful_calls=0,
            calls=[],
            notes=["yhaplo not found"],
        )

    calls: list[HaplogroupCall] = []
    start_time = time.time()
    yhaplo_out = output_dir / "yhaplo"
    yhaplo_out.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        try:
            sample_out = yhaplo_out / sample.name
            sample_out.mkdir(exist_ok=True)

            cmd = [yhaplo_path, "-i", str(sample.vcf_path), "-o", str(sample_out)]
            subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            # yhaplo outputs to haplogroups.<vcf_basename>.txt
            haplogroup = "NA"
            ycc_haplogroup = None

            # Search for haplogroups.*.txt files
            for hg_file in sample_out.glob("haplogroups.*.txt"):
                with open(hg_file) as f:
                    for line in f:
                        if line.startswith("id") or line.startswith("#") or not line.strip():
                            continue
                        # Format: sample_id<TAB>hg_snp<TAB>hg_snp<TAB>ycc_haplogroup
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            haplogroup = parts[1]  # SNP-based haplogroup
                        if len(parts) >= 4:
                            ycc_haplogroup = parts[3]  # YCC nomenclature
                        break
                break

            extra = {}
            if ycc_haplogroup:
                extra["ycc_haplogroup"] = ycc_haplogroup

            calls.append(
                HaplogroupCall(
                    tool="yhaplo",
                    sample=sample.name,
                    haplogroup=haplogroup,
                    extra=extra,
                )
            )
        except subprocess.TimeoutExpired:
            calls.append(
                HaplogroupCall(tool="yhaplo", sample=sample.name, haplogroup="TIMEOUT")
            )
        except Exception as e:
            calls.append(
                HaplogroupCall(
                    tool="yhaplo",
                    sample=sample.name,
                    haplogroup="ERROR",
                    extra={"error": str(e)},
                )
            )

    runtime = time.time() - start_time
    successful = sum(1 for c in calls if c.haplogroup not in ("NA", "ERROR", "TIMEOUT"))

    return ValidationResult(
        dataset="ancient-local",
        tool="yhaplo",
        total_samples=len(samples),
        successful_calls=successful,
        runtime_seconds=runtime,
        calls=calls,
    )


def run_yleaf_ancient(
    samples: list[SampleInfo],
    output_dir: Path,
    ancient_mode: bool = False,  # Disabled by default - causes empty predictions
) -> ValidationResult:
    """
    Run Yleaf on ancient samples.

    Note: The -aDNA flag filters too aggressively for low-coverage ancient DNA,
    often resulting in no predictions. By default, we run without -aDNA.
    """
    yleaf_path = find_tool("Yleaf")
    if not yleaf_path:
        return ValidationResult(
            dataset="ancient-local",
            tool="Yleaf",
            total_samples=len(samples),
            successful_calls=0,
            calls=[],
            notes=["Yleaf not found"],
        )

    calls: list[HaplogroupCall] = []
    start_time = time.time()
    yleaf_out = output_dir / "yleaf"
    yleaf_out.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        try:
            sample_out = yleaf_out / sample.name
            # Use absolute path for Yleaf output
            cmd = [
                yleaf_path,
                "-vcf", str(sample.vcf_path.absolute()),
                "-rg", "hg19" if sample.reference == "grch37" else "hg38",
                "-o", str(sample_out.absolute()),
                "-r", "1",
                "-q", "10",
                "-pq", "0.5",
                "-force",
            ]
            if ancient_mode:
                cmd.append("-aDNA")

            subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            # Look for hg_prediction.hg in the output directory
            haplogroup = "NA"
            qc_score = None

            # Yleaf outputs hg_prediction.hg directly in the output directory
            pred_file = sample_out / "hg_prediction.hg"
            if not pred_file.exists():
                # Also check in subdirectories
                pred_files = list(sample_out.glob("**/hg_prediction.hg"))
                if pred_files:
                    pred_file = pred_files[0]

            if pred_file.exists():
                with open(pred_file) as f:
                    for line in f:
                        if line.startswith("Sample") or line.startswith("#"):
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            # Column 1 is haplogroup, column 6 is QC-score
                            hg = parts[1]
                            if hg and hg != "NA":
                                haplogroup = hg
                            if len(parts) >= 6:
                                try:
                                    qc_score = float(parts[5])
                                except ValueError:
                                    pass
                            break

            calls.append(
                HaplogroupCall(
                    tool="Yleaf",
                    sample=sample.name,
                    haplogroup=haplogroup,
                    confidence=qc_score,
                )
            )
        except subprocess.TimeoutExpired:
            calls.append(
                HaplogroupCall(tool="Yleaf", sample=sample.name, haplogroup="TIMEOUT")
            )
        except Exception as e:
            calls.append(
                HaplogroupCall(
                    tool="Yleaf",
                    sample=sample.name,
                    haplogroup="ERROR",
                    extra={"error": str(e)},
                )
            )

    runtime = time.time() - start_time
    successful = sum(1 for c in calls if c.haplogroup not in ("NA", "ERROR", "TIMEOUT"))
    confs = [c.confidence for c in calls if c.confidence is not None]
    mean_conf = sum(confs) / len(confs) if confs else None

    return ValidationResult(
        dataset="ancient-local",
        tool="Yleaf",
        total_samples=len(samples),
        successful_calls=successful,
        mean_confidence=mean_conf,
        runtime_seconds=runtime,
        calls=calls,
    )


def run_pathphynder_ancient(
    samples: list[SampleInfo],
    output_dir: Path,
    tree_dir: Path,
    tree_filename: str,
    tree_prefix: str,
    reference_path: Path,
    isogg_snps_path: Path | None = None,
) -> ValidationResult:
    """
    Run pathPhynder on ancient samples (requires BAM).

    NOTE: pathPhynder must be run from the directory containing:
      - tree_data/<prefix>.sites.bed (from prepare step)
      - intree_folder/ (created during run)
      - results_folder/ (created during run)

    Args:
        samples: List of samples
        output_dir: Output directory for results (for copying final outputs)
        tree_dir: Directory containing tree file and tree_data/
        tree_filename: Tree file name (relative to tree_dir)
        tree_prefix: Relative prefix for tree data (e.g., "tree_data/bigtree")
        reference_path: Path to reference genome FASTA
        isogg_snps_path: Optional path to ISOGG SNPs file
    """
    pathphynder_path = find_tool("pathPhynder")
    if not pathphynder_path:
        return ValidationResult(
            dataset="ancient-local",
            tool="pathPhynder",
            total_samples=len(samples),
            successful_calls=0,
            calls=[],
            notes=["pathPhynder not found"],
        )

    # Check required files
    tree_path = tree_dir / tree_filename
    if not tree_path.exists():
        return ValidationResult(
            dataset="ancient-local",
            tool="pathPhynder",
            total_samples=len(samples),
            successful_calls=0,
            calls=[],
            notes=[f"Tree not found: {tree_path}"],
        )

    # Check if tree_data directory exists with prepared files
    tree_data_dir = tree_dir / "tree_data"
    prefix_base = tree_prefix.split("/")[-1] if "/" in tree_prefix else tree_prefix
    sites_bed = tree_data_dir / f"{prefix_base}.sites.bed"
    if not sites_bed.exists():
        return ValidationResult(
            dataset="ancient-local",
            tool="pathPhynder",
            total_samples=len(samples),
            successful_calls=0,
            calls=[],
            notes=[
                f"pathPhynder requires 'prepare' step first",
                f"tree_data/{prefix_base}.sites.bed not found",
                f"Run: cd {tree_dir} && pathPhynder -s prepare -i {tree_filename} -f <branches.snp> -p {prefix_base}",
            ],
        )

    # Ensure intree_folder and results_folder exist
    (tree_dir / "intree_folder").mkdir(exist_ok=True)
    (tree_dir / "results_folder").mkdir(exist_ok=True)

    samples_with_bam = [s for s in samples if s.bam_path and s.bam_path.exists()]
    if not samples_with_bam:
        return ValidationResult(
            dataset="ancient-local",
            tool="pathPhynder",
            total_samples=len(samples),
            successful_calls=0,
            calls=[],
            notes=["No samples with BAM files"],
        )

    calls: list[HaplogroupCall] = []
    start_time = time.time()
    pp_out = output_dir / "pathphynder"
    pp_out.mkdir(parents=True, exist_ok=True)

    for sample in samples_with_bam:
        try:
            # pathPhynder requires running from tree_dir with relative paths
            # Output prefix should be a simple name (no path separators)
            output_name = sample.name.replace("/", "_").replace(".", "_")

            cmd = [
                pathphynder_path,
                "-s", "all",
                "-i", tree_filename,
                "-p", tree_prefix,
                "-b", str(sample.bam_path.absolute()),
                "-r", str(reference_path.absolute()),
                "-o", output_name,
                "-m", "transversions",
            ]
            if isogg_snps_path and isogg_snps_path.exists():
                cmd.extend(["-G", str(isogg_snps_path.absolute())])

            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=600, cwd=str(tree_dir)
            )

            haplogroup = "NA"
            isogg_hg = None

            # Results are in tree_dir/results_folder/<output_name>.*
            results_folder = tree_dir / "results_folder"

            # Check best_node_info.txt for the best node
            best_node_file = results_folder / f"{output_name}.best_node_info.txt"
            if best_node_file.exists():
                with open(best_node_file) as f:
                    for line in f:
                        if line.startswith("sample"):
                            continue  # Header
                        parts = line.strip().split("\t")
                        if len(parts) >= 2 and parts[1]:
                            haplogroup = f"node_{parts[1]}"  # Node ID from tree
                        break

            # Check best_path_report.txt for haplogroup (hg column)
            best_path_file = results_folder / f"{output_name}.best_path_report.txt"
            if best_path_file.exists():
                with open(best_path_file) as f:
                    last_hg = None
                    for line in f:
                        if line.startswith("Edge"):
                            continue  # Header
                        parts = line.strip().split("\t")
                        # Format: Edge Node1 Node2 support conflict snp_count hg
                        if len(parts) >= 7 and parts[6] and parts[6] != "NA":
                            last_hg = parts[6]
                    if last_hg:
                        haplogroup = last_hg

            # Check ISOGG determination if available
            # Format: lines of "hg\tDer\tAnc" then "Haplogroup:\t<haplogroup>\t"
            isogg_file = results_folder / f"{output_name}.isogg_hg_determination.txt"
            if isogg_file.exists():
                with open(isogg_file) as f:
                    for line in f:
                        line = line.strip()
                        # Look for the final "Haplogroup:" line
                        if line.startswith("Haplogroup:"):
                            parts = line.split("\t")
                            if len(parts) >= 2 and parts[1]:
                                isogg_hg = parts[1]
                                # Use ISOGG haplogroup as the main haplogroup
                                haplogroup = isogg_hg
                            break

            extra = {}
            if isogg_hg:
                extra["isogg_haplogroup"] = isogg_hg
            if result.returncode != 0:
                extra["error"] = f"Exit code {result.returncode}"
                if result.stderr:
                    extra["stderr"] = result.stderr[:200]

            calls.append(
                HaplogroupCall(
                    tool="pathPhynder",
                    sample=sample.name,
                    haplogroup=haplogroup,
                    extra=extra,
                )
            )

            # Copy result files to output directory
            sample_out = pp_out / sample.name
            sample_out.mkdir(exist_ok=True)
            for f in results_folder.glob(f"{output_name}.*"):
                shutil.copy(f, sample_out / f.name)

        except subprocess.TimeoutExpired:
            calls.append(
                HaplogroupCall(
                    tool="pathPhynder", sample=sample.name, haplogroup="TIMEOUT"
                )
            )
        except Exception as e:
            calls.append(
                HaplogroupCall(
                    tool="pathPhynder",
                    sample=sample.name,
                    haplogroup="ERROR",
                    extra={"error": str(e)},
                )
            )

    # Add NA for samples without BAM
    for sample in samples:
        if sample not in samples_with_bam:
            calls.append(
                HaplogroupCall(
                    tool="pathPhynder", sample=sample.name, haplogroup="NO_BAM"
                )
            )

    runtime = time.time() - start_time
    successful = sum(
        1 for c in calls if c.haplogroup not in ("NA", "ERROR", "TIMEOUT", "NO_BAM")
    )

    return ValidationResult(
        dataset="ancient-local",
        tool="pathPhynder",
        total_samples=len(samples),
        successful_calls=successful,
        runtime_seconds=runtime,
        calls=calls,
    )


# =============================================================================
# Sample loaders
# =============================================================================


def load_aadr_ground_truth(anno_path: Path) -> dict[str, tuple[str, str]]:
    """
    Load AADR ground truth haplogroups from annotation file.

    Returns dict of sample_id -> (terminal_format, isogg_format)
    """
    if not anno_path.exists():
        return {}

    ground_truth: dict[str, tuple[str, str]] = {}

    with open(anno_path, encoding="utf-8", errors="replace") as f:
        header = None
        for line in f:
            parts = line.strip().split("\t")
            if header is None:
                header = parts
                continue

            if len(parts) < 26:
                continue

            genetic_id = parts[0]
            terminal = parts[24] if len(parts) > 24 else ""
            isogg = parts[25] if len(parts) > 25 else ""

            # Normalize missing values
            if terminal in ("..", "n/a", "n/a (female)", ""):
                terminal = ""
            if isogg in ("..", "n/a", "n/a (female)", ""):
                isogg = ""

            if terminal or isogg:
                ground_truth[genetic_id] = (terminal, isogg)

    return ground_truth


def load_ancient_samples(
    ancient_dir: Path,
    subsample: int | None = None,
    aadr_anno_path: Path | None = None,
) -> list[SampleInfo]:
    """Load sample information from the ancient genomes directory."""
    # Load AADR ground truth if available
    aadr_gt: dict[str, tuple[str, str]] = {}
    if aadr_anno_path and aadr_anno_path.exists():
        aadr_gt = load_aadr_ground_truth(aadr_anno_path)

    # Sample definitions with AADR ID mappings
    # (name, vcf_prefix, vcf_sample_id, desc, period, aadr_ids)
    sample_defs = [
        ("I0231", "I0231.390k.chrY", "SM", "Yamnaya", "Bronze Age (~3000 BCE)", ["I0231"]),
        ("I0357", "I0357.390k.chrY", "SM", "Yamnaya", "Bronze Age (~3000 BCE)", ["I0357"]),
        ("I0443", "I0443.390k.chrY", "SM", "Yamnaya", "Bronze Age (~3000 BCE)", ["I0443"]),
        (
            "Kennewick",
            "Kennewick_defaultMap1extr.realign.md.head.rmdup.chrY",
            "Ken19",
            "Kennewick Man",
            "Paleoamerican (~9000 BP)",
            ["kennewick_noUDG.SG", "Kennewick"],
        ),
        (
            "SB524A",
            "SB524A_lib.merged.markdup.chrY",
            "Cheddar",
            "Cheddar Man (lib 1)",
            "British Mesolithic (~10000 BP)",
            ["I6767", "I6767.SG", "I6767_noUDG.SG"],  # Cheddar Man
        ),
        (
            "SB524A2",
            "SB524A2_lib.merged.markdup.chrY",
            "Cheddar",
            "Cheddar Man (lib 2)",
            "British Mesolithic (~10000 BP)",
            ["I6767", "I6767.SG", "I6767_noUDG.SG"],  # Cheddar Man
        ),
        (
            "VK287",
            "VK287.final.chrY",
            "VK287",
            "Viking",
            "Viking Age (~800-1050 CE)",
            ["VK287_noUDG.SG", "VK287.SG", "VK287"],
        ),
        (
            "VK292",
            "VK292.final.chrY",
            "VK292",
            "Viking",
            "Viking Age (~800-1050 CE)",
            ["VK292_noUDG.SG", "VK292.SG", "VK292"],
        ),
        (
            "VK296",
            "VK296.final.chrY",
            "VK296",
            "Viking",
            "Viking Age (~800-1050 CE)",
            ["VK296_noUDG.SG", "VK296.SG", "VK296"],
        ),
        (
            "VK582",
            "VK582.final.chrY",
            "IA_PD_22",
            "Viking",
            "Viking Age (~800-1050 CE)",
            ["VK582_noUDG.SG", "VK582.SG", "VK582"],
        ),
    ]

    samples = []
    for name, prefix, sample_id, desc, period, aadr_ids in sample_defs:
        vcf_path = ancient_dir / f"{prefix}.vcf.gz"
        bam_path = ancient_dir / f"{prefix}.bam"

        if vcf_path.exists():
            # Find ground truth from AADR
            gt_terminal = ""
            gt_isogg = ""
            for aadr_id in aadr_ids:
                if aadr_id in aadr_gt:
                    gt_terminal, gt_isogg = aadr_gt[aadr_id]
                    if gt_terminal or gt_isogg:
                        break

            samples.append(
                SampleInfo(
                    name=name,
                    vcf_path=vcf_path,
                    vcf_sample_id=sample_id,
                    bam_path=bam_path if bam_path.exists() else None,
                    description=desc,
                    period=period,
                    reference="grch37",
                    ground_truth_terminal=gt_terminal,
                    ground_truth_isogg=gt_isogg,
                )
            )

    if subsample and len(samples) > subsample:
        random.seed(42)
        samples = random.sample(samples, subsample)

    return samples


# =============================================================================
# Benchmark runners
# =============================================================================


def benchmark_1kg(
    base_dir: Path,
    subsample: int | None = None,
    threads: int = 1,
) -> list[ValidationResult]:
    """Run 1000 Genomes Phase 3 benchmark."""
    validation_dir = base_dir / "data" / "validation"
    data_dir = base_dir / "data"

    ground_truth_path = validation_dir / "poznik2016_haplogroups.tsv"
    vcf_path = validation_dir / "1kg_chrY_phase3.vcf.gz"
    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = validation_dir / "ybrowse_snps_hg19.csv"

    for path in [ground_truth_path, vcf_path, tree_path, snp_db_path]:
        if not path.exists():
            print(f"    SKIP: {path.name} not found")
            return []

    ground_truth = load_ground_truth(ground_truth_path)
    samples = list(ground_truth.keys())

    if subsample and len(samples) > subsample:
        random.seed(42)
        samples = random.sample(samples, subsample)

    print(f"    {len(samples)} samples")

    results = []

    # Heuristic mode
    print("    Running heuristic mode...")
    results.append(
        run_yallhap_benchmark(
            vcf_path=vcf_path,
            tree_path=tree_path,
            snp_db_path=snp_db_path,
            ground_truth=ground_truth,
            samples=samples,
            reference="grch37",
            bayesian=False,
            threads=threads,
            dataset_name="1KG Phase 3",
        )
    )

    # Bayesian mode
    print("    Running Bayesian mode...")
    results.append(
        run_yallhap_benchmark(
            vcf_path=vcf_path,
            tree_path=tree_path,
            snp_db_path=snp_db_path,
            ground_truth=ground_truth,
            samples=samples,
            reference="grch37",
            bayesian=True,
            threads=threads,
            dataset_name="1KG Phase 3",
        )
    )

    return results


def benchmark_gnomad(
    base_dir: Path,
    subsample: int | None = None,
    threads: int = 1,
) -> list[ValidationResult]:
    """Run gnomAD high-coverage benchmark."""
    highcov_dir = base_dir / "data" / "validation_highcov"
    validation_dir = base_dir / "data" / "validation"
    data_dir = base_dir / "data"

    ground_truth_path = validation_dir / "poznik2016_haplogroups.tsv"
    vcf_path = highcov_dir / "vcf" / "gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz"
    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = data_dir / "ybrowse_snps.csv"

    # Check for pre-extracted subset
    subset_vcf = highcov_dir / "vcf" / "gnomad_1kg_shared_diagnostic.vcf.gz"
    if subset_vcf.exists():
        vcf_path = subset_vcf

    for path in [ground_truth_path, tree_path, snp_db_path]:
        if not path.exists():
            print(f"    SKIP: {path.name} not found")
            return []

    if not vcf_path.exists():
        print(f"    SKIP: gnomAD VCF not found")
        return []

    import pysam

    ground_truth = load_ground_truth(ground_truth_path)
    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)
    vcf.close()

    samples = [s for s in ground_truth.keys() if s in vcf_samples]

    if subsample and len(samples) > subsample:
        random.seed(42)
        samples = random.sample(samples, subsample)

    if not samples:
        print("    SKIP: No overlapping samples")
        return []

    print(f"    {len(samples)} samples")

    results = []

    print("    Running heuristic mode...")
    results.append(
        run_yallhap_benchmark(
            vcf_path=vcf_path,
            tree_path=tree_path,
            snp_db_path=snp_db_path,
            ground_truth=ground_truth,
            samples=samples,
            reference="grch38",
            bayesian=False,
            threads=threads,
            dataset_name="gnomAD High-Cov",
        )
    )

    print("    Running Bayesian mode...")
    results.append(
        run_yallhap_benchmark(
            vcf_path=vcf_path,
            tree_path=tree_path,
            snp_db_path=snp_db_path,
            ground_truth=ground_truth,
            samples=samples,
            reference="grch38",
            bayesian=True,
            threads=threads,
            dataset_name="gnomAD High-Cov",
        )
    )

    return results


def benchmark_aadr(
    base_dir: Path,
    subsample: int | None = None,
    threads: int = 1,
    samples_per_bin: int = 500,
) -> list[ValidationResult]:
    """
    Run AADR ancient DNA benchmark with variant density stratification.

    Tests both heuristic transversions-only and Bayesian ancient modes,
    stratified by variant density bins.
    """
    from collections import defaultdict

    ancient_dir = base_dir / "data" / "ancient"
    validation_dir = base_dir / "data" / "validation"
    data_dir = base_dir / "data"

    ground_truth_path = ancient_dir / "aadr_1240k_ground_truth.tsv"
    vcf_path = ancient_dir / "aadr_chrY_v2.vcf.gz"
    tree_path = data_dir / "yfull_tree.json"
    snp_db_path = validation_dir / "ybrowse_snps_hg19.csv"
    cache_path = ancient_dir / "aadr_variant_density.json"

    for path in [ground_truth_path, vcf_path, tree_path, snp_db_path]:
        if not path.exists():
            print(f"    SKIP: {path.name} not found")
            return []

    import pysam

    raw_ground_truth = load_ground_truth(ground_truth_path, "haplogroup_terminal")
    ground_truth = {
        k: v for k, v in raw_ground_truth.items() if v and "-" in v and v[0].isalpha()
    }

    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)
    vcf.close()

    valid_samples = [s for s in ground_truth.keys() if s in vcf_samples]
    print(f"    {len(valid_samples)} samples with ground truth in VCF")

    # Calculate variant density
    density_map = calculate_variant_density(vcf_path, valid_samples, cache_path)

    # Group samples by density bin
    samples_by_bin: dict[str, list[str]] = defaultdict(list)
    for sample_id in valid_samples:
        density = density_map.get(sample_id, 0.0)
        bin_name = get_density_bin(density)
        samples_by_bin[bin_name].append(sample_id)

    print("    Samples by variant density bin:")
    for bin_name in DENSITY_BINS:
        count = len(samples_by_bin.get(bin_name, []))
        print(f"      {bin_name}: {count}")

    # Select samples from each bin
    selected_samples: list[str] = []
    for bin_name in DENSITY_BINS:
        bin_samples = samples_by_bin.get(bin_name, [])
        if len(bin_samples) > samples_per_bin:
            random.seed(42)
            bin_samples = random.sample(bin_samples, samples_per_bin)
        selected_samples.extend(bin_samples)

    # Apply overall subsample if requested
    if subsample and len(selected_samples) > subsample:
        random.seed(42)
        selected_samples = random.sample(selected_samples, subsample)

    if not selected_samples:
        print("    SKIP: No overlapping samples")
        return []

    print(f"    Running validation on {len(selected_samples)} samples...")

    results = []

    print("    Running heuristic transversions-only mode...")
    results.append(
        run_yallhap_benchmark(
            vcf_path=vcf_path,
            tree_path=tree_path,
            snp_db_path=snp_db_path,
            ground_truth=ground_truth,
            samples=selected_samples,
            reference="grch37",
            bayesian=False,
            transversions_only=True,
            threads=threads,
            dataset_name="AADR v54",
        )
    )

    print("    Running Bayesian ancient mode...")
    results.append(
        run_yallhap_benchmark(
            vcf_path=vcf_path,
            tree_path=tree_path,
            snp_db_path=snp_db_path,
            ground_truth=ground_truth,
            samples=selected_samples,
            reference="grch37",
            bayesian=True,
            ancient_mode=True,
            transversions_only=False,  # Bayesian ancient mode, not transversions
            threads=threads,
            dataset_name="AADR v54",
        )
    )

    return results


# =============================================================================
# Report generation
# =============================================================================


def generate_comparison_table(
    results: list[ValidationResult],
) -> tuple[list[str], list[list[str]]]:
    """Generate comparison table from results."""
    all_samples = set()
    for r in results:
        for c in r.calls:
            all_samples.add(c.sample)

    tool_calls: dict[str, dict[str, HaplogroupCall]] = {}
    for r in results:
        if r.tool not in tool_calls:
            tool_calls[r.tool] = {}
        for c in r.calls:
            tool_calls[r.tool][c.sample] = c

    tools = sorted(tool_calls.keys())
    headers = ["Sample"] + tools

    rows = []
    for sample in sorted(all_samples):
        row = [sample]
        for tool in tools:
            call = tool_calls.get(tool, {}).get(sample)
            if call:
                hg = call.haplogroup
                # Add ISOGG haplogroup if available
                if call.extra.get("isogg_haplogroup"):
                    hg = f"{hg} (ISOGG: {call.extra['isogg_haplogroup']})"
                # Add YCC haplogroup for yhaplo
                elif call.extra.get("ycc_haplogroup"):
                    hg = f"{hg} (YCC: {call.extra['ycc_haplogroup']})"
                row.append(hg)
            else:
                row.append("N/A")
        rows.append(row)

    return headers, rows


def generate_comparison_table_with_ground_truth(
    results: list[ValidationResult],
    samples: list[SampleInfo],
) -> tuple[list[str], list[list[str]]]:
    """Generate comparison table with ground truth for ancient samples."""
    # Build sample lookup
    sample_gt: dict[str, tuple[str, str]] = {
        s.name: (s.ground_truth_terminal, s.ground_truth_isogg) for s in samples
    }

    tool_calls: dict[str, dict[str, HaplogroupCall]] = {}
    for r in results:
        if r.tool not in tool_calls:
            tool_calls[r.tool] = {}
        for c in r.calls:
            tool_calls[r.tool][c.sample] = c

    tools = sorted(tool_calls.keys())

    # Headers: Sample, AADR Terminal, AADR ISOGG, then tools
    headers = ["Sample", "AADR Terminal", "AADR ISOGG"] + tools

    rows = []
    for sample in sorted(s.name for s in samples):
        gt_terminal, gt_isogg = sample_gt.get(sample, ("", ""))
        row = [sample, gt_terminal or "-", gt_isogg or "-"]

        for tool in tools:
            call = tool_calls.get(tool, {}).get(sample)
            if call:
                hg = call.haplogroup
                if call.extra.get("isogg_haplogroup"):
                    hg = f"{hg} (ISOGG: {call.extra['isogg_haplogroup']})"
                elif call.extra.get("ycc_haplogroup"):
                    hg = f"{hg} (YCC: {call.extra['ycc_haplogroup']})"
                row.append(hg)
            else:
                row.append("N/A")
        rows.append(row)

    return headers, rows


def generate_markdown_report(
    results: list[ValidationResult],
    output_path: Path,
    ancient_samples: list[SampleInfo],
) -> None:
    """Generate comprehensive markdown report."""
    with open(output_path, "w") as f:
        f.write("# Y-Chromosome Haplogroup Validation Report\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Executive Summary
        f.write("## Executive Summary\n\n")
        f.write("| Dataset | Tool | Samples | Accuracy | Runtime |\n")
        f.write("|---------|------|---------|----------|--------|\n")
        for r in results:
            acc = f"{r.accuracy:.1f}%" if r.correct_major > 0 else "N/A"
            f.write(
                f"| {r.dataset} | {r.tool} | {r.total_samples} | {acc} | {r.runtime_seconds:.1f}s |\n"
            )
        f.write("\n")

        # Group results by dataset
        datasets: dict[str, list[ValidationResult]] = {}
        for r in results:
            if r.dataset not in datasets:
                datasets[r.dataset] = []
            datasets[r.dataset].append(r)

        # Per-dataset sections
        for dataset, ds_results in datasets.items():
            f.write(f"## {dataset}\n\n")

            if dataset == "ancient-local":
                # Sample info with ground truth
                f.write("### Sample Information\n\n")
                f.write("| Sample | Description | Period | AADR Terminal | AADR ISOGG | Has BAM |\n")
                f.write("|--------|-------------|--------|---------------|------------|--------|\n")
                for s in ancient_samples:
                    has_bam = "✓" if s.bam_path and s.bam_path.exists() else "✗"
                    gt_t = s.ground_truth_terminal or "-"
                    gt_i = s.ground_truth_isogg or "-"
                    f.write(f"| {s.name} | {s.description} | {s.period} | {gt_t} | {gt_i} | {has_bam} |\n")
                f.write("\n")

                # Comparison table with ground truth
                headers, rows = generate_comparison_table_with_ground_truth(ds_results, ancient_samples)
                f.write("### Tool Comparison (with AADR Ground Truth)\n\n")
                f.write("| " + " | ".join(headers) + " |\n")
                f.write("|" + "|".join(["---"] * len(headers)) + "|\n")
                for row in rows:
                    f.write("| " + " | ".join(row) + " |\n")
                f.write("\n")

            else:
                # Benchmark results
                f.write("| Tool | Samples | Major Match | Accuracy | Confidence | Runtime |\n")
                f.write("|------|---------|-------------|----------|------------|--------|\n")
                for r in ds_results:
                    conf = f"{r.mean_confidence:.3f}" if r.mean_confidence else "N/A"
                    f.write(
                        f"| {r.tool} | {r.total_samples} | {r.correct_major}/{r.total_samples} | "
                        f"{r.accuracy:.1f}% | {conf} | {r.runtime_seconds:.1f}s |\n"
                    )
                f.write("\n")

        # Methods
        f.write("## Methods\n\n")
        f.write("### yallHap Modes\n\n")
        f.write("- **Heuristic**: Maximum-likelihood path traversal\n")
        f.write("- **Bayesian**: Probabilistic classification with allelic depth\n")
        f.write("- **Ancient**: Filters damage-like transitions (C>T, G>A)\n")
        f.write("- **Transversions-only**: Strictest mode for degraded DNA\n")
        f.write("- **ISOGG**: Maps YFull haplogroups to ISOGG nomenclature\n\n")

        f.write("### External Tools\n\n")
        f.write("- **yhaplo**: ISOGG-based classification\n")
        f.write("- **Yleaf**: YFull tree with `-aDNA` flag\n")
        f.write("- **pathPhynder**: Bayesian phylogenetic placement (BAM required)\n\n")

        f.write("### Ground Truth\n\n")
        f.write("- **AADR Terminal**: Y haplogroup from AADR v54 (terminal mutation format)\n")
        f.write("- **AADR ISOGG**: Y haplogroup from AADR v54 (ISOGG format)\n\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Comprehensive Y-chromosome haplogroup validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Output markdown file"
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Parallel threads (default: 1)"
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=None,
        help="Limit samples per dataset for quick testing",
    )
    parser.add_argument("--skip-1kg", action="store_true", help="Skip 1000 Genomes")
    parser.add_argument("--skip-gnomad", action="store_true", help="Skip gnomAD")
    parser.add_argument("--skip-aadr", action="store_true", help="Skip AADR")
    parser.add_argument("--skip-ancient", action="store_true", help="Skip local ancient")
    parser.add_argument("--skip-pathphynder", action="store_true", help="Skip pathPhynder")
    parser.add_argument("--skip-yhaplo", action="store_true", help="Skip yhaplo")
    parser.add_argument("--skip-yleaf", action="store_true", help="Skip Yleaf")
    parser.add_argument(
        "--ancient-dir",
        type=Path,
        default=None,
        help="Directory with ancient VCF/BAM files",
    )

    args = parser.parse_args()

    project_root = Path(__file__).parent.parent
    ancient_dir = args.ancient_dir or (project_root / "data" / "ancient_genomes")
    tree_path = project_root / "data" / "yfull_tree.json"
    snp_db_path = project_root / "data" / "validation" / "ybrowse_snps_hg19.csv"

    # ISOGG databases (GRCh37 and GRCh38)
    pathphynder_data = project_root.parent / "pathPhynder" / "data"
    isogg_db_grch37 = pathphynder_data / "211108.snps_isogg_curated.txt"
    isogg_db_grch38 = project_root / "data" / "isogg_snps_grch38.txt"

    # pathPhynder resources
    # Note: pathPhynder requires a "prepare" step to generate tree data files
    # Run from tree directory: pathPhynder -s prepare -i tree.nwk -f branches.snp.gz -p bigtree
    pathphynder_tree_dir = pathphynder_data / "BigTree_Y" / "older_versions"
    pathphynder_tree_filename = "bigtree_annotated_V1.nwk"
    pathphynder_tree = pathphynder_tree_dir / pathphynder_tree_filename
    # The prefix is relative to the tree directory (e.g., "tree_data/bigtree")
    pathphynder_prefix = "tree_data/bigtree"
    pathphynder_ref = pathphynder_data / "reference_sequences" / "hs37d5_Y.fa"
    pathphynder_isogg = isogg_db_grch37

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = args.output.parent / "validation_temp"
    temp_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Y-Chromosome Haplogroup Validation")
    print("=" * 70)
    if args.subsample:
        print(f"  Subsampling: {args.subsample} samples per dataset")
    if args.threads > 1:
        print(f"  Threads: {args.threads}")
    print()

    results: list[ValidationResult] = []

    # 1. 1000 Genomes
    if not args.skip_1kg:
        print("[1/4] 1000 Genomes Phase 3...")
        results.extend(benchmark_1kg(project_root, args.subsample, args.threads))
        print()

    # 2. gnomAD
    if not args.skip_gnomad:
        print("[2/4] gnomAD High-Coverage...")
        results.extend(benchmark_gnomad(project_root, args.subsample, args.threads))
        print()

    # 3. AADR
    if not args.skip_aadr:
        print("[3/4] AADR Ancient DNA...")
        results.extend(benchmark_aadr(project_root, args.subsample, args.threads))
        print()

    # AADR annotation file for ground truth
    aadr_anno_path = project_root / "data" / "ancient" / "v54.1.p1_HO_public.anno"

    # 4. Local ancient samples
    ancient_samples: list[SampleInfo] = []
    if not args.skip_ancient:
        print("[4/4] Local Ancient Samples...")
        if ancient_dir.exists():
            ancient_samples = load_ancient_samples(ancient_dir, args.subsample, aadr_anno_path)
            print(f"    {len(ancient_samples)} samples")

            if ancient_samples:
                # yallHap modes
                print("    yallHap regular...")
                results.append(
                    run_yallhap_ancient(
                        ancient_samples, tree_path, snp_db_path, "regular", args.threads
                    )
                )
                print("    yallHap bayesian...")
                results.append(
                    run_yallhap_ancient(
                        ancient_samples, tree_path, snp_db_path, "bayesian", args.threads
                    )
                )
                print("    yallHap transversions...")
                results.append(
                    run_yallhap_ancient(
                        ancient_samples, tree_path, snp_db_path, "transversions", args.threads
                    )
                )
                # ISOGG mode
                if isogg_db_grch37.exists():
                    print("    yallHap isogg...")
                    try:
                        results.append(
                            run_yallhap_ancient(
                                ancient_samples,
                                tree_path,
                                snp_db_path,
                                "isogg",
                                args.threads,
                                isogg_db_path=isogg_db_grch37,
                            )
                        )
                    except Exception as e:
                        print(f"    ERROR: yallHap isogg failed: {e}")
                        raise
                else:
                    print(f"    ERROR: ISOGG database not found: {isogg_db_grch37}")
                    print("           Expected pathPhynder ISOGG database for GRCh37")
                    raise FileNotFoundError(f"ISOGG database not found: {isogg_db_grch37}")

                # External tools
                if not args.skip_yhaplo:
                    print("    yhaplo...")
                    results.append(run_yhaplo_ancient(ancient_samples, temp_dir))

                if not args.skip_yleaf:
                    print("    Yleaf...")
                    results.append(run_yleaf_ancient(ancient_samples, temp_dir))

                if not args.skip_pathphynder:
                    missing_files = []
                    if not pathphynder_tree.exists():
                        missing_files.append(f"tree: {pathphynder_tree}")
                    if not pathphynder_ref.exists():
                        missing_files.append(f"reference: {pathphynder_ref}")

                    if missing_files:
                        print(f"    pathPhynder: missing files, skipping")
                        for mf in missing_files:
                            print(f"      - {mf}")
                    else:
                        print("    pathPhynder...")
                        results.append(
                            run_pathphynder_ancient(
                                ancient_samples,
                                temp_dir,
                                pathphynder_tree_dir,
                                pathphynder_tree_filename,
                                pathphynder_prefix,
                                pathphynder_ref,
                                pathphynder_isogg if pathphynder_isogg.exists() else None,
                            )
                        )
        else:
            print(f"    SKIP: {ancient_dir} not found")
        print()

    if not results:
        print("No benchmarks completed!")
        return 1

    # Generate report
    print("Generating report...")
    generate_markdown_report(results, args.output, ancient_samples)
    print(f"  Saved to: {args.output}")
    print()

    # Print summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Dataset':<20} {'Tool':<18} {'Samples':>8} {'Accuracy':>10} {'Time':>8}")
    print("-" * 70)
    for r in results:
        acc = f"{r.accuracy:.1f}%" if r.correct_major > 0 else "N/A"
        print(
            f"{r.dataset:<20} {r.tool:<18} {r.total_samples:>8} {acc:>10} {r.runtime_seconds:>7.1f}s"
        )
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
