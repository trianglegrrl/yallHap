"""
Command-line interface for yclade.

Provides pipeline-friendly CLI with proper exit codes and output formats.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import TextIO

import click

from yclade import __version__
from yclade.classifier import HaplogroupCall, HaplogroupClassifier
from yclade.snps import SNPDatabase
from yclade.tree import Tree


@click.group()
@click.version_option(version=__version__, prog_name="yclade")
def main() -> None:
    """
    YClade: Modern Y-chromosome haplogroup inference.

    A pipeline-friendly tool for Y-chromosome haplogroup classification
    supporting modern and ancient DNA with probabilistic confidence scoring.
    """
    pass


@main.command()
@click.argument("vcf", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--tree",
    "-t",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to YFull tree JSON file",
)
@click.option(
    "--snp-db",
    "-s",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to SNP database CSV file",
)
@click.option(
    "--reference",
    "-r",
    type=click.Choice(["grch37", "grch38", "t2t"]),
    default="grch38",
    help="Reference genome [default: grch38]",
)
@click.option(
    "--sample",
    type=str,
    default=None,
    help="Sample name for multi-sample VCF [default: first sample]",
)
@click.option(
    "--ancient",
    is_flag=True,
    help="Enable ancient DNA mode (filter damage-like transitions)",
)
@click.option(
    "--min-depth",
    type=int,
    default=10,
    help="Minimum read depth [default: 10]",
)
@click.option(
    "--min-quality",
    type=int,
    default=20,
    help="Minimum genotype quality [default: 20]",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output file [default: stdout]",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["json", "tsv"]),
    default="json",
    help="Output format [default: json]",
)
def classify(
    vcf: Path,
    tree: Path,
    snp_db: Path,
    reference: str,
    sample: str | None,
    ancient: bool,
    min_depth: int,
    min_quality: int,
    output: Path | None,
    output_format: str,
) -> None:
    """
    Classify Y-chromosome haplogroup from VCF.

    VCF must contain Y-chromosome variants and be indexed (.tbi or .csi).
    """
    try:
        # Load tree
        click.echo(f"Loading tree from {tree}...", err=True)
        tree_obj = Tree.from_json(tree)

        # Load SNP database
        click.echo(f"Loading SNP database from {snp_db}...", err=True)
        snp_db_obj = SNPDatabase.from_csv(snp_db)

        # Create classifier
        classifier = HaplogroupClassifier(
            tree=tree_obj,
            snp_db=snp_db_obj,
            reference=reference,  # type: ignore
            min_depth=min_depth,
            min_quality=min_quality,
            ancient_mode=ancient,
        )

        # Run classification
        click.echo(f"Classifying {vcf}...", err=True)
        result = classifier.classify(vcf, sample)

        # Output result
        if output:
            with open(output, "w") as f:
                _write_result(result, f, output_format)
            click.echo(f"Result written to {output}", err=True)
        else:
            _write_result(result, sys.stdout, output_format)

        # Exit with appropriate code
        if result.haplogroup == "NA":
            sys.exit(1)  # Classification failed
        elif result.confidence < 0.95:
            sys.exit(2)  # Low confidence
        else:
            sys.exit(0)  # Success

    except FileNotFoundError as e:
        click.echo(f"Error: File not found: {e}", err=True)
        sys.exit(10)
    except ValueError as e:
        click.echo(f"Error: Invalid input: {e}", err=True)
        sys.exit(11)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(99)


def _write_result(result: HaplogroupCall, file: TextIO, format: str) -> None:
    """Write classification result to file."""
    if format == "json":
        json.dump(result.to_dict(), file, indent=2)
        file.write("\n")
    elif format == "tsv":
        # Header
        header = [
            "sample",
            "haplogroup",
            "confidence",
            "qc1",
            "qc2",
            "qc3",
            "qc4",
            "derived",
            "ancestral",
            "missing",
        ]
        file.write("\t".join(header) + "\n")

        # Data
        row = [
            result.sample,
            result.haplogroup,
            f"{result.confidence:.4f}",
            f"{result.qc_scores.qc1_backbone:.4f}",
            f"{result.qc_scores.qc2_terminal:.4f}",
            f"{result.qc_scores.qc3_path:.4f}",
            f"{result.qc_scores.qc4_posterior:.4f}",
            str(result.snp_stats.derived),
            str(result.snp_stats.ancestral),
            str(result.snp_stats.missing),
        ]
        file.write("\t".join(row) + "\n")


@main.command()
@click.argument("vcf_files", nargs=-1, type=click.Path(exists=True, path_type=Path))
@click.option(
    "--tree",
    "-t",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to YFull tree JSON file",
)
@click.option(
    "--snp-db",
    "-s",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to SNP database CSV file",
)
@click.option(
    "--reference",
    "-r",
    type=click.Choice(["grch37", "grch38", "t2t"]),
    default="grch38",
    help="Reference genome [default: grch38]",
)
@click.option(
    "--ancient",
    is_flag=True,
    help="Enable ancient DNA mode",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Output TSV file",
)
@click.option(
    "--threads",
    type=int,
    default=1,
    help="Number of parallel threads [default: 1]",
)
def batch(
    vcf_files: tuple[Path, ...],
    tree: Path,
    snp_db: Path,
    reference: str,
    ancient: bool,
    output: Path,
    threads: int,
) -> None:
    """
    Batch classify multiple VCF files.

    Writes results to a single TSV file with one row per sample.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    try:
        # Load tree and SNP database
        tree_obj = Tree.from_json(tree)
        snp_db_obj = SNPDatabase.from_csv(snp_db)

        classifier = HaplogroupClassifier(
            tree=tree_obj,
            snp_db=snp_db_obj,
            reference=reference,  # type: ignore
            ancient_mode=ancient,
        )

        results: list[HaplogroupCall] = []

        # Process files
        click.echo(f"Processing {len(vcf_files)} files...", err=True)

        if threads == 1:
            for vcf in vcf_files:
                result = classifier.classify(vcf)
                results.append(result)
                click.echo(f"  {vcf.name}: {result.haplogroup}", err=True)
        else:
            # Parallel processing not yet implemented
            click.echo("Warning: Parallel processing not yet implemented, using single thread", err=True)
            for vcf in vcf_files:
                result = classifier.classify(vcf)
                results.append(result)

        # Write output
        with open(output, "w") as f:
            header = [
                "sample",
                "haplogroup",
                "confidence",
                "qc1",
                "qc2",
                "qc3",
                "qc4",
                "derived",
                "ancestral",
                "missing",
            ]
            f.write("\t".join(header) + "\n")

            for result in results:
                row = [
                    result.sample,
                    result.haplogroup,
                    f"{result.confidence:.4f}",
                    f"{result.qc_scores.qc1_backbone:.4f}",
                    f"{result.qc_scores.qc2_terminal:.4f}",
                    f"{result.qc_scores.qc3_path:.4f}",
                    f"{result.qc_scores.qc4_posterior:.4f}",
                    str(result.snp_stats.derived),
                    str(result.snp_stats.ancestral),
                    str(result.snp_stats.missing),
                ]
                f.write("\t".join(row) + "\n")

        click.echo(f"Results written to {output}", err=True)

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(99)


@main.command()
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("data"),
    help="Output directory for downloaded files [default: ./data]",
)
def download(output_dir: Path) -> None:
    """
    Download YFull tree and SNP database.

    Downloads the latest YFull tree from GitHub and YBrowse SNP database.
    """
    import requests

    output_dir.mkdir(parents=True, exist_ok=True)

    # Download YFull tree
    click.echo("Downloading YFull tree...", err=True)
    tree_url = "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json"

    try:
        response = requests.get(tree_url, timeout=60)
        response.raise_for_status()
        tree_path = output_dir / "yfull_tree.json"
        with open(tree_path, "wb") as f:
            f.write(response.content)
        click.echo(f"  Saved to {tree_path}", err=True)
    except requests.RequestException as e:
        click.echo(f"  Failed to download YFull tree: {e}", err=True)
        sys.exit(1)

    # Download YBrowse SNP database
    click.echo("Downloading YBrowse SNP database...", err=True)
    snp_url = "http://ybrowse.org/gbrowse2/gff/snps_hg38.csv"

    try:
        response = requests.get(snp_url, timeout=120)
        response.raise_for_status()
        snp_path = output_dir / "ybrowse_snps.csv"
        with open(snp_path, "wb") as f:
            f.write(response.content)
        click.echo(f"  Saved to {snp_path}", err=True)
    except requests.RequestException as e:
        click.echo(f"  Failed to download SNP database: {e}", err=True)
        sys.exit(1)

    click.echo("Download complete!", err=True)


if __name__ == "__main__":
    main()
