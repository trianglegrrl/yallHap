"""
Command-line interface for yallhap.

Provides pipeline-friendly CLI with proper exit codes and output formats.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import TextIO

import click

from yallhap import __version__
from yallhap.classifier import HaplogroupCall, HaplogroupClassifier
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


@click.group()
@click.version_option(version=__version__, prog_name="yallhap")
def main() -> None:
    """
    yallHap: Modern Y-chromosome haplogroup inference.

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
    help="Enable ancient DNA mode (filter C>T and G>A damage-like transitions)",
)
@click.option(
    "--transversions-only",
    is_flag=True,
    help="Only use transversions (strictest ancient DNA mode, ignores all transitions)",
)
@click.option(
    "--damage-rescale",
    type=click.Choice(["none", "moderate", "aggressive"]),
    default="none",
    help="Rescale quality scores for potentially damaged variants [default: none]",
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
    transversions_only: bool,
    damage_rescale: str,
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

        # Handle T2T reference: need to liftover positions
        if reference == "t2t":
            snp_db_obj = _prepare_t2t_database(snp_db_obj)

        # Create classifier
        classifier = HaplogroupClassifier(
            tree=tree_obj,
            snp_db=snp_db_obj,
            reference=reference,  # type: ignore
            min_depth=min_depth,
            min_quality=min_quality,
            ancient_mode=ancient,
            transversions_only=transversions_only,
            damage_rescale=damage_rescale,  # type: ignore
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


def _prepare_t2t_database(snp_db: SNPDatabase) -> SNPDatabase:
    """
    Prepare SNP database for T2T reference by computing liftover positions.

    Looks for chain files in data/liftover/ and computes T2T positions.

    Args:
        snp_db: SNP database to enhance with T2T positions

    Returns:
        Same database with T2T positions populated

    Raises:
        FileNotFoundError: If chain files are not found
    """
    # Look for chain files in standard locations
    chain_dirs = [
        Path("data/liftover"),
        Path(__file__).parent.parent.parent / "data" / "liftover",
    ]

    grch38_chain = None
    grch37_chain = None

    for chain_dir in chain_dirs:
        if (chain_dir / "grch38-chm13v2.chain").exists():
            grch38_chain = chain_dir / "grch38-chm13v2.chain"
        if (chain_dir / "hg19-chm13v2.chain").exists():
            grch37_chain = chain_dir / "hg19-chm13v2.chain"

    if grch38_chain is None and grch37_chain is None:
        raise FileNotFoundError(
            "T2T liftover chain files not found. "
            "Run 'python scripts/download_liftover_chains.py' to download them."
        )

    click.echo("Computing T2T positions via liftover...", err=True)
    stats = snp_db.compute_all_t2t_positions(grch38_chain, grch37_chain)
    click.echo(
        f"  Lifted {stats['total_with_t2t']} SNPs to T2T "
        f"(GRCh38: {stats['lifted_from_grch38']}, GRCh37: {stats['lifted_from_grch37']})",
        err=True,
    )

    return snp_db


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
    help="Enable ancient DNA mode (filter C>T and G>A damage-like transitions)",
)
@click.option(
    "--transversions-only",
    is_flag=True,
    help="Only use transversions (strictest ancient DNA mode)",
)
@click.option(
    "--damage-rescale",
    type=click.Choice(["none", "moderate", "aggressive"]),
    default="none",
    help="Rescale quality scores for potentially damaged variants [default: none]",
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
    transversions_only: bool,
    damage_rescale: str,
    output: Path,
    threads: int,
) -> None:
    """
    Batch classify multiple VCF files.

    Writes results to a single TSV file with one row per sample.
    """

    try:
        # Load tree and SNP database
        tree_obj = Tree.from_json(tree)
        snp_db_obj = SNPDatabase.from_csv(snp_db)

        # Handle T2T reference: need to liftover positions
        if reference == "t2t":
            snp_db_obj = _prepare_t2t_database(snp_db_obj)

        classifier = HaplogroupClassifier(
            tree=tree_obj,
            snp_db=snp_db_obj,
            reference=reference,  # type: ignore
            ancient_mode=ancient,
            transversions_only=transversions_only,
            damage_rescale=damage_rescale,  # type: ignore
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
            click.echo(
                "Warning: Parallel processing not yet implemented, using single thread", err=True
            )
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


def _download_with_progress(
    url: str,
    output_path: Path,
    label: str,
    timeout: int = 300,
) -> None:
    """
    Download a file with streaming and progress bar.

    Args:
        url: URL to download
        output_path: Path to save the file
        label: Label to show in progress bar
        timeout: Request timeout in seconds

    Raises:
        requests.RequestException: If download fails
    """
    import requests

    response = requests.get(url, stream=True, timeout=timeout)
    response.raise_for_status()

    # Get file size if available
    total_size = response.headers.get("content-length")
    total_bytes = int(total_size) if total_size else None

    with (
        open(output_path, "wb") as f,
        click.progressbar(  # type: ignore[var-annotated]
            length=total_bytes,
            label=label,
            show_eta=True,
            show_percent=True,
            file=sys.stderr,
        ) as bar,
    ):
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                bar.update(len(chunk))


@main.command()
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("data"),
    help="Output directory for downloaded files [default: ./data]",
)
@click.option(
    "--force",
    "-f",
    is_flag=True,
    help="Force re-download even if files exist",
)
def download(output_dir: Path, force: bool) -> None:
    """
    Download YFull tree and SNP database.

    Downloads the latest YFull tree from GitHub and YBrowse SNP database.
    Skips files that already exist unless --force is specified.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define files to download
    tree_url = "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json"
    tree_path = output_dir / "yfull_tree.json"
    snp_url = "http://ybrowse.org/gbrowse2/gff/snps_hg38.csv"
    snp_path = output_dir / "ybrowse_snps.csv"

    # Check existing files
    tree_exists = tree_path.exists()
    snp_exists = snp_path.exists()

    if tree_exists and snp_exists and not force:
        click.echo(f"Files already exist in {output_dir}/", err=True)
        click.echo(f"  - {tree_path.name}", err=True)
        click.echo(f"  - {snp_path.name}", err=True)
        click.echo("Use --force to re-download.", err=True)
        sys.exit(0)

    # Download YFull tree
    if tree_exists and not force:
        click.echo(f"Skipping {tree_path.name} (exists)", err=True)
    else:
        try:
            _download_with_progress(tree_url, tree_path, "YFull tree")
            click.echo(f"  Saved to {tree_path}", err=True)
        except Exception as e:
            click.echo(f"  Failed to download YFull tree: {e}", err=True)
            sys.exit(1)

    # Download YBrowse SNP database
    if snp_exists and not force:
        click.echo(f"Skipping {snp_path.name} (exists)", err=True)
    else:
        try:
            _download_with_progress(snp_url, snp_path, "YBrowse SNPs")
            click.echo(f"  Saved to {snp_path}", err=True)
        except Exception as e:
            click.echo(f"  Failed to download SNP database: {e}", err=True)
            sys.exit(1)

    click.echo("Download complete!", err=True)


if __name__ == "__main__":
    main()
