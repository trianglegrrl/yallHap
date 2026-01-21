"""
Command-line interface for yallhap.

Provides pipeline-friendly CLI with proper exit codes and output formats.
"""

from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, TextIO, TypedDict

import click

from yallhap import __version__
from yallhap.classifier import HaplogroupCall, HaplogroupClassifier
from yallhap.snps import SNPDatabase
from yallhap.tree import Tree


class _DownloadFile(TypedDict):
    """Type definition for download file entries."""

    url: str
    path: Path
    label: str


# Module-level classifier for worker processes (initialized by _init_worker)
_worker_classifier: HaplogroupClassifier | None = None


def _is_bam_file(path: Path) -> bool:
    """
    Check if a file is a BAM file based on extension.

    Args:
        path: Path to check

    Returns:
        True if file has .bam extension (case-insensitive)
    """
    return path.suffix.lower() == ".bam"


def _load_snp_database(path: Path) -> SNPDatabase:
    """
    Load SNP database, auto-detecting format.

    Supports:
    - YBrowse GFF-style CSV (columns: seqid, source, type, start, ...)
    - Simple CSV format (columns: name, aliases, grch37_pos, grch38_pos, ...)

    Args:
        path: Path to SNP database file

    Returns:
        Populated SNPDatabase instance
    """
    import csv

    # Peek at first line to detect format
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader, [])

    # YBrowse GFF format has 'seqid' as first column
    if header and header[0].strip('"') == "seqid":
        return SNPDatabase.from_ybrowse_gff_csv(path)
    else:
        return SNPDatabase.from_csv(path)


def _init_worker(
    tree_path: Path,
    snp_db_path: Path,
    classifier_config: dict[str, Any],
) -> None:
    """
    Initialize classifier in worker process.

    Called once per worker process to load tree and SNP database,
    avoiding repeated loading for each file.

    Args:
        tree_path: Path to YFull tree JSON file
        snp_db_path: Path to SNP database CSV file
        classifier_config: Configuration dict for HaplogroupClassifier
    """
    global _worker_classifier
    tree = Tree.from_json(tree_path)
    snp_db = _load_snp_database(snp_db_path)
    _worker_classifier = HaplogroupClassifier(tree=tree, snp_db=snp_db, **classifier_config)


def _classify_file(vcf_path: Path) -> HaplogroupCall:
    """
    Classify a single VCF file using the worker's classifier.

    Must be called after _init_worker has initialized the classifier.

    Args:
        vcf_path: Path to VCF file

    Returns:
        HaplogroupCall with classification result

    Raises:
        RuntimeError: If worker classifier not initialized
    """
    if _worker_classifier is None:
        raise RuntimeError("Worker classifier not initialized")
    return _worker_classifier.classify(vcf_path)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name="yallhap")
def main() -> None:
    """
    yallHap: Modern Y-chromosome haplogroup inference.

    A pipeline-friendly tool for Y-chromosome haplogroup classification
    supporting modern and ancient DNA with probabilistic confidence scoring.
    """
    pass


@main.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
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
    "--min-base-quality",
    type=int,
    default=20,
    help="Minimum base quality for BAM pileup [default: 20]",
)
@click.option(
    "--min-mapping-quality",
    type=int,
    default=20,
    help="Minimum mapping quality for BAM pileup [default: 20]",
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
    "--bayesian",
    is_flag=True,
    help="Use Bayesian posterior calculation with allelic depth (AD) support",
)
@click.option(
    "--error-rate",
    type=float,
    default=0.001,
    help="Sequencing error rate for Bayesian mode [default: 0.001]",
)
@click.option(
    "--damage-rate",
    type=float,
    default=0.1,
    help="Ancient DNA damage rate for Bayesian mode [default: 0.1]",
)
@click.option(
    "--estimate-contamination",
    is_flag=True,
    help="Estimate Y-chromosome contamination from allele depths",
)
@click.option(
    "--max-tolerance",
    type=int,
    default=3,
    help="Maximum ancestral calls before stopping traversal in ancient mode [default: 3]",
)
@click.option(
    "--isogg",
    is_flag=True,
    help="Include ISOGG haplogroup name in output",
)
@click.option(
    "--isogg-db",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Path to ISOGG SNP database file (pathPhynder format)",
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
    input_file: Path,
    tree: Path,
    snp_db: Path,
    reference: str,
    sample: str | None,
    min_base_quality: int,
    min_mapping_quality: int,
    ancient: bool,
    transversions_only: bool,
    damage_rescale: str,
    min_depth: int,
    min_quality: int,
    bayesian: bool,
    error_rate: float,
    damage_rate: float,
    estimate_contamination: bool,
    max_tolerance: int,
    isogg: bool,
    isogg_db: Path | None,
    output: Path | None,
    output_format: str,
) -> None:
    """
    Classify Y-chromosome haplogroup from VCF or BAM.

    INPUT_FILE can be a VCF (indexed with .tbi/.csi) or BAM file (indexed with .bai).
    For BAM files, pileup is performed directly at known SNP positions.
    """
    try:
        # Auto-adjust filtering for ancient DNA with Bayesian mode
        # The Bayesian classifier uses allelic depth to model uncertainty,
        # so aggressive pre-filtering defeats its purpose
        if ancient and bayesian:
            if min_depth == 10:  # User didn't override default
                min_depth = 1
                click.echo(
                    "Note: Using min-depth=1 for ancient+bayesian mode "
                    "(override with --min-depth)",
                    err=True,
                )
            if min_quality == 20:  # User didn't override default
                min_quality = 0
                click.echo(
                    "Note: Using min-quality=0 for ancient+bayesian mode "
                    "(override with --min-quality)",
                    err=True,
                )

        # Load tree
        click.echo(f"Loading tree from {tree}...", err=True)
        tree_obj = Tree.from_json(tree)

        # Load SNP database (auto-detect format)
        click.echo(f"Loading SNP database from {snp_db}...", err=True)
        snp_db_obj = _load_snp_database(snp_db)

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
            bayesian=bayesian,
            error_rate=error_rate,
            damage_rate=damage_rate,
        )

        # Load ISOGG database if requested
        isogg_mapper = None
        if isogg:
            from yallhap.isogg import ISOGGDatabase, ISOGGMapper

            if isogg_db:
                click.echo(f"Loading ISOGG database from {isogg_db}...", err=True)
                isogg_db_obj = ISOGGDatabase.from_file(isogg_db)
            else:
                # Try default location
                default_isogg = Path("data/isogg_snps_grch38.txt")
                if default_isogg.exists():
                    click.echo(f"Loading ISOGG database from {default_isogg}...", err=True)
                    isogg_db_obj = ISOGGDatabase.from_file(default_isogg)
                else:
                    click.echo(
                        "Warning: --isogg requested but no ISOGG database found. "
                        "Use --isogg-db to specify path.",
                        err=True,
                    )
                    isogg_db_obj = None

            if isogg_db_obj:
                isogg_mapper = ISOGGMapper(tree_obj, isogg_db_obj)

        # Run classification - detect input type
        is_bam = _is_bam_file(input_file)
        click.echo(f"Classifying {input_file}{'(BAM)' if is_bam else ''}...", err=True)

        if is_bam:
            result = classifier.classify_from_bam(
                input_file,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
            )
        else:
            result = classifier.classify(input_file, sample)

        # Estimate contamination if requested
        contamination_result = None
        if estimate_contamination:
            from yallhap.contamination import ContaminationResult, estimate_contamination_with_snpdb

            if is_bam:
                click.echo(
                    "Warning: Contamination estimation from BAM not yet supported, skipping.",
                    err=True,
                )
            else:
                # Get variants from classifier
                variants = classifier._get_variants(input_file, sample)  # type: ignore
                rate, n_sites = estimate_contamination_with_snpdb(
                    variants, snp_db_obj, reference=reference, min_depth=min_depth  # type: ignore
                )
                contamination_result = ContaminationResult(rate=rate, n_sites=n_sites)
                click.echo(f"Contamination estimate: {rate:.1%} ({n_sites} sites)", err=True)

        # Add ISOGG haplogroup if available
        isogg_haplogroup = None
        if isogg_mapper and result.haplogroup != "NA":
            isogg_haplogroup = isogg_mapper.to_isogg(result.haplogroup)

        # Store extra results in result dict for output
        # (These will be included when writing result)
        extra_data = {
            "isogg_haplogroup": isogg_haplogroup,
            "contamination_rate": contamination_result.rate if contamination_result else None,
            "contamination_sites": contamination_result.n_sites if contamination_result else None,
            "max_tolerance": max_tolerance,
        }

        # Output result
        if output:
            with open(output, "w") as f:
                _write_result(result, f, output_format, extra_data)
            click.echo(f"Result written to {output}", err=True)
        else:
            _write_result(result, sys.stdout, output_format, extra_data)

        # Exit with appropriate code
        # Only use non-zero exit for actual failures, not informational states
        # Low confidence is informational - users should check the output
        if result.haplogroup == "NA":
            sys.exit(1)  # Classification failed
        else:
            sys.exit(0)  # Success (confidence level is in the output)

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


def _write_result(
    result: HaplogroupCall,
    file: TextIO,
    format: str,
    extra_data: dict[str, Any] | None = None,
) -> None:
    """Write classification result to file."""
    if format == "json":
        output_dict = result.to_dict()
        # Add extra data if provided
        if extra_data:
            for key, value in extra_data.items():
                if value is not None:
                    output_dict[key] = value
        json.dump(output_dict, file, indent=2)
        file.write("\n")
    elif format == "tsv":
        # Header - include Bayesian fields if present
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

        # Add Bayesian columns if present
        if result.posterior is not None:
            header.extend(["posterior", "credible_set_95"])

        # Add ISOGG column if present
        isogg_hg = extra_data.get("isogg_haplogroup") if extra_data else None
        if isogg_hg is not None:
            header.append("isogg_haplogroup")

        # Add contamination columns if present
        if extra_data and extra_data.get("contamination_rate") is not None:
            header.extend(["contamination_rate", "contamination_sites"])

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

        # Add Bayesian values if present
        if result.posterior is not None:
            row.extend(
                [
                    f"{result.posterior:.4f}",
                    ";".join(result.credible_set_95),
                ]
            )

        # Add ISOGG value if present
        if isogg_hg is not None:
            row.append(isogg_hg)

        # Add contamination values if present
        if extra_data and extra_data.get("contamination_rate") is not None:
            row.extend(
                [
                    f"{extra_data['contamination_rate']:.4f}",
                    str(extra_data["contamination_sites"]),
                ]
            )

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
    "--bayesian",
    is_flag=True,
    help="Use Bayesian posterior calculation with allelic depth (AD) support",
)
@click.option(
    "--error-rate",
    type=float,
    default=0.001,
    help="Sequencing error rate for Bayesian mode [default: 0.001]",
)
@click.option(
    "--damage-rate",
    type=float,
    default=0.1,
    help="Ancient DNA damage rate for Bayesian mode [default: 0.1]",
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
    bayesian: bool,
    error_rate: float,
    damage_rate: float,
    output: Path,
    threads: int,
) -> None:
    """
    Batch classify multiple VCF files.

    Writes results to a single TSV file with one row per sample.
    """

    try:
        results: list[HaplogroupCall] = []

        # Process files
        click.echo(f"Processing {len(vcf_files)} files...", err=True)

        # Build classifier configuration dict for workers
        classifier_config: dict[str, Any] = {
            "reference": reference,
            "ancient_mode": ancient,
            "transversions_only": transversions_only,
            "damage_rescale": damage_rescale,
            "bayesian": bayesian,
            "error_rate": error_rate,
            "damage_rate": damage_rate,
        }

        if threads == 1:
            # Sequential processing - load tree/snp_db once in main process
            tree_obj = Tree.from_json(tree)
            snp_db_obj = SNPDatabase.from_csv(snp_db)

            # Handle T2T reference: need to liftover positions
            if reference == "t2t":
                snp_db_obj = _prepare_t2t_database(snp_db_obj)
                classifier_config["reference"] = "t2t"

            classifier = HaplogroupClassifier(
                tree=tree_obj,
                snp_db=snp_db_obj,
                **classifier_config,
            )

            for vcf in vcf_files:
                result = classifier.classify(vcf)
                results.append(result)
                click.echo(f"  {vcf.name}: {result.haplogroup}", err=True)
        else:
            # Parallel processing using ProcessPoolExecutor
            # Note: T2T liftover happens in each worker since we can't pickle the db
            if reference == "t2t":
                click.echo(
                    "Warning: T2T reference with parallel processing requires liftover "
                    "in each worker. Consider using --threads 1 for T2T.",
                    err=True,
                )

            # Cap threads to number of files
            actual_threads = min(threads, len(vcf_files))
            click.echo(f"Using {actual_threads} parallel workers...", err=True)

            with ProcessPoolExecutor(
                max_workers=actual_threads,
                initializer=_init_worker,
                initargs=(tree, snp_db, classifier_config),
            ) as executor:
                # Submit all tasks
                future_to_vcf = {executor.submit(_classify_file, vcf): vcf for vcf in vcf_files}

                # Collect results as they complete
                for future in as_completed(future_to_vcf):
                    vcf = future_to_vcf[future]
                    result = future.result()
                    results.append(result)
                    click.echo(f"  {vcf.name}: {result.haplogroup}", err=True)

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
    Download YFull tree and SNP databases for all reference genomes.

    Downloads:
    - YFull tree from GitHub
    - YBrowse SNP database (GRCh38/hg38 positions)
    - YBrowse SNP database (GRCh37/hg19 positions)
    - ISOGG SNP database (GRCh37 - for --isogg mode)
    - ISOGG SNP database (GRCh38 - for --isogg mode)

    Use the appropriate SNP database for your VCF reference:
    - GRCh37/hg19 VCFs: -s ybrowse_snps_grch37.csv -r grch37
    - GRCh38/hg38 VCFs: -s ybrowse_snps_grch38.csv -r grch38
    - T2T VCFs:         -s ybrowse_snps_grch38.csv -r t2t

    For ISOGG mode, add --isogg and the database is auto-detected.

    Skips files that already exist unless --force is specified.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define files to download
    files_to_download: list[_DownloadFile] = [
        {
            "url": "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json",
            "path": output_dir / "yfull_tree.json",
            "label": "YFull tree",
        },
        {
            "url": "http://ybrowse.org/gbrowse2/gff/snps_hg38.csv",
            "path": output_dir / "ybrowse_snps_grch38.csv",
            "label": "YBrowse SNPs (GRCh38)",
        },
        {
            "url": "http://ybrowse.org/gbrowse2/gff/snps_hg19.csv",
            "path": output_dir / "ybrowse_snps_grch37.csv",
            "label": "YBrowse SNPs (GRCh37)",
        },
        # ISOGG SNP databases from pathPhynder (for --isogg mode)
        {
            "url": "https://raw.githubusercontent.com/ruidlpm/pathPhynder/master/data/211108.snps_isogg_curated.txt",
            "path": output_dir / "isogg_snps_grch37.txt",
            "label": "ISOGG SNPs (GRCh37)",
        },
        {
            "url": "https://raw.githubusercontent.com/ruidlpm/pathPhynder/master/data/211108.snps_isogg_curated.b38.txt",
            "path": output_dir / "isogg_snps_grch38.txt",
            "label": "ISOGG SNPs (GRCh38)",
        },
    ]

    # Check which files already exist
    all_exist = all(f["path"].exists() for f in files_to_download)

    if all_exist and not force:
        click.echo(f"All files already exist in {output_dir}/", err=True)
        for f in files_to_download:
            click.echo(f"  - {f['path'].name}", err=True)
        click.echo("Use --force to re-download.", err=True)
        sys.exit(0)

    # Download each file
    failed = False
    for f in files_to_download:
        if f["path"].exists() and not force:
            click.echo(f"Skipping {f['path'].name} (exists)", err=True)
        else:
            try:
                _download_with_progress(f["url"], f["path"], f["label"])
                click.echo(f"  Saved to {f['path']}", err=True)
            except Exception as e:
                click.echo(f"  Failed to download {f['label']}: {e}", err=True)
                failed = True

    if failed:
        click.echo("\nSome downloads failed. Re-run with --force to retry.", err=True)
        sys.exit(1)

    click.echo("\nDownload complete!", err=True)
    click.echo("\nUsage:", err=True)
    click.echo("  GRCh37/hg19 VCFs: -s ybrowse_snps_grch37.csv -r grch37", err=True)
    click.echo("  GRCh38/hg38 VCFs: -s ybrowse_snps_grch38.csv -r grch38", err=True)
    click.echo("  T2T VCFs:         -s ybrowse_snps_grch38.csv -r t2t", err=True)
    click.echo("\nISOGG mode (for standardized haplogroup nomenclature):", err=True)
    click.echo("  Add --isogg to classification commands", err=True)
    click.echo("  ISOGG databases: isogg_snps_grch37.txt, isogg_snps_grch38.txt", err=True)


if __name__ == "__main__":
    main()
