#!/usr/bin/env python3
"""
Download Allen Ancient DNA Resource (AADR) data for validation.

The AADR is curated by the Reich Lab and contains genome-wide data from
thousands of ancient individuals with published Y-chromosome haplogroups.

Data is hosted on Harvard Dataverse:
https://dataverse.harvard.edu/dataverse/reich_lab

This script downloads:
1. AADR annotation file (contains haplogroup assignments)
2. Optionally: genotype data in EIGENSTRAT or PLINK format

Note: The genotype files are large (~2GB). For Y-chromosome validation,
we primarily need the annotation file which contains haplogroup calls.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import click
import requests

# AADR v54.1 (December 2023) - latest release
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW
AADR_FILES = {
    # Annotation file with haplogroup assignments
    "anno": {
        "url": "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1/SHARE/public.dir/v54.1_1240K_public.anno",
        "filename": "aadr_v54.1_anno.tsv",
        "description": "Annotation file with Y-haplogroup assignments",
        "size_mb": 15,
    },
    # Genotype data in EIGENSTRAT format (optional, large)
    "geno": {
        "url": "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1/SHARE/public.dir/v54.1_1240K_public.geno",
        "filename": "aadr_v54.1.geno",
        "description": "Genotype file (EIGENSTRAT format)",
        "size_mb": 1800,
    },
    "snp": {
        "url": "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1/SHARE/public.dir/v54.1_1240K_public.snp",
        "filename": "aadr_v54.1.snp",
        "description": "SNP file (EIGENSTRAT format)",
        "size_mb": 50,
    },
    "ind": {
        "url": "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1/SHARE/public.dir/v54.1_1240K_public.ind",
        "filename": "aadr_v54.1.ind",
        "description": "Individual file (EIGENSTRAT format)",
        "size_mb": 1,
    },
}


def download_file(url: str, output_path: Path, label: str) -> bool:
    """Download a file with progress indication."""
    click.echo(f"Downloading {label}...", err=True)
    click.echo(f"  URL: {url}", err=True)

    try:
        response = requests.get(url, stream=True, timeout=300)
        response.raise_for_status()

        total_size = response.headers.get("content-length")
        total_bytes = int(total_size) if total_size else None

        with open(output_path, "wb") as f:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_bytes:
                        pct = (downloaded / total_bytes) * 100
                        mb = downloaded / (1024 * 1024)
                        click.echo(
                            f"\r  {pct:.1f}% ({mb:.1f} MB)", nl=False, err=True
                        )

        click.echo(f"\n  Saved to {output_path}", err=True)
        return True

    except requests.RequestException as e:
        click.echo(f"\n  Error: {e}", err=True)
        return False


def extract_y_haplogroups(anno_path: Path, output_path: Path) -> int:
    """
    Extract Y-chromosome haplogroup data from AADR annotation file.

    Returns number of samples with Y-haplogroup assignments.
    """
    click.echo(f"Extracting Y-haplogroups from {anno_path}...", err=True)

    count = 0
    with open(anno_path, "r") as f_in, open(output_path, "w") as f_out:
        # Write header
        f_out.write("sample_id\thaplogroup\tpublication\tdate_bp\tlocation\tsex\n")

        for line_num, line in enumerate(f_in):
            if line_num == 0:
                # Parse header to find column indices
                headers = line.strip().split("\t")
                try:
                    idx_id = headers.index("Genetic_ID") if "Genetic_ID" in headers else headers.index("Instance_ID")
                    idx_yhg = headers.index("Y_haplogroup") if "Y_haplogroup" in headers else -1
                    idx_pub = headers.index("Publication") if "Publication" in headers else -1
                    idx_date = headers.index("Date_mean") if "Date_mean" in headers else -1
                    idx_loc = headers.index("Locality") if "Locality" in headers else -1
                    idx_sex = headers.index("Sex") if "Sex" in headers else -1
                except ValueError as e:
                    click.echo(f"  Warning: Could not find column: {e}", err=True)
                    # Try alternative column names
                    idx_id = 0
                    idx_yhg = -1
                    for i, h in enumerate(headers):
                        if "Y" in h.upper() and "HAPLO" in h.upper():
                            idx_yhg = i
                            break
                continue

            if idx_yhg == -1:
                click.echo("  Error: No Y-haplogroup column found", err=True)
                break

            cols = line.strip().split("\t")
            if len(cols) <= max(idx_id, idx_yhg):
                continue

            sample_id = cols[idx_id] if idx_id >= 0 else ""
            y_hg = cols[idx_yhg] if idx_yhg >= 0 and idx_yhg < len(cols) else ""
            pub = cols[idx_pub] if idx_pub >= 0 and idx_pub < len(cols) else ""
            date = cols[idx_date] if idx_date >= 0 and idx_date < len(cols) else ""
            loc = cols[idx_loc] if idx_loc >= 0 and idx_loc < len(cols) else ""
            sex = cols[idx_sex] if idx_sex >= 0 and idx_sex < len(cols) else ""

            # Only include samples with Y-haplogroup (males)
            if y_hg and y_hg not in ("n/a", "NA", ".."):
                f_out.write(f"{sample_id}\t{y_hg}\t{pub}\t{date}\t{loc}\t{sex}\n")
                count += 1

    click.echo(f"  Extracted {count} samples with Y-haplogroups", err=True)
    return count


@click.command()
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("data/ancient"),
    help="Output directory [default: ./data/ancient]",
)
@click.option(
    "--full",
    is_flag=True,
    help="Download full genotype data (large, ~2GB)",
)
@click.option(
    "--force",
    "-f",
    is_flag=True,
    help="Force re-download even if files exist",
)
def main(output_dir: Path, full: bool, force: bool) -> None:
    """
    Download AADR (Allen Ancient DNA Resource) data.

    By default, downloads only the annotation file (~15MB) which contains
    Y-chromosome haplogroup assignments for thousands of ancient samples.

    Use --full to also download the genotype data (~2GB) for running
    haplogroup classification directly.

    AADR v54.1 (December 2023):
    https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Always download annotation file
    files_to_download = ["anno"]

    if full:
        files_to_download.extend(["geno", "snp", "ind"])
        click.echo("Downloading full AADR dataset (this will take a while)...", err=True)
    else:
        click.echo("Downloading AADR annotation file...", err=True)
        click.echo("Use --full to also download genotype data (~2GB)", err=True)

    success = True
    for file_key in files_to_download:
        file_info = AADR_FILES[file_key]
        output_path = output_dir / file_info["filename"]

        if output_path.exists() and not force:
            click.echo(f"Skipping {file_info['filename']} (exists)", err=True)
            continue

        click.echo(f"\n{file_info['description']} (~{file_info['size_mb']} MB)", err=True)
        if not download_file(file_info["url"], output_path, file_info["filename"]):
            success = False
            continue

    # Extract Y-haplogroups from annotation file
    anno_path = output_dir / AADR_FILES["anno"]["filename"]
    if anno_path.exists():
        y_hg_path = output_dir / "aadr_y_haplogroups.tsv"
        extract_y_haplogroups(anno_path, y_hg_path)

    if success:
        click.echo("\nDownload complete!", err=True)
        click.echo(f"\nFiles saved to: {output_dir}/", err=True)
        click.echo("\nTo validate yallHap against AADR:", err=True)
        click.echo("  1. The annotation file contains Y-haplogroup ground truth", err=True)
        click.echo("  2. For classification, you need VCF files which require", err=True)
        click.echo("     converting EIGENSTRAT format or accessing BAM files", err=True)
        click.echo("     from the European Nucleotide Archive (ENA)", err=True)
    else:
        click.echo("\nSome downloads failed. Check errors above.", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

