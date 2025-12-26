#!/usr/bin/env python3
"""
Download liftover chain files for T2T coordinate conversion.

Downloads official T2T consortium chain files from their public S3 bucket.
These are used to convert GRCh37/GRCh38 coordinates to T2T-CHM13v2.0.
"""

from __future__ import annotations

import gzip
import shutil
import sys
from pathlib import Path

import click
import requests


# T2T consortium chain file URLs (from s3://human-pangenomics)
CHAIN_FILES = {
    "grch38_to_t2t": {
        "url": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain",
        "filename": "grch38-chm13v2.chain",
    },
    "hg19_to_t2t": {
        "url": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-hg19.chain",
        "filename": "hg19-chm13v2.chain",
    },
}


def download_file(url: str, output_path: Path, label: str) -> None:
    """Download a file with progress indication."""
    click.echo(f"Downloading {label}...", err=True)

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
                    click.echo(f"\r  {pct:.1f}% ({downloaded}/{total_bytes} bytes)", nl=False, err=True)

    click.echo(f"\n  Saved to {output_path}", err=True)


@click.command()
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("data/liftover"),
    help="Output directory for chain files [default: ./data/liftover]",
)
@click.option(
    "--force",
    "-f",
    is_flag=True,
    help="Force re-download even if files exist",
)
def main(output_dir: Path, force: bool) -> None:
    """
    Download T2T liftover chain files.

    Downloads official T2T consortium chain files for coordinate conversion:
    - grch38-chm13v2.chain: GRCh38 to T2T-CHM13v2.0
    - hg19-chm13v2.chain: hg19/GRCh37 to T2T-CHM13v2.0
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    for chain_type, info in CHAIN_FILES.items():
        output_path = output_dir / info["filename"]

        if output_path.exists() and not force:
            click.echo(f"Skipping {info['filename']} (exists)", err=True)
            continue

        try:
            download_file(info["url"], output_path, info["filename"])
        except requests.RequestException as e:
            click.echo(f"Error downloading {chain_type}: {e}", err=True)
            sys.exit(1)

    click.echo("Download complete!", err=True)


if __name__ == "__main__":
    main()

