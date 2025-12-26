#!/usr/bin/env python3
"""
Download public ancient DNA samples with known Y-chromosome haplogroups.

Ancient DNA samples with published haplogroup calls are used to validate
yallHap's ancient DNA mode against established ground truth.

Data Sources:
- Allen Ancient DNA Resource (AADR): https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
- European Nucleotide Archive (ENA) for raw sequence data
- Mathieson et al. papers for haplogroup annotations

Note: Ancient DNA VCF data requires special access or processing from BAM files.
This script provides the framework for validation once data is obtained.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import click


@dataclass
class AncientSample:
    """An ancient DNA sample with known haplogroup."""

    sample_id: str
    haplogroup: str
    publication: str
    date_bp: int | None = None  # Years before present
    location: str = ""
    notes: str = ""


# Curated list of well-known ancient samples with Y-chromosome haplogroups
# From various publications using ancient DNA
KNOWN_ANCIENT_SAMPLES = [
    # Mathieson et al. 2015 - Nature
    AncientSample(
        "I0585", "R1b1a2a1a2", "Mathieson 2015", 4950, "Hungary", "Bell Beaker"
    ),
    AncientSample(
        "I0807", "R1b1a2a1a2", "Mathieson 2015", 4150, "Germany", "Bell Beaker"
    ),
    # Haak et al. 2015 - Nature
    AncientSample(
        "RISE509", "R1a1a1", "Haak 2015", 4000, "Russia", "Sintashta culture"
    ),
    AncientSample(
        "RISE511", "R1a1a1", "Haak 2015", 4000, "Russia", "Sintashta culture"
    ),
    # Olalde et al. 2018 - Nature
    AncientSample(
        "I2457", "R1b1a1a2a1a2", "Olalde 2018", 4400, "Britain", "Bell Beaker"
    ),
    AncientSample(
        "I2565", "R1b1a1a2a1a2", "Olalde 2018", 4200, "Britain", "Bell Beaker"
    ),
    # Fu et al. 2016 - Nature (early modern humans)
    AncientSample(
        "Ust_Ishim", "K2a", "Fu 2016", 45000, "Siberia", "Oldest sequenced modern human"
    ),
    AncientSample(
        "Oase1", "K", "Fu 2016", 40000, "Romania", "With Neanderthal ancestry"
    ),
    # Lazaridis et al. 2014 - Nature
    AncientSample(
        "Loschbour", "I2a1b", "Lazaridis 2014", 8000, "Luxembourg", "Mesolithic HG"
    ),
    AncientSample(
        "Stuttgart", "G2a2a", "Lazaridis 2014", 7000, "Germany", "Early Neolithic"
    ),
]


def create_ground_truth_file(output_path: Path) -> None:
    """Create ground truth TSV file from known samples."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            ["sample_id", "haplogroup", "publication", "date_bp", "location", "notes"]
        )
        for sample in KNOWN_ANCIENT_SAMPLES:
            writer.writerow(
                [
                    sample.sample_id,
                    sample.haplogroup,
                    sample.publication,
                    sample.date_bp or "",
                    sample.location,
                    sample.notes,
                ]
            )
    click.echo(f"Created ground truth file: {output_path}", err=True)


@click.command()
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("data/ancient"),
    help="Output directory [default: ./data/ancient]",
)
def main(output_dir: Path) -> None:
    """
    Prepare ancient DNA validation data.

    Creates a ground truth file with known ancient DNA Y-chromosome haplogroups.
    Actual VCF files must be obtained separately from AADR or processed from BAM files.

    Data sources:
    - Allen Ancient DNA Resource: https://reich.hms.harvard.edu/downloadable-genotypes-present-day-and-ancient-dna-data-compiled-published-papers
    - Reich Lab VCF downloads require registration
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create ground truth file
    ground_truth_path = output_dir / "ancient_ground_truth.tsv"
    create_ground_truth_file(ground_truth_path)

    click.echo("\nAncient DNA validation setup:", err=True)
    click.echo(f"  Ground truth: {ground_truth_path}", err=True)
    click.echo("", err=True)
    click.echo("To complete validation, you need to:", err=True)
    click.echo("  1. Download AADR VCF data from Reich Lab website", err=True)
    click.echo("  2. Extract Y-chromosome data for samples in ground truth", err=True)
    click.echo("  3. Run validation:", err=True)
    click.echo("     python scripts/validate_ancient.py \\", err=True)
    click.echo("       --vcf data/ancient/ancient_chrY.vcf.gz \\", err=True)
    click.echo("       --ground-truth data/ancient/ancient_ground_truth.tsv", err=True)
    click.echo("", err=True)
    click.echo("Note: Ancient DNA validation is expected to have lower", err=True)
    click.echo("accuracy than modern samples due to DNA damage and missing data.", err=True)


if __name__ == "__main__":
    main()

