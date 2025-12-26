#!/usr/bin/env python3
"""
Convert EIGENSTRAT format to VCF.

EIGENSTRAT format consists of:
- .geno file: genotypes - can be ASCII (0,1,2,9) or packed binary
- .snp file: SNP info (ID, chrom, genetic_pos, physical_pos, ref, alt)
- .ind file: individual info (ID, sex, population)

This script converts to VCF format for use with yallHap.
"""

from __future__ import annotations

import gzip
import sys
from dataclasses import dataclass
from pathlib import Path

import click


@dataclass
class SNPInfo:
    """SNP information from .snp file."""
    name: str
    chrom: str
    genetic_pos: float
    physical_pos: int
    ref: str
    alt: str


def read_ind_file(ind_path: Path) -> list[str]:
    """Read individual IDs from .ind file."""
    individuals = []
    with open(ind_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 1:
                individuals.append(parts[0])
    return individuals


def read_snp_file(snp_path: Path) -> list[SNPInfo]:
    """Read SNP information from .snp file."""
    snps = []
    with open(snp_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                snps.append(SNPInfo(
                    name=parts[0],
                    chrom=parts[1],
                    genetic_pos=float(parts[2]),
                    physical_pos=int(parts[3]),
                    ref=parts[4],
                    alt=parts[5],
                ))
    return snps


def is_packed_ancestrymap(geno_path: Path, num_individuals: int) -> bool:
    """
    Detect if geno file is in packed binary format.

    Packed format has bytes_per_snp = ceil(num_individuals / 4)
    ASCII format has one character per individual per SNP.
    """
    # Read first few bytes
    with open(geno_path, "rb") as f:
        first_bytes = f.read(100)

    # ASCII format should only contain 0, 1, 2, 9, newlines
    ascii_chars = set(b"0129\n\r")
    if all(b in ascii_chars for b in first_bytes):
        return False

    return True


def unpack_genotypes(packed_bytes: bytes, num_individuals: int) -> list[int]:
    """
    Unpack genotypes from packed binary format.

    Each byte contains 4 genotypes, packed as 2 bits each.
    Bits represent COUNT of allele1:
      00 = 0 copies of allele1 (hom allele2)
      01 = 1 copy of allele1 (het)
      10 = 2 copies of allele1 (hom allele1)
      11 = missing
    """
    genotypes = []
    for byte in packed_bytes:
        for shift in [0, 2, 4, 6]:  # Extract 4 genotypes per byte
            if len(genotypes) >= num_individuals:
                break
            geno = (byte >> shift) & 0x03
            if geno == 3:
                geno = 9  # Missing
            genotypes.append(geno)
    return genotypes[:num_individuals]


def write_vcf_header(f, individuals: list[str], reference: str = "GRCh37") -> None:
    """Write VCF header."""
    f.write("##fileformat=VCFv4.2\n")
    f.write(f"##source=eigenstrat_to_vcf.py\n")
    f.write(f"##reference={reference}\n")
    f.write("##INFO=<ID=.,Number=.,Type=String,Description=\"No INFO\">\n")
    f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
    f.write("\t".join(individuals) + "\n")


@click.command()
@click.option(
    "--geno",
    "-g",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to .geno file",
)
@click.option(
    "--snp",
    "-s",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to .snp file",
)
@click.option(
    "--ind",
    "-i",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Path to .ind file",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Output VCF file (.vcf or .vcf.gz)",
)
@click.option(
    "--chrom",
    "-c",
    type=str,
    default=None,
    help="Only output variants on this chromosome (e.g., 'Y' or '24')",
)
@click.option(
    "--samples",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="File with sample IDs to include (one per line)",
)
def main(
    geno: Path,
    snp: Path,
    ind: Path,
    output: Path,
    chrom: str | None,
    samples: Path | None,
) -> None:
    """
    Convert EIGENSTRAT format to VCF.

    Example:
        python eigenstrat_to_vcf.py \\
            --geno data.geno --snp data.snp --ind data.ind \\
            --output chrY.vcf.gz --chrom Y
    """
    click.echo("Reading individual file...", err=True)
    individuals = read_ind_file(ind)
    num_individuals = len(individuals)
    click.echo(f"  Found {num_individuals} individuals", err=True)

    # Filter samples if specified
    sample_set = None
    if samples:
        with open(samples) as f:
            sample_set = set(line.strip() for line in f if line.strip())
        click.echo(f"  Filtering to {len(sample_set)} samples", err=True)

    # Get indices of samples to include
    if sample_set:
        sample_indices = [i for i, ind_id in enumerate(individuals) if ind_id in sample_set]
        filtered_individuals = [individuals[i] for i in sample_indices]
    else:
        sample_indices = list(range(num_individuals))
        filtered_individuals = individuals

    click.echo("Reading SNP file...", err=True)
    snps = read_snp_file(snp)
    num_snps = len(snps)
    click.echo(f"  Found {num_snps} SNPs", err=True)

    # Determine target chromosomes
    chrom_variants = None
    if chrom:
        if chrom.upper() in ["Y", "24", "CHRY"]:
            chrom_variants = {"Y", "24", "chrY"}
        else:
            chrom_variants = {chrom}
        target_snps = [s for s in snps if s.chrom in chrom_variants]
        click.echo(f"  {len(target_snps)} SNPs on chromosome {chrom}", err=True)

    # Detect file format
    is_packed = is_packed_ancestrymap(geno, num_individuals)
    click.echo(f"  Geno file format: {'packed binary' if is_packed else 'ASCII'}", err=True)

    # Calculate bytes per SNP for packed format
    bytes_per_snp = (num_individuals + 3) // 4 if is_packed else num_individuals + 1  # +1 for newline

    # Open output file
    if str(output).endswith(".gz"):
        out_file = gzip.open(output, "wt")
    else:
        out_file = open(output, "w")

    try:
        # Write header
        write_vcf_header(out_file, filtered_individuals)

        click.echo("Converting genotypes...", err=True)

        variants_written = 0

        with open(geno, "rb") as geno_file:
            for snp_idx, snp_info in enumerate(snps):
                # Read genotype data for this SNP
                if is_packed:
                    packed_data = geno_file.read(bytes_per_snp)
                    if len(packed_data) < bytes_per_snp:
                        click.echo(f"Warning: unexpected end of file at SNP {snp_idx}", err=True)
                        break
                    genotypes = unpack_genotypes(packed_data, num_individuals)
                else:
                    line = geno_file.readline().decode("ascii", errors="ignore").strip()
                    genotypes = [int(c) if c in "0129" else 9 for c in line]

                # Filter by chromosome
                if chrom_variants and snp_info.chrom not in chrom_variants:
                    continue

                # Extract genotypes for selected samples
                # EIGENSTRAT genotype encodes COUNT of allele1 (first allele in .snp file):
                # 0 = 0 copies of allele1 = homozygous allele2 = VCF 1/1 (hom ALT)
                # 1 = 1 copy of allele1 = heterozygous = VCF 0/1
                # 2 = 2 copies of allele1 = homozygous allele1 = VCF 0/0 (hom REF)
                # 9 = missing = VCF ./.
                sample_genotypes = []
                for idx in sample_indices:
                    if idx < len(genotypes):
                        gt = genotypes[idx]
                        if gt == 0:
                            sample_genotypes.append("1/1")  # hom ALT (allele2)
                        elif gt == 1:
                            sample_genotypes.append("0/1")  # het
                        elif gt == 2:
                            sample_genotypes.append("0/0")  # hom REF (allele1)
                        else:
                            sample_genotypes.append("./.")  # missing
                    else:
                        sample_genotypes.append("./.")

                # Skip if all missing
                if all(gt == "./." for gt in sample_genotypes):
                    continue

                # Normalize chromosome name
                vcf_chrom = snp_info.chrom
                if vcf_chrom == "24":
                    vcf_chrom = "Y"

                # Write VCF line
                out_file.write(
                    f"{vcf_chrom}\t{snp_info.physical_pos}\t{snp_info.name}\t"
                    f"{snp_info.ref}\t{snp_info.alt}\t.\tPASS\t.\tGT\t"
                )
                out_file.write("\t".join(sample_genotypes) + "\n")
                variants_written += 1

                if variants_written % 5000 == 0:
                    click.echo(f"  Processed {variants_written} variants...", err=True)

        click.echo(f"Wrote {variants_written} variants to {output}", err=True)

    finally:
        out_file.close()

    click.echo("Done!", err=True)


if __name__ == "__main__":
    main()
