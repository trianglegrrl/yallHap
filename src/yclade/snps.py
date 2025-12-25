"""
SNP database for Y-chromosome haplogroup-defining markers.

Loads and manages the YBrowse SNP database with position lookups
across multiple reference genomes (GRCh37, GRCh38, T2T).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

ReferenceGenome = Literal["grch37", "grch38", "t2t"]


@dataclass
class SNP:
    """
    A Y-chromosome SNP marker.

    Attributes:
        name: Primary SNP name
        aliases: Alternative names (L-series, M-series, etc.)
        position_grch37: Position in GRCh37 coordinates
        position_grch38: Position in GRCh38 coordinates
        position_t2t: Position in T2T-CHM13v2.0 coordinates
        ancestral: Ancestral allele
        derived: Derived allele
        haplogroup: Associated haplogroup
    """

    name: str
    aliases: list[str] = field(default_factory=list)
    position_grch37: int | None = None
    position_grch38: int | None = None
    position_t2t: int | None = None
    ancestral: str = ""
    derived: str = ""
    haplogroup: str = ""

    def get_position(self, reference: ReferenceGenome) -> int | None:
        """Get position for specified reference genome."""
        if reference == "grch37":
            return self.position_grch37
        elif reference == "grch38":
            return self.position_grch38
        elif reference == "t2t":
            return self.position_t2t
        else:
            raise ValueError(f"Unknown reference: {reference}")

    @property
    def all_names(self) -> list[str]:
        """Return primary name plus all aliases."""
        return [self.name] + self.aliases


class SNPDatabase:
    """
    Database of Y-chromosome SNP markers.

    Provides efficient lookup by position and by name.
    """

    def __init__(self) -> None:
        self._snps: dict[str, SNP] = {}  # name -> SNP
        self._by_position_grch37: dict[int, list[SNP]] = {}
        self._by_position_grch38: dict[int, list[SNP]] = {}
        self._by_position_t2t: dict[int, list[SNP]] = {}
        self._alias_map: dict[str, str] = {}  # alias -> primary name

    @classmethod
    def from_csv(cls, path: Path | str) -> SNPDatabase:
        """
        Load database from YBrowse-format CSV.

        Expected columns: name, aliases, grch37_pos, grch38_pos, ancestral, derived, haplogroup

        Args:
            path: Path to CSV file

        Returns:
            Populated SNPDatabase instance
        """
        db = cls()
        path = Path(path)

        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                snp = cls._parse_row(row)
                db._add_snp(snp)

        return db

    @classmethod
    def from_ybrowse_vcf(cls, path: Path | str) -> SNPDatabase:
        """
        Load database from YBrowse VCF format.

        Args:
            path: Path to VCF file (can be gzipped)

        Returns:
            Populated SNPDatabase instance
        """
        import gzip

        db = cls()
        path = Path(path)

        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                snp = cls._parse_vcf_line(line)
                if snp:
                    db._add_snp(snp)

        return db

    @staticmethod
    def _parse_row(row: dict[str, str]) -> SNP:
        """Parse a CSV row into an SNP object."""
        aliases = []
        if "aliases" in row and row["aliases"]:
            aliases = [a.strip() for a in row["aliases"].split(",")]

        return SNP(
            name=row.get("name", ""),
            aliases=aliases,
            position_grch37=int(row["grch37_pos"]) if row.get("grch37_pos") else None,
            position_grch38=int(row["grch38_pos"]) if row.get("grch38_pos") else None,
            position_t2t=int(row["t2t_pos"]) if row.get("t2t_pos") else None,
            ancestral=row.get("ancestral", ""),
            derived=row.get("derived", ""),
            haplogroup=row.get("haplogroup", ""),
        )

    @staticmethod
    def _parse_vcf_line(line: str) -> SNP | None:
        """Parse a VCF line into an SNP object."""
        parts = line.strip().split("\t")
        if len(parts) < 5:
            return None

        chrom, pos, snp_id, ref, alt = parts[:5]

        # Only Y chromosome
        if chrom.lower() not in ("y", "chry"):
            return None

        # Parse INFO field for additional data
        info = {}
        if len(parts) > 7:
            for item in parts[7].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info[key] = value

        return SNP(
            name=snp_id,
            position_grch38=int(pos),  # VCF position
            ancestral=ref,
            derived=alt,
            haplogroup=info.get("HG", ""),
        )

    def _add_snp(self, snp: SNP) -> None:
        """Add SNP to database with all indexes."""
        self._snps[snp.name] = snp

        # Index by position
        if snp.position_grch37:
            self._by_position_grch37.setdefault(snp.position_grch37, []).append(snp)
        if snp.position_grch38:
            self._by_position_grch38.setdefault(snp.position_grch38, []).append(snp)
        if snp.position_t2t:
            self._by_position_t2t.setdefault(snp.position_t2t, []).append(snp)

        # Index aliases
        for alias in snp.aliases:
            self._alias_map[alias] = snp.name

    def get_by_name(self, name: str) -> SNP:
        """
        Get SNP by name (including aliases).

        Args:
            name: SNP name or alias

        Returns:
            SNP object

        Raises:
            KeyError: If SNP not found
        """
        # Check primary name first
        if name in self._snps:
            return self._snps[name]

        # Check aliases
        if name in self._alias_map:
            return self._snps[self._alias_map[name]]

        raise KeyError(f"SNP not found: {name}")

    def get_by_position(
        self, position: int, reference: ReferenceGenome = "grch38"
    ) -> list[SNP]:
        """
        Get SNPs at a genomic position.

        Args:
            position: Genomic position
            reference: Reference genome

        Returns:
            List of SNPs at position (may be empty)
        """
        if reference == "grch37":
            return self._by_position_grch37.get(position, [])
        elif reference == "grch38":
            return self._by_position_grch38.get(position, [])
        elif reference == "t2t":
            return self._by_position_t2t.get(position, [])
        else:
            raise ValueError(f"Unknown reference: {reference}")

    def __contains__(self, name: str) -> bool:
        """Check if SNP exists by name or alias."""
        return name in self._snps or name in self._alias_map

    def __len__(self) -> int:
        """Return number of SNPs in database."""
        return len(self._snps)

    def __iter__(self):
        """Iterate over all SNPs."""
        return iter(self._snps.values())

    @property
    def positions(self) -> dict[ReferenceGenome, set[int]]:
        """Return set of all positions for each reference."""
        return {
            "grch37": set(self._by_position_grch37.keys()),
            "grch38": set(self._by_position_grch38.keys()),
            "t2t": set(self._by_position_t2t.keys()),
        }
