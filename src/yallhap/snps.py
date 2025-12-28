"""
SNP database for Y-chromosome haplogroup-defining markers.

Loads and manages the YBrowse SNP database with position lookups
across multiple reference genomes (GRCh37, GRCh38, T2T).
"""

from __future__ import annotations

import csv
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    pass

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

    @classmethod
    def from_ybrowse_gff_csv(cls, path: Path | str) -> SNPDatabase:
        """
        Load database from YBrowse GFF-style CSV format.

        This is the actual format of the YBrowse snps_hg38.csv file, which uses
        GFF-style columns rather than the simpler internal format.

        Expected columns:
            seqid, source, type, start, end, score, strand, phase,
            Name, ID, allele_anc, allele_der, YCC_haplogroup, ISOGG_haplogroup,
            mutation, count_tested, count_derived, ref, comment

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
                snp = cls._parse_ybrowse_gff_row(row)
                if snp:
                    db._add_snp(snp)

        return db

    @staticmethod
    def _parse_ybrowse_gff_row(
        row: dict[str, str], detected_reference: ReferenceGenome | None = None
    ) -> SNP | None:
        """
        Parse a YBrowse GFF-style CSV row into an SNP object.

        Skips indels (rows where alleles are not single nucleotides).

        Args:
            row: CSV row as dict
            detected_reference: Reference genome detected from seqid column
        """
        # Get alleles (handle None values from malformed CSV rows)
        ancestral = (row.get("allele_anc") or "").strip()
        derived = (row.get("allele_der") or "").strip()

        # Skip indels - alleles should be single nucleotides
        valid_alleles = {"A", "C", "G", "T"}
        if ancestral not in valid_alleles or derived not in valid_alleles:
            return None

        # Get position (start column is 1-based position)
        position_str = (row.get("start") or "").strip()
        if not position_str:
            return None

        try:
            position = int(position_str)
        except ValueError:
            return None

        # Get SNP name
        name = (row.get("Name") or "").strip() or (row.get("ID") or "").strip()
        if not name:
            return None

        # Get haplogroup (prefer YCC_haplogroup)
        haplogroup = (row.get("YCC_haplogroup") or "").strip()

        # Detect reference from seqid if not provided
        if detected_reference is None:
            seqid = (row.get("seqid") or "").lower()
            if "hg19" in seqid or "grch37" in seqid:
                detected_reference = "grch37"
            elif "hg38" in seqid or "grch38" in seqid:
                detected_reference = "grch38"
            else:
                # Default to grch38
                detected_reference = "grch38"

        # Create SNP with position in correct reference
        if detected_reference == "grch37":
            return SNP(
                name=name,
                position_grch37=position,
                ancestral=ancestral,
                derived=derived,
                haplogroup=haplogroup,
            )
        else:
            return SNP(
                name=name,
                position_grch38=position,
                ancestral=ancestral,
                derived=derived,
                haplogroup=haplogroup,
            )

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

    def __iter__(self) -> Iterator[SNP]:
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

    def lift_to_t2t(
        self,
        chain_file: Path | str,
        source_reference: ReferenceGenome = "grch38",
    ) -> int:
        """
        Compute T2T positions for all SNPs using liftover.

        Uses pyliftover to convert positions from source reference to T2T-CHM13v2.0.
        Updates SNP objects in place with the computed T2T positions.

        Args:
            chain_file: Path to chain file (e.g., grch38-chm13v2.chain)
            source_reference: Source reference genome (grch37 or grch38)

        Returns:
            Number of SNPs successfully lifted over

        Raises:
            FileNotFoundError: If chain file doesn't exist
            ValueError: If source_reference is t2t (can't liftover to itself)
        """
        from pyliftover import LiftOver

        if source_reference == "t2t":
            raise ValueError("Cannot liftover from T2T to T2T")

        chain_path = Path(chain_file)
        if not chain_path.exists():
            raise FileNotFoundError(f"Chain file not found: {chain_path}")

        lo = LiftOver(str(chain_path))
        lifted_count = 0

        for snp in self._snps.values():
            # Get source position
            if source_reference == "grch37":
                source_pos = snp.position_grch37
                # pyliftover uses UCSC chromosome names
                chrom = "chrY"
            else:
                source_pos = snp.position_grch38
                chrom = "chrY"

            if source_pos is None:
                continue

            # pyliftover uses 0-based coordinates
            result = lo.convert_coordinate(chrom, source_pos - 1)

            if result and len(result) > 0:
                # Result is list of (chrom, pos, strand, score) tuples
                # Take the first (best) result and convert back to 1-based
                new_chrom, new_pos, strand, score = result[0]
                snp.position_t2t = new_pos + 1
                lifted_count += 1

                # Update T2T position index
                self._by_position_t2t.setdefault(snp.position_t2t, []).append(snp)

        return lifted_count

    def compute_all_t2t_positions(
        self,
        grch38_chain: Path | str | None = None,
        grch37_chain: Path | str | None = None,
    ) -> dict[str, int]:
        """
        Compute T2T positions from available source references.

        Tries GRCh38 first if available, then falls back to GRCh37.
        Returns statistics about the liftover process.

        Args:
            grch38_chain: Path to GRCh38-to-T2T chain file
            grch37_chain: Path to GRCh37-to-T2T chain file

        Returns:
            Dict with counts: lifted_from_grch38, lifted_from_grch37, total_with_t2t
        """
        stats = {
            "lifted_from_grch38": 0,
            "lifted_from_grch37": 0,
            "total_with_t2t": 0,
        }

        # Try GRCh38 first
        if grch38_chain and Path(grch38_chain).exists():
            stats["lifted_from_grch38"] = self.lift_to_t2t(grch38_chain, "grch38")

        # Then try GRCh37 for any remaining SNPs without T2T positions
        if grch37_chain and Path(grch37_chain).exists():
            # Only lift SNPs that don't have T2T positions yet
            snps_needing_lift = [s for s in self._snps.values() if s.position_t2t is None]
            if snps_needing_lift:
                stats["lifted_from_grch37"] = self.lift_to_t2t(grch37_chain, "grch37")

        stats["total_with_t2t"] = sum(1 for s in self._snps.values() if s.position_t2t)

        return stats
