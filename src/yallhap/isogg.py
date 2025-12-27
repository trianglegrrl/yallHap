"""
ISOGG SNP database and haplogroup mapping.

Provides ISOGG (International Society of Genetic Genealogy) haplogroup
nomenclature support for YFull haplogroup calls.

Based on pathPhynder's ISOGG SNP lists.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from yallhap.tree import Tree


@dataclass
class ISOGGSNP:
    """
    An ISOGG SNP marker.

    Attributes:
        name: SNP name (e.g., "M269", "L21")
        haplogroup: ISOGG haplogroup (e.g., "R1b1a1b")
        position: GRCh38 position
        ancestral: Ancestral allele
        derived: Derived allele
    """

    name: str
    haplogroup: str
    position: int
    ancestral: str
    derived: str


class ISOGGDatabase:
    """
    Database of ISOGG SNP markers.

    Supports loading from pathPhynder format files and lookup
    by SNP name or position.
    """

    def __init__(self) -> None:
        """Initialize empty database."""
        self._snps_by_name: dict[str, ISOGGSNP] = {}
        self._snps_by_position: dict[int, list[ISOGGSNP]] = {}
        self._haplogroups: set[str] = set()

    @classmethod
    def from_file(cls, path: Path | str) -> ISOGGDatabase:
        """
        Load ISOGG SNPs from pathPhynder format file.

        File format (tab-separated):
        SNP_name    ISOGG_haplogroup    position    ancestral    derived

        Args:
            path: Path to ISOGG SNP file

        Returns:
            Populated ISOGGDatabase
        """
        db = cls()
        path = Path(path)

        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) < 5:
                    continue

                name, haplogroup, pos_str, ancestral, derived = parts[:5]

                try:
                    position = int(pos_str)
                except ValueError:
                    continue

                snp = ISOGGSNP(
                    name=name,
                    haplogroup=haplogroup,
                    position=position,
                    ancestral=ancestral,
                    derived=derived,
                )
                db._add_snp(snp)

        return db

    def _add_snp(self, snp: ISOGGSNP) -> None:
        """Add a SNP to the database."""
        self._snps_by_name[snp.name] = snp
        if snp.position not in self._snps_by_position:
            self._snps_by_position[snp.position] = []
        self._snps_by_position[snp.position].append(snp)
        self._haplogroups.add(snp.haplogroup)

    def get_by_name(self, name: str) -> ISOGGSNP | None:
        """Get SNP by name."""
        return self._snps_by_name.get(name)

    def get_by_position(self, position: int) -> ISOGGSNP | None:
        """Get first SNP at position."""
        snps = self._snps_by_position.get(position, [])
        return snps[0] if snps else None

    def get_all_at_position(self, position: int) -> list[ISOGGSNP]:
        """Get all SNPs at position."""
        return self._snps_by_position.get(position, [])

    def get_haplogroups_at_position(self, position: int) -> list[str]:
        """Get all haplogroups with SNPs at position."""
        snps = self._snps_by_position.get(position, [])
        return [snp.haplogroup for snp in snps]

    def get_all_haplogroups(self) -> list[str]:
        """Get list of all unique haplogroups."""
        return sorted(self._haplogroups)

    def get_snp_counts_by_haplogroup(self) -> dict[str, int]:
        """Get count of SNPs per haplogroup."""
        counts: dict[str, int] = {}
        for snp in self._snps_by_name.values():
            counts[snp.haplogroup] = counts.get(snp.haplogroup, 0) + 1
        return counts

    def __len__(self) -> int:
        """Return number of SNPs in database."""
        return len(self._snps_by_name)

    def __iter__(self) -> Iterator[ISOGGSNP]:
        """Iterate over all SNPs."""
        return iter(self._snps_by_name.values())

    def __contains__(self, name: str) -> bool:
        """Check if SNP name is in database."""
        return name in self._snps_by_name


class ISOGGMapper:
    """
    Maps YFull haplogroup names to ISOGG nomenclature.

    Uses SNP database to find ISOGG equivalents for YFull haplogroups.
    """

    def __init__(self, tree: Tree, isogg_db: ISOGGDatabase) -> None:
        """
        Initialize mapper.

        Args:
            tree: YFull phylogenetic tree
            isogg_db: ISOGG SNP database
        """
        self.tree = tree
        self.isogg_db = isogg_db

        # Build mapping from YFull to ISOGG haplogroups
        # This is a simplified approach - could be enhanced with
        # more sophisticated name matching
        self._yfull_to_isogg: dict[str, str] = {}
        self._build_mapping()

    def _build_mapping(self) -> None:
        """Build YFull to ISOGG haplogroup mapping."""
        # Get all ISOGG haplogroups
        isogg_hgs = set(self.isogg_db.get_all_haplogroups())

        # For each YFull haplogroup, try to find an ISOGG match
        for node in self.tree.iter_depth_first():
            yfull_name = node.name
            isogg_match = self._find_isogg_match(yfull_name, isogg_hgs)
            if isogg_match:
                self._yfull_to_isogg[yfull_name] = isogg_match

    def _find_isogg_match(self, yfull_name: str, isogg_hgs: set[str]) -> str | None:
        """Find best ISOGG match for a YFull haplogroup name."""
        # Direct match
        if yfull_name in isogg_hgs:
            return yfull_name

        # YFull uses format like "R-L21" while ISOGG uses "R1b1a1b1a1a2c"
        # Try stripping the SNP suffix
        if "-" in yfull_name:
            base = yfull_name.split("-")[0]
            if base in isogg_hgs:
                return base

        # Look for ISOGG haplogroups that start with YFull base
        # E.g., "R1b" might match "R1b1a1b..."
        for isogg_hg in isogg_hgs:
            if isogg_hg.startswith(yfull_name):
                return isogg_hg

        return None

    def to_isogg(self, yfull_haplogroup: str) -> str:
        """
        Convert YFull haplogroup to ISOGG nomenclature.

        Args:
            yfull_haplogroup: YFull haplogroup name

        Returns:
            ISOGG haplogroup name, or original if no mapping exists
        """
        return self._yfull_to_isogg.get(yfull_haplogroup, yfull_haplogroup)

    def from_isogg(self, isogg_haplogroup: str) -> str | None:
        """
        Find YFull haplogroup for an ISOGG haplogroup.

        Args:
            isogg_haplogroup: ISOGG haplogroup name

        Returns:
            YFull haplogroup name, or None if no mapping exists
        """
        for yfull, isogg in self._yfull_to_isogg.items():
            if isogg == isogg_haplogroup:
                return yfull
        return None
