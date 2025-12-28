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

    Uses SNP database to find the most specific ISOGG equivalent for
    YFull haplogroups based on the SNPs that define them.
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

        # Build mapping from SNP name to ISOGG haplogroup
        self._snp_to_isogg: dict[str, str] = {}
        for snp in isogg_db:
            self._snp_to_isogg[snp.name] = snp.haplogroup

        # Build sorted list of ISOGG haplogroups (longer = more specific)
        self._isogg_hgs_by_length = sorted(
            isogg_db.get_all_haplogroups(),
            key=len,
            reverse=True,
        )

        # Pre-build YFull to ISOGG mapping using tree structure
        self._yfull_to_isogg: dict[str, str] = {}
        self._build_mapping()

    def _build_mapping(self) -> None:
        """Build YFull to ISOGG haplogroup mapping using SNP names and tree."""
        isogg_hgs = set(self.isogg_db.get_all_haplogroups())

        for node in self.tree.iter_depth_first():
            yfull_name = node.name
            isogg_match = self._find_isogg_for_node(node, isogg_hgs)
            if isogg_match:
                self._yfull_to_isogg[yfull_name] = isogg_match

    def _find_isogg_for_node(self, node: Node, isogg_hgs: set[str]) -> str | None:  # noqa: F821
        """
        Find the most specific ISOGG haplogroup for a tree node.

        Strategy:
        1. Check if the node's defining SNPs are in the ISOGG database
        2. Walk up the tree and collect all defining SNPs
        3. Find the most specific (longest) ISOGG haplogroup from those SNPs
        """
        # Check if this node or its ancestors have defining SNPs in ISOGG
        best_isogg: str | None = None
        best_length = 0

        # Walk from this node to root, collecting ISOGG haplogroups
        current_name: str | None = node.name
        while current_name is not None:
            try:
                current = self.tree.get(current_name)
            except KeyError:
                break

            # Check the node's defining SNPs
            for snp_name in current.snps:
                if snp_name in self._snp_to_isogg:
                    isogg_hg = self._snp_to_isogg[snp_name]
                    if len(isogg_hg) > best_length:
                        best_isogg = isogg_hg
                        best_length = len(isogg_hg)

            # Extract SNP name from YFull format (e.g., "R-L21" -> "L21")
            if "-" in current.name:
                snp_name = current.name.split("-", 1)[1]
                if snp_name in self._snp_to_isogg:
                    isogg_hg = self._snp_to_isogg[snp_name]
                    if len(isogg_hg) > best_length:
                        best_isogg = isogg_hg
                        best_length = len(isogg_hg)

            # Also check if the full name or base matches
            if current.name in isogg_hgs:
                if len(current.name) > best_length:
                    best_isogg = current.name
                    best_length = len(current.name)

            # Check base haplogroup (before hyphen)
            base = current.name.split("-")[0] if "-" in current.name else current.name
            if base in isogg_hgs:
                if len(base) > best_length:
                    best_isogg = base
                    best_length = len(base)

            # Move to parent
            current_name = current.parent_name

        return best_isogg

    def to_isogg(self, yfull_haplogroup: str) -> str:
        """
        Convert YFull haplogroup to ISOGG nomenclature.

        Args:
            yfull_haplogroup: YFull haplogroup name

        Returns:
            ISOGG haplogroup name, or original if no mapping exists
        """
        # Priority 1: Direct SNP lookup (most reliable)
        # YFull format is "LETTER-SNP" (e.g., "R-L21", "I-M438")
        if "-" in yfull_haplogroup:
            snp_name = yfull_haplogroup.split("-", 1)[1]
            if snp_name in self._snp_to_isogg:
                return self._snp_to_isogg[snp_name]

        # Priority 2: Pre-built mapping from tree traversal
        if yfull_haplogroup in self._yfull_to_isogg:
            return self._yfull_to_isogg[yfull_haplogroup]

        # Priority 3: Base haplogroup letter if it's in ISOGG
        base = yfull_haplogroup.split("-")[0] if "-" in yfull_haplogroup else yfull_haplogroup
        isogg_hgs = set(self.isogg_db.get_all_haplogroups())
        if base in isogg_hgs:
            return base

        return yfull_haplogroup

    def to_isogg_from_snps(self, derived_snp_names: list[str]) -> str | None:
        """
        Find the most specific ISOGG haplogroup from a list of derived SNPs.

        This is more accurate than name-based mapping as it uses actual
        observed derived alleles.

        Args:
            derived_snp_names: List of SNP names that are derived

        Returns:
            Most specific ISOGG haplogroup, or None if no matches
        """
        best_isogg: str | None = None
        best_length = 0

        for snp_name in derived_snp_names:
            if snp_name in self._snp_to_isogg:
                isogg_hg = self._snp_to_isogg[snp_name]
                if len(isogg_hg) > best_length:
                    best_isogg = isogg_hg
                    best_length = len(isogg_hg)

        return best_isogg

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
