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
        # Store all SNPs with same name (some SNPs map to multiple haplogroups)
        self._snps_by_name: dict[str, list[ISOGGSNP]] = {}
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
        if snp.name not in self._snps_by_name:
            self._snps_by_name[snp.name] = []
        self._snps_by_name[snp.name].append(snp)
        if snp.position not in self._snps_by_position:
            self._snps_by_position[snp.position] = []
        self._snps_by_position[snp.position].append(snp)
        self._haplogroups.add(snp.haplogroup)

    def get_by_name(self, name: str) -> ISOGGSNP | None:
        """Get first SNP by name."""
        snps = self._snps_by_name.get(name, [])
        return snps[0] if snps else None

    def get_all_by_name(self, name: str) -> list[ISOGGSNP]:
        """Get all SNPs with given name (for recurrent mutations)."""
        return self._snps_by_name.get(name, [])

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
        for snp_list in self._snps_by_name.values():
            for snp in snp_list:
                counts[snp.haplogroup] = counts.get(snp.haplogroup, 0) + 1
        return counts

    def __len__(self) -> int:
        """Return total number of SNP entries in database."""
        return sum(len(snps) for snps in self._snps_by_name.values())

    def __iter__(self) -> Iterator[ISOGGSNP]:
        """Iterate over all SNPs."""
        for snp_list in self._snps_by_name.values():
            yield from snp_list

    def __contains__(self, name: str) -> bool:
        """Check if SNP name is in database."""
        return name in self._snps_by_name


class ISOGGMapper:
    """
    Maps YFull haplogroup names to ISOGG nomenclature.

    The join between YFull and ISOGG is on SNP name:
    - YFull uses format "X-SNP" (e.g., "R-L21")
    - ISOGG uses hierarchical names but each is defined by a SNP

    Algorithm:
    1. Parse the defining SNP from YFull name
    2. Look up that SNP in ISOGG database
    3. If not found, walk UP the YFull tree trying each ancestor's SNP
    4. Return the first ISOGG match (rollup to nearest ancestor with ISOGG entry)
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
        # Note: same SNP may appear multiple times with different haplogroups
        # (recurrent mutations or database errors)
        # We keep ALL mappings and choose based on context when looking up
        self._snp_to_isogg_all: dict[str, list[str]] = {}
        for snp in isogg_db:
            if snp.name not in self._snp_to_isogg_all:
                self._snp_to_isogg_all[snp.name] = []
            self._snp_to_isogg_all[snp.name].append(snp.haplogroup)

        # Set of valid ISOGG haplogroup names (for direct name matching)
        self._isogg_hgs = set(isogg_db.get_all_haplogroups())

        # Pre-build YFull to ISOGG mapping
        self._yfull_to_isogg: dict[str, str] = {}
        self._build_mapping()

    def _parse_snps_from_node(self, node_name: str, node_snps: list[str]) -> list[str]:
        """
        Extract all individual SNP names from a node.

        YFull tree has SNPs in multiple formats:
        - In the name: "R-L21" -> extract "L21"
        - In snps list: may be comma-separated like "FGC79428, Z2665"

        Args:
            node_name: The node's haplogroup name
            node_snps: The node's snps list

        Returns:
            List of individual SNP names
        """
        snps: list[str] = []

        # Extract SNP from name (e.g., "R-L21" -> "L21")
        if "-" in node_name:
            snp_from_name = node_name.split("-", 1)[1]
            snps.append(snp_from_name)

        # Parse snps list, splitting on commas
        for snp_entry in node_snps:
            # Split on comma and clean up
            for part in snp_entry.split(","):
                part = part.strip()
                if part:
                    # Also split on "/" (YFull uses this too)
                    for subpart in part.split("/"):
                        subpart = subpart.strip()
                        if subpart:
                            snps.append(subpart)

        return snps

    def _get_isogg_for_snp(self, snp: str, major_clade_hint: str | None = None) -> str | None:
        """
        Get ISOGG haplogroup for a single SNP.

        When a SNP maps to multiple haplogroups (recurrent mutation),
        prefer the one matching the major clade hint, or the shorter one.

        Args:
            snp: SNP name to look up
            major_clade_hint: Preferred major clade (e.g., "I" for I-L596)

        Returns:
            ISOGG haplogroup, or None if not found
        """
        if snp not in self._snp_to_isogg_all:
            return None

        candidates = self._snp_to_isogg_all[snp]
        if len(candidates) == 1:
            return candidates[0]

        # Multiple candidates - need to choose
        if major_clade_hint:
            # Prefer candidate matching the major clade
            for hg in candidates:
                if hg and hg[0].upper() == major_clade_hint[0].upper():
                    return hg

        # Fall back to shortest (most ancestral)
        return min(candidates, key=len)

    def _find_isogg_for_node_snps(
        self, snps: list[str], major_clade_hint: str | None = None
    ) -> str | None:
        """
        Find ISOGG haplogroup from a list of SNPs.

        Returns the most specific (longest) ISOGG haplogroup that matches
        any of the provided SNPs, with preference for matching major clade.

        Args:
            snps: List of SNP names to check
            major_clade_hint: Preferred major clade for disambiguation

        Returns:
            Most specific ISOGG haplogroup, or None if no match
        """
        best_isogg: str | None = None
        best_length = 0

        for snp in snps:
            isogg_hg = self._get_isogg_for_snp(snp, major_clade_hint)
            if isogg_hg and len(isogg_hg) > best_length:
                best_isogg = isogg_hg
                best_length = len(isogg_hg)

        return best_isogg

    def _build_mapping(self) -> None:
        """
        Build YFull to ISOGG haplogroup mapping.

        For each YFull node:
        1. Try to find ISOGG match from the node's own SNPs
        2. If not found, walk up to ancestors until a match is found
        """
        for node in self.tree.iter_depth_first():
            isogg_match = self._find_isogg_for_yfull(node.name)
            if isogg_match:
                self._yfull_to_isogg[node.name] = isogg_match

    def _find_isogg_for_yfull(self, yfull_name: str) -> str | None:
        """
        Find ISOGG haplogroup for a YFull haplogroup name.

        Strategy:
        1. Check if the node's name is directly in ISOGG (e.g., "I2")
        2. Parse SNPs from the node and look them up
        3. If no match, walk up to parent and repeat
        4. Return the first match found (rollup)

        Args:
            yfull_name: YFull haplogroup name

        Returns:
            ISOGG haplogroup, or None if no mapping exists
        """
        # Extract major clade hint from the original name
        # E.g., "I-L596" -> "I", "R1b" -> "R"
        major_clade = yfull_name.split("-")[0] if "-" in yfull_name else yfull_name
        if major_clade:
            major_clade = major_clade[0].upper()

        current_name: str | None = yfull_name

        while current_name is not None:
            try:
                current = self.tree.get(current_name)
            except KeyError:
                break

            # Check if the base name is directly an ISOGG haplogroup
            base = current.name.split("-")[0] if "-" in current.name else current.name
            if base in self._isogg_hgs:
                return base

            # Parse all SNPs from this node
            snps = self._parse_snps_from_node(current.name, current.snps)

            # Look up SNPs in ISOGG with major clade hint for disambiguation
            isogg_match = self._find_isogg_for_node_snps(snps, major_clade)
            if isogg_match:
                return isogg_match

            # Walk up to parent
            current_name = current.parent_name

        return None

    def to_isogg(self, yfull_haplogroup: str) -> str:
        """
        Convert YFull haplogroup to ISOGG nomenclature.

        Args:
            yfull_haplogroup: YFull haplogroup name (e.g., "R-L21", "I2")

        Returns:
            ISOGG haplogroup name (e.g., "R1b1a1b1a1a2c", "I2")
            Returns original if no mapping exists
        """
        # Extract major clade hint for disambiguation
        major_clade = yfull_haplogroup.split("-")[0][0].upper() if yfull_haplogroup else None

        # Priority 1: Direct SNP lookup from name
        # YFull format is "LETTER-SNP" (e.g., "R-L21", "I-M438")
        if "-" in yfull_haplogroup:
            snp_name = yfull_haplogroup.split("-", 1)[1]
            isogg = self._get_isogg_for_snp(snp_name, major_clade)
            if isogg:
                return isogg

        # Priority 2: Check if base name is directly an ISOGG haplogroup
        base = yfull_haplogroup.split("-")[0] if "-" in yfull_haplogroup else yfull_haplogroup
        if base in self._isogg_hgs:
            return base

        # Priority 3: Pre-built mapping (includes rollup to ancestors)
        if yfull_haplogroup in self._yfull_to_isogg:
            return self._yfull_to_isogg[yfull_haplogroup]

        # No mapping found, return original with tilde to indicate uncertainty
        return f"{yfull_haplogroup}~"

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
        return self._find_isogg_for_node_snps(derived_snp_names)

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
