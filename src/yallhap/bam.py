"""
BAM file reading for Y-chromosome haplogroup classification.

Provides direct BAM reading using pysam pileup, inspired by pathPhynder's
approach. This enables Bayesian classification on ancient DNA samples
without requiring pre-called VCF files.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from types import TracebackType

import pysam

from yallhap.vcf import Variant


@dataclass
class PileupResult:
    """
    Result of querying a single position in a BAM file.

    Attributes:
        position: 1-based genomic position
        ref_allele: Reference allele at this position
        allele_counts: Dict of base -> count (A, T, C, G)
        total_depth: Total read depth at this position
    """

    position: int
    ref_allele: str
    allele_counts: dict[str, int]
    total_depth: int

    def get_allele_depth(self, ref: str, alt: str) -> tuple[int, int]:
        """
        Get allele depths for ref and alt alleles.

        Args:
            ref: Reference allele
            alt: Alternative allele

        Returns:
            Tuple of (ref_depth, alt_depth)
        """
        ref_depth = self.allele_counts.get(ref.upper(), 0)
        alt_depth = self.allele_counts.get(alt.upper(), 0)
        return (ref_depth, alt_depth)


class BAMReader:
    """
    Reader for Y-chromosome data directly from BAM files.

    Uses pysam pileup to query allele depths at specific positions,
    enabling haplogroup classification without pre-called variants.
    This is particularly useful for ancient DNA where variant calling
    may miss low-coverage positions.
    """

    # Possible Y chromosome names
    Y_CHROMS = {"Y", "chrY", "y", "chry", "24"}

    def __init__(
        self,
        path: Path | str,
        reference_path: Path | str | None = None,
        min_base_quality: int = 20,
        min_mapping_quality: int = 20,
    ):
        """
        Initialize BAM reader.

        Args:
            path: Path to BAM file (must be indexed with .bai)
            reference_path: Optional path to reference FASTA (for MD tag reconstruction)
            min_base_quality: Minimum base quality to count (default 20)
            min_mapping_quality: Minimum mapping quality to count (default 20)
        """
        self.path = Path(path)
        self.reference_path = Path(reference_path) if reference_path else None
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality
        self._bam: pysam.AlignmentFile | None = None
        self._y_chrom: str | None = None
        self._sample: str | None = None

    def __enter__(self) -> BAMReader:
        self._bam = pysam.AlignmentFile(str(self.path), "rb")
        self._detect_y_chrom()
        self._detect_sample()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        if self._bam:
            self._bam.close()

    def _detect_y_chrom(self) -> None:
        """Detect Y chromosome naming convention in BAM."""
        if self._bam is None:
            raise RuntimeError("BAM not opened")

        references = set(self._bam.references)
        for name in self.Y_CHROMS:
            if name in references:
                self._y_chrom = name
                return

        raise ValueError(f"No Y chromosome found in BAM. References: {references}")

    def _detect_sample(self) -> None:
        """Detect sample name from BAM read groups."""
        if self._bam is None:
            raise RuntimeError("BAM not opened")

        # Try to get sample name from read groups
        header_dict = self._bam.header.to_dict()
        if "RG" in header_dict:
            read_groups = header_dict["RG"]
            if read_groups and len(read_groups) > 0:
                # Use first read group's sample name
                self._sample = read_groups[0].get("SM", self.path.stem)
                return

        # Fall back to file name
        self._sample = self.path.stem

    @property
    def sample(self) -> str:
        """Return the sample name."""
        if self._sample is None:
            raise RuntimeError("Sample not detected (call __enter__ first)")
        return self._sample

    @property
    def y_chrom(self) -> str:
        """Return the Y chromosome name used in this BAM."""
        if self._y_chrom is None:
            raise RuntimeError("Y chromosome not detected (call __enter__ first)")
        return self._y_chrom

    def pileup_at_position(self, position: int) -> PileupResult | None:
        """
        Get pileup data at a specific position.

        Args:
            position: 1-based genomic position

        Returns:
            PileupResult with allele counts, or None if no coverage
        """
        if self._bam is None:
            raise RuntimeError("BAM not opened")
        if self._y_chrom is None:
            raise RuntimeError("Y chromosome not detected")

        # pysam uses 0-based coordinates
        start = position - 1
        end = position

        allele_counts: Counter[str] = Counter()
        ref_allele: str | None = None

        # Use pileup to get base calls at this position
        for pileup_column in self._bam.pileup(
            self._y_chrom,
            start,
            end,
            truncate=True,
            min_base_quality=self.min_base_quality,
            min_mapping_quality=self.min_mapping_quality,
            stepper="samtools",
        ):
            if pileup_column.reference_pos != start:
                continue

            for pileup_read in pileup_column.pileups:
                # Skip deletions and ref skips
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue

                query_pos = pileup_read.query_position
                if query_pos is None:
                    continue

                # Get the base at this position
                query_seq = pileup_read.alignment.query_sequence
                if query_seq is None:
                    continue
                base = query_seq[query_pos]
                if base:
                    allele_counts[base.upper()] += 1

        total_depth = sum(allele_counts.values())
        if total_depth == 0:
            return None

        # If we don't have ref allele, use most common as proxy (not ideal)
        if ref_allele is None:
            ref_allele = allele_counts.most_common(1)[0][0] if allele_counts else "N"

        return PileupResult(
            position=position,
            ref_allele=ref_allele,
            allele_counts=dict(allele_counts),
            total_depth=total_depth,
        )

    def pileup_at_positions(self, positions: set[int]) -> Iterator[PileupResult]:
        """
        Get pileup data at multiple positions efficiently.

        Args:
            positions: Set of 1-based genomic positions

        Yields:
            PileupResult for each position with coverage
        """
        if self._bam is None:
            raise RuntimeError("BAM not opened")
        if self._y_chrom is None:
            raise RuntimeError("Y chromosome not detected")

        if not positions:
            return

        # Sort positions for efficient iteration
        sorted_positions = sorted(positions)
        min_pos = sorted_positions[0] - 1  # 0-based
        max_pos = sorted_positions[-1]

        positions_set = set(positions)

        # Use pileup over the range
        for pileup_column in self._bam.pileup(
            self._y_chrom,
            min_pos,
            max_pos,
            truncate=True,
            min_base_quality=self.min_base_quality,
            min_mapping_quality=self.min_mapping_quality,
            stepper="samtools",
        ):
            pos_1based = pileup_column.reference_pos + 1

            if pos_1based not in positions_set:
                continue

            allele_counts: Counter[str] = Counter()

            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue

                query_pos = pileup_read.query_position
                if query_pos is None:
                    continue

                query_seq = pileup_read.alignment.query_sequence
                if query_seq is None:
                    continue
                base = query_seq[query_pos]
                if base:
                    allele_counts[base.upper()] += 1

            total_depth = sum(allele_counts.values())
            if total_depth == 0:
                continue

            ref_allele = allele_counts.most_common(1)[0][0] if allele_counts else "N"

            yield PileupResult(
                position=pos_1based,
                ref_allele=ref_allele,
                allele_counts=dict(allele_counts),
                total_depth=total_depth,
            )

    def get_variant_at_position(
        self,
        position: int,
        ref: str,
        alt: str,
    ) -> Variant | None:
        """
        Get a Variant object at a specific position.

        This creates a Variant compatible with the existing classifier,
        using BAM pileup data instead of VCF calls.

        Args:
            position: 1-based genomic position
            ref: Expected reference allele
            alt: Expected alternative allele

        Returns:
            Variant object with allele depths from BAM, or None if no coverage
        """
        pileup = self.pileup_at_position(position)
        if pileup is None:
            return None

        ref_depth, alt_depth = pileup.get_allele_depth(ref, alt)
        total_depth = pileup.total_depth

        # Determine genotype based on allele depths
        # Use majority vote with threshold
        if total_depth == 0:
            genotype = None
        elif alt_depth == 0:
            genotype = 0  # Reference
        elif ref_depth == 0:
            genotype = 1  # Alt
        else:
            # Mixed - use majority
            genotype = 1 if alt_depth > ref_depth else 0

        return Variant(
            chrom="Y",
            position=position,
            ref=ref,
            alt=(alt,),
            genotype=genotype,
            depth=total_depth,
            quality=None,  # No quality score from pileup
            allele_depth=(ref_depth, alt_depth),
        )

    def get_variants_at_snp_positions(
        self,
        snp_positions: dict[int, tuple[str, str]],
    ) -> dict[int, Variant]:
        """
        Get Variants for a set of SNP positions.

        This is the main method for classification - it queries the BAM
        at all known SNP positions and returns Variant objects compatible
        with the Bayesian classifier.

        Args:
            snp_positions: Dict of position -> (ref, alt) for each SNP

        Returns:
            Dict of position -> Variant for positions with coverage
        """
        if not snp_positions:
            return {}

        variants: dict[int, Variant] = {}
        positions = set(snp_positions.keys())

        for pileup in self.pileup_at_positions(positions):
            ref, alt = snp_positions[pileup.position]
            ref_depth, alt_depth = pileup.get_allele_depth(ref, alt)
            total_depth = pileup.total_depth

            # Determine genotype
            if total_depth == 0:
                continue
            elif alt_depth == 0:
                genotype = 0
            elif ref_depth == 0:
                genotype = 1
            else:
                genotype = 1 if alt_depth > ref_depth else 0

            variants[pileup.position] = Variant(
                chrom="Y",
                position=pileup.position,
                ref=ref,
                alt=(alt,),
                genotype=genotype,
                depth=total_depth,
                quality=None,
                allele_depth=(ref_depth, alt_depth),
            )

        return variants


def read_bam_variants(
    bam_path: Path | str,
    snp_positions: dict[int, tuple[str, str]],
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
) -> tuple[str, dict[int, Variant]]:
    """
    Convenience function to read variants from BAM at SNP positions.

    Args:
        bam_path: Path to BAM file
        snp_positions: Dict of position -> (ref, alt)
        min_base_quality: Minimum base quality
        min_mapping_quality: Minimum mapping quality

    Returns:
        Tuple of (sample_name, variants_dict)
    """
    with BAMReader(
        bam_path,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
    ) as reader:
        variants = reader.get_variants_at_snp_positions(snp_positions)
        return (reader.sample, variants)
