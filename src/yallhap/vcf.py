"""
VCF file parsing for Y-chromosome variants.

Extracts Y-chromosome variants from VCF files with support for
both single-sample and multi-sample VCFs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import pysam


@dataclass
class Variant:
    """
    A variant call from a VCF file.

    Attributes:
        chrom: Chromosome (normalized to 'Y')
        position: 1-based genomic position
        ref: Reference allele
        alt: Alternative allele(s)
        genotype: Called genotype (0=ref, 1=alt, None=missing)
        depth: Read depth (DP)
        quality: Genotype quality (GQ)
        allele_depth: Per-allele depth (AD)
    """

    chrom: str
    position: int
    ref: str
    alt: tuple[str, ...]
    genotype: int | None = None
    depth: int | None = None
    quality: int | None = None
    allele_depth: tuple[int, ...] | None = None

    @property
    def called_allele(self) -> str | None:
        """Return the called allele based on genotype."""
        if self.genotype is None:
            return None
        if self.genotype == 0:
            return self.ref
        if self.genotype > 0 and self.genotype <= len(self.alt):
            return self.alt[self.genotype - 1]
        return None

    @property
    def is_snp(self) -> bool:
        """Check if variant is a SNP (not indel)."""
        if len(self.ref) != 1:
            return False
        return all(len(a) == 1 for a in self.alt)

    @property
    def total_allelic_depth(self) -> int | None:
        """Return total read depth from allelic depths (sum of AD field)."""
        if self.allele_depth is None:
            return None
        return sum(self.allele_depth)

    @property
    def support_ratio(self) -> float | None:
        """
        Calculate the fraction of reads supporting the called allele.

        Returns:
            Ratio of reads supporting the called genotype (0.0-1.0),
            or None if AD is missing, genotype is missing, or total depth is 0.
        """
        if self.allele_depth is None:
            return None
        if self.genotype is None:
            return None

        total = sum(self.allele_depth)
        if total == 0:
            return None

        # Get depth for the called allele
        # genotype 0 = ref (index 0), genotype 1 = first alt (index 1), etc.
        if self.genotype < len(self.allele_depth):
            called_depth = self.allele_depth[self.genotype]
            return called_depth / total

        return None

    @property
    def is_ambiguous(self) -> bool:
        """
        Check if variant call is ambiguous due to low read support.

        A call is ambiguous if the support ratio is below 0.7 (70%).
        Missing AD data is NOT considered ambiguous (we just don't know).

        Returns:
            True if support ratio is < 0.7, False otherwise
        """
        ratio = self.support_ratio
        if ratio is None:
            # Can't determine ambiguity without AD data
            return False
        return ratio < 0.7


class VCFReader:
    """
    Reader for Y-chromosome variants from VCF files.

    Handles both single-sample and multi-sample VCFs,
    with automatic detection of Y chromosome naming convention.
    """

    # Possible Y chromosome names in VCF
    Y_CHROMS = {"Y", "chrY", "y", "chry", "24"}

    def __init__(self, path: Path | str, sample: str | None = None):
        """
        Initialize VCF reader.

        Args:
            path: Path to VCF file (can be gzipped, must be indexed)
            sample: Sample name for multi-sample VCF (uses first if None)
        """
        self.path = Path(path)
        self._vcf: pysam.VariantFile | None = None
        self._sample: str | None = sample
        self._y_chrom: str | None = None

    def __enter__(self) -> VCFReader:
        self._vcf = pysam.VariantFile(str(self.path))
        self._detect_y_chrom()
        self._resolve_sample()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if self._vcf:
            self._vcf.close()

    def _detect_y_chrom(self) -> None:
        """Detect Y chromosome naming convention in VCF."""
        if self._vcf is None:
            raise RuntimeError("VCF not opened")

        contigs = set(self._vcf.header.contigs.keys())
        for name in self.Y_CHROMS:
            if name in contigs:
                self._y_chrom = name
                return

        raise ValueError(f"No Y chromosome found in VCF. Contigs: {contigs}")

    def _resolve_sample(self) -> None:
        """Resolve sample name for variant extraction."""
        if self._vcf is None:
            raise RuntimeError("VCF not opened")

        samples = list(self._vcf.header.samples)
        if not samples:
            raise ValueError("No samples found in VCF")

        if self._sample is None:
            self._sample = samples[0]
        elif self._sample not in samples:
            raise ValueError(f"Sample '{self._sample}' not found. Available: {samples}")

    @property
    def sample(self) -> str:
        """Return the sample name being read."""
        if self._sample is None:
            raise RuntimeError("Sample not resolved (call __enter__ first)")
        return self._sample

    @property
    def samples(self) -> list[str]:
        """Return all sample names in VCF."""
        if self._vcf is None:
            raise RuntimeError("VCF not opened")
        return list(self._vcf.header.samples)

    def iter_variants(self) -> Iterator[Variant]:
        """
        Iterate over Y chromosome variants.

        Yields:
            Variant objects for each Y chromosome variant
        """
        if self._vcf is None:
            raise RuntimeError("VCF not opened")
        if self._y_chrom is None:
            raise RuntimeError("Y chromosome not detected")

        for record in self._vcf.fetch(self._y_chrom):
            variant = self._record_to_variant(record)
            if variant:
                yield variant

    def iter_variants_at_positions(self, positions: set[int]) -> Iterator[Variant]:
        """
        Iterate over variants at specific positions.

        Args:
            positions: Set of positions to extract

        Yields:
            Variant objects at specified positions
        """
        for variant in self.iter_variants():
            if variant.position in positions:
                yield variant

    def _record_to_variant(self, record: pysam.VariantRecord) -> Variant | None:
        """Convert pysam record to Variant object."""
        if self._sample is None:
            return None

        sample_data = record.samples[self._sample]

        # Extract genotype
        gt = sample_data.get("GT")
        genotype: int | None = None
        if gt is not None and gt[0] is not None:
            # Y chromosome is haploid, take first allele
            genotype = gt[0]

        # Extract depth
        depth = sample_data.get("DP")

        # Extract quality
        quality = sample_data.get("GQ")

        # Extract allele depths (handle scalar or tuple)
        ad = sample_data.get("AD")
        if ad is None:
            allele_depth = None
        elif isinstance(ad, int):
            allele_depth = (ad,)
        else:
            allele_depth = tuple(ad)

        return Variant(
            chrom="Y",  # Normalize
            position=record.pos,
            ref=record.ref,
            alt=tuple(str(a) for a in record.alts) if record.alts else (),
            genotype=genotype,
            depth=depth,
            quality=quality,
            allele_depth=allele_depth,
        )

    def get_variant_at_position(self, position: int) -> Variant | None:
        """
        Get variant at specific position.

        Args:
            position: 1-based genomic position

        Returns:
            Variant at position or None if not found
        """
        if self._vcf is None:
            raise RuntimeError("VCF not opened")
        if self._y_chrom is None:
            raise RuntimeError("Y chromosome not detected")

        for record in self._vcf.fetch(self._y_chrom, position - 1, position):
            if record.pos == position:
                return self._record_to_variant(record)
        return None


def read_vcf(path: Path | str, sample: str | None = None) -> list[Variant]:
    """
    Convenience function to read all Y variants from VCF.

    Args:
        path: Path to VCF file
        sample: Sample name (optional)

    Returns:
        List of all Y chromosome variants
    """
    with VCFReader(path, sample) as reader:
        return list(reader.iter_variants())
