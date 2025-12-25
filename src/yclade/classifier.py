"""
Y-chromosome haplogroup classification.

Implements likelihood-based haplogroup inference with QC scoring
inspired by Yleaf and pathPhynder approaches.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from yclade.snps import ReferenceGenome, SNP, SNPDatabase
from yclade.tree import Tree
from yclade.vcf import Variant, VCFReader


@dataclass
class SNPStats:
    """Statistics about SNP calls used in classification."""

    informative_tested: int = 0
    derived: int = 0
    ancestral: int = 0
    missing: int = 0
    filtered_damage: int = 0

    @property
    def total_called(self) -> int:
        return self.derived + self.ancestral


@dataclass
class QCScores:
    """
    Quality control scores for haplogroup classification.

    qc1_backbone: Backbone consistency (intermediate markers match expected)
    qc2_terminal: Terminal marker consistency (defining markers for haplogroup)
    qc3_path: Within-haplogroup consistency (path from major to terminal)
    qc4_posterior: Posterior probability from likelihood calculation
    """

    qc1_backbone: float = 0.0
    qc2_terminal: float = 0.0
    qc3_path: float = 0.0
    qc4_posterior: float = 0.0

    @property
    def combined(self) -> float:
        """Combined QC score (geometric mean)."""
        scores = [self.qc1_backbone, self.qc2_terminal, self.qc3_path, self.qc4_posterior]
        if any(s <= 0 for s in scores):
            return 0.0
        product = 1.0
        for s in scores:
            product *= s
        return product ** (1 / len(scores))


@dataclass
class HaplogroupCall:
    """
    Result of haplogroup classification.

    Attributes:
        sample: Sample identifier
        haplogroup: Called haplogroup name
        confidence: Overall confidence score [0-1]
        qc_scores: Detailed QC scores
        path: Path from root to called haplogroup
        defining_snps: SNPs that define the called haplogroup
        alternatives: Alternative calls with posterior probabilities
        snp_stats: Statistics about SNPs used
        reference: Reference genome used
        tree_version: Version of phylogenetic tree used
    """

    sample: str
    haplogroup: str
    confidence: float
    qc_scores: QCScores = field(default_factory=QCScores)
    path: list[str] = field(default_factory=list)
    defining_snps: list[str] = field(default_factory=list)
    alternatives: list[tuple[str, float]] = field(default_factory=list)
    snp_stats: SNPStats = field(default_factory=SNPStats)
    reference: str = ""
    tree_version: str = ""

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON output."""
        return {
            "sample": self.sample,
            "haplogroup": self.haplogroup,
            "confidence": self.confidence,
            "reference": self.reference,
            "tree_version": self.tree_version,
            "snp_stats": {
                "informative_tested": self.snp_stats.informative_tested,
                "derived": self.snp_stats.derived,
                "ancestral": self.snp_stats.ancestral,
                "missing": self.snp_stats.missing,
                "filtered_damage": self.snp_stats.filtered_damage,
            },
            "quality_scores": {
                "qc1_backbone": self.qc_scores.qc1_backbone,
                "qc2_terminal": self.qc_scores.qc2_terminal,
                "qc3_path": self.qc_scores.qc3_path,
                "qc4_posterior": self.qc_scores.qc4_posterior,
            },
            "path": self.path,
            "defining_snps": self.defining_snps,
            "alternative_calls": [
                {"haplogroup": hg, "posterior": p} for hg, p in self.alternatives
            ],
        }


class HaplogroupClassifier:
    """
    Classifier for Y-chromosome haplogroups.

    Uses a likelihood-based approach with tree traversal to find
    the most probable haplogroup assignment.
    """

    def __init__(
        self,
        tree: Tree,
        snp_db: SNPDatabase,
        reference: ReferenceGenome = "grch38",
        min_depth: int = 1,
        min_quality: int = 20,
        ancient_mode: bool = False,
    ):
        """
        Initialize classifier.

        Args:
            tree: YFull phylogenetic tree
            snp_db: SNP database with position mappings
            reference: Reference genome for position lookup
            min_depth: Minimum read depth to use variant
            min_quality: Minimum genotype quality to use variant
            ancient_mode: Enable ancient DNA damage filtering
        """
        self.tree = tree
        self.snp_db = snp_db
        self.reference = reference
        self.min_depth = min_depth
        self.min_quality = min_quality
        self.ancient_mode = ancient_mode

        # Build position -> haplogroup mapping
        self._build_position_index()

    def _build_position_index(self) -> None:
        """Build index of positions to haplogroups."""
        self._pos_to_haplogroup: dict[int, list[str]] = {}

        for snp in self.snp_db:
            pos = snp.get_position(self.reference)
            if pos is not None and snp.haplogroup:
                self._pos_to_haplogroup.setdefault(pos, []).append(snp.haplogroup)

    def classify(
        self,
        vcf_path: Path | str,
        sample: str | None = None,
    ) -> HaplogroupCall:
        """
        Classify a sample's Y-chromosome haplogroup.

        Args:
            vcf_path: Path to VCF file
            sample: Sample name (optional, uses first if not specified)

        Returns:
            HaplogroupCall with classification result
        """
        # Read variants
        variants = self._read_variants(vcf_path, sample)
        actual_sample = sample or "unknown"

        # Get sample name from VCF if not specified
        with VCFReader(vcf_path, sample) as reader:
            actual_sample = reader.sample

        # Score haplogroups
        haplogroup_scores = self._score_haplogroups(variants)

        # Find best haplogroup
        best_hg, confidence, qc_scores, stats = self._find_best_haplogroup(
            haplogroup_scores, variants
        )

        # Build result
        path = self.tree.path_from_root(best_hg) if best_hg in self.tree else []

        # Get defining SNPs
        defining_snps: list[str] = []
        if best_hg in self.tree:
            defining_snps = self.tree.get(best_hg).snps

        # Get alternatives
        alternatives = self._get_alternatives(haplogroup_scores, best_hg)

        return HaplogroupCall(
            sample=actual_sample,
            haplogroup=best_hg,
            confidence=confidence,
            qc_scores=qc_scores,
            path=path,
            defining_snps=defining_snps,
            alternatives=alternatives,
            snp_stats=stats,
            reference=self.reference,
            tree_version="YFull",  # TODO: get actual version
        )

    def _read_variants(
        self, vcf_path: Path | str, sample: str | None
    ) -> dict[int, Variant]:
        """Read Y-chromosome variants, indexed by position."""
        variants: dict[int, Variant] = {}

        with VCFReader(vcf_path, sample) as reader:
            for variant in reader.iter_variants():
                # Apply quality filters
                if variant.depth is not None and variant.depth < self.min_depth:
                    continue
                if variant.quality is not None and variant.quality < self.min_quality:
                    continue
                if not variant.is_snp:
                    continue

                variants[variant.position] = variant

        return variants

    def _score_haplogroups(
        self, variants: dict[int, Variant]
    ) -> dict[str, dict[str, int]]:
        """
        Score each haplogroup based on observed variants.

        Returns dict of {haplogroup: {"derived": n, "ancestral": n, "missing": n}}
        """
        scores: dict[str, dict[str, int]] = {}

        # Get all informative positions
        informative_positions = self.snp_db.positions[self.reference]

        for pos in informative_positions:
            snps_at_pos = self.snp_db.get_by_position(pos, self.reference)

            for snp in snps_at_pos:
                hg = snp.haplogroup
                if not hg or hg not in self.tree:
                    continue

                if hg not in scores:
                    scores[hg] = {"derived": 0, "ancestral": 0, "missing": 0}

                if pos not in variants:
                    scores[hg]["missing"] += 1
                    continue

                variant = variants[pos]
                called = variant.called_allele

                if called is None:
                    scores[hg]["missing"] += 1
                elif called == snp.derived:
                    # Ancient DNA filter: skip C>T and G>A transitions
                    if self.ancient_mode and self._is_damage_like(snp, called):
                        scores[hg]["missing"] += 1  # Treat as missing
                        continue
                    scores[hg]["derived"] += 1
                elif called == snp.ancestral:
                    scores[hg]["ancestral"] += 1
                else:
                    # Discordant genotype - neither ancestral nor derived
                    scores[hg]["missing"] += 1

        return scores

    def _is_damage_like(self, snp: SNP, called: str) -> bool:
        """Check if variant looks like ancient DNA damage."""
        # C>T damage (on reference strand)
        if snp.ancestral == "C" and called == "T":
            return True
        # G>A damage (reverse complement of C>T)
        if snp.ancestral == "G" and called == "A":
            return True
        return False

    def _find_best_haplogroup(
        self,
        haplogroup_scores: dict[str, dict[str, int]],
        variants: dict[int, Variant],
    ) -> tuple[str, float, QCScores, SNPStats]:
        """
        Find the best haplogroup assignment using tree traversal.

        Returns (haplogroup, confidence, qc_scores, snp_stats)
        """
        # Sort haplogroups by depth (most specific first)
        sorted_hgs = sorted(
            haplogroup_scores.keys(),
            key=lambda hg: self.tree.get(hg).depth if hg in self.tree else 0,
            reverse=True,
        )

        best_hg = "NA"
        best_confidence = 0.0
        best_qc = QCScores()
        best_stats = SNPStats()

        covered_paths: set[str] = set()

        for hg in sorted_hgs:
            if hg not in self.tree:
                continue

            # Skip if we've already covered this path (more specific descendant found)
            path = self.tree.path_to_root(hg)
            if any(p in covered_paths for p in path[1:]):  # Skip self
                continue

            scores = haplogroup_scores[hg]
            total = scores["derived"] + scores["ancestral"]

            if total == 0:
                continue

            # Calculate QC2: terminal marker consistency
            qc2 = scores["derived"] / total if total > 0 else 0

            # Calculate QC3: path consistency
            qc3 = self._calculate_path_score(hg, haplogroup_scores)

            # Calculate QC1: backbone consistency
            qc1 = self._calculate_backbone_score(hg, haplogroup_scores)

            # Combined confidence
            if qc1 > 0 and qc2 > 0 and qc3 > 0:
                confidence = (qc1 * qc2 * qc3) ** (1 / 3)
            else:
                confidence = 0

            if confidence > best_confidence:
                best_hg = hg
                best_confidence = confidence
                best_qc = QCScores(
                    qc1_backbone=qc1,
                    qc2_terminal=qc2,
                    qc3_path=qc3,
                    qc4_posterior=confidence,  # Simplified for now
                )
                best_stats = SNPStats(
                    informative_tested=total + scores["missing"],
                    derived=scores["derived"],
                    ancestral=scores["ancestral"],
                    missing=scores["missing"],
                )

            # Mark path as covered
            for p in path:
                covered_paths.add(p)

        return best_hg, best_confidence, best_qc, best_stats

    def _calculate_path_score(
        self, haplogroup: str, scores: dict[str, dict[str, int]]
    ) -> float:
        """Calculate QC3: within-haplogroup path consistency."""
        path = self.tree.path_to_root(haplogroup)
        if len(path) <= 1:
            return 1.0

        matching = 0
        total = 0

        # Check each ancestor (skip self)
        for ancestor in path[1:]:
            if ancestor not in scores:
                continue

            # Only check ancestors in same major haplogroup
            major = haplogroup[0] if haplogroup else ""
            if major and major in ancestor:
                ancestor_scores = scores[ancestor]
                derived = ancestor_scores["derived"]
                ancestral = ancestor_scores["ancestral"]

                if derived + ancestral > 0:
                    # Ancestor should be derived
                    if derived >= ancestral:
                        matching += 1
                    total += 1

        return matching / total if total > 0 else 1.0

    def _calculate_backbone_score(
        self, haplogroup: str, scores: dict[str, dict[str, int]]
    ) -> float:
        """Calculate QC1: backbone consistency."""
        # Simplified backbone check - verify major haplogroup markers
        # In full implementation, this would check intermediate markers
        # like A0-T, BT, CT, CF, etc.

        # For now, return a simplified score based on path depth
        path = self.tree.path_to_root(haplogroup)

        matching = 0
        total = 0

        for node_name in path:
            if node_name not in scores:
                continue

            node_scores = scores[node_name]
            derived = node_scores["derived"]
            ancestral = node_scores["ancestral"]

            if derived + ancestral > 0:
                total += 1
                if derived > 0:
                    matching += 1

        return matching / total if total > 0 else 0.5

    def _get_alternatives(
        self, scores: dict[str, dict[str, int]], best_hg: str
    ) -> list[tuple[str, float]]:
        """Get alternative haplogroup calls with posterior probabilities."""
        alternatives: list[tuple[str, float]] = []

        for hg, hg_scores in scores.items():
            if hg == best_hg:
                continue

            total = hg_scores["derived"] + hg_scores["ancestral"]
            if total == 0:
                continue

            # Simple posterior approximation
            posterior = hg_scores["derived"] / total
            if posterior > 0.5:  # Only include reasonable alternatives
                alternatives.append((hg, round(posterior, 3)))

        # Sort by posterior probability
        alternatives.sort(key=lambda x: x[1], reverse=True)

        return alternatives[:5]  # Return top 5 alternatives


def classify(
    vcf_path: Path | str,
    tree: Tree | None = None,
    snp_db: SNPDatabase | None = None,
    reference: ReferenceGenome = "grch38",
    sample: str | None = None,
    ancient: bool = False,
    min_depth: int = 10,
    min_quality: int = 20,
) -> HaplogroupCall:
    """
    Convenience function for haplogroup classification.

    Args:
        vcf_path: Path to VCF file
        tree: YFull tree (loads default if None)
        snp_db: SNP database (loads default if None)
        reference: Reference genome
        sample: Sample name
        ancient: Enable ancient DNA mode
        min_depth: Minimum read depth
        min_quality: Minimum genotype quality

    Returns:
        HaplogroupCall with classification result
    """
    if tree is None:
        raise ValueError("Tree must be provided (default loading not yet implemented)")
    if snp_db is None:
        raise ValueError("SNP database must be provided (default loading not yet implemented)")

    classifier = HaplogroupClassifier(
        tree=tree,
        snp_db=snp_db,
        reference=reference,
        min_depth=min_depth,
        min_quality=min_quality,
        ancient_mode=ancient,
    )

    return classifier.classify(vcf_path, sample)
