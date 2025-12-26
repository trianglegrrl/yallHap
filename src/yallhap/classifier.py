"""
Y-chromosome haplogroup classification.

Implements likelihood-based haplogroup inference with QC scoring
inspired by Yleaf and pathPhynder approaches.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from yallhap.ancient import DamageRescaleMode, apply_damage_rescale
from yallhap.snps import SNP, ReferenceGenome, SNPDatabase
from yallhap.tree import Tree
from yallhap.vcf import Variant, VCFReader


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
        transversions_only: bool = False,
        damage_rescale: DamageRescaleMode = "none",
    ):
        """
        Initialize classifier.

        Args:
            tree: YFull phylogenetic tree
            snp_db: SNP database with position mappings
            reference: Reference genome for position lookup
            min_depth: Minimum read depth to use variant
            min_quality: Minimum genotype quality to use variant
            ancient_mode: Enable ancient DNA damage filtering (C>T, G>A)
            transversions_only: Only use transversions (strictest ancient mode)
            damage_rescale: Quality score rescaling for potentially damaged variants
        """
        self.tree = tree
        self.snp_db = snp_db
        self.reference = reference
        self.min_depth = min_depth
        self.min_quality = min_quality
        self.ancient_mode = ancient_mode
        self.transversions_only = transversions_only
        self.damage_rescale = damage_rescale

        # Build position -> haplogroup mapping
        self._build_position_index()

    def _build_position_index(self) -> None:
        """
        Build indices for SNP lookups.

        Creates:
        - _snp_name_to_pos: Maps SNP name → position
        - _pos_to_snp: Maps position → SNP object
        - _node_to_positions: Maps tree node name → list of positions for its SNPs
        """
        self._snp_name_to_pos: dict[str, int] = {}
        self._pos_to_snp: dict[int, SNP] = {}
        self._node_to_positions: dict[str, list[int]] = {}

        # Build SNP name → position index from SNP database
        for snp in self.snp_db:
            pos = snp.get_position(self.reference)
            if pos is not None:
                # Index by primary name and aliases
                self._snp_name_to_pos[snp.name] = pos
                self._pos_to_snp[pos] = snp
                for alias in snp.aliases:
                    self._snp_name_to_pos[alias] = pos

        # Build node → positions index from tree
        for node in self.tree.iter_depth_first():
            positions: list[int] = []
            for snp_group in node.snps:
                # SNPs in tree can be comma-separated
                for snp_name in snp_group.split(","):
                    snp_name = snp_name.strip()
                    if snp_name and snp_name in self._snp_name_to_pos:
                        positions.append(self._snp_name_to_pos[snp_name])
            if positions:
                self._node_to_positions[node.name] = positions

    def classify_batch(
        self,
        vcf_path: Path | str,
        samples: list[str],
    ) -> list[HaplogroupCall]:
        """
        Classify multiple samples from a single VCF file.

        Opens the VCF once and reads all samples in a single pass,
        which is much faster for multi-sample VCFs.

        Args:
            vcf_path: Path to VCF file
            samples: List of sample names to classify

        Returns:
            List of HaplogroupCall results, one per sample
        """
        import pysam

        vcf_path = Path(vcf_path)
        vcf = pysam.VariantFile(str(vcf_path))

        # Get sample indices
        vcf_samples = list(vcf.header.samples)
        sample_indices = {}
        for s in samples:
            if s in vcf_samples:
                sample_indices[s] = vcf_samples.index(s)
            else:
                raise ValueError(f"Sample {s} not found in VCF")

        # Detect Y chromosome name
        y_chrom = None
        for name in ["Y", "chrY", "y", "24"]:
            if name in vcf.header.contigs:
                y_chrom = name
                break

        if y_chrom is None:
            raise ValueError("No Y chromosome found in VCF")

        # Read all variants for all samples in one pass
        sample_variants: dict[str, dict[int, Variant]] = {s: {} for s in samples}

        for record in vcf.fetch(y_chrom):
            pos = record.pos
            ref = record.ref
            alts = record.alts or ()

            # Skip non-SNPs
            if len(ref) != 1 or not all(len(a) == 1 for a in alts):
                continue

            for sample_name, idx in sample_indices.items():
                sample_data = record.samples[idx]
                gt = sample_data.get("GT", (None,))

                # Determine genotype (haploid Y, take first allele)
                genotype = None if gt is None or gt[0] is None else gt[0]

                variant = Variant(
                    chrom=y_chrom,
                    position=pos,
                    ref=ref,
                    alt=alts,
                    genotype=genotype,
                    depth=sample_data.get("DP"),
                    quality=sample_data.get("GQ"),
                )
                sample_variants[sample_name][pos] = variant

        vcf.close()

        # Now classify each sample
        results = []
        for sample_name in samples:
            variants = sample_variants[sample_name]
            haplogroup_scores = self._score_haplogroups(variants)
            best_hg, confidence, qc_scores, stats = self._find_best_haplogroup(
                haplogroup_scores, variants
            )

            path = self.tree.path_from_root(best_hg) if best_hg in self.tree else []
            defining_snps: list[str] = []
            if best_hg in self.tree:
                defining_snps = self.tree.get(best_hg).snps
            alternatives = self._get_alternatives(haplogroup_scores, best_hg)

            results.append(
                HaplogroupCall(
                    sample=sample_name,
                    haplogroup=best_hg,
                    confidence=confidence,
                    qc_scores=qc_scores,
                    path=path,
                    defining_snps=defining_snps,
                    alternatives=alternatives,
                    snp_stats=stats,
                    reference=self.reference,
                    tree_version="YFull",
                )
            )

        return results

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

        Uses the tree structure to map nodes to their defining SNP positions,
        then checks whether variants show derived or ancestral alleles.

        Returns dict of {haplogroup: {"derived": n, "ancestral": n, "missing": n}}
        """
        from yallhap.ancient import is_transversion

        scores: dict[str, dict[str, int]] = {}

        # Score each tree node based on its defining SNPs
        for node_name, positions in self._node_to_positions.items():
            if node_name not in scores:
                scores[node_name] = {"derived": 0, "ancestral": 0, "missing": 0}

            for pos in positions:
                snp = self._pos_to_snp.get(pos)
                if snp is None:
                    continue

                if pos not in variants:
                    scores[node_name]["missing"] += 1
                    continue

                variant = variants[pos]
                called = variant.called_allele

                if called is None:
                    scores[node_name]["missing"] += 1
                    continue

                # Transversions-only mode: skip all transitions
                if self.transversions_only and not is_transversion(
                    snp.ancestral, snp.derived
                ):
                    scores[node_name]["missing"] += 1
                    continue

                # Apply damage rescaling if enabled
                if self.damage_rescale != "none" and variant.quality is not None:
                    rescaled_quality = apply_damage_rescale(
                        variant.quality,
                        snp.ancestral,
                        called,
                        mode=self.damage_rescale,
                    )
                    # If rescaled quality is too low, treat as missing
                    if rescaled_quality < self.min_quality:
                        scores[node_name]["missing"] += 1
                        continue

                if called == snp.derived:
                    # Ancient DNA filter: skip C>T and G>A transitions
                    if self.ancient_mode and self._is_damage_like(snp, called):
                        scores[node_name]["missing"] += 1  # Treat as missing
                        continue
                    scores[node_name]["derived"] += 1
                elif called == snp.ancestral:
                    scores[node_name]["ancestral"] += 1
                else:
                    # Discordant genotype - neither ancestral nor derived
                    scores[node_name]["missing"] += 1

        return scores

    def _is_damage_like(self, snp: SNP, called: str) -> bool:
        """Check if variant looks like ancient DNA damage."""
        # C>T damage (on reference strand) or G>A (reverse complement)
        return (snp.ancestral == "C" and called == "T") or (
            snp.ancestral == "G" and called == "A"
        )

    def _find_best_haplogroup(
        self,
        haplogroup_scores: dict[str, dict[str, int]],
        _variants: dict[int, Variant],  # Reserved for future QC calculations
    ) -> tuple[str, float, QCScores, SNPStats]:
        """
        Find the best haplogroup assignment.

        Strategy:
        1. Evaluate all haplogroups that have derived calls
        2. Score each by: derived ratio, path consistency, and specificity
        3. Pick the most specific haplogroup with high confidence

        Returns (haplogroup, confidence, qc_scores, snp_stats)
        """
        # Find all candidate haplogroups with derived calls
        candidates: list[tuple[str, float, QCScores, SNPStats]] = []

        for hg, scores in haplogroup_scores.items():
            if hg not in self.tree:
                continue

            # Must have at least some derived calls
            if scores["derived"] == 0:
                continue

            total = scores["derived"] + scores["ancestral"]
            if total == 0:
                continue

            # Calculate QC2: terminal marker consistency (derived ratio)
            qc2 = scores["derived"] / total

            # Calculate QC3: path consistency
            qc3 = self._calculate_path_score(hg, haplogroup_scores)

            # Calculate QC1: backbone consistency
            qc1 = self._calculate_backbone_score(hg, haplogroup_scores)

            # Combined confidence (geometric mean)
            confidence = (
                (qc1 * qc2 * qc3) ** (1 / 3) if qc1 > 0 and qc2 > 0 and qc3 > 0 else 0
            )

            if confidence > 0:
                qc = QCScores(
                    qc1_backbone=qc1,
                    qc2_terminal=qc2,
                    qc3_path=qc3,
                    qc4_posterior=confidence,
                )
                stats = SNPStats(
                    informative_tested=total + scores["missing"],
                    derived=scores["derived"],
                    ancestral=scores["ancestral"],
                    missing=scores["missing"],
                )
                candidates.append((hg, confidence, qc, stats))

        if not candidates:
            return "NA", 0.0, QCScores(), SNPStats()

        # Strategy: Balance specificity (depth) with evidence AND path consistency.
        # A haplogroup with ancestral calls in its path shouldn't be picked over
        # one with a clean path, even if the former is deeper.
        #
        # Key insight: We want the DEEPEST haplogroup that has GOOD path consistency.
        # Path consistency (backbone score) filters out wrong lineages.

        # Filter: Balance backbone score (path consistency) with evidence
        # Try progressively looser thresholds until we find good candidates

        def filter_candidates(
            cands: list, min_backbone: float, min_derived: int
        ) -> list:
            return [
                c
                for c in cands
                if c[2].qc1_backbone >= min_backbone and c[3].derived >= min_derived
            ]

        # Try strict thresholds first, then relax
        selected = (
            filter_candidates(candidates, 0.9, 5)
            or filter_candidates(candidates, 0.9, 3)
            or filter_candidates(candidates, 0.8, 3)
            or filter_candidates(candidates, 0.8, 2)
            or filter_candidates(candidates, 0.7, 2)
            or [c for c in candidates if c[1] >= 0.5]
            or candidates
        )

        # Strategy: Pick the deepest haplogroup that has good evidence AND
        # whose ancestors also have good evidence.
        #
        # We score each candidate by checking if its parent haplogroups have derived calls.
        # A deep haplogroup with no evidence in its parents is likely a false positive.

        def score_lineage_support(hg: str) -> int:
            """Count how many ancestors have derived calls."""
            path = list(self.tree.path_to_root(hg))
            support = 0
            for ancestor in path[1:]:  # Skip the haplogroup itself
                if (
                    ancestor in haplogroup_scores
                    and haplogroup_scores[ancestor]["derived"] >= 3
                ):
                    support += 1
            return support

        # Add lineage support to selection criteria
        selected.sort(
            key=lambda c: (
                score_lineage_support(c[0]),  # lineage support - PRIMARY
                self.tree.get(c[0]).depth,  # depth - SECONDARY
                c[3].derived,  # derived count - TERTIARY
            ),
            reverse=True,
        )

        # Return the best candidate
        best_hg, best_confidence, best_qc, best_stats = selected[0]
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
