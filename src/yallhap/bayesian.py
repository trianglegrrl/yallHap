"""
Bayesian haplogroup classification.

Provides true probabilistic classification with:
- Branch likelihood calculation
- Path likelihood accumulation
- Posterior probability computation
- 95% credible set calculation
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from yallhap.scoring import score_variant_with_ad
from yallhap.snps import SNP, ReferenceGenome, SNPDatabase
from yallhap.tree import Tree
from yallhap.vcf import Variant

# Type alias for prior types
PriorType = Literal["uniform", "coalescent"]


def compute_coalescent_prior(
    tree: Tree,
    prior_type: PriorType = "uniform",
) -> dict[str, float]:
    """
    Compute prior probabilities for all haplogroups in tree.

    Based on pathPhynder's calculate.prior function. For coalescent prior,
    edge lengths (number of SNPs) are used to weight the prior probability.
    Longer edges (more mutations) have higher prior probability as they
    represent more "sampling space" in the coalescent.

    Args:
        tree: Y-chromosome phylogenetic tree
        prior_type: Type of prior to compute:
            - "uniform": Equal weight to all haplogroups
            - "coalescent": Weight proportional to edge length (SNP count)

    Returns:
        Dict of haplogroup -> prior probability (sums to 1.0)
    """
    all_nodes = list(tree.nodes.keys())
    n = len(all_nodes)

    if n == 0:
        return {}

    if prior_type == "uniform":
        # Equal weight to all haplogroups
        prior = 1.0 / n
        return dict.fromkeys(all_nodes, prior)

    # Coalescent prior: weight by edge length (SNP count)
    # If no SNP counts available, fall back to uniform
    edge_lengths: dict[str, float] = {}
    has_snp_counts = False

    for node_name in all_nodes:
        node = tree.get(node_name)
        # Try to get SNP count (edge length)
        snp_count = getattr(node, "_snp_count", None)
        if snp_count is not None and snp_count > 0:
            has_snp_counts = True
            edge_lengths[node_name] = float(snp_count)
        else:
            # Default edge length of 1
            edge_lengths[node_name] = 1.0

    if not has_snp_counts:
        # No edge length info - fall back to uniform
        prior = 1.0 / n
        return dict.fromkeys(all_nodes, prior)

    # Normalize edge lengths to probabilities
    total_length = sum(edge_lengths.values())
    if total_length == 0:
        total_length = 1.0

    return {hg: length / total_length for hg, length in edge_lengths.items()}


@dataclass
class BranchLikelihood:
    """
    Likelihood contribution from a single branch.

    Attributes:
        branch_id: Haplogroup name for this branch
        log_likelihood: Log-likelihood contribution
        derived_count: Number of derived calls on this branch
        ancestral_count: Number of ancestral calls on this branch
    """

    branch_id: str
    log_likelihood: float
    derived_count: int
    ancestral_count: int


def compute_branch_likelihood(
    branch_id: str,
    snps: list[SNP],
    variants: dict[int, Variant],
    error_rate: float = 0.001,
    damage_rate: float = 0.1,
    is_ancient: bool = False,
    reference: ReferenceGenome = "grch38",
) -> BranchLikelihood:
    """
    Compute log-likelihood RATIO for a single branch.

    Uses a Bayes factor approach comparing:
    - P(data | descended through this branch) vs
    - P(data | NOT descended through this branch)

    The log-likelihood ratio:
    - Derived match: log((1-e)/e) ≈ +7 (strong support for being on branch)
    - Ancestral mismatch: log(e/(1-e)) ≈ -7 (strong evidence against)

    This ensures that derived matches ADD positive evidence,
    allowing deep haplogroups to outscore shallow ones.

    Args:
        branch_id: Haplogroup name
        snps: List of SNPs defining this branch
        variants: Dict of position -> Variant for observed data
        error_rate: Expected sequencing error rate
        damage_rate: Expected ancient DNA damage rate
        is_ancient: Whether to apply damage modeling
        reference: Reference genome for position lookup

    Returns:
        BranchLikelihood with log-likelihood ratio and counts
    """
    if not snps:
        # No SNPs for this branch = neutral contribution
        return BranchLikelihood(
            branch_id=branch_id,
            log_likelihood=0.0,
            derived_count=0,
            ancestral_count=0,
        )

    total_log_lik = 0.0
    derived_count = 0
    ancestral_count = 0

    # For ancient samples, use higher error rate to account for damage
    effective_error_rate = error_rate
    if is_ancient:
        # Ancient samples have higher error rate due to damage
        # Be more lenient with mismatches
        effective_error_rate = max(error_rate, 0.05)

    # Precompute log-likelihood ratios
    # On branch: expect derived with prob (1-e), ancestral with prob e
    # Off branch: expect ancestral with prob (1-e), derived with prob e
    # Log-ratio for derived: log((1-e)/e)
    # Log-ratio for ancestral: log(e/(1-e)) = -log((1-e)/e)
    log_ratio_derived = math.log((1.0 - effective_error_rate) / effective_error_rate)
    log_ratio_ancestral = -log_ratio_derived

    for snp in snps:
        position = snp.get_position(reference)
        if position is None:
            continue

        variant = variants.get(position)
        if variant is None:
            # No data at this position - neutral (ratio = 1, log = 0)
            continue

        # Score this variant
        score = score_variant_with_ad(
            variant=variant,
            snp=snp,
            error_rate=error_rate,
            damage_rate=damage_rate,
            is_ancient=is_ancient,
        )

        # Weight the evidence by call quality
        weight = max(0.1, score.weight)

        if score.is_derived:
            derived_count += 1
            # Derived call SUPPORTS being on this branch
            # Add positive log-likelihood ratio, weighted by quality
            total_log_lik += log_ratio_derived * weight
        else:
            ancestral_count += 1
            # Ancestral call CONTRADICTS being on this branch
            # Add negative log-likelihood ratio, weighted by quality
            total_log_lik += log_ratio_ancestral * weight

    return BranchLikelihood(
        branch_id=branch_id,
        log_likelihood=total_log_lik,
        derived_count=derived_count,
        ancestral_count=ancestral_count,
    )


def compute_path_likelihood(
    tree: Tree,
    haplogroup: str,
    branch_likelihoods: dict[str, BranchLikelihood],
) -> float:
    """
    Compute total log-likelihood for a path from root to haplogroup.

    The path likelihood is the sum of branch log-likelihoods along
    the path from root to the target haplogroup.

    Args:
        tree: The phylogenetic tree
        haplogroup: Target haplogroup name
        branch_likelihoods: Pre-computed branch likelihoods

    Returns:
        Total log-likelihood for the path
    """
    path = tree.path_to_root(haplogroup)
    total_log_lik = 0.0

    for node_name in path:
        if node_name in branch_likelihoods:
            total_log_lik += branch_likelihoods[node_name].log_likelihood

    return total_log_lik


class BayesianClassifier:
    """
    Bayesian haplogroup classifier using tree-aware likelihood calculation.

    Computes true posterior probabilities for haplogroup assignments
    by evaluating path likelihoods through the phylogenetic tree.
    """

    def __init__(
        self,
        tree: Tree,
        snp_db: SNPDatabase,
        error_rate: float = 0.001,
        damage_rate: float = 0.1,
        reference: ReferenceGenome = "grch38",
        use_haplogroup_prior: bool = False,
        prior_type: PriorType = "uniform",
    ):
        """
        Initialize Bayesian classifier.

        Args:
            tree: Y-chromosome phylogenetic tree
            snp_db: SNP database with haplogroup markers
            error_rate: Expected sequencing error rate
            damage_rate: Expected ancient DNA damage rate (for is_ancient mode)
            reference: Reference genome for position lookup
            use_haplogroup_prior: Whether to use population-based priors (deprecated)
            prior_type: Type of prior to use ("uniform" or "coalescent")
        """
        self.tree = tree
        self.snp_db = snp_db
        self.error_rate = error_rate
        self.damage_rate = damage_rate
        self.reference = reference
        self.use_haplogroup_prior = use_haplogroup_prior
        self.prior_type = prior_type

        # Compute prior based on prior_type
        self._prior = compute_coalescent_prior(tree, prior_type)

        # Build mapping from haplogroup -> SNPs
        # Strategy:
        # 1. First try tree's SNP assignments (node.snps) - this is authoritative
        # 2. Fall back to SNP database haplogroup field for test fixtures
        self._haplogroup_snps: dict[str, list[SNP]] = {}

        # Build SNP name to SNP lookup from database
        snp_name_to_snp: dict[str, SNP] = {}
        for snp in snp_db:
            snp_name_to_snp[snp.name] = snp
            # Handle comma-separated aliases
            for alias in snp.aliases:
                snp_name_to_snp[alias] = snp

        # Map each tree node's SNPs to actual SNP objects
        tree_has_snps = False
        for node in tree.iter_depth_first():
            snps_for_node: list[SNP] = []
            for snp_name in node.snps:
                tree_has_snps = True
                # Handle comma-separated SNP names (e.g., "MF48436, PF6419")
                for name_part in snp_name.split(","):
                    name_part = name_part.strip()
                    if name_part in snp_name_to_snp:
                        snps_for_node.append(snp_name_to_snp[name_part])
            if snps_for_node:
                self._haplogroup_snps[node.name] = snps_for_node

        # Fallback: if tree has no SNPs defined, use SNP database haplogroup field
        # This supports test fixtures and simple databases
        if not tree_has_snps:
            for snp in snp_db:
                hg = snp.haplogroup
                if hg and hg in tree:
                    if hg not in self._haplogroup_snps:
                        self._haplogroup_snps[hg] = []
                    self._haplogroup_snps[hg].append(snp)

    def _get_leaf_haplogroups(self) -> list[str]:
        """Get all leaf haplogroups (terminal nodes) in the tree."""
        leaves = []
        for node in self.tree.iter_depth_first():
            if not node.children_names:
                leaves.append(node.name)
        return leaves

    def _get_all_haplogroups(self) -> list[str]:
        """Get all haplogroups in the tree."""
        return list(self.tree.nodes.keys())

    def _compute_all_branch_likelihoods(
        self,
        variants: dict[int, Variant],
        is_ancient: bool = False,
    ) -> dict[str, BranchLikelihood]:
        """
        Compute branch likelihoods for all haplogroups.

        Args:
            variants: Dict of position -> Variant
            is_ancient: Whether to apply damage modeling

        Returns:
            Dict of haplogroup -> BranchLikelihood
        """
        branch_likelihoods: dict[str, BranchLikelihood] = {}

        for haplogroup in self.tree.nodes:
            snps = self._haplogroup_snps.get(haplogroup, [])
            bl = compute_branch_likelihood(
                branch_id=haplogroup,
                snps=snps,
                variants=variants,
                error_rate=self.error_rate,
                damage_rate=self.damage_rate if is_ancient else 0.0,
                is_ancient=is_ancient,
                reference=self.reference,
            )
            branch_likelihoods[haplogroup] = bl

        return branch_likelihoods

    def compute_posteriors(
        self,
        variants: dict[int, Variant],
        is_ancient: bool = False,
    ) -> list[tuple[str, float]]:
        """
        Compute posterior probabilities for all haplogroups.

        Args:
            variants: Dict of position -> Variant for observed data
            is_ancient: Whether to apply ancient DNA damage modeling

        Returns:
            List of (haplogroup, posterior) tuples, sorted by posterior descending
        """
        # Compute all branch likelihoods
        branch_likelihoods = self._compute_all_branch_likelihoods(variants, is_ancient)

        # Compute path likelihoods for all haplogroups
        path_log_likelihoods: dict[str, float] = {}

        for haplogroup in self.tree.nodes:
            path_ll = compute_path_likelihood(self.tree, haplogroup, branch_likelihoods)
            path_log_likelihoods[haplogroup] = path_ll

        # Convert log-likelihoods to posteriors using log-sum-exp trick
        # for numerical stability
        log_liks = list(path_log_likelihoods.values())
        if not log_liks:
            return []

        max_log_lik = max(log_liks)

        # Compute unnormalized posteriors with prior
        unnormalized: dict[str, float] = {}
        for hg, log_lik in path_log_likelihoods.items():
            # Get prior for this haplogroup
            prior = self._prior.get(hg, 1.0 / len(path_log_likelihoods))

            # Posterior ∝ likelihood × prior
            # In log space: log(posterior) ∝ log(likelihood) + log(prior)
            log_prior = math.log(prior) if prior > 0 else -100.0

            unnormalized[hg] = math.exp(log_lik - max_log_lik + log_prior)

        # Normalize to sum to 1
        total = sum(unnormalized.values())
        if total == 0:
            # All zero = uniform distribution
            n = len(unnormalized)
            posteriors = [(hg, 1.0 / n) for hg in unnormalized]
        else:
            posteriors = [(hg, prob / total) for hg, prob in unnormalized.items()]

        # Sort by posterior descending
        posteriors.sort(key=lambda x: x[1], reverse=True)

        return posteriors

    def get_credible_set(
        self,
        posteriors: list[tuple[str, float]],
        threshold: float = 0.95,
    ) -> list[str]:
        """
        Get the smallest set of haplogroups containing given probability mass.

        Args:
            posteriors: List of (haplogroup, posterior) tuples, sorted descending
            threshold: Cumulative probability threshold (default 0.95)

        Returns:
            List of haplogroup names in the credible set
        """
        credible_set: list[str] = []
        cumulative = 0.0

        for hg, prob in posteriors:
            credible_set.append(hg)
            cumulative += prob
            if cumulative >= threshold:
                break

        return credible_set

    def classify(
        self,
        variants: dict[int, Variant],
        is_ancient: bool = False,
    ) -> tuple[str, float, list[str]]:
        """
        Classify sample to best haplogroup with posterior and credible set.

        Args:
            variants: Dict of position -> Variant
            is_ancient: Whether to apply ancient DNA damage modeling

        Returns:
            Tuple of (best_haplogroup, posterior_probability, credible_set_95)
        """
        posteriors = self.compute_posteriors(variants, is_ancient)

        if not posteriors:
            return (self.tree.root.name, 0.0, [self.tree.root.name])

        best_hg, best_prob = posteriors[0]
        credible_set = self.get_credible_set(posteriors, threshold=0.95)

        return (best_hg, best_prob, credible_set)
