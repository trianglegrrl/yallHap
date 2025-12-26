"""
Validation framework for Y-chromosome haplogroup classification.

Provides tools for loading ground truth data, comparing haplogroup calls,
and computing accuracy metrics against reference datasets like Poznik 2016.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from yallhap.classifier import HaplogroupClassifier
from yallhap.snps import ReferenceGenome, SNPDatabase
from yallhap.tree import Tree


def load_ground_truth(path: Path | str) -> dict[str, str]:
    """
    Load ground truth haplogroup assignments from TSV file.

    Expected format: sample_id, haplogroup, population (tab-separated)

    Args:
        path: Path to TSV file with ground truth assignments

    Returns:
        Dictionary mapping sample_id to haplogroup

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is empty or has no samples
    """
    path = Path(path)

    ground_truth: dict[str, str] = {}

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_id = row.get("sample_id", "").strip()
            haplogroup = row.get("haplogroup", "").strip()
            if sample_id and haplogroup:
                ground_truth[sample_id] = haplogroup

    if not ground_truth:
        raise ValueError(f"No samples found in ground truth file: {path}")

    return ground_truth


@dataclass
class ComparisonResult:
    """
    Result of comparing two haplogroup assignments.

    Attributes:
        exact_match: True if haplogroups are identical
        major_match: True if major haplogroup (first letter) matches
        depth_difference: Difference in haplogroup depth/specificity
    """

    exact_match: bool
    major_match: bool
    depth_difference: int


def compare_haplogroups(expected: str, called: str) -> ComparisonResult:
    """
    Compare expected and called haplogroup assignments.

    Handles various naming conventions:
    - ISOGG style: R1b1a2a1a2c1
    - YFull style: R-L21
    - Major haplogroup: R, J, E, etc.

    Args:
        expected: Expected/ground truth haplogroup
        called: Called/predicted haplogroup

    Returns:
        ComparisonResult with match information
    """
    # Handle NA or empty values
    if not expected or not called or called == "NA" or expected == "NA":
        return ComparisonResult(exact_match=False, major_match=False, depth_difference=-1)

    # Exact match check
    exact_match = expected == called

    # Extract major haplogroup (first letter, ignoring any prefix)
    expected_major = _extract_major_haplogroup(expected)
    called_major = _extract_major_haplogroup(called)

    major_match = expected_major == called_major

    # Calculate depth difference (simplified)
    expected_depth = _estimate_depth(expected)
    called_depth = _estimate_depth(called)
    depth_difference = abs(expected_depth - called_depth)

    return ComparisonResult(
        exact_match=exact_match,
        major_match=major_match,
        depth_difference=depth_difference,
    )


def _extract_major_haplogroup(haplogroup: str) -> str:
    """
    Extract major haplogroup letter from haplogroup name.

    Handles both ISOGG (R1b1a2) and YFull (R-L21) formats.
    """
    if not haplogroup:
        return ""

    # Remove any leading quotes or whitespace
    hg = haplogroup.strip().strip('"')

    # Handle YFull format like "R-L21" - take letter before dash
    if "-" in hg:
        hg = hg.split("-")[0]

    # Handle ROOT specially
    if hg.upper().startswith("ROOT"):
        return "ROOT"

    # Major haplogroup is the first letter
    for char in hg:
        if char.isalpha():
            return char.upper()

    return ""


def _estimate_depth(haplogroup: str) -> int:
    """
    Estimate haplogroup depth/specificity.

    Simple heuristic based on length and structure.
    """
    if not haplogroup:
        return 0

    # For ISOGG format like R1b1a2a1a2c1, count alphanumeric chars
    # For YFull format like R-L21, count parts

    if "-" in haplogroup:
        # YFull format - count number of dashes + 1
        return haplogroup.count("-") + 1
    else:
        # ISOGG format - count transitions between letters and numbers
        depth = 0
        prev_is_digit = False
        for char in haplogroup:
            if char.isdigit() != prev_is_digit:
                depth += 1
                prev_is_digit = char.isdigit()
        return max(depth, 1)


@dataclass
class ValidationResult:
    """
    Result of validating a single sample.

    Attributes:
        sample_id: Sample identifier
        expected_haplogroup: Ground truth haplogroup
        called_haplogroup: Haplogroup called by classifier
        confidence: Confidence score of the call
        comparison: Comparison result between expected and called
    """

    sample_id: str
    expected_haplogroup: str
    called_haplogroup: str
    confidence: float
    comparison: ComparisonResult


@dataclass
class ValidationMetrics:
    """
    Aggregate validation metrics.

    Attributes:
        total_samples: Number of samples validated
        exact_match_rate: Proportion of exact haplogroup matches
        major_match_rate: Proportion of major haplogroup matches
        mean_confidence: Mean confidence score across all samples
        mean_depth_difference: Mean depth difference for major matches
    """

    total_samples: int
    exact_match_rate: float
    major_match_rate: float
    mean_confidence: float = 0.0
    mean_depth_difference: float = 0.0


def compute_metrics(results: list[ValidationResult]) -> ValidationMetrics:
    """
    Compute aggregate validation metrics from individual results.

    Args:
        results: List of ValidationResult objects

    Returns:
        ValidationMetrics with aggregate statistics
    """
    if not results:
        return ValidationMetrics(
            total_samples=0,
            exact_match_rate=0.0,
            major_match_rate=0.0,
        )

    total = len(results)
    exact_matches = sum(1 for r in results if r.comparison.exact_match)
    major_matches = sum(1 for r in results if r.comparison.major_match)

    confidences = [r.confidence for r in results]
    mean_confidence = sum(confidences) / len(confidences) if confidences else 0.0

    # Calculate mean depth difference for major matches only
    major_match_depths = [
        r.comparison.depth_difference
        for r in results
        if r.comparison.major_match and r.comparison.depth_difference >= 0
    ]
    mean_depth_diff = (
        sum(major_match_depths) / len(major_match_depths) if major_match_depths else 0.0
    )

    return ValidationMetrics(
        total_samples=total,
        exact_match_rate=exact_matches / total,
        major_match_rate=major_matches / total,
        mean_confidence=mean_confidence,
        mean_depth_difference=mean_depth_diff,
    )


class ValidationRunner:
    """
    Runner for validating haplogroup classification against ground truth.

    Coordinates loading of tree, SNP database, ground truth, and
    running classification on samples with comparison to expected values.
    """

    def __init__(
        self,
        tree: Tree,
        snp_db: SNPDatabase,
        ground_truth: dict[str, str],
        reference: ReferenceGenome = "grch38",
    ):
        """
        Initialize validation runner.

        Args:
            tree: YFull phylogenetic tree
            snp_db: SNP database with position mappings
            ground_truth: Dictionary mapping sample_id to expected haplogroup
            reference: Reference genome for position lookup
        """
        self.tree = tree
        self.snp_db = snp_db
        self.ground_truth = ground_truth
        self.reference = reference

        self._classifier = HaplogroupClassifier(
            tree=tree,
            snp_db=snp_db,
            reference=reference,
        )

    def get_expected_haplogroup(self, sample_id: str) -> str:
        """
        Get expected haplogroup for a sample.

        Args:
            sample_id: Sample identifier

        Returns:
            Expected haplogroup from ground truth

        Raises:
            KeyError: If sample not in ground truth
        """
        if sample_id not in self.ground_truth:
            raise KeyError(f"Sample not in ground truth: {sample_id}")
        return self.ground_truth[sample_id]

    def validate_sample(self, sample_id: str, vcf_path: Path | str) -> ValidationResult:
        """
        Validate classification for a single sample.

        Args:
            sample_id: Sample identifier
            vcf_path: Path to VCF file containing the sample

        Returns:
            ValidationResult with comparison to ground truth
        """
        expected = self.get_expected_haplogroup(sample_id)

        # Run classification
        call = self._classifier.classify(vcf_path, sample=sample_id)

        # Compare results
        comparison = compare_haplogroups(expected, call.haplogroup)

        return ValidationResult(
            sample_id=sample_id,
            expected_haplogroup=expected,
            called_haplogroup=call.haplogroup,
            confidence=call.confidence,
            comparison=comparison,
        )

    def validate_all(
        self,
        vcf_path: Path | str,
        sample_ids: list[str] | None = None,
    ) -> list[ValidationResult]:
        """
        Validate classification for multiple samples.

        Args:
            vcf_path: Path to multi-sample VCF file
            sample_ids: List of sample IDs to validate (default: all in ground truth)

        Returns:
            List of ValidationResult objects
        """
        if sample_ids is None:
            sample_ids = list(self.ground_truth.keys())

        results: list[ValidationResult] = []
        for sample_id in sample_ids:
            try:
                result = self.validate_sample(sample_id, vcf_path)
                results.append(result)
            except Exception as e:
                # Create a failed result for samples that error
                results.append(
                    ValidationResult(
                        sample_id=sample_id,
                        expected_haplogroup=self.ground_truth.get(sample_id, ""),
                        called_haplogroup=f"ERROR: {e}",
                        confidence=0.0,
                        comparison=ComparisonResult(
                            exact_match=False,
                            major_match=False,
                            depth_difference=-1,
                        ),
                    )
                )

        return results
