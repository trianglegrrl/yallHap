"""
Tests for path-aware classification with tolerance stopping.

Following TDD - these tests are written first, before implementation.
"""

from __future__ import annotations

from yallhap.tree import Tree


class TestPathTraversalWithTolerance:
    """Tests for path-aware classification with tolerance stopping."""

    def test_stops_at_tolerance_threshold(self) -> None:
        """Traversal stops when ancestral calls exceed tolerance."""
        from yallhap.classifier import traverse_with_tolerance

        # Create a simple tree for testing
        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
            "R1b": ["R-L21"],
        }
        tree = Tree.from_dict(tree_dict)

        # Path with too many ancestral calls at R1b should stop at R1
        scores = {
            "R": {"derived": 5, "ancestral": 0, "missing": 0},
            "R1": {"derived": 3, "ancestral": 1, "missing": 0},
            "R1b": {"derived": 1, "ancestral": 4, "missing": 0},  # 4 > tolerance=3
            "R-L21": {"derived": 2, "ancestral": 0, "missing": 0},
        }
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=3)
        assert hg == "R1"  # Stops at R1, doesn't go to R1b
        assert "R1b" not in path

    def test_continues_below_tolerance(self) -> None:
        """Traversal continues when ancestral calls within tolerance."""
        from yallhap.classifier import traverse_with_tolerance

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
        }
        tree = Tree.from_dict(tree_dict)

        # All nodes within tolerance
        scores = {
            "R": {"derived": 5, "ancestral": 0, "missing": 0},
            "R1": {"derived": 3, "ancestral": 2, "missing": 0},  # 2 <= tolerance=3
            "R1b": {"derived": 4, "ancestral": 1, "missing": 0},  # 1 <= tolerance=3
        }
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=3)
        assert hg == "R1b"
        assert "R1b" in path

    def test_path_includes_all_traversed_nodes(self) -> None:
        """Returned path includes all nodes traversed."""
        from yallhap.classifier import traverse_with_tolerance

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
            "R1b": ["R-L21"],
        }
        tree = Tree.from_dict(tree_dict)

        scores = {
            "R": {"derived": 5, "ancestral": 0, "missing": 0},
            "R1": {"derived": 3, "ancestral": 0, "missing": 0},
            "R1b": {"derived": 2, "ancestral": 0, "missing": 0},
            "R-L21": {"derived": 1, "ancestral": 0, "missing": 0},
        }
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=3)
        assert path == ["R", "R1", "R1b", "R-L21"]
        assert hg == "R-L21"

    def test_cumulative_ancestral_tracking(self) -> None:
        """Stats include cumulative derived and ancestral counts."""
        from yallhap.classifier import traverse_with_tolerance

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
        }
        tree = Tree.from_dict(tree_dict)

        scores = {
            "R": {"derived": 5, "ancestral": 1, "missing": 0},
            "R1": {"derived": 3, "ancestral": 2, "missing": 0},
            "R1b": {"derived": 2, "ancestral": 1, "missing": 0},
        }
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=3)
        # Stats should have cumulative counts
        assert stats["total_derived"] == 10  # 5 + 3 + 2
        assert stats["total_ancestral"] == 4  # 1 + 2 + 1

    def test_no_derived_calls_returns_root(self) -> None:
        """When no derived calls, return root of tree."""
        from yallhap.classifier import traverse_with_tolerance

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
        }
        tree = Tree.from_dict(tree_dict)

        # No derived calls anywhere
        scores = {
            "R": {"derived": 0, "ancestral": 3, "missing": 0},
            "R1": {"derived": 0, "ancestral": 2, "missing": 0},
        }
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=3)
        # Should return root or first node with no evidence
        assert hg == 'ROOT (Y-Chromosome "Adam")' or len(path) == 0

    def test_multiple_paths_chooses_best(self) -> None:
        """With multiple paths, chooses the one with most derived support."""
        from yallhap.classifier import traverse_with_tolerance

        # Tree with branching
        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["I", "R"],
            "I": ["I1", "I2"],
            "R": ["R1"],
            "R1": ["R1b"],
        }
        tree = Tree.from_dict(tree_dict)

        scores = {
            "I": {"derived": 1, "ancestral": 5, "missing": 0},  # Wrong path
            "I1": {"derived": 0, "ancestral": 3, "missing": 0},
            "I2": {"derived": 0, "ancestral": 2, "missing": 0},
            "R": {"derived": 5, "ancestral": 0, "missing": 0},  # Right path
            "R1": {"derived": 4, "ancestral": 0, "missing": 0},
            "R1b": {"derived": 3, "ancestral": 0, "missing": 0},
        }
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=3)
        assert hg == "R1b"
        assert "R" in path
        assert "I" not in path

    def test_empty_scores_returns_root(self) -> None:
        """Empty scores dict returns root."""
        from yallhap.classifier import traverse_with_tolerance

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
        }
        tree = Tree.from_dict(tree_dict)

        hg, path, stats = traverse_with_tolerance(tree, {}, max_tolerance=3)
        assert hg == 'ROOT (Y-Chromosome "Adam")'
        assert stats["total_derived"] == 0


class TestTraversalEdgeCases:
    """Edge cases for path traversal."""

    def test_single_node_tree(self) -> None:
        """Tree with only root node."""
        from yallhap.classifier import traverse_with_tolerance

        tree_dict: dict[str, list[str]] = {'ROOT (Y-Chromosome "Adam")': []}
        tree = Tree.from_dict(tree_dict)

        scores: dict[str, dict[str, int]] = {}
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=3)
        assert hg == 'ROOT (Y-Chromosome "Adam")'

    def test_tolerance_zero(self) -> None:
        """With tolerance=0, any ancestral call stops traversal."""
        from yallhap.classifier import traverse_with_tolerance

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
        }
        tree = Tree.from_dict(tree_dict)

        scores = {
            "R": {"derived": 5, "ancestral": 0, "missing": 0},
            "R1": {"derived": 3, "ancestral": 1, "missing": 0},  # 1 > 0
            "R1b": {"derived": 2, "ancestral": 0, "missing": 0},
        }
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=0)
        assert hg == "R"  # Stops before R1

    def test_high_tolerance_allows_all(self) -> None:
        """With very high tolerance, allows high ancestral counts."""
        from yallhap.classifier import traverse_with_tolerance

        tree_dict = {
            'ROOT (Y-Chromosome "Adam")': ["R"],
            "R": ["R1"],
            "R1": ["R1b"],
        }
        tree = Tree.from_dict(tree_dict)

        scores = {
            "R": {"derived": 5, "ancestral": 10, "missing": 0},
            "R1": {"derived": 3, "ancestral": 20, "missing": 0},
            "R1b": {"derived": 2, "ancestral": 15, "missing": 0},
        }
        hg, path, stats = traverse_with_tolerance(tree, scores, max_tolerance=100)
        assert hg == "R1b"  # Goes all the way despite high ancestral
