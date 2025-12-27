"""
Unit tests for yallhap.tree module.
"""

import json
from pathlib import Path

import pytest

from yallhap.tree import Node, Tree


class TestNode:
    """Tests for Node dataclass."""

    def test_node_creation(self) -> None:
        """Test basic node creation."""
        node = Node(name="R-L21", snps=["L21", "M529"])
        assert node.name == "R-L21"
        assert node.snps == ["L21", "M529"]
        assert node.parent_name is None
        assert node.children_names == []
        assert node.depth == 0

    def test_node_equality(self) -> None:
        """Nodes with same name are equal."""
        node1 = Node(name="R-L21")
        node2 = Node(name="R-L21", depth=5)  # Different depth
        node3 = Node(name="R1b")

        assert node1 == node2
        assert node1 != node3

    def test_node_hashable(self) -> None:
        """Nodes can be used in sets."""
        node1 = Node(name="R-L21")
        node2 = Node(name="R-L21")
        node3 = Node(name="R1b")

        nodes = {node1, node2, node3}
        assert len(nodes) == 2  # node1 and node2 hash the same


class TestTree:
    """Tests for Tree class."""

    def test_tree_from_dict(self, sample_tree_dict: dict) -> None:
        """Test tree construction from dictionary."""
        tree = Tree.from_dict(sample_tree_dict)

        assert len(tree) > 0
        assert tree.root.name == 'ROOT (Y-Chromosome "Adam")'
        assert "A00" in tree
        assert "R-L21" in tree

    def test_tree_from_json(self, tmp_path: Path, sample_tree_dict: dict) -> None:
        """Test tree construction from JSON file."""
        json_path = tmp_path / "tree.json"
        with open(json_path, "w") as f:
            json.dump(sample_tree_dict, f)

        tree = Tree.from_json(json_path)
        assert tree.root.name == 'ROOT (Y-Chromosome "Adam")'

    def test_tree_empty_raises(self) -> None:
        """Empty tree data raises ValueError."""
        with pytest.raises(ValueError, match="Empty tree data"):
            Tree.from_dict({})

    def test_tree_get_node(self, sample_tree_dict: dict) -> None:
        """Test getting node by name."""
        tree = Tree.from_dict(sample_tree_dict)

        node = tree.get("R-L21")
        assert node.name == "R-L21"

    def test_tree_get_missing_raises(self, sample_tree_dict: dict) -> None:
        """Getting missing node raises KeyError."""
        tree = Tree.from_dict(sample_tree_dict)

        with pytest.raises(KeyError):
            tree.get("NOT_EXISTS")

    def test_tree_contains(self, sample_tree_dict: dict) -> None:
        """Test __contains__ operator."""
        tree = Tree.from_dict(sample_tree_dict)

        assert "R-L21" in tree
        assert "NOT_EXISTS" not in tree

    def test_tree_depth_calculation(self, sample_tree_dict: dict) -> None:
        """Test that depth is correctly calculated."""
        tree = Tree.from_dict(sample_tree_dict)

        root = tree.root
        assert root.depth == 0

        a00 = tree.get("A00")
        assert a00.depth == 1

        r_l21 = tree.get("R-L21")
        assert r_l21.depth > 0
        assert r_l21.depth > tree.get("R").depth

    def test_path_to_root(self, sample_tree_dict: dict) -> None:
        """Test path from node to root."""
        tree = Tree.from_dict(sample_tree_dict)

        path = tree.path_to_root("R-L21")
        assert path[0] == "R-L21"
        assert path[-1] == 'ROOT (Y-Chromosome "Adam")'
        assert "R1b" in path
        assert "R" in path

    def test_path_from_root(self, sample_tree_dict: dict) -> None:
        """Test path from root to node."""
        tree = Tree.from_dict(sample_tree_dict)

        path = tree.path_from_root("R-L21")
        assert path[0] == 'ROOT (Y-Chromosome "Adam")'
        assert path[-1] == "R-L21"

    def test_get_parent(self, sample_tree_dict: dict) -> None:
        """Test getting parent node."""
        tree = Tree.from_dict(sample_tree_dict)

        parent = tree.get_parent("R-L21")
        assert parent is not None
        assert parent.name == "R1b"

        root_parent = tree.get_parent('ROOT (Y-Chromosome "Adam")')
        assert root_parent is None

    def test_get_children(self, sample_tree_dict: dict) -> None:
        """Test getting child nodes."""
        tree = Tree.from_dict(sample_tree_dict)

        children = tree.get_children("A00")
        names = [c.name for c in children]
        assert "A00a" in names
        assert "A00b" in names

    def test_iter_depth_first(self, sample_tree_dict: dict) -> None:
        """Test depth-first iteration."""
        tree = Tree.from_dict(sample_tree_dict)

        nodes = list(tree.iter_depth_first())
        assert nodes[0].name == 'ROOT (Y-Chromosome "Adam")'
        # Should visit all nodes
        assert len(nodes) == len(tree)

    def test_iter_breadth_first(self, sample_tree_dict: dict) -> None:
        """Test breadth-first iteration."""
        tree = Tree.from_dict(sample_tree_dict)

        nodes = list(tree.iter_breadth_first())
        assert nodes[0].name == 'ROOT (Y-Chromosome "Adam")'
        # Level 1 nodes should come before deeper nodes
        a00_idx = next(i for i, n in enumerate(nodes) if n.name == "A00")
        r_l21_idx = next(i for i, n in enumerate(nodes) if n.name == "R-L21")
        assert a00_idx < r_l21_idx

    def test_common_ancestor(self, sample_tree_dict: dict) -> None:
        """Test finding common ancestor."""
        tree = Tree.from_dict(sample_tree_dict)

        # A00a and A00b share A00 as common ancestor
        ancestor = tree.common_ancestor("A00a", "A00b")
        assert ancestor == "A00"

        # R-L21 and A00 share ROOT as common ancestor
        ancestor = tree.common_ancestor("R-L21", "A00")
        assert ancestor == 'ROOT (Y-Chromosome "Adam")'


class TestTreeVersionInfo:
    """Tests for tree version metadata."""

    def test_version_info_from_json(self, tmp_path: Path, sample_tree_dict: dict) -> None:
        """Test version_info property for tree loaded from JSON."""
        json_path = tmp_path / "tree.json"
        with open(json_path, "w") as f:
            json.dump(sample_tree_dict, f)

        tree = Tree.from_json(json_path)
        info = tree.version_info

        assert info["source"] == "YFull"
        assert info["node_count"] == len(tree)
        assert isinstance(info["snp_count"], int)
        assert info["file_hash"] is not None
        assert len(info["file_hash"]) == 8  # First 8 chars of SHA256

    def test_version_info_from_dict(self, sample_tree_dict: dict) -> None:
        """Test version_info for tree built from dict (no file hash)."""
        tree = Tree.from_dict(sample_tree_dict)
        info = tree.version_info

        assert info["source"] == "YFull"
        assert info["node_count"] == len(tree)
        assert info["file_hash"] is None  # No file loaded

    def test_version_string_with_hash(self, tmp_path: Path, sample_tree_dict: dict) -> None:
        """Test version_string includes hash when loaded from file."""
        json_path = tmp_path / "tree.json"
        with open(json_path, "w") as f:
            json.dump(sample_tree_dict, f)

        tree = Tree.from_json(json_path)
        version_str = tree.version_string

        assert "YFull" in version_str
        assert "SNPs" in version_str
        assert "hash:" in version_str

    def test_version_string_without_hash(self, sample_tree_dict: dict) -> None:
        """Test version_string without hash when built from dict."""
        tree = Tree.from_dict(sample_tree_dict)
        version_str = tree.version_string

        assert "YFull" in version_str
        assert "SNPs" in version_str
        assert "hash:" not in version_str

    def test_hash_consistent(self, tmp_path: Path, sample_tree_dict: dict) -> None:
        """Test hash is consistent for same content."""
        json_path = tmp_path / "tree.json"
        with open(json_path, "w") as f:
            json.dump(sample_tree_dict, f)

        tree1 = Tree.from_json(json_path)
        tree2 = Tree.from_json(json_path)

        assert tree1.version_info["file_hash"] == tree2.version_info["file_hash"]

    def test_snp_count_accurate(self) -> None:
        """Test SNP count is accurate for tree with known SNPs."""
        # Create tree with nodes that have SNPs
        tree = Tree()
        tree._nodes = {
            "A": Node(name="A", snps=["M91", "P97"]),
            "B": Node(name="B", snps=["M168"]),
            "C": Node(name="C", snps=[]),  # No SNPs
        }
        tree._root = tree._nodes["A"]

        assert tree.snp_count == 3  # 2 + 1 + 0


class TestTreeEdgeCases:
    """Edge case tests for Tree class."""

    def test_single_node_tree(self) -> None:
        """Tree with single node."""
        tree = Tree.from_dict({"ROOT": []})
        assert len(tree) == 1
        assert tree.root.name == "ROOT"

    def test_linear_tree(self) -> None:
        """Linear tree (no branching)."""
        tree = Tree.from_dict({
            "A": ["B"],
            "B": ["C"],
            "C": ["D"],
        })

        path = tree.path_to_root("D")
        assert len(path) == 4
        assert path == ["D", "C", "B", "A"]

    def test_wide_tree(self) -> None:
        """Wide tree (many children)."""
        tree = Tree.from_dict({
            "ROOT": ["A", "B", "C", "D", "E"],
        })

        children = tree.get_children("ROOT")
        assert len(children) == 5
