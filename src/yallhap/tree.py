"""
YFull tree parser and data structures.

Parses the YFull Y-chromosome phylogenetic tree (JSON format)
and provides efficient lookup and traversal operations.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator


@dataclass
class Node:
    """
    A node in the Y-chromosome phylogenetic tree.

    Attributes:
        name: Haplogroup name (e.g., "R-L21")
        snps: List of defining SNP names
        parent_name: Name of parent node (None for root)
        children_names: Names of child nodes
        depth: Distance from root (0 for root)
        formed: Estimated years before present when haplogroup formed
        tmrca: Time to most recent common ancestor
    """

    name: str
    snps: list[str] = field(default_factory=list)
    parent_name: str | None = None
    children_names: list[str] = field(default_factory=list)
    depth: int = 0
    formed: int | None = None
    tmrca: int | None = None

    def __hash__(self) -> int:
        return hash(self.name)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Node):
            return NotImplemented
        return self.name == other.name


class Tree:
    """
    Y-chromosome phylogenetic tree based on YFull data.

    Provides efficient lookup by haplogroup name and traversal operations.
    The tree is immutable after construction.
    """

    ROOT_NAME: str = 'ROOT (Y-Chromosome "Adam")'

    def __init__(self) -> None:
        self._nodes: dict[str, Node] = {}
        self._root: Node | None = None

    @classmethod
    def from_json(cls, path: Path | str) -> Tree:
        """
        Load tree from YFull JSON file.

        Args:
            path: Path to YFull tree JSON file

        Returns:
            Populated Tree instance

        Raises:
            FileNotFoundError: If file doesn't exist
            json.JSONDecodeError: If file is not valid JSON
            ValueError: If JSON structure is invalid
        """
        tree = cls()
        path = Path(path)

        with open(path) as f:
            data = json.load(f)

        tree._build_from_dict(data)
        return tree

    @classmethod
    def from_dict(cls, data: dict) -> Tree:
        """
        Build tree from dictionary (for testing or in-memory construction).

        Args:
            data: Dictionary with haplogroup names as keys and child lists as values

        Returns:
            Populated Tree instance
        """
        tree = cls()
        tree._build_from_dict(data)
        return tree

    def _build_from_dict(self, data: dict) -> None:
        """
        Internal method to construct tree from dictionary.

        Handles two formats:
        1. Flat format: {"parent_name": ["child1", "child2", ...], ...}
        2. YFull nested format: {"id": "name", "snps": [...], "children": [...], ...}
        """
        if not data:
            raise ValueError("Empty tree data")

        # Detect format: YFull nested format has 'children' key with list of dicts
        if "children" in data and isinstance(data.get("children"), list):
            self._build_from_yfull_nested(data)
        else:
            self._build_from_flat_dict(data)

    def _build_from_yfull_nested(
        self, data: dict, parent_name: str | None = None
    ) -> None:
        """
        Build tree from YFull's nested JSON format.

        YFull format:
        {
            "id": "A00",
            "snps": ["A00-M91", ...],
            "formed": 235900,
            "tmrca": 212800,
            "children": [{...}, {...}]
        }
        """
        # Get node name - root has empty id
        name = data.get("id", "")
        if not name:
            name = self.ROOT_NAME

        # Parse SNPs (can be string with "/" separators)
        snps_raw = data.get("snps", "")
        if isinstance(snps_raw, str):
            snps = [s.strip() for s in snps_raw.split("/") if s.strip()]
        elif isinstance(snps_raw, list):
            snps = snps_raw
        else:
            snps = []

        # Create node
        node = Node(
            name=name,
            snps=snps,
            parent_name=parent_name,
            formed=data.get("formed"),
            tmrca=data.get("tmrca"),
        )
        self._nodes[name] = node

        # Set as root if no parent
        if parent_name is None:
            self._root = node

        # Recursively process children
        children = data.get("children", [])
        if isinstance(children, list):
            for child_data in children:
                if isinstance(child_data, dict):
                    child_name = child_data.get("id", "")
                    if child_name:
                        node.children_names.append(child_name)
                        self._build_from_yfull_nested(child_data, parent_name=name)

        # Calculate depths after building (only from root call)
        if parent_name is None:
            self._calculate_depths(self._root, 0)

    def _build_from_flat_dict(self, data: dict) -> None:
        """
        Build tree from flat format: {"parent_name": ["child1", "child2", ...], ...}

        Used for testing and simpler tree representations.
        """
        # First pass: create all nodes
        all_names: set[str] = set()
        for parent, children in data.items():
            all_names.add(parent)
            if isinstance(children, list):
                all_names.update(children)

        for name in all_names:
            self._nodes[name] = Node(name=name)

        # Second pass: establish relationships
        for parent_name, children_names in data.items():
            if not isinstance(children_names, list):
                continue
            parent_node = self._nodes[parent_name]
            parent_node.children_names = list(children_names)

            for child_name in children_names:
                child_node = self._nodes[child_name]
                child_node.parent_name = parent_name

        # Find root (node with no parent)
        roots = [n for n in self._nodes.values() if n.parent_name is None]
        if len(roots) != 1:
            # YFull tree has ROOT as parent of A00 and A0-T
            # If no explicit root, create one
            if self.ROOT_NAME in self._nodes:
                self._root = self._nodes[self.ROOT_NAME]
            else:
                raise ValueError(
                    f"Expected 1 root node, found {len(roots)}: {[r.name for r in roots]}"
                )
        else:
            self._root = roots[0]

        # Third pass: calculate depths
        self._calculate_depths(self._root, 0)

    def _calculate_depths(self, node: Node, depth: int) -> None:
        """Recursively calculate depth for each node."""
        node.depth = depth
        for child_name in node.children_names:
            self._calculate_depths(self._nodes[child_name], depth + 1)

    @property
    def root(self) -> Node:
        """Return the root node of the tree."""
        if self._root is None:
            raise ValueError("Tree not initialized")
        return self._root

    @property
    def nodes(self) -> dict[str, Node]:
        """Return dictionary of all nodes keyed by name."""
        return self._nodes

    def get(self, name: str) -> Node:
        """
        Get node by haplogroup name.

        Args:
            name: Haplogroup name (e.g., "R-L21")

        Returns:
            Node with matching name

        Raises:
            KeyError: If haplogroup not found
        """
        return self._nodes[name]

    def __contains__(self, name: str) -> bool:
        """Check if haplogroup exists in tree."""
        return name in self._nodes

    def __len__(self) -> int:
        """Return number of nodes in tree."""
        return len(self._nodes)

    def path_to_root(self, name: str) -> list[str]:
        """
        Get path from node to root.

        Args:
            name: Starting haplogroup name

        Returns:
            List of haplogroup names from node to root (inclusive)

        Raises:
            KeyError: If haplogroup not found
        """
        path: list[str] = []
        node = self.get(name)

        while node is not None:
            path.append(node.name)
            if node.parent_name is None:
                break
            node = self._nodes[node.parent_name]

        return path

    def path_from_root(self, name: str) -> list[str]:
        """
        Get path from root to node.

        Args:
            name: Target haplogroup name

        Returns:
            List of haplogroup names from root to node (inclusive)
        """
        return list(reversed(self.path_to_root(name)))

    def get_parent(self, name: str) -> Node | None:
        """
        Get parent node.

        Args:
            name: Haplogroup name

        Returns:
            Parent Node or None if root
        """
        node = self.get(name)
        if node.parent_name is None:
            return None
        return self._nodes[node.parent_name]

    def get_children(self, name: str) -> list[Node]:
        """
        Get child nodes.

        Args:
            name: Haplogroup name

        Returns:
            List of child Nodes
        """
        node = self.get(name)
        return [self._nodes[child_name] for child_name in node.children_names]

    def iter_depth_first(self, start: str | None = None) -> Iterator[Node]:
        """
        Iterate through tree in depth-first order.

        Args:
            start: Starting node name (defaults to root)

        Yields:
            Nodes in depth-first order
        """
        if start is None:
            start_node = self.root
        else:
            start_node = self.get(start)

        stack = [start_node]
        while stack:
            node = stack.pop()
            yield node
            # Add children in reverse order so leftmost is processed first
            for child_name in reversed(node.children_names):
                stack.append(self._nodes[child_name])

    def iter_breadth_first(self, start: str | None = None) -> Iterator[Node]:
        """
        Iterate through tree in breadth-first order.

        Args:
            start: Starting node name (defaults to root)

        Yields:
            Nodes in breadth-first order
        """
        from collections import deque

        if start is None:
            start_node = self.root
        else:
            start_node = self.get(start)

        queue: deque[Node] = deque([start_node])
        while queue:
            node = queue.popleft()
            yield node
            for child_name in node.children_names:
                queue.append(self._nodes[child_name])

    def common_ancestor(self, name1: str, name2: str) -> str:
        """
        Find most recent common ancestor of two haplogroups.

        Args:
            name1: First haplogroup name
            name2: Second haplogroup name

        Returns:
            Name of most recent common ancestor
        """
        path1 = set(self.path_to_root(name1))
        for ancestor in self.path_to_root(name2):
            if ancestor in path1:
                return ancestor
        return self.root.name  # Should never reach here if tree is valid
