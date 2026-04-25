"""Newick tree parsing and tree structure utilities for phylogenetic analysis.

Implements standards-compliant Newick format parsing with error handling and
validation.
"""

from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class TreeNode:
    """Single node in a phylogenetic tree."""

    label: Optional[str] = None  # Tip label or internal node ID
    branch_length: Optional[float] = None
    children: list[TreeNode] = field(default_factory=list)
    is_leaf: bool = False

    def get_all_leaves(self) -> list[str]:
        """Get all tip labels from this node and its descendants."""
        leaves = []
        if self.is_leaf and self.label:
            leaves.append(self.label)
        for child in self.children:
            leaves.extend(child.get_all_leaves())
        return leaves

    def get_clade_size(self) -> int:
        """Count tips in this clade."""
        if self.is_leaf:
            return 1
        return sum(child.get_clade_size() for child in self.children)

    def to_newick(self) -> str:
        """Serialize node and descendants to Newick format."""
        if self.is_leaf:
            result = self.label or ""
        else:
            subtrees = [child.to_newick() for child in self.children]
            result = f"({','.join(subtrees)})"
            if self.label:
                result += self.label

        if self.branch_length is not None:
            result += f":{self.branch_length}"

        return result


@dataclass
class Tree:
    """Complete phylogenetic tree structure."""

    root: TreeNode
    num_tips: int = 0

    @classmethod
    def from_newick(cls, newick_string: str) -> Tree:
        """Parse Newick format tree string into Tree object.

        Args:
            newick_string: Newick format string (may omit trailing semicolon)

        Returns:
            Tree object with parsed structure

        Raises:
            ValueError: If Newick format is invalid
        """
        # Clean input
        newick_string = newick_string.strip()
        if newick_string.endswith(";"):
            newick_string = newick_string[:-1]

        if not newick_string:
            raise ValueError("Empty Newick string")

        if not (newick_string.startswith("(") or newick_string[0].isalnum() or newick_string[0] == "'"):
            raise ValueError(f"Invalid Newick format: must start with '(' or label, got: {newick_string[0]}")

        # Parse
        root, pos = _parse_newick_node(newick_string, 0)

        # Validate consumed entire string
        if pos < len(newick_string):
            remaining = newick_string[pos:].strip()
            if remaining:
                raise ValueError(f"Unexpected characters after tree: {remaining}")

        # Count tips
        num_tips = len(root.get_all_leaves())

        return cls(root=root, num_tips=num_tips)

    def get_all_tips(self) -> list[str]:
        """Get all tip labels in tree."""
        return self.root.get_all_leaves()

    def to_newick(self) -> str:
        """Serialize tree back to Newick format with trailing semicolon."""
        return self.root.to_newick() + ";"

    def validate_tips(self, accessions: list[str]) -> tuple[list[str], list[str]]:
        """Check which accessions can be mapped to tree tips.

        Args:
            accessions: List of accession identifiers

        Returns:
            Tuple of (mappable_accessions, unmappable_accessions)
        """
        tips = set(self.get_all_tips())
        mappable = [a for a in accessions if a in tips]
        unmappable = [a for a in accessions if a not in tips]
        return mappable, unmappable


def _parse_newick_node(s: str, pos: int) -> tuple[TreeNode, int]:
    """Recursively parse Newick node starting at position pos."""
    node = TreeNode()

    # Check if internal node (starts with '(')
    if pos < len(s) and s[pos] == "(":
        pos += 1  # Skip '('

        # Parse children
        while pos < len(s) and s[pos] != ")":
            if s[pos] == ",":
                pos += 1  # Skip comma
                continue
            child, pos = _parse_newick_node(s, pos)
            node.children.append(child)

        if pos >= len(s) or s[pos] != ")":
            raise ValueError(f"Expected ')' at position {pos}")
        pos += 1  # Skip ')'
    else:
        # Leaf node
        node.is_leaf = True

    # Parse label (if present)
    label_match = re.match(r"(['\"]?)([^,:;'\"()]*)\1", s[pos:])
    if label_match:
        label = label_match.group(2)
        if label:
            node.label = label
        pos += label_match.end()

    # Parse branch length (if present)
    if pos < len(s) and s[pos] == ":":
        pos += 1  # Skip ':'
        branch_match = re.match(r"([\d.eE+-]+)", s[pos:])
        if branch_match:
            node.branch_length = float(branch_match.group(1))
            pos += branch_match.end()

    return node, pos


def parse_newick_file(file_path: str | Path) -> tuple[Tree, str]:
    """Load and parse Newick tree file.

    Args:
        file_path: Path to Newick file

    Returns:
        Tuple of (Tree object, file content hash)

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If Newick format is invalid
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Tree file not found: {file_path}")

    with open(file_path, "r") as f:
        content = f.read()

    # Compute hash
    file_hash = hashlib.sha256(content.encode()).hexdigest()

    # Parse (handle multiple trees - use first)
    trees = content.split(";")
    newick_string = trees[0].strip()

    if not newick_string:
        raise ValueError("Empty Newick content in file")

    tree = Tree.from_newick(newick_string)
    logger.info(f"Parsed tree from {file_path}: {tree.num_tips} tips, hash={file_hash[:8]}")

    return tree, file_hash


def validate_tree_structure(tree: Tree) -> list[str]:
    """Validate tree structure and return list of warnings/errors.

    Args:
        tree: Tree object to validate

    Returns:
        List of validation messages (empty if valid)
    """
    messages = []

    if tree.num_tips < 1:
        messages.append("Error: Tree has no tips")

    if tree.num_tips < 2:
        messages.append("Warning: Tree has fewer than 2 tips")

    tips = tree.get_all_tips()
    if len(tips) != len(set(tips)):
        duplicates = [t for t in set(tips) if tips.count(t) > 1]
        messages.append(f"Error: Duplicate tip labels: {duplicates}")

    # Check for empty labels
    all_tips = tree.root.get_all_leaves()
    if None in all_tips or "" in all_tips:
        messages.append("Warning: Some tips have empty labels")

    return messages
