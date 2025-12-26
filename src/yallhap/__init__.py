"""
yallHap: Modern Y-chromosome haplogroup inference.

A pipeline-friendly tool for Y-chromosome haplogroup classification
supporting modern and ancient DNA with probabilistic confidence scoring.
"""

__version__ = "0.1.0"

from yallhap.tree import Tree, Node
from yallhap.classifier import classify, HaplogroupCall

__all__ = ["Tree", "Node", "classify", "HaplogroupCall", "__version__"]
