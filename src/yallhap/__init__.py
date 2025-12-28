"""
yallHap: Modern Y-chromosome haplogroup inference.

A pipeline-friendly tool for Y-chromosome haplogroup classification
supporting modern and ancient DNA with probabilistic confidence scoring.
"""

__version__ = "0.2.0"

from yallhap.classifier import HaplogroupCall, classify
from yallhap.tree import Node, Tree

__all__ = ["Tree", "Node", "classify", "HaplogroupCall", "__version__"]
