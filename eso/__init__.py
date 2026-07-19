"""ESO - Evolutionary Stability Optimizer.

Detects hypermutable sites (recombination, replication-slippage, and
methylation-motif hotspots) in engineered DNA sequences and optimizes them
away using DNAChisel, while preserving the amino-acid translation.
"""

from eso.pipeline import main, suspect_site_extractor
from eso.optimize import optimization_engine

__all__ = ["main", "suspect_site_extractor", "optimization_engine"]
