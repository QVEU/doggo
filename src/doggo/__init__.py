"""doggo: resampling-based gene set enrichment for arbitrary per-gene scores."""

from doggo.genesets import filter_gene_sets, read_gmt
from doggo.enrichment import EnrichmentConfig, run_enrichment, resample_test

__version__ = "0.4.0"

__all__ = [
    "EnrichmentConfig",
    "filter_gene_sets",
    "read_gmt",
    "run_enrichment",
    "resample_test",
    "__version__",
]
