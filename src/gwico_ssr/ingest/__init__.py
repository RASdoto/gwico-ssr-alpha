"""Input ingestion and metadata normalization for GWICO-SSR."""

from gwico_ssr.ingest.csv_loader import (
    IngestSummary,
    RowError,
    load_csv,
    register_local_files,
    validate_header,
)
from gwico_ssr.ingest.normalizers import (
    CompositeNormalizationResult,
    normalize_accession,
    normalize_completeness,
    normalize_date,
    normalize_geo_location,
    normalize_length,
    split_composite_fasta,
    split_composite_genbank,
)

__all__ = [
    "IngestSummary",
    "RowError",
    "load_csv",
    "register_local_files",
    "validate_header",
    "normalize_accession",
    "normalize_completeness",
    "normalize_date",
    "normalize_geo_location",
    "normalize_length",
    "CompositeNormalizationResult",
    "split_composite_fasta",
    "split_composite_genbank",
]
