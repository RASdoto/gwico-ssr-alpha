# GWICO-SSR Alpha Release Notes

**Version:** 0.1.0-alpha1  
**Date:** April 2026

## Summary

GWICO-SSR alpha is the first internal release of a unified, database-backed
SSR analysis platform. It replaces 15+ standalone scripts with a single CLI
application that can process complete genome datasets from raw input to
publication-ready outputs.

## What's Included

### Core Pipeline (9 stages)
- **Ingest** — CSV metadata loading with normalization and deduplication
- **Download** — NCBI Entrez acquisition with retry, rate-limiting, and manifests
- **Parse** — FASTA, GenBank, and GFF3 parsing with structured error handling
- **Detect** — Native perfect SSR detection (motif sizes 1–6), replacing PERF
- **Annotate** — Interval-tree gene mapping (O(n log n + k) vs former O(n*m))
- **Metrics** — Per-accession RA, RD, motif counts, GC content, dominant motif
- **Analyze** — Chi-square, Kruskal-Wallis, correlation, Shannon entropy with BH-FDR
- **Visualize** — 9 figure types: bar charts, heatmaps, boxplots, choropleth, sunburst
- **Export** — CSV, JSON, BED, GFF3, publication tables, run manifest with provenance

### Infrastructure
- 12 CLI commands with `--help` documentation
- SQLite database with WAL mode (PostgreSQL optional)
- TOML + environment variable configuration system
- Structured JSON/text logging with run-scoped context
- Pipeline orchestration with checkpointing and resume
- Failed accession tracking and retry support

### Quality
- 490 automated tests (unit, integration, benchmark, validation)
- Gold-standard validation fixtures with deterministic SSR positions
- Benchmark harness for detection and annotation performance
- Reproducibility verification (repeated runs produce identical output)

## Performance Characteristics

| Operation | Scale | Time |
|-----------|-------|------|
| SSR detection | 30 Kbp genome (SARS-CoV-2) | < 500 ms |
| SSR detection | 100 Kbp genome | < 5 s |
| Detection scaling | Linear with genome size | Verified |
| Annotation mapping | 1K SSRs × 100 features | < 100 ms |

## Known Deviations from PERF

| Aspect | GWICO-SSR | PERF |
|--------|-----------|------|
| Coordinates | 0-based half-open `[start, end)` | 0-based inclusive |
| Motif canonical | Lexicographic-minimum rotation (both strands) | Raw motif as found |
| Sub-repeats | Filtered (AAGAAG→AAG×2) | May report both |
| Overlap resolution | Shorter motif preferred | Not specified |

## What's Deferred

See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for the full list.

Key items not included in alpha:
- Web dashboard / REST API
- Imperfect and compound SSR detection
- Lineage-aware analysis (Pango/Nextstrain)
- Phylogenetic integration
- ML-based clustering
- Cloud-scale distributed execution
- Temporal SSR tracking
- Genome browser integration
