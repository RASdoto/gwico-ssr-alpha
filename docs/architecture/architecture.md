# GWICO-SSR Architecture

## Overview

GWICO-SSR is a CLI-first, database-backed SSR analysis platform organized as a
single Python package with 12 functional modules.

## Module Dependency Graph

```
cli/
 ├── config/         (settings loading)
 ├── logging/        (structured logging)
 ├── db/             (engine, sessions, repository)
 │    └── models/    (ORM schema)
 ├── ingest/         (CSV + file ingestion)
 ├── parsers/        (FASTA, GenBank, GFF3)
 ├── ssr/            (detection engine)
 ├── annotation/     (interval-tree mapping)
 ├── metrics/        (per-accession + cohort)
 ├── analysis/       (statistics + FDR)
 ├── visualization/  (figures)
 ├── export/         (CSV, JSON, BED, GFF3)
 └── orchestration/  (pipeline + checkpointing)
```

## Database Schema (11 Tables)

| Table | Purpose | Key Columns |
|-------|---------|-------------|
| `datasets` | Input dataset metadata | name, organism, source_type |
| `runs` | Pipeline execution records | dataset_id, status, stage, config_hash |
| `accessions` | Genome accession metadata | accession, country, species, gc_content |
| `sequence_records` | Parsed sequence metadata | accession, sequence_hash, fasta_path |
| `feature_records` | Gene annotations | accession, gene_name, start, end, strand |
| `ssr_records` | Detected SSR loci | accession, motif_canonical, start, end, repeat_units |
| `ssr_annotations` | SSR-to-feature mapping | ssr_id, feature_id, overlap_bp |
| `accession_metrics` | Per-genome summary metrics | accession, ra, rd, ssr_count_total |
| `statistical_results` | Analysis outputs | test_name, p_value, p_value_corrected, effect_size |
| `stage_checkpoints` | Pipeline stage tracking | run_id, stage, status |
| `failed_accessions` | Error tracking | accession, stage, error_message, retry_count |

## Configuration System

```
CLI flags  →  Environment variables  →  TOML file  →  Built-in defaults
(highest)                                                (lowest)
```

- Pydantic Settings models with `GWICO_SSR_` env prefix
- Nested sections: `[database]`, `[ncbi]`, `[ssr]`, `[logging]`, `[output]`

## SSR Detection Algorithm

For each motif size 1–6:
1. Slide a window of `motif_size` along the sequence
2. Count consecutive identical motif occurrences
3. If count ≥ threshold, emit SSR record
4. Filter sub-repeats (e.g., AAGAAG as hexa when AAG is the fundamental unit)
5. Resolve overlaps: prefer shorter motif size

Motif canonicalization: lexicographic minimum rotation across forward and
reverse complement strands.

## Annotation Mapping

1. Build interval tree from gene features: O(n log n)
2. Query each SSR position: O(log n + k) per query
3. Classify non-overlapping SSRs as intergenic
4. Support one-to-many (SSR spanning multiple features)

## Pipeline Orchestration

```
DETECT → ANNOTATE → METRICS → ANALYZE
```

Each stage has a checkpoint. Resume skips completed stages.
Failed accessions are tracked per-stage with retry support.
