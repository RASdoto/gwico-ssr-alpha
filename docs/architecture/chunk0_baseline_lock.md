# Chunk 0 Baseline Lock

This document records the explicit implementation baseline for the current alpha repository before any batch-first beta work begins.

## Scope

Chunk 0 is a baseline lock only.

- No runtime behavior changes are introduced here.
- No schema, configuration, CLI, downloader, parser, or orchestration semantics are changed.
- The purpose of this document is to tie later batch-first work to the actual repository state in this workspace.

## Repository Audit

### Package and CLI Baseline

- Package metadata is still alpha: `gwico-ssr` version `0.1.0a1`.
- The command surface is still the 12-command alpha CLI: `info`, `init-db`, `ingest`, `download`, `parse`, `detect`, `annotate`, `metrics`, `analyze`, `visualize`, `export`, `run`.
- No beta-only command exists yet for composite-file normalization, lineage, phylogeny, browser tracks, clustering, or temporal analysis.

### Configuration Baseline

- Configuration is still limited to the alpha sections: `[database]`, `[ncbi]`, `[ssr]`, `[logging]`, `[output]`.
- The NCBI configuration exposes a single `batch_size` field with default `500` requests per Entrez batch.
- The current config model does not distinguish API request batch size from persisted artifact batch size or downstream parse chunk size.
- No config exists yet for raw composite artifact retention, normalization manifests, duplicate handling policy, or parse chunking.

### Downloader Baseline

- The downloader writes deterministic accession-level outputs only: `sequences/fasta/{accession}.fasta` and `sequences/genbank/{accession}.gb`.
- Download orchestration iterates accession by accession and updates a single `SequenceRecord` per accession.
- Retry handling is accession-granular through a JSON retry manifest.
- The downloader does not persist composite FASTA or GenBank artifacts as first-class objects.
- The current alpha repository still assumes individual accession downloads as the primary workflow contract.

### Parser Baseline

- The FASTA parser can read multi-record files, but the persistence path still resolves one accession at a time and persists only the accession-matched record into the current accession-level contract.
- The GenBank parser reads all records from a file and then explicitly takes `records[0]`, which means composite GenBank is not supported as a first-class workflow input.
- The parse orchestration expects files to exist at accession-derived paths under `sequences/fasta` and `sequences/genbank`.
- The standard parse workflow therefore remains one-file-per-accession even though raw FASTA parsing is multi-record capable.

### ORM Schema Baseline

- The alpha ORM contains 11 tables: `datasets`, `runs`, `accessions`, `sequence_records`, `feature_records`, `ssr_records`, `ssr_annotations`, `accession_metrics`, `statistical_results`, `stage_checkpoints`, and `failed_accessions`.
- Existing provenance fields are limited to dataset manifest hash, run config hash, run metadata, per-accession sequence file paths, and retry/checkpoint tables.
- There is no explicit schema concept yet for raw composite artifacts, normalized split outputs, batch membership, artifact checksum, record index inside a composite file, or composite-to-accession provenance.

### Orchestration Baseline

- Pipeline ordering exists for `INGEST`, `DOWNLOAD`, `PARSE`, `DETECT`, `ANNOTATE`, `METRICS`, `ANALYZE`, `VISUALIZE`, and `EXPORT`.
- The default orchestrator only executes the alpha core stages `DETECT`, `ANNOTATE`, `METRICS`, and `ANALYZE`.
- Failed-accession recovery and stage checkpointing already exist, but they are accession-oriented rather than composite-artifact-aware.
- The current orchestrator does not model normalization, composite-artifact recovery, artifact manifests, or batch-derived parse boundaries.

### Analysis and Visualization Baseline

- Statistical analysis is alpha-scoped: chi-square, Kruskal-Wallis, correlation, Shannon entropy, and FDR correction on persisted accession metrics.
- Visualization is built on dataset- and accession-level persisted results and does not yet include batch provenance, repeat classes beyond perfect SSRs, or batch-derived workflow reporting.
- Export and figure generation are downstream of the accession-level parse and detection contract.

### Test Baseline

- The repository includes unit, regression, integration, validation, and benchmark test directories.
- Chunk 0 adds a CLI regression test that locks the 12-command alpha surface and rejects known beta-only commands.
- Full validation for this chunk is the existing alpha test suite plus the added CLI surface regression.

## Stable-State Report

The following commands are the stable-state validation commands for Chunk 0 and are recorded exactly as run from the workspace venv:

```powershell
& "f:\Gama build\.venv\Scripts\python.exe" -m gwico_ssr --version
& "f:\Gama build\.venv\Scripts\python.exe" -m gwico_ssr --help
& "f:\Gama build\.venv\Scripts\python.exe" -m pytest tests\unit\test_cli.py -q
& "f:\Gama build\.venv\Scripts\python.exe" -m pytest -q
```

Observed results from the actual Chunk 0 validation run:

- `& "f:\Gama build\.venv\Scripts\python.exe" -m gwico_ssr --version` reported `gwico-ssr, version 0.1.0a1`.
- `& "f:\Gama build\.venv\Scripts\python.exe" -m gwico_ssr --help` exposed the 12-command alpha CLI surface and no beta-only commands.
- `& "f:\Gama build\.venv\Scripts\python.exe" -m pytest tests\unit\test_cli.py -q` passed.
- `& "f:\Gama build\.venv\Scripts\python.exe" -m pytest -q` completed with `491 passed, 52 warnings in 38.60s` and `EXIT:0`.

The pre-Chunk-0 alpha baseline remained the verified `490 passed, 53 warnings` state from the prompt pack. The current `491 passed` result reflects the additional CLI baseline regression test added in this chunk.

## Explicit Gap List: Alpha State vs. Batch-First Beta Target

### Batch Acquisition Gaps

- The current downloader does not separate API request batch size from persisted composite artifact batch size.
- The current system does not persist large composite FASTA or GenBank artifacts.
- The current retry model is accession-level only and has no artifact-level provenance.

### Batch Parsing Gaps

- Composite FASTA is not a first-class operator workflow even though the low-level parser can read multi-record files.
- Composite GenBank is not supported because only the first record is used.
- No normalization boundary exists yet for composite artifacts to deterministic downstream accession outputs.
- No parse chunk size or stream-oriented artifact processing contract is defined yet.

### Provenance Gaps

- No manifest exists linking accession manifests to raw batch artifacts and then to normalized per-accession outputs.
- No schema exists for batch identity, composite checksum, accession membership in a batch artifact, or normalized-output provenance.

### Workflow Gaps

- The large-cohort operator story still assumes one file per accession as the practical workflow.
- `run` and `parse` do not yet expose a batch-first large-cohort path.
- The current docs do not yet define batch terminology or the required distinction between request batch size, artifact batch size, and parse chunk size.

### Detector and Downstream Gaps

- Imperfect SSR detection is not implemented.
- Compound SSR detection is not implemented.
- Downstream metrics, analysis, visualization, and export are not yet batch-provenance-aware or repeat-class-aware.

## Alpha Compatibility Rules

The following rules are locked as the compatibility baseline for subsequent chunks:

1. Perfect-only alpha SSR behavior must remain reproducible.
2. The existing 12-command alpha CLI surface must remain stable unless a later chunk explicitly documents an approved public change.
3. Accession-level workflows remain supported even after batch-first infrastructure is added.
4. Existing alpha tests must continue to pass unless an intentional interface change is explicitly documented.
5. IMEx is a semantic reference model only and must not be invoked as the runtime engine.
6. Batch-first infrastructure must precede imperfect or compound detector expansion for the intended large-cohort beta workflow.
7. Request batch size, persisted artifact batch size, and downstream parse chunk size must remain separate concepts in later chunks.

## Implementation Consequence for Later Chunks

Batch-first infrastructure must be implemented before detector expansion.

The reason is architectural rather than optional: the intended beta workflow requires large-cohort download, normalization, parse, provenance, and recovery behavior to exist before imperfect and compound detection can become the default large-cohort path.

## Chunk 0 Exit Condition

Chunk 0 is complete only when:

- the audit is tied to the actual repository files,
- the alpha compatibility rules are documented,
- the batch-first target gaps are explicit,
- and the alpha validation commands complete successfully.