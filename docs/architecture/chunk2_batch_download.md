# Chunk 2 Configurable Large-Batch Download Support

Chunk 2 integrates batch-aware acquisition into runtime download behavior while preserving accession-level compatibility.

## What Changed

1. Request batching for remote retrieval

- Downloader now supports request-batched retrieval using `ncbi.request_batch_size`.
- FASTA uses batch retrieval calls when request chunks contain multiple accessions.
- GenBank now also supports batch retrieval with deterministic record splitting and fallback to per-accession retrieval.

2. Persisted composite artifact batching

- Downloader groups pending accessions by `batch.artifact_batch_size`.
- For each artifact chunk and file type, deterministic composite artifacts are written to:
  - `outputs/batches/raw/fasta/batch_000001.fasta`
  - `outputs/batches/raw/genbank/batch_000001.gb`
- Artifact naming is stable and sequential within each file type.

3. Accession-level recovery remains explicit

- Retry manifest behavior remains accession-granular.
- Batch failures are still recorded per accession in summary and retry manifest inputs.

4. Batch provenance scaffolding is populated

- Download now persists batch provenance rows to chunk1 scaffold tables:
  - `batch_artifacts`
  - `batch_normalized_records`

## Operator Sizing Guidance

- Recommended API request size: `500`
- Recommended persisted artifact size for million-accession scale: `10000`
- Parse chunk sizing remains separate and stream-oriented (`batch.parse_chunk_size`)

These are separate controls and are intentionally not conflated.

## Compatibility

- Existing accession-level output files remain generated as before.
- Existing alpha download command remains valid.
- Small runs still work in accession-compatible mode.
