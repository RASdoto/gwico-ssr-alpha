# Chunk 1 Batch-First Foundation and Configuration Uplift

This document records the Chunk 1 foundation changes for batch-first workflows.

## Scope

Chunk 1 adds only configuration, schema scaffolding, and orchestration extension points.

- No imperfect or compound detection behavior is introduced in this chunk.
- Existing accession-level alpha behavior remains valid and supported.
- The additions are foundation-only for later batch acquisition and normalization chunks.

## Configuration Additions

The configuration model now separates request sizing from artifact sizing and parse chunk sizing.

### NCBI Section

- `ncbi.batch_size` remains present for alpha compatibility.
- `ncbi.request_batch_size` is now explicit and defaults to `500`.

### Batch Section

New `[batch]` section:

- `artifact_batch_size` (default `10000`)
- `parse_chunk_size` (default `1000`)
- `retain_raw_artifacts` (default `true`)
- `manifest_policy` (default `required`)
- `duplicate_handling` (default `error`)

This encodes the required distinction between:

1. API request batch size
2. Persisted composite artifact batch size
3. Downstream parse chunk size

## Schema Scaffolding

Two new tables were added to represent batch provenance without altering existing alpha record contracts:

1. `batch_artifacts`

- Raw composite artifact metadata, file type/source, sizing metadata, checksum, status, and manifest path.

2. `batch_normalized_records`

- Mapping from source artifact and record index to accession-level normalized output paths and parse status.

These tables are scaffolding only in Chunk 1 and are not yet consumed by download/parse execution paths.

## Orchestration Extension Points

Chunk 1 adds foundation helpers in orchestration for later integration:

- resolved batch plan generation from settings,
- deterministic artifact base path resolution.

No default pipeline stage ordering or runtime stage execution behavior changed in this chunk.

## Compatibility

- Existing alpha accession-level workflows remain supported.
- Existing alpha CLI commands remain unchanged.
- Perfect-only alpha SSR behavior remains unchanged.
