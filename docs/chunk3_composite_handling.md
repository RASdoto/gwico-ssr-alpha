# Chunk 3: Composite FASTA and GenBank Handling

**Status:** ✅ Complete  
**Date:** 2026-04-23  
**Version:** 0.1.0-alpha1 + Chunk 3 changes

## Overview

Chunk 3 implements first-class support for multi-record (composite) FASTA and GenBank files. This enables batch-first workflows where large genomic datasets can be ingested as composite files without requiring manual pre-splitting.

## What Changed

### 1. New Composite File Utilities (`src/gwico_ssr/ingest/normalizers.py`)

Added high-level functions for splitting composite files into normalized single-record outputs:

#### `split_composite_fasta()`
Splits multi-record FASTA files into individual normalized files.

```python
from gwico_ssr.ingest.normalizers import split_composite_fasta

result = split_composite_fasta(
    composite_path="sequences/batch_001.fasta",
    output_dir="sequences/normalized",
    overwrite=False,
    duplicate_policy="first",  # or "error", "skip"
)

print(f"Records: {result.record_count}")
print(f"Unique accessions: {result.accession_count}")
for rec in result.normalized_records:
    print(f"  {rec['accession']} -> {rec['output_path']}")
```

**Features:**
- Handles multi-record FASTA files transparently
- Extracts accession IDs from FASTA headers (first whitespace-delimited token)
- Computes SHA-256 checksums for integrity verification
- Configurable duplicate handling:
  - `"first"` (default): Keep first occurrence, silently skip duplicates
  - `"error"`: Fail on duplicates with error message
  - `"skip"`: Skip all duplicates, ignore extras
- Creates output directory automatically
- Returns `CompositeNormalizationResult` with detailed provenance

#### `split_composite_genbank()`
Splits multi-record GenBank files into individual normalized files.

```python
from gwico_ssr.ingest.normalizers import split_composite_genbank

result = split_composite_genbank(
    composite_path="annotations/batch_001.gb",
    output_dir="annotations/normalized",
    duplicate_policy="first",
)
```

**Features:**
- Same interface as `split_composite_fasta()` but for GenBank files
- Preserves all GenBank annotations and features
- Extracts accession from LOCUS line or ID field

### 2. Enhanced GenBank Parser (`src/gwico_ssr/parsers/genbank_parser.py`)

Added `parse_genbank_composite()` function for parsing multi-record GenBank files:

```python
from gwico_ssr.parsers.genbank_parser import parse_genbank_composite

results = parse_genbank_composite("batch_annotations.gb")
for result in results:
    if result.success:
        print(f"Parsed: {result.sequence_info.accession}")
        print(f"  Features: {result.feature_count}")
        print(f"  CDS: {result.cds_count}")
    else:
        print(f"Error: {result.errors}")
```

**Key Changes:**
- Original `parse_genbank()` now emits a warning when multi-record files are detected
- Backward compatibility preserved: still processes first record only
- New `parse_genbank_composite()` returns `list[GenBankParseResult]`, one per record

### 3. FASTA Parser Already Supported Multi-Record

The existing `parse_fasta()` function already supported multi-record FASTA files (via `SeqIO.parse()` loop). Chunk 3 validates this behavior and confirms it works correctly.

## API Reference

### `CompositeNormalizationResult`

Dataclass returned by `split_composite_*()` functions:

```python
@dataclass
class CompositeNormalizationResult:
    source_path: str                          # Path to input file
    file_type: str                            # "fasta" or "genbank"
    checksum_sha256: str                      # SHA-256 hash of input
    record_count: int                         # Total records in file
    accession_count: int                      # Unique accessions
    normalized_records: list[dict]            # [{accession, output_path}, ...]
    errors: list[dict]                        # [{accession, message}, ...]
    duplicates: list[dict]                    # [{accession, record_index, policy}, ...]
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
```

## Workflow Examples

### Example 1: Split and Reparse Composite FASTA

```python
from pathlib import Path
from gwico_ssr.ingest.normalizers import split_composite_fasta
from gwico_ssr.parsers.fasta_parser import parse_fasta

# Split composite file
split_result = split_composite_fasta(
    "large_batch.fasta",
    "output/normalized_fasta"
)

print(f"Split {split_result.record_count} records into "
      f"{split_result.accession_count} accessions")

# Process each normalized file
for rec in split_result.normalized_records:
    accession = rec["accession"]
    normalized_path = rec["output_path"]
    
    # Parse individual file
    parse_result = parse_fasta(normalized_path)
    if parse_result.success:
        seq = parse_result.records[0]
        print(f"{accession}: {seq.sequence_length} bp, "
              f"GC={seq.gc_content:.2%}")
```

### Example 2: Batch GenBank Handling with Feature Extraction

```python
from gwico_ssr.ingest.normalizers import split_composite_genbank
from gwico_ssr.parsers.genbank_parser import parse_genbank

# Split composite GenBank
split_result = split_composite_genbank(
    "ncbi_batch_001.gb",
    "output/normalized_genbank"
)

# Persist mappings (for downstream annotation)
for rec in split_result.normalized_records:
    accession = rec["accession"]
    normalized_path = rec["output_path"]
    
    # Parse and extract features
    parse_result = parse_genbank(normalized_path)
    if parse_result.success:
        print(f"{accession}:")
        print(f"  Length: {parse_result.sequence_info.sequence_length} bp")
        print(f"  CDS features: {parse_result.cds_count}")
        for feat in parse_result.features[:3]:  # First 3 features
            print(f"    {feat.gene_name}: {feat.start}-{feat.end}")
```

### Example 3: Duplicate Handling

```python
# Option A: Keep first occurrence (default)
result = split_composite_fasta(
    "data.fasta",
    "output",
    duplicate_policy="first"
)
# Silently skips additional occurrences; logs in duplicates list

# Option B: Fail on duplicates
result = split_composite_fasta(
    "data.fasta",
    "output",
    duplicate_policy="error"
)
if result.duplicates:
    print(f"ERROR: Found {len(result.duplicates)} duplicate accessions")
    for dup in result.duplicates:
        print(f"  {dup['accession']} (record {dup['record_index']})")

# Option C: Skip all duplicates entirely
result = split_composite_fasta(
    "data.fasta",
    "output",
    duplicate_policy="skip"
)
print(f"Processed {result.accession_count} unique accessions from "
      f"{result.record_count} records")
```

## Testing

All composite file functionality is covered by tests:

```bash
# Run composite handling tests only
pytest tests/unit/test_composite_handling.py -v

# Run full test suite
pytest tests/
```

**Test Coverage:**
- ✅ Composite FASTA splitting (basic, accession extraction, output validation)
- ✅ Composite GenBank splitting (basic, accession extraction, output validation)
- ✅ Duplicate handling (all three policies: first, error, skip)
- ✅ Error handling (missing files, empty files)
- ✅ Multi-record GenBank parsing
- ✅ Round-trip consistency (split → reparse → verify)

**Test Fixtures:**
- `tests/fixtures/multi_record.fasta` - 2-record FASTA fixture
- `tests/fixtures/multi_record.gb` - 2-record GenBank fixture

## Backward Compatibility

✅ **Fully backward compatible.**

- Existing single-record workflows unchanged
- FASTA parser still works exactly as before
- GenBank parser warns but processes first record if multi-record file provided
- No schema changes required
- All existing tests pass

## Integration with Batch-First Pipeline (Chunks 0-2)

Chunk 3 enables the normalization boundary that Chunks 0-2 designed for:

1. **Chunk 0** established baseline and gap analysis
2. **Chunk 1** added batch config and schema (`request_batch_size`, `artifact_batch_size`)
3. **Chunk 2** implemented configurable batch download
4. **Chunk 3** (this) adds composite file normalization ← **You are here**

Together, these enable:
```
Raw composite files (10,000 accessions) 
  ↓ split_composite_fasta()
Normalized single-record files (10,000 × {accession}.fasta)
  ↓ existing parse, detect, annotate pipeline (unchanged)
Processed SSR data, metrics, exports
```

## Known Limitations

1. **Large files**: For very large composite files (>1 GB), memory usage grows linearly during split. Use streaming approach in future versions if needed.
2. **Malformed records**: If a record is malformed within a composite file, that record fails but processing continues for others.
3. **Version parsing**: GenBank version numbers (e.g., "NC_045512.2") are preserved in accession field. Normalize separately if needed for exact matching.

## Next Steps (Chunk 4)

Chunk 4 will integrate composite normalization into the batch-aware parse and orchestration layer, making composite files a first-class default input for large cohorts.

## Files Modified

- `src/gwico_ssr/ingest/normalizers.py` - Added composite splitting functions
- `src/gwico_ssr/parsers/genbank_parser.py` - Added multi-record parsing
- `tests/unit/test_composite_handling.py` - New comprehensive test suite (18 tests)
- `tests/fixtures/multi_record.gb` - New GenBank composite fixture
- `tests/fixtures/multi_record.fasta` - Already existed, validated

## Questions?

For usage questions or issues, refer to:
- [Architecture guide](architecture/architecture.md)
- [Validation guide](validation/validation.md)
- Test examples: `tests/unit/test_composite_handling.py`
