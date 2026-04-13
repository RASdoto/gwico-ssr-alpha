# GWICO-SSR Validation and Benchmarking

## Gold-Standard Validation

### Fixtures

Two synthetic sequences with SSRs at known positions:

| Sequence | Length | SSR Hits | Motif Sizes Covered |
|----------|--------|----------|---------------------|
| GOLD_001 | 500 bp | 12 | 1, 2, 3, 4, 5, 6 |
| GOLD_002 | 300 bp | 4 | 1, 2, 3 |

Files:
- `tests/fixtures/gold_standard/gold_standard.fasta`
- `tests/fixtures/gold_standard/gold_standard_expected.json`
- `tests/fixtures/gold_standard/generate_fixtures.py` (reproducible generator)

### Correctness Verification

All 16 expected SSR hits are verified with exact match on:
- Start/end coordinates
- Raw motif and canonical motif
- Motif size and repeat unit count
- Strand assignment
- Actual repeat sequence

### Invariant Properties Verified

- `end - start == repeat_length_bp` (half-open coordinates)
- `repeat_length_bp == motif_size × repeat_units`
- `actual_repeat == sequence[start:end]`
- All coordinates within `[0, sequence_length)`
- Hits sorted by start position
- Motif sizes in range 1–6
- Strand is "+" or "-"
- Repeat units ≥ threshold for motif size

## Reproducibility

- **5× repeated detection runs** on both gold-standard sequences produce
  identical output (positions, motifs, counts).
- **Synthetic genome reproducibility**: same seed → same genome → same SSRs
  across 3 runs on 10K and 30K genomes.
- **Accession-independent**: same sequence with different accession IDs
  produces identical SSR hits.

## Performance Benchmarks

### SSR Detection

| Genome Size | Detection Time | SSRs Found |
|-------------|---------------|------------|
| 1 Kbp | < 10 ms | varies |
| 10 Kbp | < 50 ms | varies |
| 30 Kbp | < 500 ms | varies |
| 100 Kbp | < 5,000 ms | varies |

**Acceptance threshold:** 30 Kbp genome (SARS-CoV-2 scale) < 500 ms.

**Scaling:** Approximately linear with genome size.

### Annotation Mapping

Interval-tree and naive approaches produce identical results at all scales.
Python `intervaltree` has higher constant factors; the tree approach becomes
advantageous at very large feature counts (> 10K features).

## Known Deviations from PERF

| Aspect | GWICO-SSR | PERF |
|--------|-----------|------|
| Coordinates | 0-based half-open `[start, end)` | 0-based inclusive |
| Motif canonical | Lexicographic-min rotation (both strands) | Raw as found |
| Sub-repeats | Filtered (AAGAAG → AAG×2) | May report both |
| Overlap resolution | Shorter motif preferred | Not specified |

These deviations are intentional design choices for consistency and
correctness. They are documented in the benchmark test suite docstrings.

## Test Coverage

| Test File | Tests | Coverage Area |
|-----------|-------|---------------|
| test_validation.py | 28 | Gold-standard, reproducibility, invariants |
| test_benchmarks.py | 22 | Performance, scaling, report generation |
| test_ssr.py | 41 | Detector unit tests, motif utilities |
| test_annotation.py | 40 | Mapper, interval tree, naive comparison |
| Total benchmark+validation | 50 | |
| Total all tests | 490 | Full pipeline |
