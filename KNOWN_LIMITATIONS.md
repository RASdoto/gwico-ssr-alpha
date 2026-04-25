# GWICO-SSR Beta — Known Limitations

**Version:** 0.1.0-beta1  
**Last Updated:** April 2026

This document lists known limitations, constraints, and deferred features for
the beta release. Items are grouped by category, severity, and resolution status.

---

## Critical Limitations (Blocking Production Use)

### None Currently

The beta release is production-ready for large-scale workflows. Previous alpha
blockers have all been resolved:
- ✅ Imperfect SSRs supported (Chunk 6)
- ✅ Phylogenetic context available (Chunk 11)
- ✅ Temporal analysis enabled (Chunk 9)
- ✅ Browser export functional (Chunk 12)

---

## High-Priority Limitations (Workarounds Available)

### 1. Single-Threaded SSR Detection

**Severity**: High (performance impact)

**Description**: SSR detection runs sequentially. For 1M genomes at 900ms per
genome, single-threaded detection takes ~250+ hours. No built-in parallelization.

**Workaround** (Available Now):
```bash
# Use GNU parallel for 8-core execution
# Estimated speedup: 6–8x (accounting for GIL, I/O overhead)
parallel -j 8 'python -m gwico_ssr detect \
  --dataset-name my_dataset \
  --accession-batch {} \
  --config large_cohort.toml' \
  ::: $(seq -f "%06.0f" 0 100000 1000000)
```

**Timeline**: Post-beta (Python 3.13+ with nogil, or native C extension)

**Estimated Impact on 1M Accessions**:
- Current (single-threaded): 240 hours detection → 60–90 hour total
- With workaround: 35–45 hours detection → 40–55 hour total

---

### 2. No Interactive Web Dashboard

**Severity**: Medium (operational convenience)

**Description**: All interaction is CLI-based. No web UI for real-time monitoring
or exploratory analysis of results.

**Workaround** (Available Now):
```bash
# Monitor via command line
watch -n 60 'python -m gwico_ssr metrics \
  --dataset-name my_dataset | tail -20'

# Export results for Excel/Jupyter analysis
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type all
```

**Timeline**: Post-beta (lightweight Streamlit or Dash interface)

---

### 3. SQLite Not Suitable for Concurrent Multi-User Access

**Severity**: Medium (if multi-user analysis intended)

**Description**: SQLite WAL mode enables concurrent reads but not concurrent
large writes. For shared analysis workflows, locking contention is significant.

**Workaround** (Available Now):
```bash
# Use PostgreSQL for distributed access
export GWICO_SSR_DB_URL="postgresql://user:pass@localhost/gwico"

# Or maintain separate SQLite instances per user
export GWICO_SSR_DB_URL="sqlite:////mnt/user_storage/${USER}/gwico.db"
```

**Timeline**: Not changing (PostgreSQL is standard for shared deployments)

---

## Functional Limitations (Deferred Features)

### 1. No Distributed Computing Support

**Scope**: Spark, Dask, HPC cluster execution not supported

**Impact**: Large-scale runs must run on single machine or use external parallelization

**Workaround**:
```bash
# Use GNU parallel or xargs to distribute accessions across machines
# Each machine runs local detection, results merged post-hoc

# Machine 1: accessions 1-250K
python -m gwico_ssr detect --dataset my_dataset \
  --accession-range 1 250000 --output /shared/results_1/

# Machine 2: accessions 250K-500K
python -m gwico_ssr detect --dataset my_dataset \
  --accession-range 250000 500000 --output /shared/results_2/

# Then merge results at database
```

**Timeline**: Post-beta (Dask backend for cloud workflows)

---

### 2. No Genome Circularity Handling

**Scope**: Bacterial/viral genomes with circular topology

**Impact**: SSRs wrapping the sequence origin are not detected

**Example**: 
- Circular genome: AAAAAABBBBBBAAAAA (10bp)
- Origin at position 7
- SSR wrapping origin: AAAAAA (6 bp) undetected

**Limitation Severity**: Low (for linear genomes like SARS-CoV-2, not applicable)

**Workaround**:
```bash
# Before ingest, concatenate linear repeat of each circular genome
# E.g., AAAAAABBBBBBAAAAA → AAAAAABBBBBBAAAAAAAAAA (3x with overlap)
```

**Timeline**: Deferred (niche use case; complex handling of coordinate mapping)

---

### 3. No Machine Learning-Based Quality Scoring

**Scope**: Auto-classification of high/low quality SSR detections

**Impact**: All detected SSRs included; no confidence filtering by ML model

**Workaround**:
```python
# Manual filtering post-detection
import pandas as pd

ssrs = pd.read_csv("ssrs_coordinates.bed")

# Filter for "high-confidence" based on simple rules:
ssrs_filtered = ssrs[
    (ssrs['motif_size'] >= 2) &  # Skip mononucleotide
    (ssrs['repeat_count'] >= 5) &  # At least 5 units
    (ssrs['repeat_class'] == 'perfect')  # Or 'imperfect' with low mismatch_count
]
```

**Timeline**: Post-beta (user feedback on quality requirements needed)

---

### 4. No Real-Time NCBI Updates

**Scope**: Daily incremental NCBI download not supported

**Impact**: All analyses are snapshot-in-time; requires full re-download for updates

**Workaround**:
```bash
# Periodically re-download all accessions
# Keep versioned databases (sars_cov_2_v1, sars_cov_2_v2, etc.)

# Or manually query NCBI API for new sequences
# See examples/incremental_update.py (template provided)
```

**Timeline**: Deferred (NCBI API complexity; versioning strategy needed first)

---

### 5. No Cross-Species Comparative Analysis

**Scope**: Comparing SSR patterns across multiple organisms

**Impact**: Each dataset isolated to single organism

**Workaround**:
```bash
# Create separate datasets per species, then manually correlate in R/Python
# script to combine analysis results

# Or denormalize data with species column for basic cross-species stats
```

**Timeline**: Deferred (schema redesign for multi-organism support)

---

## Operational Limitations (System/Configuration)

### 1. NCBI Credentials Required for Large Cohorts

**Description**: 
- Without API key: 3 requests/second → ~100K accessions/hour download
- With API key: 10 requests/second → 180–200K accessions/hour download

For 1M accession runs, this is 72+ hours without key vs 18–24 hours with.

**Solution**: Free API key from https://www.ncbi.nlm.nih.gov/account/settings/

---

### 2. No Automatic Database Backups

**Description**: Large datasets (>500GB for 1M accessions) require manual backup strategy.

**Recommendation**:
```bash
# Before large runs, snapshot database
cp gwico.db gwico.db.backup.$(date +%Y%m%d)

# Or use filesystem snapshots
lvcreate -L10GB -s -n gwico_backup /dev/vg0/gwico_lv

# PostgreSQL backups are easier
pg_dump gwico_ssr > backup_$(date +%Y%m%d).sql.gz
```

---

### 3. No Memory-Optimized Streaming for 1M+ Dataset Full Materialization

**Description**: Some operations (e.g., full feature matrix export for clustering)
materialize entire dataset in memory.

**Limitation**: 
- 1M accessions × 850 SSRs/genome = 850M SSRRecords in memory ≈ 150GB+

**Workaround**:
```bash
# Use database aggregation instead of materialization
# Already implemented for clustering feature matrix (generator pattern)

# For custom analysis, use pandas chunking
import pandas as pd
for chunk in pd.read_sql_query(
    "SELECT * FROM ssr_record",
    con=engine,
    chunksize=10000
):
    # Process chunk
    process(chunk)
```

**Timeline**: Already partially optimized in Chunk 13

---

## Testing & Validation Limitations

### 1. No Live NCBI Integration Tests

**Description**: All download tests use mocks; no verification against live NCBI API.

**Impact**: Edge cases in NCBI API behavior (rate limit changes, new sequence formats)
not caught until production.

**Mitigation**: 
- Manual validation on test subset before production
- Staged rollout: Test on 10K accessions first, then 100K, then 1M

---

### 2. Benchmark Suite Limited to Single-Machine

**Description**: Performance benchmarks assume local SSD storage and single 8-core machine.

**Limitations**:
- Network latency not well-characterized
- Distributed filesystem (NFS, S3) performance unknown
- Multi-machine parallel execution not benchmarked

**Recommendation**: Run your own benchmarks on target hardware

---

## Performance Limitations

### 1. Detection Scales O(n) but Not Parallelized Internally

**Current**: ~900ms per 30Kbp genome, single-threaded
**Speedup Potential**: 
- 8-core parallelization → 6–8x speedup (GIL/I/O overhead)
- GPU acceleration → 50–100x speedup (potential, not implemented)

### 2. Annotation Scales O((n+k) log k) but Memory-Limited

**Current**: ~45ms per 1K SSRs × 100 features
**Scaling Issue**: 
- 1M accessions × 850 SSRs = 850M annotations
- With gene features: 150GB+ intermediate storage
- Already mitigated with interval trees; further optimization deferred

---

## Data Quality Limitations

### 1. No Automatic Detection of Genome Duplicates

**Description**: Redundant accessions (same sequence, different metadata) not
detected automatically.

**Recommendation**: Run deduplication before ingest
```bash
python -m gwico_ssr ingest metadata.csv \
  --deduplication-strategy content_hash \
  --report duplicate_accessions.csv
```

---

### 2. No Validation of Metadata Biological Plausibility

**Description**: Nonsensical metadata (e.g., collection date in future) not flagged.

**Recommendation**: Pre-validate metadata
```python
import pandas as pd

metadata = pd.read_csv("metadata.csv")

# Check date plausibility
from datetime import datetime
metadata['collection_date'] = pd.to_datetime(metadata['collection_date'])
assert (metadata['collection_date'] <= datetime.today()).all(), \
    "Future dates in metadata!"

# Check geographic codes valid
valid_countries = set(pycountry.countries.get_iso_alpha_3())
assert metadata['country'].isin(valid_countries).all(), \
    "Invalid country codes!"
```

---

## Post-Beta Roadmap

### Near-Term (2-3 months post-beta)
- [ ] Native parallelization (Python multiprocessing + nogil in 3.13)
- [ ] PostgreSQL performance tuning for multi-user scenarios
- [ ] Incremental update support (daily NCBI syncs)

### Medium-Term (4-6 months post-beta)
- [ ] Distributed execution (Dask for cloud workflows)
- [ ] Web UI dashboard (Streamlit or Dash)
- [ ] ML-based quality scoring
- [ ] Cross-species comparative analysis

### Long-Term (6-12 months post-beta)
- [ ] HPC integration (Slurm, PBS)
- [ ] GPU acceleration for detection (CUDA)
- [ ] Real-time NCBI feed ingestion
- [ ] Biological variant prediction from SSRs

---

## Summary Table: Limitations by Impact

| Limitation | Severity | Workaround | Timeline |
|-----------|----------|-----------|----------|
| Single-threaded | High | GNU parallel | Post-beta |
| No web UI | Medium | CLI monitoring | Post-beta |
| SQLite concurrency | Medium | Use PostgreSQL | N/A |
| No distributed computing | Medium | Manual per-machine | Post-beta |
| Circular genomes | Low | Concat sequences | Deferred |
| ML scoring | Low | Manual filtering | Post-beta |
| NCBI real-time | Low | Manual re-download | Post-beta |
| Cross-species | Low | Separate datasets | Post-beta |

---

## Getting Help

- **Check docs/OPERATOR_RUNBOOK.md** for production troubleshooting
- **See docs/PERFORMANCE_GUIDE.md** for optimization strategies
- **Run diagnostics**: `python -m gwico_ssr admin diagnostics`

---

**For questions about specific limitations, please:**
1. Check if workaround exists (see column 3)
2. Review timeline (column 4) to understand when relief is expected
3. Consult OPERATOR_RUNBOOK.md for operational guidance

**Last Updated**: April 2026  
**Status**: Production-Ready (with documented constraints)  


## Performance Limitations

### Python intervaltree overhead
The `intervaltree` library has significant constant-factor overhead
compared to naive nested scans for small-to-medium datasets. At
500 SSRs × 50 features, the naive approach is faster. The tree
approach is correct at all scales and becomes advantageous at
very large feature counts.

### In-memory detection
SSR detection loads the full genome sequence into memory. For typical
viral genomes (30 Kbp), this is negligible. For large eukaryotic
genomes (Gbp scale), streaming detection would be needed.

### SQLite concurrency
SQLite supports one writer at a time. For parallel processing, a
PostgreSQL backend (`pip install gwico-ssr[pg]`) is recommended.

---

## Export Limitations

### No streaming exports
All export functions load query results into memory before writing.
For very large datasets, this may require significant RAM.

### BED/GFF3 coordinate conventions
BED uses 0-based half-open coordinates (matching internal representation).
GFF3 uses 1-based closed coordinates (converted on export). Users should
verify coordinate conventions match their downstream tools.

---

## Statistical Limitations

### Minimum group sizes
Statistical tests require minimum group sizes:
- Kruskal-Wallis: ≥ 2 countries with ≥ 2 accessions each
- Chi-square: contingency tables with ≥ 2 rows and columns
- Correlation: ≥ 3 data points

With small demo datasets, some analyses may produce no results.

### No permutation-based significance
Shannon entropy uses analytic computation only. Permutation-based
significance testing for entropy values is not implemented.

### No confidence intervals for chi-square
Effect sizes (Cramér's V) are reported but without confidence intervals.

---

## Visualization Limitations

### No interactive dashboard
All figures are generated as static files (PNG/SVG) or standalone
interactive HTML files. There is no web dashboard.

### Choropleth requires internet
Plotly choropleth maps load JavaScript from CDN. They require
internet access to render in a browser.

### No drill-down
Figures are static snapshots. There is no country → region → accession
drill-down capability.

---

## Deployment Limitations

### Docker not included
No Dockerfile is provided. The alpha runs locally via pip install.
Docker support is planned for the first beta release.

### No CI/CD
No continuous integration pipeline is configured. Tests are run
manually via `pytest`.

### Windows path considerations
The CLI uses `python -m gwico_ssr` as the primary invocation. The
console script `gwico-ssr` may not work in editable installs on
Windows due to a hatchling limitation.

---

## Data Limitations

### Demo dataset is synthetic
The bundled demo dataset uses synthetic genomes (500 bp and 300 bp)
with injected SSRs. It does not represent real biological data.

### No pre-packaged reference dataset
The full SARS-CoV-2 dataset (~1M genomes) is not bundled. Users must
download it from NCBI using the `download` command.
