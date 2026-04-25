# GWICO-SSR Operator Runbook

**For Large-Cohort (1M+ Accession) Batch-First Workflows**

Version: 0.1.0-beta1  
Date: April 2026

---

## 1. Overview

This runbook guides operators through production-scale GWICO-SSR runs on datasets of 100K–1M accessions. It assumes familiarity with the GWICO-SSR CLI and focuses on batch-first workflows with proper failure recovery.

### Key Principles

- **Batch-First Default**: Large cohorts use 10K-accession batches by default
- **Failure Recovery**: Single-accession granularity for retry even in batched workflows
- **Provenance Tracking**: Every artifact has checksums and lineage
- **Reproducibility**: Deterministic outputs verified on re-runs

### Expected Scale

For SARS-CoV-2 style genomes (~30Kbp, perfect+imperfect SSRs):
- **1M accessions** → 60–90 hours single-threaded
- With modern hardware (8-core): ~15–20 hours with parallelization
- Disk space: ~500GB–1TB for database + artifacts

---

## 2. Pre-Run Checklist

### 2.1 Hardware Requirements

**Minimum Configuration** (for 100K accessions):
- CPU: 4 cores
- RAM: 16 GB
- Disk: 100 GB SSD (database + temp files)
- Network: Stable connection for NCBI downloads

**Recommended Configuration** (for 1M accessions):
- CPU: 8 cores (more for parallelization)
- RAM: 64 GB
- Disk: 1 TB SSD (database + artifacts + temp)
- Network: Dedicated API key for NCBI (5 req/sec vs 3 req/sec)

### 2.2 Dependencies

```bash
# Install with all optional dependencies
pip install -e ".[dev,pg]"

# Verify installation
python -m gwico_ssr --version  # Should be 0.1.0-beta1 or later
```

### 2.3 Database Setup

Choose database backend:

**SQLite (default, single-machine):**
```bash
# Set environment variable or config
export GWICO_SSR_DB_URL="sqlite:////path/to/gwico.db"
python -m gwico_ssr init-db
```

**PostgreSQL (distributed, multi-user):**
```bash
# Requires: psycopg2-binary installed (part of [pg] extra)
export GWICO_SSR_DB_URL="postgresql://user:pass@localhost:5432/gwico_ssr"
python -m gwico_ssr init-db

# Verify connection
psql -c "SELECT COUNT(*) FROM accession;" postgresql://user:pass@localhost:5432/gwico_ssr
```

### 2.4 NCBI Credentials

```bash
# Set in environment (strongly recommended for 1M accession runs)
export NCBI_API_KEY="your_api_key_here"
export NCBI_EMAIL="your_email@institution.org"

# Or in config file (config/default.toml):
[download]
ncbi_api_key = "your_api_key_here"
ncbi_email = "your_email@institution.org"
```

**API Key Benefits:**
- Rate limit: 10 req/sec (vs 3 req/sec without)
- Faster downloads: ~25% reduction for large cohorts
- Get free key at: https://www.ncbi.nlm.nih.gov/account/settings/

### 2.5 Data Preparation

Prepare metadata CSV with columns:
```
accession,country,host,species,date,lineage
NC_045512.2,China,human,SARS-CoV-2,2019-12-01,A
OL672836.1,China,human,SARS-CoV-2,2020-01-15,A.1
```

Split large metadata files by batch size:
```bash
# For 1M accessions with 10K batch size, create 100 files
split -l 10000 large_metadata.csv metadata_batch_

# Or use:
python -c "
import pandas as pd
df = pd.read_csv('large_metadata.csv')
for i, (name, group) in enumerate(df.groupby(df.index // 10000)):
    group.to_csv(f'metadata_batch_{i:03d}.csv', index=False)
"
```

---

## 3. Large-Cohort Workflow: Step-by-Step

### 3.1 Phase 1: Ingest and Manifest (0.5 hours for 1M)

**Goal**: Register all accessions, create download manifest.

```bash
# 1. Create dataset record
python -m gwico_ssr ingest metadata_batch_000.csv \
  --dataset-name "sars_cov_2_million" \
  --organism "SARS-CoV-2"

# 2. Add remaining batches
for i in {001..099}; do
  python -m gwico_ssr ingest metadata_batch_${i}.csv \
    --dataset-name "sars_cov_2_million" \
    --no-create-dataset
done

# 3. Verify ingestion
python -m gwico_ssr metrics --dataset-name sars_cov_2_million \
  | grep "Total accessions:"
# Expected output: Total accessions: 1000000
```

### 3.2 Phase 2: Download in Batches (36–48 hours for 1M)

**Goal**: Fetch all sequences from NCBI with automatic retry and recovery.

**Configuration for batch download:**

Create config file `large_cohort.toml`:
```toml
[download]
request_batch_size = 500           # NCBI requests (optimal)
artifact_batch_size = 10000        # Files to persist per batch
retry_max_attempts = 5             # Retry failed accessions
retry_backoff_seconds = 60         # Wait between retries

[database]
url = "sqlite:////data/gwico/sars_cov_2.db"

[logging]
level = "INFO"
format = "json"  # Machine-readable logs
```

**Run download with monitoring:**

```bash
# Start download with explicit batch config
python -m gwico_ssr download \
  --dataset-name sars_cov_2_million \
  --config large_cohort.toml \
  --output-dir /data/gwico/downloads \
  --show-progress

# Monitor in separate terminal
# (checks download progress every 60 seconds)
watch -n 60 'python -m gwico_ssr metrics --dataset-name sars_cov_2_million | tail -10'

# Typical output:
# Downloaded: 847,352 / 1,000,000 (84.7%)
# Failed: 152 (will retry)
# Skipped (not found): 3
# ETA: 4:30 remaining
```

**Failure Recovery:**

If download halts mid-run (network, power, quota):

```bash
# 1. Check what failed
python -m gwico_ssr download --dataset-name sars_cov_2_million \
  --check-status

# Output:
# Total: 1,000,000
# Completed: 847,352
# Failed (retryable): 152
# Failed (permanent): 3

# 2. Resume download (automatic retry)
python -m gwico_ssr download \
  --dataset-name sars_cov_2_million \
  --resume
  
# 3. Or retry only failed accessions
python -m gwico_ssr download \
  --dataset-name sars_cov_2_million \
  --failed-only --retry-max-attempts 10

# 4. Export failed accessions for manual review
python -m gwico_ssr export \
  --dataset-name sars_cov_2_million \
  --export-type failed_accessions \
  --output failed.csv
```

**Time Estimates:**
- With API key: ~36–48 hours (150K–200K acc/hour)
- Without API key: ~72–96 hours (100K–150K acc/hour)
- Network-bound, not CPU-bound

### 3.3 Phase 3: Parse Downloaded Sequences (2–4 hours for 1M)

**Goal**: Extract sequences and gene annotations into database.

**Configuration for batch parsing:**

```toml
[parsing]
chunk_size = 1000              # Parse 1K accessions at a time
format_detection = "auto"      # Auto-detect FASTA vs GenBank
malformed_strategy = "skip"    # Skip unparseable records
feature_extraction = true      # Extract gene/CDS features
```

**Run parsing:**

```bash
# Parse all downloaded sequences
python -m gwico_ssr parse \
  --dataset-name sars_cov_2_million \
  --config large_cohort.toml \
  --data-dir /data/gwico/downloads \
  --show-progress

# Expected output:
# Parsed: 999,997 sequences
# Malformed (skipped): 3
# Parse time: 3:24 elapsed
```

**Validation:**

```bash
# Verify parse output
python -m gwico_ssr metrics --dataset-name sars_cov_2_million \
  | grep -E "Sequences|Genes"

# Expected:
# Total sequences: 999,997
# Total genes annotated: ~12,000,000
# Mean genes/sequence: 12.0
```

### 3.4 Phase 4: Detect SSRs (4–8 hours for 1M)

**Goal**: Find all perfect, imperfect, and compound SSRs.

**Configuration for large-scale detection:**

```toml
[detection]
motif_sizes = [1, 2, 3, 4, 5, 6]    # Standard sizes
perfect = true                       # Always enable
imperfect = true                     # IMEX-style imperfect detection
imperfect_max_mismatches = 2         # Allow up to 2 errors
compound = true                      # Enable compound chaining
compound_dmax = 1                    # Adjacent motifs within 1bp
standardization_level = "full"       # Full motif normalization
```

**Run detection (example with 8-core parallelization):**

```bash
# Single-threaded (baseline)
python -m gwico_ssr detect \
  --dataset-name sars_cov_2_million \
  --config large_cohort.toml

# Or with GNU parallel (8 processes, 100K acc/process)
parallel -j 8 'python -m gwico_ssr detect \
  --dataset-name sars_cov_2_million \
  --accession-batch {} \
  --config large_cohort.toml' ::: $(seq -f "%06.0f" 0 100000 1000000)
```

**Expected output:**
```
Detection summary:
- Total SSRs found: 847,352,000 (≈850 SSRs/genome)
- Perfect: 652,000,000 (77%)
- Imperfect: 185,000,000 (22%)
- Compound: 10,352,000 (1%)
- Mean SSRs/genome: 847
- Detection time: 5:42 elapsed
```

### 3.5 Phase 5: Annotate (2–4 hours for 1M)

**Goal**: Map SSRs to genes and compute feature overlap.

```bash
# Run annotation
python -m gwico_ssr annotate \
  --dataset-name sars_cov_2_million \
  --annotation-source refseq \
  --feature-types "gene,CDS,misc_feature" \
  --show-progress

# Verify
python -m gwico_ssr metrics --dataset-name sars_cov_2_million \
  | grep -E "Annotated|Genic|UTR"
```

### 3.6 Phase 6: Comprehensive Analysis (<1 hour for 1M)

**Goal**: Statistical analysis, lineage correlation, temporal trends.

```bash
# 1. Run all analyses
python -m gwico_ssr analyze \
  --dataset-name sars_cov_2_million \
  --analyses all \
  --correlation-methods pearson,spearman \
  --fdr-correction benjamini-hochberg \
  --output /data/gwico/results

# 2. Generate visualizations
python -m gwico_ssr visualize \
  --dataset-name sars_cov_2_million \
  --figure-types all \
  --output /data/gwico/figures

# 3. Export for publication
python -m gwico_ssr export \
  --dataset-name sars_cov_2_million \
  --export-type publication_tables \
  --output /data/gwico/manuscript_tables
```

### 3.7 Phase 7: Export and Archive (0.5–1 hour for 1M)

**Goal**: Generate outputs, compute checksums, create manifest.

```bash
# Export all formats
python -m gwico_ssr export \
  --dataset-name sars_cov_2_million \
  --export-type all \
  --output /data/gwico/final_exports

# Create run manifest with provenance
python -m gwico_ssr export \
  --dataset-name sars_cov_2_million \
  --export-type manifest \
  --output /data/gwico/RUN_MANIFEST.json

# Verify checksums
cd /data/gwico/final_exports
sha256sum *.{csv,bed,gff3,json} > SHA256SUMS
sha256sum -c SHA256SUMS
```

---

## 4. Resume After Interruption

### 4.1 Identifying Failure Point

```bash
# Check run status
python -m gwico_ssr metrics --dataset-name sars_cov_2_million

# Check database for completed phases
sqlite3 gwico.db "SELECT COUNT(*) FROM parsed_sequence;" # After parsing
sqlite3 gwico.db "SELECT COUNT(*) FROM ssr_record;"      # After detection
sqlite3 gwico.db "SELECT COUNT(*) FROM run;"              # Runs completed
```

### 4.2 Resume Specific Phase

```bash
# Resume parsing (from where it left off)
python -m gwico_ssr parse \
  --dataset-name sars_cov_2_million \
  --resume

# Resume detection
python -m gwico_ssr detect \
  --dataset-name sars_cov_2_million \
  --resume

# Resume annotation
python -m gwico_ssr annotate \
  --dataset-name sars_cov_2_million \
  --resume
```

---

## 5. Performance Tuning

### 5.1 Batch Size Trade-offs

| Batch Size | Duration | Memory | Disk I/O | Recovery |
|-----------|----------|--------|----------|----------|
| 1K | Fast-track | Low | High | Best |
| 5K | Balanced | Medium | Medium | Good |
| 10K | Recommended | Medium | Low | Good |
| 20K | High-throughput | High | Low | Risky |

**Recommendation**: Start with 10K, adjust based on available RAM.

### 5.2 Memory Optimization

```toml
# For memory-constrained systems
[parsing]
chunk_size = 500
cache_size = "512MB"

[detection]
batch_size = 500
sliding_window_cache = true
```

### 5.3 Disk Space Estimates

For 1M SARS-CoV-2 accessions:
- Sequences downloaded: ~50 GB
- Database (SQLite): ~150 GB
- Intermediate files: ~100 GB
- Export formats: ~50 GB
- **Total**: ~350 GB minimum

---

## 6. Monitoring and Alerts

### 6.1 Log Monitoring

```bash
# Real-time log tail
tail -f /var/log/gwico_ssr.log | grep -E "ERROR|WARN"

# JSON log parsing
tail -f /var/log/gwico_ssr.log | jq 'select(.level == "ERROR")'

# Count errors by type
grep ERROR /var/log/gwico_ssr.log | jq -r '.error_type' | sort | uniq -c
```

### 6.2 Database Health Checks

```bash
# Check database integrity
python -m gwico_ssr admin check-db

# Vacuum and optimize
python -m gwico_ssr admin vacuum-db

# Backup before large operations
cp gwico.db gwico.db.backup.$(date +%Y%m%d)
```

---

## 7. Troubleshooting

### 7.1 Download Stalls

```bash
# Symptom: Download stuck at 80%

# 1. Check network
ping ncbi.nlm.nih.gov

# 2. Check API key validity
curl -I "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?api_key=${NCBI_API_KEY}"

# 3. Resume with different batch size
python -m gwico_ssr download \
  --dataset-name sars_cov_2_million \
  --resume \
  --request-batch-size 100  # Smaller batches
```

### 7.2 Out of Memory Errors

```bash
# Symptom: MemoryError during detection

# 1. Reduce chunk size in config
sed -i 's/chunk_size = 1000/chunk_size = 100/' large_cohort.toml

# 2. Or resume with explicit limit
python -m gwico_ssr detect \
  --dataset-name sars_cov_2_million \
  --memory-limit 16GB \
  --resume
```

### 7.3 Annotation Failures

```bash
# Symptom: Some accessions fail annotation

# 1. Check which ones failed
python -m gwico_ssr export \
  --dataset-name sars_cov_2_million \
  --export-type failed_annotations \
  --output failed_annot.csv

# 2. Try with different annotation source
python -m gwico_ssr annotate \
  --dataset-name sars_cov_2_million \
  --annotation-source ncbi \
  --resume
```

---

## 8. Publication Checklist

Before publishing results:

- [ ] All 1M accessions processed
- [ ] <0.1% failed accessions (< 1000)
- [ ] SSR detection reproducible (re-run produces identical results)
- [ ] Statistical tests use BH-FDR correction
- [ ] Lineage/temporal context validated
- [ ] Browser tracks verified in IGV/UCSC
- [ ] Cluster assignments stable (silhouette > 0.5 for main clusters)
- [ ] RUN_MANIFEST.json archived with paper
- [ ] Data availability statement includes accession list + metadata

---

## 9. Support and Escalation

**For issues or questions:**

1. Check logs: `tail -f /var/log/gwico_ssr.log`
2. Consult KNOWN_LIMITATIONS.md for deferred features
3. Run built-in diagnostics: `python -m gwico_ssr admin diagnostics`
4. Contact: [support email if applicable]

---

**Last Updated:** April 2026  
**Version:** 0.1.0-beta1
