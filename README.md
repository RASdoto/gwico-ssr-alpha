# GWICO-SSR

**In-silico SSR identification and characterization platform for complete genome datasets.**

GWICO-SSR automates the detection, annotation, and statistical analysis of Simple Sequence Repeats (SSRs / microsatellites) across large genomic datasets. It replaces a collection of standalone scripts with a single, database-backed, reproducible CLI pipeline.

**Version:** 0.1.0-alpha1

## Quick Start

### Installation

```bash
# Clone and install in development mode
cd "GWICO-SSR Alpha"
pip install -e ".[dev]"

# Verify installation
python -m gwico_ssr --version
python -m gwico_ssr --help
```

### Demo Run

Run the bundled end-to-end demo to verify your installation:

```bash
python examples/demo_run.py
```

This runs the full pipeline on two synthetic genomes and produces outputs in `outputs/demo_run/`.

### Step-by-Step Workflow

```bash
# 1. Initialize the database
python -m gwico_ssr init-db

# 2. Ingest accession metadata from CSV
python -m gwico_ssr ingest data/demo/demo_metadata.csv --dataset-name demo

# 3. Parse sequences and gene annotations
python -m gwico_ssr parse demo --data-dir data/demo

# 4. Detect SSRs (motif sizes 1-6)
python -m gwico_ssr detect demo --data-dir data/demo

# 5. Annotate SSRs against gene features
python -m gwico_ssr annotate demo

# 6. Compute per-accession metrics (RA, RD, motif counts)
python -m gwico_ssr metrics demo

# 7. Run statistical analyses with FDR correction
python -m gwico_ssr analyze demo

# 8. Generate publication-ready figures
python -m gwico_ssr visualize demo --output-dir outputs/figures --format png

# 9. Export in multiple formats
python -m gwico_ssr export demo --output-dir outputs/exports --formats csv,json,bed,gff3,tables
```

### One-Command Pipeline

The `run` command executes the full pipeline with checkpointing:

```bash
# Run all stages
python -m gwico_ssr run demo --data-dir data/demo

# Resume an interrupted run
python -m gwico_ssr run demo --resume <RUN_ID>

# Dry run to preview stages
python -m gwico_ssr run demo --data-dir data/demo --dry-run
```

## Configuration

Configuration is loaded in priority order:

1. CLI flags (`--log-level`, `--log-format`, `--config`)
2. Environment variables (`GWICO_SSR_DB_URL`, `GWICO_SSR_LOG_LEVEL`, etc.)
3. TOML config file (explicit path → `GWICO_SSR_CONFIG` env var → `config/default.toml`)
4. Built-in defaults

### Config Profiles

| Profile | File | Purpose |
|---------|------|---------|
| Default | `config/default.toml` | Standard local execution |
| Development | `config/development.toml` | Verbose logging, SQL echo |
| Batch | `config/batch.toml` | Large dataset processing |

```bash
# Use a specific config profile
python -m gwico_ssr --config config/batch.toml run demo --data-dir data/

# Override config with environment variables
export GWICO_SSR_DB_URL=postgresql://user:pass@host/gwico
export GWICO_SSR_NCBI_API_KEY=your_key_here
```

See [.env.example](.env.example) for all available environment variables.

### SSR Detection Thresholds

Default thresholds (minimum repeat units before an SSR is reported):

| Motif Size | Name | Default |
|-----------|------|---------|
| 1 bp | Mononucleotide | ≥ 10 |
| 2 bp | Dinucleotide | ≥ 5 |
| 3 bp | Trinucleotide | ≥ 3 |
| 4 bp | Tetranucleotide | ≥ 3 |
| 5 bp | Pentanucleotide | ≥ 3 |
| 6 bp | Hexanucleotide | ≥ 2 |

These are configurable in the TOML file under `[ssr]` or via `GWICO_SSR_SSR_MIN_REPEATS_*` env vars.

## CLI Reference

| Command | Description |
|---------|-------------|
| `info` | Show current configuration and SSR thresholds |
| `init-db` | Create database tables (use `--drop` to reset) |
| `ingest` | Load accession metadata from CSV |
| `download` | Download sequences from NCBI Entrez |
| `parse` | Parse FASTA/GenBank/GFF3 files |
| `detect` | Run SSR detection on parsed sequences |
| `annotate` | Map SSRs to gene features via interval trees |
| `metrics` | Compute per-accession SSR metrics |
| `analyze` | Run statistical analyses with FDR correction |
| `visualize` | Generate static and interactive figures |
| `export` | Export CSV, JSON, BED, GFF3, publication tables |
| `run` | Execute full pipeline with checkpointing |

All data-processing commands support `--json-summary` for machine-readable output and `--force` to reprocess existing data.

## Architecture

```
gwico-ssr/
  pyproject.toml          # Package definition and dependencies
  config/
    default.toml          # Default configuration
    batch.toml            # Large-dataset configuration
    development.toml      # Development configuration
  examples/
    demo_run.py           # End-to-end demo script
  data/
    demo/                 # Demo dataset (2 synthetic genomes)
  src/gwico_ssr/
    cli/                  # CLI commands (Click)
    config/               # TOML + env config loading (Pydantic)
    logging/              # Structured JSON/text logging
    db/                   # Database engine, sessions, repositories
    models/               # ORM models (11 tables)
    ingest/               # CSV ingestion, NCBI download
    parsers/              # FASTA, GenBank, GFF3 parsers
    ssr/                  # SSR detection engine + motif utilities
    annotation/           # Interval-tree gene mapping
    metrics/              # Per-genome and cohort metrics
    analysis/             # Statistical analysis with FDR
    visualization/        # Publication-ready figures
    export/               # Multi-format export + provenance
    orchestration/        # Pipeline orchestration + checkpointing
    utils/                # Shared helpers
  tests/                  # 490 tests (unit, integration, benchmark)
  docs/                   # Architecture, methods, validation docs
```

### Pipeline Flow

```
CSV metadata → ingest → DB
                         ↓
NCBI / local → download → parse → sequences + features in DB
                                          ↓
                                    detect SSRs
                                          ↓
                                    annotate (interval tree)
                                          ↓
                                    compute metrics
                                          ↓
                                    statistical analysis
                                          ↓
                              visualize + export (CSV/JSON/BED/GFF3)
```

### Database

- **Default**: SQLite with WAL mode and foreign-key pragmas
- **Optional**: PostgreSQL via `pip install gwico-ssr[pg]`
- **11 tables**: Dataset, Run, Accession, SequenceRecord, FeatureRecord, SSRRecord, SSRAnnotation, AccessionMetrics, StatisticalResult, StageCheckpoint, FailedAccession
- **Incremental**: Skip-if-exists logic at every stage; `--force` to reprocess

### SSR Detection Engine

Native Python implementation replacing PERF:
- Perfect SSR detection for motif sizes 1–6 bp
- Lexicographic-minimum-rotation canonical motifs across both strands
- Sub-repeat filtering (e.g., AAGAAG as 6-mer suppressed when AAG×2 is the fundamental unit)
- 0-based half-open coordinates `[start, end)`
- Deterministic, reproducible output (validated against gold-standard fixtures)

### Statistical Methods

- **Kruskal-Wallis** with eta-squared effect size (SSR burden by country)
- **Chi-square** with Cramér's V (gene×country, motif×country)
- **Pearson/Spearman** correlation with Fisher-z confidence intervals
- **Shannon entropy** with Pielou's evenness index
- **Base composition** chi-square goodness-of-fit test
- **Benjamini-Hochberg FDR** correction applied across all tests

## Development

```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run tests with coverage
pytest --cov=gwico_ssr --cov-report=term-missing

# Run only benchmark tests
pytest tests/benchmarks/

# Run validation tests
pytest tests/unit/test_validation.py -v
```

### Technology Stack

| Component | Library |
|-----------|---------|
| Language | Python 3.10+ |
| CLI | Click |
| Config | Pydantic + TOML |
| Database | SQLAlchemy 2.0 (SQLite / PostgreSQL) |
| Bio | BioPython |
| Stats | SciPy, statsmodels, scikit-learn |
| Visualization | Plotly, Matplotlib, Seaborn |
| Testing | pytest + pytest-cov + pytest-mock |

## License

MIT
