# GWICO-SSR Alpha — Known Limitations

**Version:** 0.1.0-alpha1

This document lists known limitations, constraints, and deferred features for
the alpha release. Items are grouped by category and ranked by severity.

---

## Critical Limitations

### No parallel processing
SSR detection and annotation run sequentially. For the full SARS-CoV-2 dataset
(~1M genomes), processing would take significant time single-threaded.

**Mitigation:** The pipeline is designed for per-accession parallelism.
The orchestration layer supports this architecturally; multiprocessing
integration is planned for the first post-alpha release.

### NCBI download requires real credentials
The `download` command requires a valid NCBI API key and email.
Without an API key, rate limiting is 3 requests/second. All download
tests use mocks — no integration tests against live NCBI.

---

## Functional Limitations

### Perfect SSRs only
Only perfect (exact) tandem repeats are detected. Imperfect SSRs
(with mismatches) and compound SSRs (adjacent different motifs)
are not supported.

### No lineage-aware analysis
Pango lineage, Nextstrain clade, and other variant annotations are
not integrated. SSR analysis is not correlated with viral evolution.

### No temporal analysis
Collection dates are stored but not used for trend analysis or
pandemic-wave correlation.

### Single-organism datasets
Each dataset is assumed to contain one organism. Cross-species
comparative analysis is not implemented.

### No genome circularity handling
SSRs that wrap around the origin of circular genomes are not detected.
For SARS-CoV-2 (linear), this is not relevant.

---

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
