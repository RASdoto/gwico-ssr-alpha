# GWICO-SSR Statistical Methods

## Overview

GWICO-SSR implements six statistical analyses with Benjamini-Hochberg
FDR correction applied across all tests in a suite run.

## Analyses

### 1. Kruskal-Wallis Test (SSR Burden by Country)

**Purpose:** Test whether SSR count per genome differs significantly
across countries.

**Method:** Non-parametric one-way ANOVA on ranks.

**Effect size:** Eta-squared (η² = H / (N - 1)).

**Input:** Per-accession SSR counts grouped by country.

**Minimum requirement:** ≥ 2 countries with ≥ 2 accessions each.

### 2. Chi-Square Test (Gene × Country)

**Purpose:** Test whether SSR distribution across genes is
independent of geographic origin.

**Method:** Pearson chi-square test on contingency table.

**Effect size:** Cramér's V = √(χ² / (N × min(r-1, c-1))).

### 3. Chi-Square Test (Motif × Country)

**Purpose:** Test whether canonical motif frequencies vary
by country.

**Method:** Same as gene × country but for motif counts.

### 4. Correlation: Genome Length vs SSR Count

**Purpose:** Quantify the relationship between genome size
and SSR burden.

**Methods:**
- Pearson r (linear association)
- Spearman ρ (monotonic association)

**Confidence intervals:** Fisher z-transformation for Pearson r.

### 5. Correlation: GC Content vs SSR Count

**Purpose:** Test whether GC composition predicts SSR density.

**Methods:** Same as length correlation.

### 6. Shannon Entropy of Motif Distributions

**Purpose:** Quantify motif diversity.

**Method:** H = -Σ(pᵢ × log₂(pᵢ)) for each accession or country.

**Normalization:** Pielou's evenness J = H / log₂(S) where S is the
number of unique motifs.

### 7. Base Composition Test

**Purpose:** Test whether AT/GC composition deviates from
uniform (25% each base).

**Method:** Chi-square goodness-of-fit test.

## Multiple Testing Correction

**Method:** Benjamini-Hochberg False Discovery Rate (FDR).

All p-values from a suite run are collected and corrected together.
Both raw and corrected p-values are stored.

## Output Schema

Each result is stored as a `StatisticalResult` record with:
- `analysis_name`: identifier for the analysis
- `grouping`: variable used for grouping (e.g., "country")
- `metric`: response variable (e.g., "ssr_count_total")
- `test_name`: statistical test (e.g., "kruskal_wallis")
- `statistic`: test statistic value
- `p_value`: raw p-value
- `p_value_corrected`: BH-FDR corrected p-value
- `effect_size`: effect size measure
- `ci_lower`, `ci_upper`: confidence interval bounds
- `n`: sample size
- `metadata_json`: additional method-specific details

## Metric Definitions

| Metric | Formula | Units |
|--------|---------|-------|
| RA (Relative Abundance) | SSR_count / (genome_length / 1000) | SSRs/kb |
| RD (Relative Density) | SSR_bp_total / (genome_length / 1000000) | bp/Mb |
| GC Content | (G + C) / total_bases | fraction (0–1) |
| Dominant Motif | canonical motif with highest count | string |
