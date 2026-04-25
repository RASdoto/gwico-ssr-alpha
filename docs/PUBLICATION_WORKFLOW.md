# GWICO-SSR Publication Workflow

**Step-by-Step Guide to Publication-Ready Results**

Version: 0.1.0-beta1  
Date: April 2026

---

## 1. Overview

This guide walks researchers through generating publication-quality outputs from GWICO-SSR analysis. It covers:

- Data preparation and validation
- Statistical rigor and FDR correction
- Figure generation and formatting
- Supplementary table preparation
- Reproducibility documentation
- Data availability statements

**Audience**: Biologists and bioinformaticians preparing manuscripts.

---

## 2. Pre-Analysis Checklist

### 2.1 Data Quality Validation

Before beginning analysis, verify:

```bash
# 1. Check metadata completeness
python -m gwico_ssr admin check-metadata \
  --dataset-name my_dataset \
  --output quality_report.txt

# Expected output:
# Total accessions: 1,000,000
# Metadata completeness: 99.8% (1,998 missing country)
# Lineage coverage: 98.5% (15,000 unknown)
# Collection date range: 2019-12-01 to 2023-12-31

# 2. Verify SSR detection consistency
python -c "
import pandas as pd
from sqlalchemy import create_engine, func, select
from gwico_ssr.models import Accession, SSRRecord

engine = create_engine('sqlite:///gwico.db')
with engine.connect() as conn:
    # Check SSR count distribution
    stmt = select(func.count(SSRRecord.id)).group_by(SSRRecord.accession_id)
    result = conn.execute(stmt).scalars().all()
    print(f'Mean SSRs/genome: {sum(result)/len(result):.1f}')
    print(f'Min: {min(result)}, Max: {max(result)}')
"

# 3. Verify run manifest
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type manifest \
  --output RUN_MANIFEST.json

# Verify SHA256 checksums
cd /output_dir
sha256sum -c SHA256SUMS.txt
```

### 2.2 Metadata Annotation

Ensure metadata includes:
- Accession (strain ID)
- Collection date (YYYY-MM-DD)
- Country of origin
- Host species
- Lineage/clade (if applicable)
- Geographic coordinates (optional, for choropleth)

Example CSV format:
```csv
accession,collection_date,country,host,species,lineage
NC_045512.2,2019-12-01,China,Human,SARS-CoV-2,A
OL672836.1,2020-01-15,China,Human,SARS-CoV-2,A.1
```

---

## 3. Full Publication Workflow

### 3.1 Phase 1: Run Analysis with FDR Correction

```bash
# Run comprehensive analysis with BH-FDR correction
python -m gwico_ssr analyze \
  --dataset-name my_dataset \
  --analyses all \
  --fdr-correction benjamini-hochberg \
  --fdr-alpha 0.05 \
  --correlation-methods pearson,spearman \
  --output analysis_output/

# Expected outputs:
# analysis_output/
#   ├── chi_square_results.csv
#   ├── kruskal_wallis_results.csv
#   ├── correlation_matrix.csv
#   ├── entropy_results.csv
#   └── lineage_comparison.csv
```

### 3.2 Phase 2: Generate Publication Figures

```bash
# Generate all figure types
python -m gwico_ssr visualize \
  --dataset-name my_dataset \
  --figure-types all \
  --output figures/manuscript/ \
  --dpi 300 \
  --format {pdf,png}

# Creates publication-ready figures:
# figures/manuscript/
#   ├── figure_1_ssr_distribution.pdf
#   ├── figure_2_repeat_class_heatmap.pdf
#   ├── figure_3_geographic_choropleth.pdf
#   ├── figure_4_phylogenetic_trees.pdf
#   ├── figure_5_clustering_pca.pdf
#   ├── supplementary_boxplots_by_lineage.pdf
#   └── supplementary_temporal_trends.pdf
```

### 3.3 Phase 3: Create Supplementary Tables

```bash
# Export publication-ready tables
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type publication_tables \
  --output tables/

# Creates formatted tables:
# tables/
#   ├── Table_S1_Accession_Summary.xlsx
#   ├── Table_S2_SSR_Characteristics.xlsx
#   ├── Table_S3_Statistical_Analysis.xlsx
#   ├── Table_S4_Lineage_Enrichment.xlsx
#   └── Table_S5_Clustering_Results.xlsx
```

### 3.4 Phase 4: Generate Manuscript Supplementary Files

```bash
# Create BED tracks for supplementary
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type browser_tracks \
  --output supp_tracks/

# GFF3 for genome browsers
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type gff3 \
  --feature-types perfect,imperfect,compound \
  --output supp_gff/
```

---

## 4. Statistical Rigor Checklist

### 4.1 Multiple Testing Correction

✅ **Implemented in GWICO-SSR**:
- Benjamini-Hochberg FDR correction (default α=0.05)
- Bonferroni correction (alternative)
- Adjusted p-values reported in all results

**Verification**:
```python
import pandas as pd

results = pd.read_csv('chi_square_results.csv')

# Check FDR correction
print(results[['p_value', 'p_adjusted', 'significant']].head(10))

# All significant results should have p_adjusted ≤ 0.05
assert (results[results['significant']] ['p_adjusted'] <= 0.05).all()
```

### 4.2 Effect Size Reporting

**Always report alongside p-values:**

```bash
# Create summary table with effect sizes
python -c "
import pandas as pd

stats = pd.read_csv('kruskal_wallis_results.csv')
stats['effect_size'] = stats['eta_squared']  # Effect size from test

# Filter significant results (p_adj < 0.05) with meaningful effect (eta > 0.1)
significant = stats[(stats['p_adjusted'] < 0.05) & (stats['effect_size'] > 0.1)]

# Create publication table
significant[['feature', 'statistic', 'p_value', 'p_adjusted', 
             'effect_size']].to_csv('Table_Main_Statistics.csv', index=False)

print(f'Reported {len(significant)} significant findings')
"
```

### 4.3 Reproducibility Verification

```bash
# Re-run identical analysis, verify outputs match exactly
python -m gwico_ssr analyze \
  --dataset-name my_dataset \
  --analyses all \
  --output analysis_rerun/

# Verify SHA256 checksums match
sha256sum analysis_output/*.csv > run1.sha256
sha256sum analysis_rerun/*.csv > run2.sha256

# Should be identical (bit-identical, not approximate)
diff run1.sha256 run2.sha256 || echo "REPRODUCIBILITY CHECK FAILED"
```

---

## 5. Figure Guidelines

### 5.1 Figure 1: SSR Overview

**Components**:
- Distribution of SSRs per genome (histogram)
- Repeat class breakdown (stacked bar)
- Motif size distribution (line graph)

**Code**:
```python
from gwico_ssr.visualization import SSRDistributionPlotter

plotter = SSRDistributionPlotter(dataset_name='my_dataset')
fig = plotter.create_figure_1()
fig.savefig('figures/Figure_1_SSR_Overview.pdf', dpi=300, bbox_inches='tight')
```

### 5.2 Figure 2: Geographic Distribution

**Components**:
- World choropleth of SSR density by country
- Top 10 countries by SSR count
- Geographic clustering

**Reproducibility Note**: Include exact GeoJSON source in Methods.

### 5.3 Figure 3: Temporal Trends

**Components**:
- SSR characteristics over time (line plots)
- Lineage frequency over time (stacked area)
- Detected lineage turnover points

### 5.4 Figure 4: Clustering Results

**Components**:
- PCA scatter plot with cluster labels
- Cluster silhouette plot
- Cluster composition (by lineage)

**Include in figure legend**:
- Algorithm: K-means with k=10
- Seed: 42 (for reproducibility)
- Features: [list 7 features used]
- Silhouette score: 0.X

### 5.5 Figure 5: Phylogenetic Context

**Components**:
- Phylogenetic tree colored by dominant repeat class
- SSR count mapped to branch length
- Clade-level enrichment heatmap

---

## 6. Supplementary Material Preparation

### 6.1 Supplementary Data Files

**File S1: Accession Metadata**
- All 1,000,000+ accessions with original metadata

```bash
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type accession_metadata \
  --output File_S1_Metadata.csv
```

**File S2: SSR Coordinates (BED Format)**
- All detected SSRs with coordinates, repeat class, motif

```bash
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type bed \
  --feature-types perfect,imperfect,compound \
  --output File_S2_SSRs.bed.gz
```

**File S3: Genome Browser Tracks**
- UCSC-compatible BED12 files, one per repeat class

```bash
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type browser_tracks \
  --output File_S3_Tracks/
```

### 6.2 Supplementary Methods Section

```markdown
## Supplementary Methods

### Data Processing and SSR Detection

GWICO-SSR pipeline version 0.1.0-beta1 was used as follows:

1. **Ingest**: {accession_count} accessions loaded from NCBI
2. **Download**: Sequences retrieved via Entrez with rate-limiting
3. **Parse**: FASTA/GenBank sequences parsed with feature extraction
4. **Detect**: 
   - Perfect SSRs: Motif sizes 1–6
   - Imperfect SSRs: Up to 2 mismatches allowed
   - Compound SSRs: Adjacent motifs within 1 bp
5. **Annotate**: Mapped to {gene_count} reference genes
6. **Analysis**: 
   - Chi-square test for repeat class distributions
   - Kruskal-Wallis test for continuous features
   - Spearman correlation for trait associations
   - Benjamini-Hochberg FDR correction (α=0.05)

### Reproducibility

The analysis is fully reproducible:
- Random seed: 42 (clustering, PCA)
- Database manifest: RUN_MANIFEST.json
- SHA256 checksums: SHA256SUMS.txt
- Full command history: command_history.txt

### Data Availability

All input and output files are available at: [ZENODO/DRYAD/FigShare URL]
```

---

## 7. Data Availability Statement

Use this template:

```markdown
## Data Availability

The SSR coordinates, metadata, and analysis results are available at:

**Repository**: Zenodo (DOI: 10.5281/zenodo.XXXXXXX)
- Accession list: accessions.csv.gz
- SSR coordinates: ssrs.bed.gz
- Metadata: metadata.csv.gz
- Analysis results: analysis_results.xlsx
- Reproducibility manifest: RUN_MANIFEST.json
- Checksums: SHA256SUMS.txt

**Source sequences**: NCBI GenBank (accession list in File S1)

**GWICO-SSR source code**: GitHub (https://github.com/...)

**Reproducibility**: The full pipeline can be re-run using:
```bash
python examples/publication_workflow.py --dataset my_dataset
```
```

---

## 8. Methods Section Template

```markdown
## Methods

### SSR Detection

Simple Sequence Repeats (SSRs) were detected using GWICO-SSR (version 0.1.0-beta1).

**Definition**: Perfect SSRs consist of an exact tandem repeat of a motif (1–6 nucleotides) at least twice. Imperfect SSRs allow up to 2 nucleotide mismatches. Compound SSRs are adjacent SSRs of different motifs within 1 bp.

**Implementation**: Detection uses a sliding-window algorithm with O(n) time complexity. Interval trees enable O((log n + k)) annotation mapping, where n is the genome length and k is the number of features.

### Statistical Analysis

All tests used Benjamini-Hochberg FDR correction with significance threshold α=0.05.

**Chi-square**: Test for independence of SSR repeat class and discrete traits (lineage, country).

**Kruskal-Wallis**: Non-parametric test for SSR density differences across groups.

**Spearman correlation**: Association between continuous features (SSR density, GC content) and viral traits.

### Clustering

K-means clustering (k=10, seed=42) was performed on standardized 7-feature vectors:
1. SSR count per accession
2. Total repeat length (bp)
3. Total repeat units
4. Mean repeat length
5. Perfect SSR count
6. Imperfect SSR count
7. Compound SSR count

Features were standardized using z-score normalization. Silhouette score: 0.XX.

### Phylogenetic Context

Maximum-likelihood phylogenetic trees were constructed using IQ-TREE2 with GTR+G model. SSR characteristics were mapped to branches using GARD-inferred recombination-aware alignment.

### Software and Reproducibility

All analyses were performed using GWICO-SSR v0.1.0-beta1. The complete command history, parameter configuration, and checksums are provided in RUN_MANIFEST.json. Results are reproducible to the bit (verified by re-running identical analysis).
```

---

## 9. Figure Quality Checklist

Before submitting figures to journal:

- [ ] Figure resolution ≥300 DPI for print
- [ ] Fonts ≥10pt and consistent across figures
- [ ] Color schemes accessible (colorblind-safe)
- [ ] Legend includes sample size (n=1,000,000)
- [ ] Error bars or confidence intervals shown
- [ ] Statistical test and p-value in caption
- [ ] Figure generated from **final** dataset (no revisions)
- [ ] SHA256 checksum provided for traceability

---

## 10. Revision-Proof Workflow

### 10.1 Version Control for Analysis

```bash
# 1. Tag data snapshot
git tag -a "analysis_v1.0" -m "Final for submission"

# 2. Export manifests and checksums
python -m gwico_ssr export \
  --dataset-name my_dataset \
  --export-type manifest \
  --output RUN_MANIFEST_v1.json

sha256sum analysis_output/* > SHA256SUMS_v1.txt

# 3. Create reproducibility package
tar -czf analysis_v1.tar.gz \
  RUN_MANIFEST_v1.json \
  SHA256SUMS_v1.txt \
  analysis_output/ \
  figures/ \
  tables/

# Upload as supplementary material
```

### 10.2 Handling Reviewer Requests

**Scenario**: Reviewer asks "Can you re-run with different SSR motif sizes?"

```bash
# Create new analysis variant
python -m gwico_ssr detect \
  --dataset-name my_dataset \
  --motif-sizes 1,2,3 \  # Reduced set
  --output analysis_motif_1_3/

# Generate new figures for response
python -m gwico_ssr visualize \
  --dataset-name my_dataset \
  --output figures/revision_motif_1_3/

# Document changes
echo "Reviewer response: Motif sizes reduced to 1-3. SSR counts decreased by 15%." \
  > REVISION_NOTES.txt
```

---

## 11. Quality Assurance Before Submission

### 11.1 Final Validation Script

```python
#!/usr/bin/env python
"""Validate publication-ready outputs."""

import subprocess
import pandas as pd
import hashlib
from pathlib import Path

def validate_publication():
    """Run all checks before submission."""
    
    print("🔍 PUBLICATION VALIDATION CHECKLIST")
    print("=" * 50)
    
    # 1. Check SHA256 checksums
    print("\n✓ Verifying data integrity...")
    result = subprocess.run(["sha256sum", "-c", "SHA256SUMS.txt"], 
                          capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ❌ FAILED: {result.stderr}")
        return False
    print("  ✅ All files verified")
    
    # 2. Check statistical results
    print("\n✓ Checking statistical analysis...")
    stats = pd.read_csv("analysis_output/kruskal_wallis_results.csv")
    sig_results = stats[stats['p_adjusted'] < 0.05]
    print(f"  ✅ {len(sig_results)} significant results (p_adj < 0.05)")
    
    # 3. Verify figures exist and are not empty
    print("\n✓ Verifying figures...")
    fig_dir = Path("figures/manuscript")
    for fig in fig_dir.glob("*.pdf"):
        if fig.stat().st_size < 10000:
            print(f"  ⚠ WARNING: {fig.name} is suspiciously small")
        else:
            print(f"  ✅ {fig.name}")
    
    # 4. Check supplementary files
    print("\n✓ Checking supplementary files...")
    supp_files = [
        "File_S1_Metadata.csv",
        "File_S2_SSRs.bed.gz",
        "File_S3_Tracks/browser_tracks.bed",
    ]
    for fname in supp_files:
        if Path(fname).exists():
            print(f"  ✅ {fname}")
        else:
            print(f"  ❌ {fname} MISSING")
            return False
    
    # 5. Reproducibility check
    print("\n✓ Reproducibility verification...")
    manifest = pd.read_json("RUN_MANIFEST.json")
    print(f"  ✅ Manifest created at {manifest['timestamp']}")
    print(f"  ✅ Dataset: {manifest['dataset_name']}")
    print(f"  ✅ Total accessions: {manifest['accession_count']}")
    
    print("\n" + "=" * 50)
    print("✅ ALL CHECKS PASSED - Ready for submission!")
    return True

if __name__ == "__main__":
    validate_publication()
```

---

## 12. Common Issues and Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| "p-values all < 0.05" | Multiple testing not corrected | Enable FDR, check for p-hacking |
| Figures not reproducible | Random seed not set | Use `seed=42`, export RUN_MANIFEST |
| Missing effect sizes | Reporting p-values only | Add eta², Cohen's d, OR to table |
| Browser tracks won't load | Coordinate format error | Verify BED12 is 0-based half-open |
| Conflicting SSR counts | Detection parameters differ | Document exact --motif-sizes used |

---

## 13. References

- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate. *JRSB*, 57(1), 289–300.
- Greenland, S. (2019). *Statistical tests, P values, confidence intervals, and power*. SAGE Publications.
- Reproducibility guidelines: https://www.nature.com/articles/s41596-020-00409-y

---

**Last Updated:** April 2026  
**Version:** 0.1.0-beta1
