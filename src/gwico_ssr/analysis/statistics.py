"""Statistical analysis workflows for GWICO-SSR.

Implements chi-square, Kruskal-Wallis, Pearson/Spearman correlation,
Shannon entropy, and Benjamini-Hochberg FDR correction.

All analyses persist results via the StatisticalResult model and return
structured result objects for programmatic access.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from typing import Sequence

import numpy as np
from scipy import stats as sp_stats
from statsmodels.stats.multitest import multipletests
from sqlalchemy import select, func
from sqlalchemy.orm import Session

from gwico_ssr.db.repository import insert_statistical_result
from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    SSRAnnotation,
    SSRRecord,
)


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class AnalysisResult:
    """Single statistical test output."""

    analysis_name: str
    grouping: str | None
    metric: str | None
    test_name: str
    statistic: float | None
    p_value: float | None
    p_value_corrected: float | None = None
    effect_size: float | None = None
    ci_lower: float | None = None
    ci_upper: float | None = None
    n: int | None = None
    metadata: dict | None = None

    def to_dict(self, run_id: int) -> dict:
        """Convert to kwargs suitable for insert_statistical_result."""
        return {
            "run_id": run_id,
            "analysis_name": self.analysis_name,
            "grouping": self.grouping,
            "metric": self.metric,
            "test_name": self.test_name,
            "statistic": self.statistic,
            "p_value": self.p_value,
            "p_value_corrected": self.p_value_corrected,
            "effect_size": self.effect_size,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "n": self.n,
            "metadata_json": json.dumps(self.metadata) if self.metadata else None,
        }


@dataclass
class AnalysisSuite:
    """Collection of analysis results with FDR correction support."""

    results: list[AnalysisResult] = field(default_factory=list)

    def apply_fdr(self, method: str = "fdr_bh") -> None:
        """Apply Benjamini-Hochberg (or other) FDR correction to all results."""
        p_vals = [r.p_value for r in self.results if r.p_value is not None]
        if not p_vals:
            return
        indices = [i for i, r in enumerate(self.results) if r.p_value is not None]
        _, corrected, _, _ = multipletests(p_vals, method=method)
        for idx, corr_p in zip(indices, corrected):
            self.results[idx].p_value_corrected = float(corr_p)

    def persist(self, session: Session, run_id: int) -> int:
        """Persist all results and return count."""
        count = 0
        for r in self.results:
            insert_statistical_result(session, **r.to_dict(run_id))
            count += 1
        session.flush()
        return count

    def to_list(self) -> list[dict]:
        """Serialize results for JSON output."""
        return [asdict(r) for r in self.results]


# ---------------------------------------------------------------------------
# Kruskal-Wallis: SSR burden by country
# ---------------------------------------------------------------------------

def kruskal_wallis_by_country(
    session: Session,
    run_id: int,
    metric_col: str = "ssr_count_total",
    min_group_size: int = 5,
) -> AnalysisResult | None:
    """Kruskal-Wallis H-test for SSR metric differences across countries.

    Groups with fewer than *min_group_size* accessions are excluded.
    Effect size: eta-squared (H / (N-1)).
    """
    col = getattr(AccessionMetrics, metric_col, None)
    if col is None:
        return None

    stmt = (
        select(Accession.country, col)
        .join(AccessionMetrics, AccessionMetrics.accession == Accession.accession)
        .where(AccessionMetrics.run_id == run_id)
        .where(Accession.country.is_not(None))
    )
    rows = session.execute(stmt).all()
    if not rows:
        return None

    # Group values by country
    groups: dict[str, list[float]] = {}
    for country, val in rows:
        if val is not None:
            groups.setdefault(country, []).append(float(val))

    # Filter small groups
    groups = {k: v for k, v in groups.items() if len(v) >= min_group_size}
    if len(groups) < 2:
        return None

    arrays = list(groups.values())
    n_total = sum(len(a) for a in arrays)

    h_stat, p_value = sp_stats.kruskal(*arrays)
    # Eta-squared: H / (N - 1)
    eta_sq = float(h_stat) / (n_total - 1) if n_total > 1 else None

    return AnalysisResult(
        analysis_name="kruskal_wallis_by_country",
        grouping="country",
        metric=metric_col,
        test_name="kruskal-wallis",
        statistic=float(h_stat),
        p_value=float(p_value),
        effect_size=eta_sq,
        n=n_total,
        metadata={
            "n_groups": len(groups),
            "group_sizes": {k: len(v) for k, v in groups.items()},
        },
    )


# ---------------------------------------------------------------------------
# Chi-square: gene × country contingency
# ---------------------------------------------------------------------------

def chi_square_gene_country(
    session: Session,
    dataset_id: int,
    min_cell_count: int = 5,
) -> AnalysisResult | None:
    """Chi-square test of independence for gene × country SSR distribution.

    Effect size: Cramér's V = sqrt(chi2 / (N * (min(r, c) - 1))).
    """
    stmt = (
        select(
            SSRAnnotation.gene_name,
            Accession.country,
            func.count().label("cnt"),
        )
        .join(Accession, SSRAnnotation.accession == Accession.accession)
        .where(SSRAnnotation.gene_name.is_not(None))
        .where(Accession.country.is_not(None))
    )
    # scope to dataset via accession
    stmt = stmt.where(Accession.dataset_id == dataset_id)
    stmt = stmt.group_by(SSRAnnotation.gene_name, Accession.country)
    rows = session.execute(stmt).all()

    if not rows:
        return None

    # Build contingency table
    genes: list[str] = sorted({r[0] for r in rows})
    countries: list[str] = sorted({r[1] for r in rows})
    if len(genes) < 2 or len(countries) < 2:
        return None

    gene_idx = {g: i for i, g in enumerate(genes)}
    country_idx = {c: i for i, c in enumerate(countries)}
    table = np.zeros((len(genes), len(countries)), dtype=float)
    for gene, country, cnt in rows:
        table[gene_idx[gene], country_idx[country]] = cnt

    n_total = int(table.sum())
    chi2, p_value, dof, _ = sp_stats.chi2_contingency(table)
    # Cramér's V
    k = min(len(genes), len(countries))
    cramers_v = math.sqrt(chi2 / (n_total * (k - 1))) if n_total > 0 and k > 1 else None

    return AnalysisResult(
        analysis_name="chi_square_gene_country",
        grouping="gene_x_country",
        metric="ssr_count",
        test_name="chi-square",
        statistic=float(chi2),
        p_value=float(p_value),
        effect_size=cramers_v,
        n=n_total,
        metadata={
            "dof": int(dof),
            "n_genes": len(genes),
            "n_countries": len(countries),
            "genes": genes,
        },
    )


# ---------------------------------------------------------------------------
# Chi-square: motif × country contingency
# ---------------------------------------------------------------------------

def chi_square_motif_country(
    session: Session,
    dataset_id: int,
) -> AnalysisResult | None:
    """Chi-square test for motif_size × country distribution."""
    stmt = (
        select(
            SSRRecord.motif_size,
            Accession.country,
            func.count().label("cnt"),
        )
        .join(Accession, SSRRecord.accession == Accession.accession)
        .where(Accession.country.is_not(None))
        .where(Accession.dataset_id == dataset_id)
        .group_by(SSRRecord.motif_size, Accession.country)
    )
    rows = session.execute(stmt).all()
    if not rows:
        return None

    sizes = sorted({r[0] for r in rows})
    countries = sorted({r[1] for r in rows})
    if len(sizes) < 2 or len(countries) < 2:
        return None

    size_idx = {s: i for i, s in enumerate(sizes)}
    country_idx = {c: i for i, c in enumerate(countries)}
    table = np.zeros((len(sizes), len(countries)), dtype=float)
    for size, country, cnt in rows:
        table[size_idx[size], country_idx[country]] = cnt

    n_total = int(table.sum())
    chi2, p_value, dof, _ = sp_stats.chi2_contingency(table)
    k = min(len(sizes), len(countries))
    cramers_v = math.sqrt(chi2 / (n_total * (k - 1))) if n_total > 0 and k > 1 else None

    return AnalysisResult(
        analysis_name="chi_square_motif_country",
        grouping="motif_size_x_country",
        metric="ssr_count",
        test_name="chi-square",
        statistic=float(chi2),
        p_value=float(p_value),
        effect_size=cramers_v,
        n=n_total,
        metadata={
            "dof": int(dof),
            "n_motif_sizes": len(sizes),
            "n_countries": len(countries),
        },
    )


# ---------------------------------------------------------------------------
# Pearson / Spearman correlation
# ---------------------------------------------------------------------------

def correlation_length_ssr(
    session: Session,
    run_id: int,
    method: str = "pearson",
) -> AnalysisResult | None:
    """Correlation between genome length and SSR count.

    *method*: 'pearson' or 'spearman'.
    Reports 95% CI for Pearson r via Fisher z-transform.
    """
    stmt = (
        select(Accession.genome_length, AccessionMetrics.ssr_count_total)
        .join(AccessionMetrics, AccessionMetrics.accession == Accession.accession)
        .where(AccessionMetrics.run_id == run_id)
        .where(Accession.genome_length.is_not(None))
    )
    rows = session.execute(stmt).all()
    if len(rows) < 3:
        return None

    x = np.array([float(r[0]) for r in rows])
    y = np.array([float(r[1]) for r in rows])

    if method == "spearman":
        corr, p_value = sp_stats.spearmanr(x, y)
        ci_lo, ci_hi = None, None
    else:
        corr, p_value = sp_stats.pearsonr(x, y)
        # Fisher z-transform CI
        ci_lo, ci_hi = _pearson_ci(float(corr), len(rows))

    return AnalysisResult(
        analysis_name=f"correlation_length_ssr_{method}",
        grouping=None,
        metric="genome_length_vs_ssr_count",
        test_name=method,
        statistic=float(corr),
        p_value=float(p_value),
        effect_size=float(corr),
        ci_lower=ci_lo,
        ci_upper=ci_hi,
        n=len(rows),
    )


def correlation_gc_ssr(
    session: Session,
    run_id: int,
    method: str = "pearson",
) -> AnalysisResult | None:
    """Correlation between GC content and SSR count (or RA)."""
    stmt = (
        select(Accession.gc_content, AccessionMetrics.ra)
        .join(AccessionMetrics, AccessionMetrics.accession == Accession.accession)
        .where(AccessionMetrics.run_id == run_id)
        .where(Accession.gc_content.is_not(None))
        .where(AccessionMetrics.ra.is_not(None))
    )
    rows = session.execute(stmt).all()
    if len(rows) < 3:
        return None

    x = np.array([float(r[0]) for r in rows])
    y = np.array([float(r[1]) for r in rows])

    if method == "spearman":
        corr, p_value = sp_stats.spearmanr(x, y)
        ci_lo, ci_hi = None, None
    else:
        corr, p_value = sp_stats.pearsonr(x, y)
        ci_lo, ci_hi = _pearson_ci(float(corr), len(rows))

    return AnalysisResult(
        analysis_name=f"correlation_gc_ssr_{method}",
        grouping=None,
        metric="gc_content_vs_ra",
        test_name=method,
        statistic=float(corr),
        p_value=float(p_value),
        effect_size=float(corr),
        ci_lower=ci_lo,
        ci_upper=ci_hi,
        n=len(rows),
    )


def _pearson_ci(r: float, n: int, alpha: float = 0.05) -> tuple[float | None, float | None]:
    """95% CI for Pearson r via Fisher z-transform."""
    if n < 4 or abs(r) >= 1.0:
        return None, None
    z = np.arctanh(r)
    se = 1.0 / math.sqrt(n - 3)
    z_crit = sp_stats.norm.ppf(1 - alpha / 2)
    lo = float(np.tanh(z - z_crit * se))
    hi = float(np.tanh(z + z_crit * se))
    return lo, hi


# ---------------------------------------------------------------------------
# Shannon entropy for motif diversity
# ---------------------------------------------------------------------------

def shannon_entropy_motifs(
    session: Session,
    dataset_id: int,
) -> AnalysisResult | None:
    """Shannon entropy of canonical motif distribution across the dataset.

    H = -sum(p_i * log2(p_i)) for each unique motif's proportion.
    Reports max possible entropy as metadata for normalization.
    """
    stmt = (
        select(
            SSRRecord.motif_canonical,
            func.count().label("cnt"),
        )
        .join(Accession, SSRRecord.accession == Accession.accession)
        .where(Accession.dataset_id == dataset_id)
        .group_by(SSRRecord.motif_canonical)
    )
    rows = session.execute(stmt).all()
    if not rows:
        return None

    counts = np.array([float(r[1]) for r in rows])
    total = counts.sum()
    if total == 0:
        return None

    probs = counts / total
    # Filter zero probabilities
    probs = probs[probs > 0]
    h = float(-np.sum(probs * np.log2(probs)))
    h_max = math.log2(len(counts)) if len(counts) > 1 else 0.0
    evenness = h / h_max if h_max > 0 else None

    return AnalysisResult(
        analysis_name="shannon_entropy_motifs",
        grouping=None,
        metric="motif_diversity",
        test_name="shannon-entropy",
        statistic=h,
        p_value=None,  # Entropy is descriptive, no p-value
        effect_size=evenness,  # Pielou's evenness J = H / H_max
        n=int(total),
        metadata={
            "n_unique_motifs": len(counts),
            "h_max": h_max,
            "top_motifs": [
                {"motif": rows[i][0], "count": int(rows[i][1])}
                for i in np.argsort(counts)[::-1][:10]
            ],
        },
    )


# ---------------------------------------------------------------------------
# Shannon entropy per country (geographic diversity)
# ---------------------------------------------------------------------------

def shannon_entropy_by_country(
    session: Session,
    run_id: int,
) -> list[AnalysisResult]:
    """Shannon entropy of motif-size distribution per country.

    Returns one AnalysisResult per country.
    """
    stmt = (
        select(
            Accession.country,
            AccessionMetrics.mono_count,
            AccessionMetrics.di_count,
            AccessionMetrics.tri_count,
            AccessionMetrics.tetra_count,
            AccessionMetrics.penta_count,
            AccessionMetrics.hexa_count,
        )
        .join(AccessionMetrics, AccessionMetrics.accession == Accession.accession)
        .where(AccessionMetrics.run_id == run_id)
        .where(Accession.country.is_not(None))
    )
    rows = session.execute(stmt).all()
    if not rows:
        return []

    # Aggregate motif size counts per country
    country_counts: dict[str, np.ndarray] = {}
    for row in rows:
        country = row[0]
        vals = np.array([float(v) for v in row[1:7]])
        if country in country_counts:
            country_counts[country] += vals
        else:
            country_counts[country] = vals.copy()

    results = []
    for country, counts in sorted(country_counts.items()):
        total = counts.sum()
        if total == 0:
            continue
        probs = counts / total
        probs = probs[probs > 0]
        h = float(-np.sum(probs * np.log2(probs)))
        h_max = math.log2(6)  # 6 motif sizes
        results.append(AnalysisResult(
            analysis_name="shannon_entropy_by_country",
            grouping=country,
            metric="motif_size_diversity",
            test_name="shannon-entropy",
            statistic=h,
            p_value=None,
            effect_size=h / h_max if h_max > 0 else None,
            n=int(total),
        ))

    return results


# ---------------------------------------------------------------------------
# Base composition uniformity (chi-square goodness of fit)
# ---------------------------------------------------------------------------

def base_composition_test(
    session: Session,
    dataset_id: int,
) -> AnalysisResult | None:
    """Chi-square goodness-of-fit test for motif base composition.

    Tests whether A, T, G, C are uniformly represented across all SSR
    motif sequences. Expected proportions under null: 0.25 each.
    """
    stmt = (
        select(SSRRecord.actual_repeat)
        .join(Accession, SSRRecord.accession == Accession.accession)
        .where(Accession.dataset_id == dataset_id)
        .where(SSRRecord.actual_repeat.is_not(None))
    )
    rows = session.execute(stmt).all()
    if not rows:
        return None

    counts = {"A": 0, "T": 0, "G": 0, "C": 0}
    total = 0
    for (repeat,) in rows:
        for base in repeat.upper():
            if base in counts:
                counts[base] += 1
                total += 1

    if total == 0:
        return None

    observed = np.array([counts["A"], counts["T"], counts["G"], counts["C"]], dtype=float)
    expected = np.full(4, total / 4.0)
    chi2, p_value = sp_stats.chisquare(observed, expected)

    # Effect size: Cramér's V for 1×4 table (df=3)
    cramers_v = math.sqrt(chi2 / (total * (4 - 1))) if total > 0 else None

    at_pct = (counts["A"] + counts["T"]) / total * 100 if total > 0 else 0
    gc_pct = (counts["G"] + counts["C"]) / total * 100 if total > 0 else 0

    return AnalysisResult(
        analysis_name="base_composition_uniformity",
        grouping=None,
        metric="base_composition",
        test_name="chi-square-gof",
        statistic=float(chi2),
        p_value=float(p_value),
        effect_size=cramers_v,
        n=total,
        metadata={
            "base_counts": counts,
            "at_percent": round(at_pct, 2),
            "gc_percent": round(gc_pct, 2),
        },
    )


# ---------------------------------------------------------------------------
# Full analysis suite runner
# ---------------------------------------------------------------------------

def run_all_analyses(
    session: Session,
    run_id: int,
    dataset_id: int,
) -> AnalysisSuite:
    """Run all standard analyses and return an FDR-corrected suite.

    This is the main entry point for the analysis layer.
    """
    suite = AnalysisSuite()

    # Kruskal-Wallis (SSR count by country)
    kw_count = kruskal_wallis_by_country(session, run_id, "ssr_count_total")
    if kw_count:
        suite.results.append(kw_count)

    # Kruskal-Wallis (RA by country)
    kw_ra = kruskal_wallis_by_country(session, run_id, "ra")
    if kw_ra:
        kw_ra.analysis_name = "kruskal_wallis_ra_by_country"
        suite.results.append(kw_ra)

    # Chi-square (gene × country)
    chi_gc = chi_square_gene_country(session, dataset_id)
    if chi_gc:
        suite.results.append(chi_gc)

    # Chi-square (motif × country)
    chi_mc = chi_square_motif_country(session, dataset_id)
    if chi_mc:
        suite.results.append(chi_mc)

    # Pearson correlation (length vs SSR)
    corr_len_p = correlation_length_ssr(session, run_id, "pearson")
    if corr_len_p:
        suite.results.append(corr_len_p)

    # Spearman correlation (length vs SSR)
    corr_len_s = correlation_length_ssr(session, run_id, "spearman")
    if corr_len_s:
        suite.results.append(corr_len_s)

    # Pearson correlation (GC vs RA)
    corr_gc_p = correlation_gc_ssr(session, run_id, "pearson")
    if corr_gc_p:
        suite.results.append(corr_gc_p)

    # Spearman correlation (GC vs RA)
    corr_gc_s = correlation_gc_ssr(session, run_id, "spearman")
    if corr_gc_s:
        suite.results.append(corr_gc_s)

    # Shannon entropy (motif diversity)
    entropy = shannon_entropy_motifs(session, dataset_id)
    if entropy:
        suite.results.append(entropy)

    # Base composition uniformity
    base_comp = base_composition_test(session, dataset_id)
    if base_comp:
        suite.results.append(base_comp)

    # Per-country entropy
    country_entropy = shannon_entropy_by_country(session, run_id)
    suite.results.extend(country_entropy)

    # Apply BH-FDR correction across all results with p-values
    suite.apply_fdr()

    return suite
