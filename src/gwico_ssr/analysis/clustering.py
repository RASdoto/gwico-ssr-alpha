"""Reproducible clustering and comparative analysis for SSR-derived features.

Provides:
- Feature matrix generation from accession metrics with temporal/lineage context
- Preprocessing (normalization, scaling, missing data handling)
- Multiple clustering algorithms (K-means, hierarchical, DBSCAN, PCA)
- Reproducible seeding and parameter storage
- Comparative analysis by lineage, geography, temporal bins
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from hashlib import sha256
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from scipy.spatial.distance import pdist, squareform
from sqlalchemy import select, func
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    Run,
    SSRRecord,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class ClusteringResult:
    """Single clustering run result."""

    algorithm: str  # 'kmeans', 'hierarchical', 'dbscan', 'pca'
    n_clusters: int
    n_samples: int
    labels: np.ndarray | None = None
    centroids: np.ndarray | None = None
    silhouette: float | None = None
    davies_bouldin: float | None = None
    inertia: float | None = None
    seed: int = 42
    parameters: dict = field(default_factory=dict)
    preprocessor_config: dict = field(default_factory=dict)
    feature_names: list[str] = field(default_factory=list)
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def reproducibility_hash(self) -> str:
        """Generate reproducibility hash for validation."""
        key_data = {
            "algorithm": self.algorithm,
            "n_clusters": self.n_clusters,
            "seed": self.seed,
            "parameters": json.dumps(self.parameters, sort_keys=True),
            "preprocessor_config": json.dumps(self.preprocessor_config, sort_keys=True),
        }
        key_str = json.dumps(key_data, sort_keys=True)
        return sha256(key_str.encode()).hexdigest()

    def to_dict(self) -> dict:
        """Convert to dictionary for storage."""
        return {
            "algorithm": self.algorithm,
            "n_clusters": self.n_clusters,
            "n_samples": self.n_samples,
            "silhouette_score": self.silhouette,
            "davies_bouldin_score": self.davies_bouldin,
            "inertia": self.inertia,
            "seed": self.seed,
            "parameters_json": json.dumps(self.parameters),
            "preprocessor_config_json": json.dumps(self.preprocessor_config),
            "feature_names_json": json.dumps(self.feature_names),
            "reproducibility_hash": self.reproducibility_hash(),
            "timestamp": self.timestamp,
        }


@dataclass
class ClusterAssignment:
    """Single sample cluster assignment."""

    accession: str
    cluster_label: int
    distance_to_centroid: float | None = None


# ---------------------------------------------------------------------------
# Feature Matrix Generation
# ---------------------------------------------------------------------------

class FeatureMatrixGenerator:
    """Generate feature matrices from SSR data with optional context."""

    def __init__(self, session: Session):
        """Initialize with database session.

        Args:
            session: SQLAlchemy session for queries
        """
        self.session = session

    def _get_accession_metrics(self, run_id: int | None = None) -> pd.DataFrame:
        """Query accession metrics with optional run filter.

        Returns:
            DataFrame with columns: accession, ssr_count_total, repeat_length_total, etc.
        """
        stmt = select(
            Accession.accession,
            func.count(SSRRecord.ssr_id).label("ssr_count"),
            func.sum(SSRRecord.repeat_length_bp).label("repeat_length_bp_total"),
            func.sum(SSRRecord.repeat_units).label("repeat_units_total"),
            func.avg(SSRRecord.repeat_length_bp).label("repeat_length_bp_mean"),
        ).join(SSRRecord, Accession.accession == SSRRecord.accession)

        if run_id is not None:
            stmt = stmt.where(SSRRecord.run_id == run_id)

        stmt = stmt.group_by(Accession.accession)

        results = self.session.execute(stmt).fetchall()
        df = pd.DataFrame(
            results,
            columns=[
                "accession",
                "ssr_count",
                "repeat_length_bp_total",
                "repeat_units_total",
                "repeat_length_bp_mean",
            ],
        )
        return df

    def _add_repeat_class_composition(self, df: pd.DataFrame, run_id: int | None = None) -> pd.DataFrame:
        """Add repeat class composition columns.

        Returns DataFrame with added columns for perfect/imperfect/compound counts.
        """
        stmt = select(
            SSRRecord.accession,
            SSRRecord.repeat_class,
            func.count(SSRRecord.ssr_id).label("count"),
        ).group_by(SSRRecord.accession, SSRRecord.repeat_class)

        if run_id is not None:
            stmt = stmt.where(SSRRecord.run_id == run_id)

        results = self.session.execute(stmt).fetchall()
        class_df = pd.DataFrame(results, columns=["accession", "repeat_class", "count"])

        # Pivot to get columns per repeat class
        pivot_df = class_df.pivot(index="accession", columns="repeat_class", values="count").fillna(0)

        # Ensure all classes present
        for cls in ["perfect", "imperfect", "compound_component"]:
            if cls not in pivot_df.columns:
                pivot_df[cls] = 0

        # Rename
        pivot_df = pivot_df.rename(
            columns={
                "perfect": "perfect_count",
                "imperfect": "imperfect_count",
                "compound_component": "compound_component_count",
            }
        )

        # Merge with main df
        df = df.merge(pivot_df, left_on="accession", right_index=True, how="left")
        return df

    def _add_geographic_context(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add geographic (country) context."""
        stmt = select(Accession.accession, Accession.country).distinct()
        results = self.session.execute(stmt).fetchall()
        geo_df = pd.DataFrame(results, columns=["accession", "country"])
        df = df.merge(geo_df, on="accession", how="left")
        return df

    def generate_from_accession_metrics(
        self,
        run_id: int | None = None,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Generate feature matrix from SSR accession metrics.

        Args:
            run_id: Optional filter to specific run

        Returns:
            (feature_matrix, feature_names)
        """
        logger.info("Generating feature matrix from accession metrics")

        # Start with basic metrics
        df = self._get_accession_metrics(run_id=run_id)
        logger.info(f"Retrieved metrics for {len(df)} accessions")

        # Add repeat class composition
        df = self._add_repeat_class_composition(df, run_id=run_id)

        # Add geographic context
        df = self._add_geographic_context(df)

        # Select numeric features for clustering
        feature_columns = [
            "ssr_count",
            "repeat_length_bp_total",
            "repeat_units_total",
            "repeat_length_bp_mean",
            "perfect_count",
            "imperfect_count",
            "compound_component_count",
        ]

        # Ensure all feature columns exist
        for col in feature_columns:
            if col not in df.columns:
                df[col] = 0

        # Extract feature matrix
        X = df[feature_columns].fillna(0).values
        accessions = df["accession"].values

        logger.info(f"Generated feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
        return X, feature_columns, accessions, df[["accession", "country"]]

    def generate_for_run(
        self,
        run_id: int,
    ) -> tuple[pd.DataFrame, list[str], np.ndarray]:
        """Generate feature matrix for specific run.

        Args:
            run_id: Run ID to generate features for

        Returns:
            (X_matrix, feature_names, accessions)
        """
        X, feature_names, accessions, context_df = self.generate_from_accession_metrics(run_id=run_id)
        return X, feature_names, accessions


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

class ClusteringPreprocessor:
    """Preprocessing for clustering with reproducible scaling."""

    def __init__(self, scaling_method: str = "standard"):
        """Initialize preprocessor.

        Args:
            scaling_method: 'standard' (StandardScaler) or 'minmax' (MinMaxScaler)
        """
        self.scaling_method = scaling_method
        self.scaler = None
        self.feature_means = None
        self.feature_stds = None

    def fit_and_transform(self, X: np.ndarray) -> tuple[np.ndarray, dict]:
        """Fit scaler on data and transform.

        Args:
            X: Feature matrix (n_samples × n_features)

        Returns:
            (transformed_X, config_dict)
        """
        # Handle missing values
        X_clean = np.where(np.isnan(X), np.nanmean(X, axis=0), X)

        if self.scaling_method == "standard":
            self.scaler = StandardScaler()
        elif self.scaling_method == "minmax":
            self.scaler = MinMaxScaler()
        else:
            raise ValueError(f"Unknown scaling method: {self.scaling_method}")

        X_scaled = self.scaler.fit_transform(X_clean)

        config = {
            "scaling_method": self.scaling_method,
            "n_features": X.shape[1],
            "n_samples": X.shape[0],
        }

        if self.scaling_method == "standard":
            config["mean"] = self.scaler.mean_.tolist()
            config["scale"] = self.scaler.scale_.tolist()
        elif self.scaling_method == "minmax":
            config["feature_min"] = self.scaler.data_min_.tolist()
            config["feature_max"] = self.scaler.data_max_.tolist()

        return X_scaled, config

    def transform(self, X: np.ndarray) -> np.ndarray:
        """Transform data using fitted scaler.

        Args:
            X: Feature matrix

        Returns:
            Scaled feature matrix
        """
        if self.scaler is None:
            raise ValueError("Preprocessor not fitted. Call fit_and_transform first.")
        X_clean = np.where(np.isnan(X), np.nanmean(X, axis=0), X)
        return self.scaler.transform(X_clean)


# ---------------------------------------------------------------------------
# Clustering Algorithms
# ---------------------------------------------------------------------------

class ClusteringAnalysis:
    """Multiple clustering algorithms with reproducible seeding."""

    def __init__(self, seed: int = 42):
        """Initialize clustering analyzer.

        Args:
            seed: Random seed for reproducibility
        """
        self.seed = seed
        np.random.seed(seed)

    def kmeans(
        self,
        X: np.ndarray,
        n_clusters: int,
        n_init: int = 10,
    ) -> ClusteringResult:
        """K-means clustering.

        Args:
            X: Feature matrix (n_samples × n_features)
            n_clusters: Number of clusters
            n_init: Number of initializations

        Returns:
            ClusteringResult with cluster labels and metrics
        """
        logger.info(f"Running K-means with k={n_clusters}, seed={self.seed}")

        kmeans = KMeans(
            n_clusters=n_clusters,
            random_state=self.seed,
            n_init=n_init,
        )
        labels = kmeans.fit_predict(X)

        silhouette = silhouette_score(X, labels)
        davies_bouldin = davies_bouldin_score(X, labels)

        logger.info(
            f"K-means complete: silhouette={silhouette:.3f}, "
            f"davies_bouldin={davies_bouldin:.3f}, inertia={kmeans.inertia_:.2f}"
        )

        result = ClusteringResult(
            algorithm="kmeans",
            n_clusters=n_clusters,
            n_samples=X.shape[0],
            labels=labels,
            centroids=kmeans.cluster_centers_,
            silhouette=silhouette,
            davies_bouldin=davies_bouldin,
            inertia=kmeans.inertia_,
            seed=self.seed,
            parameters={
                "n_clusters": n_clusters,
                "n_init": n_init,
            },
        )
        return result

    def hierarchical(
        self,
        X: np.ndarray,
        n_clusters: int,
        linkage: str = "ward",
    ) -> ClusteringResult:
        """Hierarchical agglomerative clustering.

        Args:
            X: Feature matrix
            n_clusters: Number of clusters
            linkage: 'ward', 'complete', 'average', 'single'

        Returns:
            ClusteringResult with cluster labels and metrics
        """
        logger.info(f"Running hierarchical clustering with k={n_clusters}, linkage={linkage}")

        hc = AgglomerativeClustering(
            n_clusters=n_clusters,
            linkage=linkage,
        )
        labels = hc.fit_predict(X)

        silhouette = silhouette_score(X, labels)
        davies_bouldin = davies_bouldin_score(X, labels)

        logger.info(
            f"Hierarchical complete: silhouette={silhouette:.3f}, "
            f"davies_bouldin={davies_bouldin:.3f}"
        )

        result = ClusteringResult(
            algorithm="hierarchical",
            n_clusters=n_clusters,
            n_samples=X.shape[0],
            labels=labels,
            silhouette=silhouette,
            davies_bouldin=davies_bouldin,
            seed=self.seed,
            parameters={
                "n_clusters": n_clusters,
                "linkage": linkage,
            },
        )
        return result

    def dbscan(
        self,
        X: np.ndarray,
        eps: float = 0.5,
        min_samples: int = 5,
    ) -> ClusteringResult:
        """DBSCAN density-based clustering.

        Args:
            X: Feature matrix
            eps: Neighborhood radius
            min_samples: Minimum samples in neighborhood

        Returns:
            ClusteringResult with cluster labels and metrics
        """
        logger.info(f"Running DBSCAN with eps={eps}, min_samples={min_samples}")

        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        labels = dbscan.fit_predict(X)

        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = list(labels).count(-1)

        logger.info(f"DBSCAN complete: {n_clusters} clusters, {n_noise} noise points")

        # Compute silhouette only if we have clusters
        if n_clusters > 1 and n_noise < len(labels):
            # Filter out noise points for silhouette calculation
            mask = labels != -1
            if mask.sum() > 0:
                silhouette = silhouette_score(X[mask], labels[mask])
            else:
                silhouette = None
        else:
            silhouette = None

        result = ClusteringResult(
            algorithm="dbscan",
            n_clusters=n_clusters,
            n_samples=X.shape[0],
            labels=labels,
            silhouette=silhouette,
            seed=self.seed,
            parameters={
                "eps": eps,
                "min_samples": min_samples,
                "n_noise_points": n_noise,
            },
        )
        return result

    def pca(
        self,
        X: np.ndarray,
        n_components: int | float = 0.95,
    ) -> ClusteringResult:
        """PCA dimensionality reduction.

        Args:
            X: Feature matrix
            n_components: Number of components or variance fraction

        Returns:
            ClusteringResult with transformed data
        """
        logger.info(f"Running PCA with n_components={n_components}")

        pca = PCA(n_components=n_components, random_state=self.seed)
        X_transformed = pca.fit_transform(X)

        logger.info(
            f"PCA complete: {pca.n_components_} components, "
            f"explained_variance={pca.explained_variance_ratio_.sum():.3f}"
        )

        result = ClusteringResult(
            algorithm="pca",
            n_clusters=pca.n_components_,
            n_samples=X.shape[0],
            labels=None,  # PCA doesn't assign labels
            silhouette=pca.explained_variance_ratio_.sum(),
            seed=self.seed,
            parameters={
                "n_components_request": str(n_components),
                "n_components_actual": pca.n_components_,
                "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
            },
        )
        return result


# ---------------------------------------------------------------------------
# Comparative Analysis
# ---------------------------------------------------------------------------

def clusters_by_lineage(
    session: Session,
    accessions: np.ndarray,
    labels: np.ndarray,
) -> pd.DataFrame:
    """Analyze cluster composition by lineage.

    Args:
        session: SQLAlchemy session
        accessions: Array of accession names
        labels: Array of cluster labels

    Returns:
        DataFrame with cluster × lineage composition
    """
    # Build mapping
    data = []
    for acc, label in zip(accessions, labels):
        data.append({"accession": acc, "cluster": label})

    df = pd.DataFrame(data)

    # Query lineage info
    stmt = select(Accession.accession, Accession.host).distinct()
    results = session.execute(stmt).fetchall()
    lineage_df = pd.DataFrame(results, columns=["accession", "lineage"])

    df = df.merge(lineage_df, on="accession", how="left")

    # Create crosstab
    crosstab = pd.crosstab(df["cluster"], df["lineage"], margins=True)
    logger.info(f"Cluster × Lineage:\n{crosstab}")

    return crosstab


def clusters_by_geography(
    session: Session,
    accessions: np.ndarray,
    labels: np.ndarray,
) -> pd.DataFrame:
    """Analyze cluster composition by geography.

    Args:
        session: SQLAlchemy session
        accessions: Array of accession names
        labels: Array of cluster labels

    Returns:
        DataFrame with cluster × country composition
    """
    data = []
    for acc, label in zip(accessions, labels):
        data.append({"accession": acc, "cluster": label})

    df = pd.DataFrame(data)

    # Query geography info
    stmt = select(Accession.accession, Accession.country).distinct()
    results = session.execute(stmt).fetchall()
    geo_df = pd.DataFrame(results, columns=["accession", "country"])

    df = df.merge(geo_df, on="accession", how="left")

    # Create crosstab
    crosstab = pd.crosstab(df["cluster"], df["country"], margins=True)
    logger.info(f"Cluster × Geography:\n{crosstab}")

    return crosstab
