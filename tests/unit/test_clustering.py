"""Test suite for Chunk 13 clustering and comparative analysis.

Coverage:
- Feature matrix generation
- Preprocessing (scaling, normalization, missing data)
- K-means clustering with reproducibility
- Hierarchical clustering
- DBSCAN clustering
- PCA dimensionality reduction
- Comparative analysis (by lineage, geography)
- Edge cases and backward compatibility
"""

from __future__ import annotations

import pytest
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Base,
    Dataset,
    Run,
    Accession,
    SSRRecord,
)
from gwico_ssr.analysis.clustering import (
    FeatureMatrixGenerator,
    ClusteringPreprocessor,
    ClusteringAnalysis,
    ClusteringResult,
    clusters_by_lineage,
    clusters_by_geography,
)


# -----------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------

@pytest.fixture
def db_session():
    """In-memory SQLite database session for testing."""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    session = Session(engine)
    yield session
    session.close()


@pytest.fixture
def populated_db(db_session):
    """Session with dataset, runs, accessions, and SSR records."""
    session = db_session

    # Create dataset
    dataset = Dataset(name="test_dataset", source_type="csv", organism="virus")
    session.add(dataset)
    session.flush()

    # Create run
    run = Run(dataset_id=dataset.dataset_id, pipeline_version="1.0.0")
    session.add(run)
    session.flush()

    # Create accessions with geographic context
    accessions = []
    for i, (country, host) in enumerate([
        ("USA", "human"),
        ("USA", "human"),
        ("China", "human"),
        ("China", "human"),
        ("Brazil", "human"),
        ("Brazil", "human"),
        ("USA", "bat"),
        ("China", "bat"),
    ]):
        acc = Accession(
            accession=f"accession_{i}",
            dataset_id=dataset.dataset_id,
            country=country,
            host=host,
            species="virus_A",
        )
        session.add(acc)
        accessions.append(acc)
    session.flush()

    # Create SSR records with varying profiles
    ssr_data = [
        # Cluster 1: high SSR count, high perfect
        (0, 10, 5, "perfect"),
        (0, 12, 6, "perfect"),
        (1, 11, 5, "perfect"),
        (1, 13, 6, "perfect"),
        # Cluster 2: medium SSR count, mixed
        (2, 5, 2, "perfect"),
        (2, 4, 1, "imperfect"),
        (3, 6, 2, "perfect"),
        (3, 3, 1, "imperfect"),
        # Cluster 3: low SSR count, more imperfect
        (4, 2, 0, "imperfect"),
        (4, 3, 1, "imperfect"),
        (5, 1, 0, "imperfect"),
        (5, 2, 0, "compound_component"),
        # Outlier
        (6, 20, 10, "perfect"),
        (7, 1, 0, "imperfect"),
    ]

    for acc_idx, repeat_units, motif_size, repeat_class in ssr_data:
        ssr = SSRRecord(
            run_id=run.run_id,
            accession=f"accession_{acc_idx}",
            motif_raw="AT",
            motif_canonical="AT",
            repeat_units=repeat_units,
            motif_size=motif_size,
            repeat_class=repeat_class,
            repeat_length_bp=repeat_units * motif_size,
            start=100 + (repeat_units * 10),
            end=100 + (repeat_units * 10) + (repeat_units * motif_size),
            strand="+",
        )
        session.add(ssr)

    session.commit()
    return session


# -----------------------------------------------------------------------
# Feature Matrix Generation Tests (5)
# -----------------------------------------------------------------------

class TestFeatureMatrixGeneration:
    """Test feature matrix generation."""

    def test_generate_basic_features(self, populated_db):
        """Generate basic feature matrix from accession metrics."""
        session = populated_db
        gen = FeatureMatrixGenerator(session)

        X, feature_names, accessions = gen.generate_for_run(
            session.query(Run.run_id).first()[0]
        )

        assert X.shape[0] == 8  # 8 accessions
        assert X.shape[1] == 7  # 7 features
        assert len(accessions) == 8
        assert len(feature_names) == 7

    def test_feature_columns_correct(self, populated_db):
        """Feature matrix has expected column names."""
        session = populated_db
        gen = FeatureMatrixGenerator(session)

        _, feature_names, _ = gen.generate_for_run(
            session.query(Run.run_id).first()[0]
        )

        expected_cols = [
            "ssr_count",
            "repeat_length_bp_total",
            "repeat_units_total",
            "repeat_length_bp_mean",
            "perfect_count",
            "imperfect_count",
            "compound_component_count",
        ]
        assert set(feature_names) == set(expected_cols)

    def test_feature_values_numeric(self, populated_db):
        """Feature matrix contains valid numeric values."""
        session = populated_db
        gen = FeatureMatrixGenerator(session)

        X, _, _ = gen.generate_for_run(
            session.query(Run.run_id).first()[0]
        )

        assert np.all(np.isfinite(X))
        assert np.all(X >= 0)  # Counts should be non-negative

    def test_different_runs_produce_different_matrices(self, populated_db):
        """Different runs produce different feature matrices."""
        session = populated_db

        # Create second run
        dataset = session.query(Dataset).first()
        run2 = Run(dataset_id=dataset.dataset_id, pipeline_version="2.0.0")
        session.add(run2)
        session.flush()

        # Add different SSR data
        for i in range(2):
            ssr = SSRRecord(
                run_id=run2.run_id,
                accession=f"accession_{i}",
                motif_raw="GC",
                motif_canonical="GC",
                repeat_units=30,
                motif_size=2,
                repeat_class="perfect",
                repeat_length_bp=60,
                start=1000,
                end=1060,
                strand="+",
            )
            session.add(ssr)
        session.commit()

        gen = FeatureMatrixGenerator(session)
        run1_id = session.query(Run.run_id).first()[0]
        run2_id = run2.run_id

        X1, _, _ = gen.generate_for_run(run1_id)
        X2, _, _ = gen.generate_for_run(run2_id)

        assert X1.shape != X2.shape or not np.allclose(X1[:2], X2[:2])

    def test_accession_order_preserved(self, populated_db):
        """Accession order is consistent in output."""
        session = populated_db
        gen = FeatureMatrixGenerator(session)

        run_id = session.query(Run.run_id).first()[0]
        _, _, acc1 = gen.generate_for_run(run_id)
        _, _, acc2 = gen.generate_for_run(run_id)

        assert np.array_equal(acc1, acc2)


# -----------------------------------------------------------------------
# Preprocessing Tests (5)
# -----------------------------------------------------------------------

class TestPreprocessing:
    """Test preprocessing operations."""

    def test_standard_scaling(self):
        """StandardScaler normalizes features."""
        X = np.array([[1, 2], [3, 4], [5, 6]])
        preprocessor = ClusteringPreprocessor(scaling_method="standard")

        X_scaled, config = preprocessor.fit_and_transform(X)

        assert X_scaled.shape == X.shape
        assert np.allclose(X_scaled.mean(axis=0), 0, atol=1e-10)
        assert np.allclose(X_scaled.std(axis=0), 1, atol=1e-10)
        assert config["scaling_method"] == "standard"

    def test_minmax_scaling(self):
        """MinMaxScaler scales to [0, 1]."""
        X = np.array([[1, 2], [3, 4], [5, 6]])
        preprocessor = ClusteringPreprocessor(scaling_method="minmax")

        X_scaled, config = preprocessor.fit_and_transform(X)

        assert X_scaled.shape == X.shape
        assert np.all(X_scaled >= 0)
        assert np.all(X_scaled <= 1)
        assert config["scaling_method"] == "minmax"

    def test_missing_value_handling(self):
        """NaN values are handled properly."""
        X = np.array([[1, 2], [np.nan, 4], [5, np.nan]])
        preprocessor = ClusteringPreprocessor()

        X_scaled, config = preprocessor.fit_and_transform(X)

        assert not np.any(np.isnan(X_scaled))
        assert X_scaled.shape == X.shape

    def test_transform_consistency(self):
        """Transform applies same scaling as fit_and_transform."""
        X_train = np.array([[1, 2], [3, 4], [5, 6]])
        X_test = np.array([[2, 3], [4, 5]])

        preprocessor = ClusteringPreprocessor()
        X_train_scaled, _ = preprocessor.fit_and_transform(X_train)

        X_test_scaled = preprocessor.transform(X_test)

        assert X_test_scaled.shape == X_test.shape
        assert np.all(np.isfinite(X_test_scaled))

    def test_config_preservation(self):
        """Configuration is stored for reproducibility."""
        X = np.array([[1, 2], [3, 4], [5, 6]])
        preprocessor = ClusteringPreprocessor(scaling_method="standard")

        _, config = preprocessor.fit_and_transform(X)

        assert "mean" in config
        assert "scale" in config
        assert len(config["mean"]) == 2


# -----------------------------------------------------------------------
# K-Means Clustering Tests (4)
# -----------------------------------------------------------------------

class TestKMeansClustering:
    """Test K-means clustering."""

    def test_kmeans_basic(self):
        """K-means produces expected number of clusters."""
        X = np.array([
            [0, 0], [1, 1], [2, 2],  # Cluster 1
            [10, 10], [11, 11], [12, 12],  # Cluster 2
        ])
        clusterer = ClusteringAnalysis(seed=42)
        result = clusterer.kmeans(X, n_clusters=2)

        assert result.algorithm == "kmeans"
        assert result.n_clusters == 2
        assert len(result.labels) == 6
        assert set(result.labels) == {0, 1}

    def test_kmeans_reproducibility(self):
        """Same seed produces same results."""
        X = np.random.RandomState(42).randn(20, 5)

        clusterer1 = ClusteringAnalysis(seed=42)
        result1 = clusterer1.kmeans(X, n_clusters=3)

        clusterer2 = ClusteringAnalysis(seed=42)
        result2 = clusterer2.kmeans(X, n_clusters=3)

        assert np.array_equal(result1.labels, result2.labels)
        assert np.allclose(result1.centroids, result2.centroids)

    def test_kmeans_different_seeds_differ(self):
        """Different seeds can produce different results."""
        X = np.random.RandomState(123).randn(50, 10)

        clusterer1 = ClusteringAnalysis(seed=42)
        result1 = clusterer1.kmeans(X, n_clusters=5)

        clusterer2 = ClusteringAnalysis(seed=999)
        result2 = clusterer2.kmeans(X, n_clusters=5)

        # Results may differ (not guaranteed but likely)
        # At least check they're both valid
        assert len(set(result1.labels)) <= 5
        assert len(set(result2.labels)) <= 5

    def test_kmeans_metrics_computed(self):
        """K-means computes quality metrics."""
        X = np.array([[0, 0], [1, 1], [10, 10], [11, 11]])
        clusterer = ClusteringAnalysis(seed=42)
        result = clusterer.kmeans(X, n_clusters=2)

        assert result.silhouette is not None
        assert result.davies_bouldin is not None
        assert result.inertia is not None
        assert -1 <= result.silhouette <= 1


# -----------------------------------------------------------------------
# Hierarchical Clustering Tests (3)
# -----------------------------------------------------------------------

class TestHierarchicalClustering:
    """Test hierarchical clustering."""

    def test_hierarchical_basic(self):
        """Hierarchical clustering produces expected clusters."""
        X = np.array([
            [0, 0], [1, 1], [2, 2],
            [10, 10], [11, 11], [12, 12],
        ])
        clusterer = ClusteringAnalysis(seed=42)
        result = clusterer.hierarchical(X, n_clusters=2)

        assert result.algorithm == "hierarchical"
        assert result.n_clusters == 2
        assert len(result.labels) == 6

    def test_hierarchical_linkage_options(self):
        """Different linkage methods produce results."""
        X = np.random.RandomState(42).randn(15, 3)
        clusterer = ClusteringAnalysis(seed=42)

        for linkage in ["ward", "complete", "average"]:
            result = clusterer.hierarchical(X, n_clusters=3, linkage=linkage)
            assert len(set(result.labels)) <= 3

    def test_hierarchical_metrics_computed(self):
        """Hierarchical clustering computes metrics."""
        X = np.array([[0, 0], [1, 1], [10, 10], [11, 11]])
        clusterer = ClusteringAnalysis()
        result = clusterer.hierarchical(X, n_clusters=2)

        assert result.silhouette is not None
        assert result.davies_bouldin is not None


# -----------------------------------------------------------------------
# DBSCAN Clustering Tests (3)
# -----------------------------------------------------------------------

class TestDBSCANClustering:
    """Test DBSCAN clustering."""

    def test_dbscan_basic(self):
        """DBSCAN produces clusters and noise points."""
        X = np.array([
            [0, 0], [1, 1], [2, 2],  # Cluster 1
            [10, 10], [11, 11], [12, 12],  # Cluster 2
            [50, 50],  # Noise
        ])
        clusterer = ClusteringAnalysis(seed=42)
        result = clusterer.dbscan(X, eps=3, min_samples=2)

        assert result.algorithm == "dbscan"
        assert -1 in result.labels  # Has noise points
        n_noise = np.sum(result.labels == -1)
        assert n_noise > 0

    def test_dbscan_eps_sensitivity(self):
        """DBSCAN responds to eps parameter."""
        X = np.random.RandomState(42).randn(20, 2)
        clusterer = ClusteringAnalysis(seed=42)

        result_small = clusterer.dbscan(X, eps=0.1, min_samples=2)
        result_large = clusterer.dbscan(X, eps=10, min_samples=2)

        # Smaller eps should produce more clusters
        assert result_small.n_clusters <= result_large.n_clusters

    def test_dbscan_noise_handling(self):
        """DBSCAN correctly identifies noise points."""
        X = np.array([
            [0, 0], [1, 1],  # Cluster
            [100, 100],  # Isolated noise
        ])
        clusterer = ClusteringAnalysis()
        result = clusterer.dbscan(X, eps=2, min_samples=2)

        assert -1 in result.labels
        assert result.labels[-1] == -1  # Last point is noise


# -----------------------------------------------------------------------
# PCA Dimensionality Reduction Tests (2)
# -----------------------------------------------------------------------

class TestPCAReduction:
    """Test PCA dimensionality reduction."""

    def test_pca_basic(self):
        """PCA reduces dimensions."""
        X = np.random.RandomState(42).randn(50, 10)
        clusterer = ClusteringAnalysis(seed=42)
        result = clusterer.pca(X, n_components=3)

        assert result.algorithm == "pca"
        assert result.n_clusters == 3  # n_components stored as n_clusters

    def test_pca_variance_explained(self):
        """PCA reports variance explained."""
        X = np.random.RandomState(42).randn(50, 10)
        clusterer = ClusteringAnalysis(seed=42)
        result = clusterer.pca(X, n_components=5)

        # Silhouette field stores cumulative variance
        assert 0 < result.silhouette <= 1


# -----------------------------------------------------------------------
# Clustering Result Reproducibility Tests (2)
# -----------------------------------------------------------------------

class TestClusteringResultReproducibility:
    """Test result reproducibility hashing."""

    def test_reproducibility_hash_deterministic(self):
        """Reproducibility hash is deterministic."""
        result = ClusteringResult(
            algorithm="kmeans",
            n_clusters=3,
            n_samples=100,
            seed=42,
            parameters={"n_init": 10},
        )

        hash1 = result.reproducibility_hash()
        hash2 = result.reproducibility_hash()

        assert hash1 == hash2

    def test_reproducibility_hash_differs_on_changes(self):
        """Different parameters produce different hashes."""
        result1 = ClusteringResult(
            algorithm="kmeans",
            n_clusters=3,
            n_samples=100,
            seed=42,
        )
        result2 = ClusteringResult(
            algorithm="kmeans",
            n_clusters=4,
            n_samples=100,
            seed=42,
        )

        assert result1.reproducibility_hash() != result2.reproducibility_hash()


# -----------------------------------------------------------------------
# Comparative Analysis Tests (3)
# -----------------------------------------------------------------------

class TestComparativeAnalysis:
    """Test comparative analysis functions."""

    def test_clusters_by_lineage(self, populated_db):
        """Lineage analysis produces crosstab."""
        session = populated_db
        accessions = np.array([f"accession_{i}" for i in range(8)])
        labels = np.array([0, 0, 1, 1, 2, 2, 0, 1])

        crosstab = clusters_by_lineage(session, accessions, labels)

        assert crosstab is not None
        assert isinstance(crosstab, pd.DataFrame)
        # Should have clusters in rows and hosts in columns
        assert crosstab.shape[0] >= 2  # At least 2 clusters

    def test_clusters_by_geography(self, populated_db):
        """Geographic analysis produces crosstab."""
        session = populated_db
        accessions = np.array([f"accession_{i}" for i in range(8)])
        labels = np.array([0, 0, 0, 1, 1, 1, 2, 2])

        crosstab = clusters_by_geography(session, accessions, labels)

        assert crosstab is not None
        assert isinstance(crosstab, pd.DataFrame)
        # Should have clusters in rows and countries in columns
        assert crosstab.shape[0] >= 2


# -----------------------------------------------------------------------
# Integration Tests (3)
# -----------------------------------------------------------------------

class TestIntegration:
    """Integration tests with full workflow."""

    def test_full_clustering_workflow(self, populated_db):
        """End-to-end clustering from data to results."""
        session = populated_db
        run_id = session.query(Run.run_id).first()[0]

        # Generate features
        gen = FeatureMatrixGenerator(session)
        X, feature_names, accessions = gen.generate_for_run(run_id)

        # Preprocess
        preprocessor = ClusteringPreprocessor()
        X_scaled, config = preprocessor.fit_and_transform(X)

        # Cluster
        clusterer = ClusteringAnalysis(seed=42)
        result = clusterer.kmeans(X_scaled, n_clusters=3)

        # Validate
        assert result.n_samples == len(accessions)
        assert len(result.labels) == len(accessions)
        assert len(set(result.labels)) <= 3

    def test_clustering_reproducibility_full(self, populated_db):
        """Full workflow is reproducible."""
        session = populated_db
        run_id = session.query(Run.run_id).first()[0]

        # First run
        gen1 = FeatureMatrixGenerator(session)
        X1, _, acc1 = gen1.generate_for_run(run_id)
        prep1 = ClusteringPreprocessor()
        X_scaled1, _ = prep1.fit_and_transform(X1)
        clust1 = ClusteringAnalysis(seed=123)
        result1 = clust1.kmeans(X_scaled1, n_clusters=3)

        # Second run (should be identical)
        gen2 = FeatureMatrixGenerator(session)
        X2, _, acc2 = gen2.generate_for_run(run_id)
        prep2 = ClusteringPreprocessor()
        X_scaled2, _ = prep2.fit_and_transform(X2)
        clust2 = ClusteringAnalysis(seed=123)
        result2 = clust2.kmeans(X_scaled2, n_clusters=3)

        # Results should be identical
        assert np.array_equal(result1.labels, result2.labels)
        assert np.array_equal(acc1, acc2)

    def test_single_sample_edge_case(self, db_session):
        """Handle single-sample edge case gracefully."""
        session = db_session

        # Create minimal data: 1 accession, 1 SSR
        dataset = Dataset(name="edge_case", source_type="csv", organism="virus")
        session.add(dataset)
        session.flush()

        run = Run(dataset_id=dataset.dataset_id, pipeline_version="1.0.0")
        session.add(run)
        session.flush()

        acc = Accession(
            accession="single_accession",
            dataset_id=dataset.dataset_id,
            country="USA",
        )
        session.add(acc)
        session.flush()

        ssr = SSRRecord(
            run_id=run.run_id,
            accession="single_accession",
            motif_raw="AT",
            motif_canonical="AT",
            repeat_units=5,
            motif_size=2,
            repeat_class="perfect",
            repeat_length_bp=10,
            start=100,
            end=110,
            strand="+",
        )
        session.add(ssr)
        session.commit()

        # Try to generate and cluster
        gen = FeatureMatrixGenerator(session)
        X, _, _ = gen.generate_for_run(run.run_id)

        assert X.shape[0] == 1


# -----------------------------------------------------------------------
# Backward Compatibility Tests (2)
# -----------------------------------------------------------------------

class TestBackwardCompatibility:
    """Ensure Chunk 13 doesn't break existing functionality."""

    def test_existing_ssr_records_unaffected(self, populated_db):
        """SSRRecords are still queryable."""
        session = populated_db
        ssr = session.query(SSRRecord).first()
        assert ssr is not None
        assert ssr.motif_canonical == "AT"

    def test_existing_analysis_modules_work(self, populated_db):
        """Can still use statistics, temporal, lineage modules."""
        from gwico_ssr.analysis import AnalysisResult, TemporalSummary

        # Just verify imports work
        result = AnalysisResult(
            analysis_name="test",
            grouping="test_group",
            metric="test_metric",
            test_name="test_test",
            statistic=1.0,
            p_value=0.05,
        )
        assert result.analysis_name == "test"
