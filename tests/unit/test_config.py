"""Tests for configuration loading."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from gwico_ssr.config import Settings, load_settings


def test_default_settings():
    """Settings have sensible defaults without any config file."""
    settings = load_settings("/nonexistent/path.toml")
    assert settings.database.url == "sqlite:///gwico_ssr.db"
    assert settings.ssr.min_repeats_mono == 10
    assert settings.ssr.min_repeats_tri == 3
    assert settings.logging.level == "INFO"


def test_load_from_toml(sample_toml: Path):
    """Settings load from a TOML file."""
    settings = load_settings(str(sample_toml))
    assert settings.database.url == "sqlite:///test.db"
    assert settings.database.echo is True
    assert settings.ncbi.email == "test@example.com"
    assert settings.ncbi.max_retries == 5
    assert settings.ssr.min_repeats_mono == 12
    assert settings.ssr.min_repeats_di == 6
    # Unset values keep defaults
    assert settings.ssr.min_repeats_tri == 3
    assert settings.logging.level == "DEBUG"
    assert settings.logging.format == "text"
    assert settings.output.dir == "test_outputs"


def test_env_override(monkeypatch: pytest.MonkeyPatch):
    """Environment variables override defaults."""
    monkeypatch.setenv("GWICO_SSR_DB_URL", "postgresql://localhost/gwico")
    monkeypatch.setenv("GWICO_SSR_LOG_LEVEL", "WARNING")
    settings = load_settings("/nonexistent/path.toml")
    assert settings.database.url == "postgresql://localhost/gwico"
    assert settings.logging.level == "WARNING"


def test_min_repeats_for_size():
    """SSR min_repeats_for_size returns correct values for each motif size."""
    settings = load_settings("/nonexistent/path.toml")
    assert settings.ssr.min_repeats_for_size(1) == 10
    assert settings.ssr.min_repeats_for_size(2) == 5
    assert settings.ssr.min_repeats_for_size(3) == 3
    assert settings.ssr.min_repeats_for_size(6) == 2


def test_config_from_env_variable(sample_toml: Path, monkeypatch: pytest.MonkeyPatch):
    """GWICO_SSR_CONFIG env var is used when no explicit path is given."""
    monkeypatch.setenv("GWICO_SSR_CONFIG", str(sample_toml))
    settings = load_settings()
    assert settings.database.url == "sqlite:///test.db"
