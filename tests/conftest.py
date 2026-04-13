"""Shared pytest fixtures for GWICO-SSR tests."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def tmp_dir(tmp_path: Path) -> Path:
    """Provide a temporary directory for test artifacts."""
    return tmp_path


@pytest.fixture
def sample_toml(tmp_path: Path) -> Path:
    """Create a sample TOML config file for testing."""
    config_content = """\
[database]
url = "sqlite:///test.db"
echo = true

[ncbi]
email = "test@example.com"
max_retries = 5

[ssr]
min_repeats_mono = 12
min_repeats_di = 6

[logging]
level = "DEBUG"
format = "text"
file = ""

[output]
dir = "test_outputs"
"""
    config_file = tmp_path / "test_config.toml"
    config_file.write_text(config_content, encoding="utf-8")
    return config_file


@pytest.fixture(autouse=True)
def clean_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """Remove GWICO_SSR_ env vars to prevent test pollution."""
    for key in list(os.environ.keys()):
        if key.startswith("GWICO_SSR_"):
            monkeypatch.delenv(key, raising=False)
