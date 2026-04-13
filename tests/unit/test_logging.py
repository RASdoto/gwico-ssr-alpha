"""Tests for structured logging setup."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from gwico_ssr.logging import setup_logging, set_run_id


def test_setup_logging_returns_logger():
    """setup_logging returns a configured logger."""
    logger = setup_logging(level="DEBUG", fmt="text")
    assert isinstance(logger, logging.Logger)
    assert logger.name == "gwico_ssr"
    assert logger.level == logging.DEBUG


def test_setup_logging_json_format(capsys):
    """JSON formatter emits valid JSON to stderr."""
    logger = setup_logging(level="INFO", fmt="json")
    logger.info("test message")
    captured = capsys.readouterr()
    # JSON goes to stderr
    log_line = captured.err.strip().split("\n")[-1]
    parsed = json.loads(log_line)
    assert parsed["message"] == "test message"
    assert parsed["level"] == "INFO"
    assert "timestamp" in parsed


def test_setup_logging_text_format(capsys):
    """Text formatter emits human-readable output."""
    logger = setup_logging(level="INFO", fmt="text")
    logger.info("hello world")
    captured = capsys.readouterr()
    assert "hello world" in captured.err


def test_setup_logging_with_run_id(capsys):
    """Run ID appears in JSON log output."""
    logger = setup_logging(level="INFO", fmt="json", run_id="test-run-001")
    logger.info("with run id")
    captured = capsys.readouterr()
    log_line = captured.err.strip().split("\n")[-1]
    parsed = json.loads(log_line)
    assert parsed["run_id"] == "test-run-001"


def test_set_run_id_updates_context(capsys):
    """set_run_id updates the run_id on subsequent log records."""
    logger = setup_logging(level="INFO", fmt="json")
    set_run_id("updated-run-42")
    logger.info("after update")
    captured = capsys.readouterr()
    log_line = captured.err.strip().split("\n")[-1]
    parsed = json.loads(log_line)
    assert parsed["run_id"] == "updated-run-42"


def test_log_file_creation(tmp_path: Path):
    """Log file is created when a path is specified."""
    log_file = tmp_path / "test.log"
    logger = setup_logging(level="INFO", fmt="text", log_file=str(log_file))
    logger.info("file log test")
    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8")
    assert "file log test" in content


def test_package_import():
    """gwico_ssr package imports successfully."""
    import gwico_ssr
    assert hasattr(gwico_ssr, "__version__")
    assert gwico_ssr.__version__
