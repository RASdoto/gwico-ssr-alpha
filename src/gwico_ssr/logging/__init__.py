"""Structured logging for GWICO-SSR.

Provides JSON or text formatted logging with run-scoped context.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Optional


class JSONFormatter(logging.Formatter):
    """Emit log records as single-line JSON objects."""

    def format(self, record: logging.LogRecord) -> str:
        log_entry = {
            "timestamp": self.formatTime(record, self.datefmt),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        if hasattr(record, "run_id"):
            log_entry["run_id"] = record.run_id
        if record.exc_info and record.exc_info[1]:
            log_entry["exception"] = self.formatException(record.exc_info)
        return json.dumps(log_entry, default=str)


class RunContextFilter(logging.Filter):
    """Inject run_id into log records when set."""

    def __init__(self, run_id: Optional[str] = None) -> None:
        super().__init__()
        self.run_id = run_id

    def filter(self, record: logging.LogRecord) -> bool:
        record.run_id = self.run_id  # type: ignore[attr-defined]
        return True

    def set_run_id(self, run_id: str) -> None:
        self.run_id = run_id


# Module-level filter reference for updating run_id after init
_run_filter: Optional[RunContextFilter] = None


def set_run_id(run_id: str) -> None:
    """Update the run_id on the global log context filter."""
    global _run_filter
    if _run_filter is not None:
        _run_filter.set_run_id(run_id)


def setup_logging(
    level: str = "INFO",
    fmt: str = "json",
    log_file: str = "",
    run_id: Optional[str] = None,
) -> logging.Logger:
    """Configure the root gwico_ssr logger.

    Args:
        level: Log level name (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        fmt: 'json' for structured JSON output, 'text' for human-readable.
        log_file: Path to a log file. Empty string means stderr only.
        run_id: Optional run identifier injected into every log record.

    Returns:
        The configured gwico_ssr root logger.
    """
    global _run_filter

    logger = logging.getLogger("gwico_ssr")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.handlers.clear()

    _run_filter = RunContextFilter(run_id=run_id)
    logger.addFilter(_run_filter)

    if fmt == "json":
        formatter = JSONFormatter()
    else:
        text_fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        if run_id:
            text_fmt = f"%(asctime)s [%(levelname)s] [run:%(run_id)s] %(name)s: %(message)s"
        formatter = logging.Formatter(text_fmt)

    # stderr handler
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    # optional file handler
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(str(log_path), encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger
