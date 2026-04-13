"""Tests for CLI entrypoint."""

from __future__ import annotations

from click.testing import CliRunner

from gwico_ssr import __version__
from gwico_ssr.cli import cli


def test_cli_help():
    """CLI --help exits cleanly and shows usage."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "GWICO-SSR" in result.output
    assert "In-silico SSR" in result.output


def test_cli_version():
    """CLI --version shows the package version."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_cli_info_command():
    """CLI info subcommand runs without error."""
    runner = CliRunner()
    result = runner.invoke(cli, ["info"])
    assert result.exit_code == 0
    assert "GWICO-SSR v" in result.output
    assert "Database:" in result.output


def test_cli_info_with_log_level():
    """CLI accepts --log-level override."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--log-level", "DEBUG", "info"])
    assert result.exit_code == 0


def test_cli_info_with_log_format():
    """CLI accepts --log-format override."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--log-format", "text", "info"])
    assert result.exit_code == 0
