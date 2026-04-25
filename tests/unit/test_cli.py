"""Tests for CLI entrypoint."""

from __future__ import annotations

from click.testing import CliRunner

from gwico_ssr import __version__
from gwico_ssr.cli import cli


EXPECTED_ALPHA_COMMANDS = {
    "analyze",
    "annotate",
    "detect",
    "download",
    "export",
    "info",
    "ingest",
    "init-db",
    "metrics",
    "parse",
    "run",
    "visualize",
}

UNEXPECTED_BETA_COMMANDS = {
    "browser",
    "cluster",
    "lineage",
    "normalize-batch",
    "phylo",
    "temporal",
}


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


def test_cli_command_surface_matches_alpha_baseline():
    """CLI help exposes the alpha command surface and no beta commands."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])

    assert result.exit_code == 0

    for command in EXPECTED_ALPHA_COMMANDS:
        assert command in result.output

    for command in UNEXPECTED_BETA_COMMANDS:
        assert command not in result.output
