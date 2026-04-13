"""Configuration loading for GWICO-SSR.

Loads settings from TOML config files with environment variable overrides.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional

from pydantic import Field
from pydantic_settings import BaseSettings

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

# Default config search path
_DEFAULT_CONFIG_NAME = "default.toml"
_CONFIG_DIR = Path(__file__).resolve().parents[3] / "config"


def _find_config_file(explicit_path: Optional[str] = None) -> Optional[Path]:
    """Locate the configuration file.

    Priority:
    1. Explicit path passed as argument
    2. GWICO_SSR_CONFIG environment variable
    3. config/default.toml relative to project root
    """
    if explicit_path:
        p = Path(explicit_path)
        if p.is_file():
            return p
        return None

    env_path = os.environ.get("GWICO_SSR_CONFIG")
    if env_path:
        p = Path(env_path)
        if p.is_file():
            return p

    default = _CONFIG_DIR / _DEFAULT_CONFIG_NAME
    if default.is_file():
        return default

    return None


class DatabaseSettings(BaseSettings):
    url: str = Field(default="sqlite:///gwico_ssr.db", description="Database connection URL")
    echo: bool = Field(default=False, description="Echo SQL statements")

    model_config = {"env_prefix": "GWICO_SSR_DB_"}


class NCBISettings(BaseSettings):
    api_key: str = Field(default="", description="NCBI Entrez API key")
    email: str = Field(default="", description="Email for NCBI Entrez")
    max_retries: int = Field(default=3, description="Maximum download retries")
    rate_limit: float = Field(default=3.0, description="Requests per second (10 with API key)")
    batch_size: int = Field(default=500, description="Accessions per Entrez batch")

    model_config = {"env_prefix": "GWICO_SSR_NCBI_"}


class SSRSettings(BaseSettings):
    min_repeats_mono: int = Field(default=10, description="Min repeats for mononucleotide SSRs")
    min_repeats_di: int = Field(default=5, description="Min repeats for dinucleotide SSRs")
    min_repeats_tri: int = Field(default=3, description="Min repeats for trinucleotide SSRs")
    min_repeats_tetra: int = Field(default=3, description="Min repeats for tetranucleotide SSRs")
    min_repeats_penta: int = Field(default=3, description="Min repeats for pentanucleotide SSRs")
    min_repeats_hexa: int = Field(default=2, description="Min repeats for hexanucleotide SSRs")

    model_config = {"env_prefix": "GWICO_SSR_SSR_"}

    def min_repeats_for_size(self, motif_size: int) -> int:
        mapping = {
            1: self.min_repeats_mono,
            2: self.min_repeats_di,
            3: self.min_repeats_tri,
            4: self.min_repeats_tetra,
            5: self.min_repeats_penta,
            6: self.min_repeats_hexa,
        }
        return mapping[motif_size]


class LoggingSettings(BaseSettings):
    level: str = Field(default="INFO", description="Log level")
    format: str = Field(default="json", description="Log format: json or text")
    file: str = Field(default="", description="Log file path (empty = stderr only)")

    model_config = {"env_prefix": "GWICO_SSR_LOG_"}


class OutputSettings(BaseSettings):
    dir: str = Field(default="outputs", description="Base output directory")

    model_config = {"env_prefix": "GWICO_SSR_OUTPUT_"}


class Settings(BaseSettings):
    database: DatabaseSettings = Field(default_factory=DatabaseSettings)
    ncbi: NCBISettings = Field(default_factory=NCBISettings)
    ssr: SSRSettings = Field(default_factory=SSRSettings)
    logging: LoggingSettings = Field(default_factory=LoggingSettings)
    output: OutputSettings = Field(default_factory=OutputSettings)

    model_config = {"env_prefix": "GWICO_SSR_"}


def load_settings(config_path: Optional[str] = None) -> Settings:
    """Load settings from TOML file, then apply environment variable overrides.

    TOML values are loaded first, then Pydantic's environment variable
    support overrides any value set in the environment.
    """
    toml_data: dict = {}
    cfg_file = _find_config_file(config_path)
    if cfg_file is not None:
        with open(cfg_file, "rb") as f:
            toml_data = tomllib.load(f)

    # Build nested settings from TOML, letting env vars override
    db_data = toml_data.get("database", {})
    ncbi_data = toml_data.get("ncbi", {})
    ssr_data = toml_data.get("ssr", {})
    log_data = toml_data.get("logging", {})
    output_data = toml_data.get("output", {})

    return Settings(
        database=DatabaseSettings(**db_data),
        ncbi=NCBISettings(**ncbi_data),
        ssr=SSRSettings(**ssr_data),
        logging=LoggingSettings(**log_data),
        output=OutputSettings(**output_data),
    )
