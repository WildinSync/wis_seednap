"""Configuration loading and validation utilities."""

import sys
from pathlib import Path
from typing import Any, Dict, Optional

import yaml
from pydantic import ValidationError

from seednap.config.models import PipelineConfig


class ConfigError(Exception):
    """Configuration error exception."""

    pass


def load_yaml(file_path: Path) -> Dict[str, Any]:
    """
    Load YAML file and return as dictionary.

    Args:
        file_path: Path to YAML file

    Returns:
        Dictionary with configuration

    Raises:
        ConfigError: If file cannot be read or parsed
    """
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            config_dict = yaml.safe_load(f)

        if config_dict is None:
            raise ConfigError(f"Config file is empty: {file_path}")

        if not isinstance(config_dict, dict):
            raise ConfigError(f"Config file must contain a YAML dictionary: {file_path}")

        return config_dict

    except FileNotFoundError as e:
        raise ConfigError(f"Config file not found: {file_path}") from e
    except yaml.YAMLError as e:
        raise ConfigError(f"Invalid YAML in config file {file_path}: {e}") from e
    except Exception as e:
        raise ConfigError(f"Error reading config file {file_path}: {e}") from e


def merge_configs(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively merge two configuration dictionaries.

    The override dictionary takes precedence. Nested dictionaries are merged
    recursively, lists and other values are replaced entirely.

    Args:
        base: Base configuration dictionary
        override: Override configuration dictionary

    Returns:
        Merged configuration dictionary
    """
    merged = base.copy()

    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            # Recursively merge nested dictionaries
            merged[key] = merge_configs(merged[key], value)
        else:
            # Override value (including lists and primitives)
            merged[key] = value

    return merged


def load_config(
    config_path: Path, defaults_path: Optional[Path] = None, validate: bool = True
) -> PipelineConfig:
    """
    Load and validate pipeline configuration.

    This function:
    1. Loads default configuration (if provided)
    2. Loads user configuration
    3. Merges them (user config takes precedence)
    4. Validates using Pydantic model
    5. Returns validated configuration object

    Args:
        config_path: Path to user configuration YAML file
        defaults_path: Optional path to default configuration file
        validate: Whether to validate the configuration (default: True)

    Returns:
        Validated PipelineConfig object

    Raises:
        ConfigError: If configuration is invalid or cannot be loaded
    """
    # Load user config
    config_dict = load_yaml(config_path)

    # Load and merge defaults if provided
    if defaults_path is not None and defaults_path.exists():
        defaults_dict = load_yaml(defaults_path)
        config_dict = merge_configs(defaults_dict, config_dict)

    # Validate using Pydantic model
    if validate:
        try:
            config = PipelineConfig(**config_dict)
            return config
        except ValidationError as e:
            # Format validation errors in a user-friendly way
            error_messages = []
            for error in e.errors():
                location = " -> ".join(str(loc) for loc in error["loc"])
                message = error["msg"]
                error_messages.append(f"  • {location}: {message}")

            error_text = "\n".join(error_messages)
            raise ConfigError(
                f"Configuration validation failed for {config_path}:\n{error_text}"
            ) from e
    else:
        # Return unvalidated dict (for debugging or special cases)
        return config_dict  # type: ignore


def validate_config_file(config_path: Path) -> tuple[bool, Optional[str]]:
    """
    Validate a configuration file without loading it fully.

    Useful for CLI validation commands.

    Args:
        config_path: Path to configuration file

    Returns:
        Tuple of (is_valid, error_message)
        - is_valid: True if valid, False otherwise
        - error_message: Error message if invalid, None if valid
    """
    try:
        load_config(config_path, validate=True)
        return True, None
    except ConfigError as e:
        return False, str(e)
    except Exception as e:
        return False, f"Unexpected error: {e}"


def print_config_error(error_message: str, exit_code: int = 1) -> None:
    """
    Print a formatted configuration error message and exit.

    Args:
        error_message: Error message to display
        exit_code: Exit code (default: 1)
    """
    print(f"\n❌ Configuration Error:\n{error_message}\n", file=sys.stderr)
    sys.exit(exit_code)


def get_default_config_path() -> Optional[Path]:
    """
    Get the default configuration file path.

    Looks for config file in the following order:
    1. ./config.yaml
    2. ./config/default.yaml
    3. ~/.config/seednap/config.yaml

    Returns:
        Path to default config file, or None if not found
    """
    candidates = [
        Path.cwd() / "config.yaml",
        Path.cwd() / "config" / "default.yaml",
        Path.home() / ".config" / "seednap" / "config.yaml",
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    return None


def create_example_config(output_path: Path, marker: str = "teleo") -> None:
    """
    Create an example configuration file.

    Args:
        output_path: Where to write the example config
        marker: Marker name for the example (default: 'teleo')

    Raises:
        ConfigError: If file cannot be written
    """
    example_config = f"""# Seednap Pipeline Configuration for {marker}
version: "0.1.0"

marker:
  name: {marker}
  description: "Example {marker} marker configuration"
  primers:
    forward: "ACACCGCCCGTCACTCT"
    reverse: "CTTCCGGTACACTTACCATG"
    name: "{marker.capitalize()}"
    target: "12S rRNA"

paths:
  raw_data: "data/raw"
  output: "outputs"
  logs: "logs"
  references: "references/{marker}"

demultiplex:
  enabled: false
  protocol: "none"

trimming:
  tool: "cutadapt"
  min_length: 20
  max_error_rate: 0.1
  cores: 4
  discard_untrimmed: true

dada2:
  filter:
    max_ee: 2.0
    trunc_q: 11
    max_n: 0
    rm_phix: true
  merge:
    min_overlap: 20
    max_mismatch: 0
  chimera:
    method: "consensus"

taxonomy:
  method: "dada2"
  databases:
    dada2:
      all: "references/{marker}/dada2_all.fasta"
      species: "references/{marker}/dada2_species.fasta"
    blast:
      fasta: "references/{marker}/blast_db.fasta"
      perc_identity: 80.0
      qcov_hsp_perc: 80.0
      evalue: 1.0e-25
      max_target_seqs: 5
      threshold_species: 98.0
      threshold_genus: 96.0
      threshold_family: 86.5
    ecotag:
      tree: "references/{marker}/taxonomy/"
      fasta: "references/{marker}/ecotag_db.fasta"
    decipher:
      trained: "references/{marker}/decipher_trained.rds"

export:
  formats:
    - "csv"
  gbif:
    enabled: true
    add_rank: true
    add_taxon: true

metrics:
  generate_plots: true
  plot_format: "png"
  metrics:
    - "read_counts"
    - "quality_scores"
    - "length_distribution"

logging:
  level: "INFO"
  format: "detailed"
  file: true
  console: true

resources:
  max_cores: 4
  memory_limit: "16G"
  parallel_samples: 2

pipeline:
  steps:
    - "trim"
    - "dada2"
    - "taxonomy"
    - "export"
  skip: []
"""

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(example_config)
    except Exception as e:
        raise ConfigError(f"Failed to write example config to {output_path}: {e}") from e
