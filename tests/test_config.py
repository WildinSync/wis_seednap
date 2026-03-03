"""Tests for configuration system."""

from pathlib import Path
from typing import Any, Dict

import pytest
from pydantic import ValidationError

from seednap.config import (
    ConfigError,
    PipelineConfig,
    PrimerConfig,
    load_config,
    load_yaml,
    merge_configs,
    validate_config_file,
)


class TestPrimerConfig:
    """Tests for PrimerConfig model."""

    def test_valid_primer_config(self, sample_primer_config: Dict[str, str]) -> None:
        """Test that valid primer config is accepted."""
        config = PrimerConfig(**sample_primer_config)
        assert config.forward == "ACACCGCCCGTCACTCT"
        assert config.reverse == "CTTCCGGTACACTTACCATG"
        assert config.name == "Teleo"

    def test_primer_sequences_converted_to_uppercase(self) -> None:
        """Test that primer sequences are converted to uppercase."""
        config = PrimerConfig(forward="acaccgcccgtcactct", reverse="cttccggtacacttaccatg")
        assert config.forward == "ACACCGCCCGTCACTCT"
        assert config.reverse == "CTTCCGGTACACTTACCATG"

    def test_invalid_dna_sequence_rejected(self) -> None:
        """Test that invalid DNA sequences are rejected."""
        with pytest.raises(ValidationError) as exc_info:
            PrimerConfig(forward="ACACGCCCGTCACTCT", reverse="CTTCCGGTACXCTTACCATG")

        assert "Invalid DNA sequence" in str(exc_info.value)

    def test_ambiguous_bases_accepted(self) -> None:
        """Test that IUPAC ambiguity codes are accepted."""
        config = PrimerConfig(forward="ATCGRYMKSWHBVDN", reverse="GCTASWKMRYN")
        assert config.forward == "ATCGRYMKSWHBVDN"

    def test_minimum_length_enforced(self) -> None:
        """Test that minimum primer length is enforced."""
        with pytest.raises(ValidationError) as exc_info:
            PrimerConfig(forward="ATCG", reverse="GCTA")

        assert "at least 10 characters" in str(exc_info.value)


class TestPipelineConfig:
    """Tests for PipelineConfig model."""

    def test_minimal_config_valid(self, minimal_config: Dict[str, Any]) -> None:
        """Test that minimal valid config is accepted."""
        config = PipelineConfig(**minimal_config)
        assert config.marker.name == "teleo"
        assert config.taxonomy.method == "dada2"

    def test_full_config_valid(self, full_config: Dict[str, Any]) -> None:
        """Test that full config with all sections is accepted."""
        config = PipelineConfig(**full_config)
        assert config.marker.name == "teleo"
        assert config.trimming.cores == 4
        assert config.dada2.filter.max_ee == 2.0

    def test_default_values_applied(self, minimal_config: Dict[str, Any]) -> None:
        """Test that default values are applied for optional fields."""
        config = PipelineConfig(**minimal_config)

        # Check defaults
        assert config.trimming.tool == "cutadapt"
        assert config.trimming.min_length == 20
        assert config.dada2.filter.max_ee == 2.0
        assert config.export.gbif.enabled is True
        assert config.logging.level == "INFO"

    def test_invalid_taxonomy_method_rejected(self, minimal_config: Dict[str, Any]) -> None:
        """Test that invalid taxonomy method is rejected."""
        minimal_config["taxonomy"]["method"] = "invalid_method"

        with pytest.raises(ValidationError) as exc_info:
            PipelineConfig(**minimal_config)

        assert "taxonomy" in str(exc_info.value).lower()

    def test_missing_database_for_method_rejected(self, minimal_config: Dict[str, Any]) -> None:
        """Test that missing database config for selected method is rejected."""
        # Set method to blast but don't provide blast database
        minimal_config["taxonomy"]["method"] = "blast"

        with pytest.raises(ValidationError) as exc_info:
            PipelineConfig(**minimal_config)

        assert "blast" in str(exc_info.value).lower()

    def test_paths_expanded_to_absolute(self, minimal_config: Dict[str, Any]) -> None:
        """Test that paths are expanded to absolute paths."""
        config = PipelineConfig(**minimal_config)

        assert config.paths.raw_data.is_absolute()
        assert config.paths.output.is_absolute()
        assert config.paths.logs.is_absolute()


class TestConfigLoader:
    """Tests for configuration loader functions."""

    def test_load_yaml_success(self, config_yaml_file: Path) -> None:
        """Test successfully loading a YAML file."""
        config_dict = load_yaml(config_yaml_file)

        assert isinstance(config_dict, dict)
        assert "marker" in config_dict
        assert "taxonomy" in config_dict

    def test_load_yaml_file_not_found(self, temp_dir: Path) -> None:
        """Test loading non-existent file raises ConfigError."""
        non_existent = temp_dir / "does_not_exist.yaml"

        with pytest.raises(ConfigError) as exc_info:
            load_yaml(non_existent)

        assert "not found" in str(exc_info.value).lower()

    def test_load_yaml_invalid_syntax(self, invalid_config_yaml_file: Path) -> None:
        """Test loading file with invalid YAML syntax raises ConfigError."""
        with pytest.raises(ConfigError) as exc_info:
            load_yaml(invalid_config_yaml_file)

        assert "invalid yaml" in str(exc_info.value).lower()

    def test_load_config_success(self, config_yaml_file: Path) -> None:
        """Test successfully loading and validating config."""
        config = load_config(config_yaml_file)

        assert isinstance(config, PipelineConfig)
        assert config.marker.name == "teleo"

    def test_load_config_validation_failure(self, temp_dir: Path) -> None:
        """Test that invalid config raises ConfigError with details."""
        invalid_config_file = temp_dir / "invalid.yaml"
        with open(invalid_config_file, "w") as f:
            f.write("marker:\n  name: test\n")  # Missing required fields

        with pytest.raises(ConfigError) as exc_info:
            load_config(invalid_config_file)

        assert "validation failed" in str(exc_info.value).lower()

    def test_validate_config_file_valid(self, config_yaml_file: Path) -> None:
        """Test validating a valid config file."""
        is_valid, error_message = validate_config_file(config_yaml_file)

        assert is_valid is True
        assert error_message is None

    def test_validate_config_file_invalid(self, temp_dir: Path) -> None:
        """Test validating an invalid config file."""
        invalid_config_file = temp_dir / "invalid.yaml"
        with open(invalid_config_file, "w") as f:
            f.write("not: a: valid: config\n")

        is_valid, error_message = validate_config_file(invalid_config_file)

        assert is_valid is False
        assert error_message is not None


class TestConfigMerging:
    """Tests for configuration merging."""

    def test_merge_configs_simple(self) -> None:
        """Test merging two simple dictionaries."""
        base = {"a": 1, "b": 2}
        override = {"b": 3, "c": 4}

        merged = merge_configs(base, override)

        assert merged == {"a": 1, "b": 3, "c": 4}

    def test_merge_configs_nested(self) -> None:
        """Test merging nested dictionaries."""
        base = {
            "section1": {"key1": "value1", "key2": "value2"},
            "section2": {"key3": "value3"},
        }
        override = {
            "section1": {"key2": "override_value"},
            "section3": {"key4": "value4"},
        }

        merged = merge_configs(base, override)

        assert merged["section1"]["key1"] == "value1"
        assert merged["section1"]["key2"] == "override_value"
        assert merged["section2"]["key3"] == "value3"
        assert merged["section3"]["key4"] == "value4"

    def test_merge_configs_lists_replaced(self) -> None:
        """Test that lists are replaced, not merged."""
        base = {"steps": ["trim", "dada2"]}
        override = {"steps": ["taxonomy", "export"]}

        merged = merge_configs(base, override)

        assert merged["steps"] == ["taxonomy", "export"]

    def test_merge_configs_with_defaults(
        self, temp_dir: Path, minimal_config: Dict[str, Any]
    ) -> None:
        """Test loading config with defaults merging."""
        # Create defaults file
        defaults_file = temp_dir / "defaults.yaml"
        defaults = {
            "logging": {"level": "DEBUG"},
            "trimming": {"cores": 8},
        }
        import yaml

        with open(defaults_file, "w") as f:
            yaml.dump(defaults, f)

        # Create user config file
        user_config_file = temp_dir / "user_config.yaml"
        with open(user_config_file, "w") as f:
            yaml.dump(minimal_config, f)

        # Load with defaults
        config = load_config(user_config_file, defaults_path=defaults_file)

        # Defaults should be applied
        assert config.logging.level == "DEBUG"
        assert config.trimming.cores == 8
