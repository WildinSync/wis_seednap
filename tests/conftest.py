"""Pytest configuration and shared fixtures."""

import tempfile
from pathlib import Path
from typing import Any, Dict

import pytest
import yaml


@pytest.fixture
def temp_dir(tmp_path: Path) -> Path:
    """Create a temporary directory for test files."""
    return tmp_path


@pytest.fixture
def sample_primer_config() -> Dict[str, str]:
    """Sample primer configuration."""
    return {
        "forward": "ACACCGCCCGTCACTCT",
        "reverse": "CTTCCGGTACACTTACCATG",
        "name": "Teleo",
        "target": "12S rRNA",
    }


@pytest.fixture
def sample_marker_config(sample_primer_config: Dict[str, str]) -> Dict[str, Any]:
    """Sample marker configuration."""
    return {
        "name": "teleo",
        "description": "Teleost fish eDNA",
        "primers": sample_primer_config,
    }


@pytest.fixture
def minimal_config(sample_marker_config: Dict[str, Any], temp_dir: Path) -> Dict[str, Any]:
    """
    Minimal valid pipeline configuration.

    Uses temporary directory for paths so tests don't interfere with real files.
    """
    return {
        "version": "0.1.0",
        "marker": sample_marker_config,
        "paths": {
            "raw_data": str(temp_dir / "raw"),
            "output": str(temp_dir / "outputs"),
            "logs": str(temp_dir / "logs"),
            "references": str(temp_dir / "references"),
        },
        "taxonomy": {
            "method": "dada2",
            "databases": {
                "dada2": {
                    "all": str(temp_dir / "references" / "dada2_all.fasta"),
                    "species": str(temp_dir / "references" / "dada2_species.fasta"),
                }
            },
        },
    }


@pytest.fixture
def full_config(minimal_config: Dict[str, Any], temp_dir: Path) -> Dict[str, Any]:
    """
    Full pipeline configuration with all sections.

    This represents a complete, production-ready configuration.
    """
    config = minimal_config.copy()

    config.update(
        {
            "demultiplex": {
                "enabled": False,
                "protocol": "none",
            },
            "trimming": {
                "tool": "cutadapt",
                "min_length": 20,
                "max_error_rate": 0.1,
                "cores": 4,
                "discard_untrimmed": True,
            },
            "dada2": {
                "filter": {
                    "max_ee": 2.0,
                    "trunc_q": 11,
                    "max_n": 0,
                    "rm_phix": True,
                },
                "merge": {
                    "min_overlap": 20,
                    "max_mismatch": 0,
                },
                "chimera": {
                    "method": "consensus",
                },
            },
            "export": {
                "formats": ["csv"],
                "gbif": {
                    "enabled": True,
                    "add_rank": True,
                    "add_taxon": True,
                },
            },
            "metrics": {
                "generate_plots": True,
                "plot_format": "png",
                "metrics": ["read_counts", "quality_scores", "length_distribution"],
            },
            "logging": {
                "level": "INFO",
                "format": "detailed",
                "file": True,
                "console": True,
            },
            "pipeline": {
                "steps": ["trim", "dada2", "taxonomy", "export"],
                "skip": [],
            },
        }
    )

    # Add database configs for all taxonomic methods
    config["taxonomy"]["databases"].update(
        {
            "blast": {
                "fasta": str(temp_dir / "references" / "blast_db.fasta"),
                "perc_identity": 80.0,
                "qcov_hsp_perc": 80.0,
                "evalue": 1e-25,
                "max_target_seqs": 5,
                "threshold_species": 98.0,
                "threshold_genus": 96.0,
                "threshold_family": 86.5,
            },
            "ecotag": {
                "tree": str(temp_dir / "references" / "taxonomy"),
                "fasta": str(temp_dir / "references" / "ecotag_db.fasta"),
            },
            "decipher": {
                "trained": str(temp_dir / "references" / "decipher_trained.rds"),
            },
        }
    )

    return config


@pytest.fixture
def config_yaml_file(temp_dir: Path, minimal_config: Dict[str, Any]) -> Path:
    """Create a temporary YAML config file."""
    config_file = temp_dir / "config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(minimal_config, f)
    return config_file


@pytest.fixture
def invalid_config_yaml_file(temp_dir: Path) -> Path:
    """Create a temporary invalid YAML config file."""
    config_file = temp_dir / "invalid_config.yaml"
    with open(config_file, "w") as f:
        f.write("this is not:\nvalid: yaml: syntax:\n  - broken")
    return config_file


@pytest.fixture
def sample_fastq_r1(temp_dir: Path) -> Path:
    """Create a minimal FASTQ R1 file for testing."""
    fastq_file = temp_dir / "sample1_R1.fastq.gz"
    # In real tests, you'd create a gzipped FASTQ
    # For now, just create the path
    fastq_file.touch()
    return fastq_file


@pytest.fixture
def sample_fastq_r2(temp_dir: Path) -> Path:
    """Create a minimal FASTQ R2 file for testing."""
    fastq_file = temp_dir / "sample1_R2.fastq.gz"
    fastq_file.touch()
    return fastq_file


@pytest.fixture
def sample_fasta(temp_dir: Path) -> Path:
    """Create a sample FASTA file for testing."""
    fasta_file = temp_dir / "sequences.fasta"
    with open(fasta_file, "w") as f:
        f.write(">seq1\n")
        f.write("ACGTACGTACGTACGT\n")
        f.write(">seq2\n")
        f.write("TGCATGCATGCATGCA\n")
    return fasta_file
