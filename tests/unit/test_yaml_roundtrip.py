"""Smoke test: every shipped marker YAML must load cleanly under strict validation.

Catches typos and missing/extra keys at the config layer before they reach the
pipeline runtime.
"""

from pathlib import Path

import pytest

from seednap.config.loader import load_config
from seednap.config.models import PipelineConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
MARKER_DIR = REPO_ROOT / "config" / "markers"


def _marker_yamls():
    return sorted(MARKER_DIR.glob("*.yaml"))


@pytest.mark.parametrize("yaml_path", _marker_yamls(), ids=lambda p: p.name)
def test_marker_yaml_loads(yaml_path: Path) -> None:
    """Each marker YAML must round-trip through the strict Pydantic loader."""
    cfg = load_config(yaml_path)
    assert isinstance(cfg, PipelineConfig)
    assert cfg.marker.name, "marker.name must be non-empty"
    assert cfg.marker.primers.forward, "marker.primers.forward required"
    assert cfg.marker.primers.reverse, "marker.primers.reverse required"


@pytest.mark.parametrize("yaml_path", _marker_yamls(), ids=lambda p: p.name)
def test_blast_db_config_roundtrip(yaml_path: Path) -> None:
    """For markers using BLAST, the DB config must include the new fields."""
    cfg = load_config(yaml_path)
    if cfg.taxonomy.method != "blast":
        pytest.skip(f"{yaml_path.name} is not BLAST-based")

    db = cfg.taxonomy.get_database_config()
    # New fields introduced by Commits A, B, C
    assert hasattr(db, "task")
    assert hasattr(db, "threshold_order")
    assert hasattr(db, "threshold_class")
    assert hasattr(db, "top_bitscore_pct")
    # Contaminants moved to TaxonomicAssignmentConfig in Commit F so all
    # methods (BLAST, DECIPHER, ecotag, DADA2) can share it.
    assert hasattr(cfg.taxonomy, "contaminants")

    # Sanity ranges
    assert 0 <= db.threshold_class <= db.threshold_order <= db.threshold_family
    assert db.threshold_family <= db.threshold_genus <= db.threshold_species <= 100
    assert 0 <= db.top_bitscore_pct <= 100
    assert db.task in ("megablast", "blastn", "dc-megablast", "blastn-short")


def test_typo_in_yaml_is_rejected(tmp_path: Path) -> None:
    """Strict validation: a typo in a taxonomy database block is rejected at LOAD time.

    taxonomy.databases is an open Dict[str, Any] (to support multiple methods), but
    validate_databases parses every present method block into its strict (extra="forbid")
    model during load, so a misspelled field errors in load_config() rather than lazily
    mid-run at get_database_config().
    """
    bad_yaml = tmp_path / "bad.yaml"
    bad_yaml.write_text(
        """
version: "0.1.0"
marker:
  name: "test"
  primers:
    forward: "ACGTACGTACGT"
    reverse: "ACGTACGTACGT"
paths:
  raw_data: "/tmp/raw"
  output: "/tmp/output"
  logs: "/tmp/logs"
taxonomy:
  method: "blast"
  databases:
    blast:
      fasta: "/tmp/ref.fasta"
      pecr_identity: 80.0
"""
    )
    with pytest.raises(Exception):
        load_config(bad_yaml)
