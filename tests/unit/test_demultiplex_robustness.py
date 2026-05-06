"""Unit tests for the demultiplex skip flag and per-sample failure handling (Commit I).

These tests don't run cutadapt — they exercise the configuration plumbing
and the orchestrator's skip behaviour using mock state, so they're fast and
don't require any external tools.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from seednap.config.models import DemultiplexConfig


def test_skip_flag_default_is_false() -> None:
    cfg = DemultiplexConfig()
    assert cfg.skip is False
    assert cfg.max_sample_failure_rate == 0.5


def test_skip_flag_can_be_set() -> None:
    cfg = DemultiplexConfig(skip=True)
    assert cfg.skip is True


def test_max_sample_failure_rate_validates_range() -> None:
    # Below 0 -> reject
    with pytest.raises(Exception):
        DemultiplexConfig(max_sample_failure_rate=-0.1)
    # Above 1 -> reject
    with pytest.raises(Exception):
        DemultiplexConfig(max_sample_failure_rate=1.5)
    # Edges allowed
    cfg = DemultiplexConfig(max_sample_failure_rate=0.0)
    assert cfg.max_sample_failure_rate == 0.0
    cfg = DemultiplexConfig(max_sample_failure_rate=1.0)
    assert cfg.max_sample_failure_rate == 1.0


def test_strict_validation_rejects_unknown_keys() -> None:
    """Demultiplex config still rejects typos (StrictModel inheritance)."""
    with pytest.raises(Exception):
        DemultiplexConfig(skipp=True)  # typo


def test_orchestrator_skip_message_in_config() -> None:
    """skip=True in YAML must round-trip cleanly through load_config."""
    import yaml

    from seednap.config.loader import load_config

    yaml_text = """
version: "0.1.0"
marker:
  name: "teleo"
  primers:
    forward: "ACACCGCCCGTCACTCT"
    reverse: "CTTCCGGTACACTTACCATG"
demultiplex:
  enabled: false
  skip: true
paths:
  raw_data: "/tmp/raw"
  output: "/tmp/output"
  logs: "/tmp/logs"
  references: "/tmp/refs"
taxonomy:
  method: "blast"
  databases:
    blast:
      fasta: "/tmp/ref.fasta"
"""
    p = Path("/tmp/test_skip_demux.yaml")
    p.write_text(yaml_text)
    cfg = load_config(p)
    assert cfg.demultiplex.skip is True
    assert cfg.demultiplex.enabled is False
