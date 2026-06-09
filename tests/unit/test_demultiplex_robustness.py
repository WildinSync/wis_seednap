"""Unit tests for DemultiplexConfig and per-sample failure handling.

These tests don't run cutadapt -- they exercise the configuration plumbing only.
Demultiplexing now runs iff "demultiplex" is listed in pipeline.steps (there is no
separate enabled/skip flag); these tests cover the remaining demultiplex params.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from seednap.config.models import DemultiplexConfig


def test_demultiplex_defaults() -> None:
    cfg = DemultiplexConfig()
    assert cfg.protocol == "none"
    assert cfg.max_sample_failure_rate == 0.5


def test_max_sample_failure_rate_validates_range() -> None:
    # Below 0 -> reject
    with pytest.raises(Exception):
        DemultiplexConfig(max_sample_failure_rate=-0.1)
    # Above 1 -> reject
    with pytest.raises(Exception):
        DemultiplexConfig(max_sample_failure_rate=1.5)
    # Edges allowed
    assert DemultiplexConfig(max_sample_failure_rate=0.0).max_sample_failure_rate == 0.0
    assert DemultiplexConfig(max_sample_failure_rate=1.0).max_sample_failure_rate == 1.0


def test_strict_validation_rejects_removed_and_unknown_keys() -> None:
    """The removed enable/skip gates are now rejected, as are typos (StrictModel)."""
    with pytest.raises(Exception):
        DemultiplexConfig(enabled=True)  # removed: opt in via pipeline.steps instead
    with pytest.raises(Exception):
        DemultiplexConfig(skip=True)  # removed: omit "demultiplex" from steps instead
    with pytest.raises(Exception):
        DemultiplexConfig(protocoll="ligation")  # typo


def test_demultiplex_section_round_trips() -> None:
    """A demultiplex section round-trips; the step runs only when listed in pipeline.steps."""
    from seednap.config.loader import load_config

    yaml_text = """
marker:
  name: "teleo"
  primers:
    forward: "ACACCGCCCGTCACTCT"
    reverse: "CTTCCGGTACACTTACCATG"
demultiplex:
  protocol: "ligation"
paths:
  raw_data: "/tmp/raw"
  output: "/tmp/output"
  logs: "/tmp/logs"
taxonomy:
  method: "blast"
  databases:
    blast:
      fasta: "/tmp/ref.fasta"
"""
    p = Path("/tmp/test_demux_roundtrip.yaml")
    p.write_text(yaml_text)
    cfg = load_config(p)
    assert cfg.demultiplex.protocol == "ligation"
    # Not run unless explicitly listed in pipeline.steps.
    assert "demultiplex" not in cfg.pipeline.steps
