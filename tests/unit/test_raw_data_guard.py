"""Regression tests for GitHub issue #2 (wrong data processed).

Two safeguards:
1. The shipped generic teleo *reference* config must NOT carry a concrete real
   dataset path in paths.raw_data -- it ships a placeholder, so an unedited run
   fails the preflight check instead of silently processing the wrong dataset.
2. Sample discovery raises loudly when paths.raw_data exists but holds no flat
   per-sample FASTQs (e.g. a per-library/subdirectory or multiplexed layout),
   rather than returning an empty list and producing an empty run.
"""

from pathlib import Path

import pytest
import yaml

from seednap.config.loader import load_config
from seednap.pipeline.orchestrator import PipelineOrchestrator

_REPO = Path(__file__).resolve().parents[2]


def test_teleo_reference_config_ships_placeholder_raw_data():
    cfg = yaml.safe_load((_REPO / "config" / "markers" / "teleo.yaml").read_text())
    raw_data = cfg["paths"]["raw_data"]
    assert raw_data.startswith("/path/to"), (
        f"teleo.yaml is a reference template and must ship a placeholder raw_data, "
        f"not a concrete dataset path (got {raw_data!r})."
    )


def test_get_sample_list_raises_when_no_flat_fastqs(tmp_path, monkeypatch):
    # Contain any relative dir-creation side effects inside tmp_path.
    monkeypatch.chdir(tmp_path)
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo_rhone.yaml"))

    raw = tmp_path / "raw"
    (raw / "VL492___MB0725C1__").mkdir(parents=True)  # a per-library subdir, no flat FASTQs
    cfg.paths.raw_data = raw
    cfg.paths.output = tmp_path / "out"
    cfg.paths.logs = tmp_path / "out" / "logs"

    orch = PipelineOrchestrator(cfg)
    with pytest.raises(FileNotFoundError, match="No forward-read FASTQ"):
        orch._get_sample_list()
