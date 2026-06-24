"""The 'darwincore' pipeline step: ordering + load-time metadata requirement.

The DarwinCore occurrence file (the GBIF-ready output) is built in-pipeline by the
'darwincore' step, which consumes the long-format 'export' output and per-sample +
per-project metadata. It must come after 'export', and because it needs metadata it is
required (validated up front) to have report.sample_metadata / report.project_metadata set.
"""

from pathlib import Path

import pytest

from seednap.config.loader import load_config
from seednap.config.models.operational import VALID_STEPS, PipelineStepsConfig
from seednap.errors.preflight import preflight_checks
from seednap.pipeline.orchestrator import PipelineOrchestrator

_REPO = Path(__file__).resolve().parents[2]


def test_darwincore_registered():
    assert "darwincore" in VALID_STEPS
    assert hasattr(PipelineOrchestrator, "run_darwincore")


def test_darwincore_requires_export_before_it():
    # valid: export precedes darwincore
    PipelineStepsConfig(steps=["trim", "dada2", "taxonomy", "export", "darwincore", "report"])
    # missing export
    with pytest.raises(ValueError, match="export"):
        PipelineStepsConfig(steps=["trim", "dada2", "taxonomy", "darwincore", "report"])
    # export after darwincore
    with pytest.raises(ValueError, match="export"):
        PipelineStepsConfig(steps=["trim", "dada2", "taxonomy", "darwincore", "export", "report"])


def test_preflight_requires_metadata_when_darwincore_enabled(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo_rhone.yaml"))
    cfg.pipeline.steps = ["trim", "swarm", "taxonomy", "export", "darwincore", "report"]
    cfg.report.sample_metadata = None
    cfg.report.project_metadata = None
    msgs = " ".join(str(p) for p in preflight_checks(cfg))
    assert "darwincore" in msgs
    assert "report.sample_metadata" in msgs
    assert "report.project_metadata" in msgs
