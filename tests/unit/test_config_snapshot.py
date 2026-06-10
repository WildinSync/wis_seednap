"""Reproducibility: the orchestrator must snapshot the effective merged config
into the output tree and reference it from the state JSON, and must stamp the
running seednap version into the state.

A run has to be reconstructable from its outputs. On start the orchestrator writes
the fully-merged config (defaults + marker YAML) to
``<paths.output>/.<marker>_config.snapshot.yaml`` and records that path plus the
running ``seednap_version`` on the pipeline state, both of which then land in the
state JSON.

The env's installed ``seednap`` may resolve to a different checkout than the one
under test, so we put this repo's ``src`` first on sys.path before importing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import yaml  # noqa: E402

from seednap.__version__ import __version__ as SEEDNAP_VERSION  # noqa: E402
from seednap.config.loader import load_config  # noqa: E402
from seednap.pipeline.orchestrator import PipelineOrchestrator  # noqa: E402


def _write_minimal_config(root: Path) -> Path:
    """Write a minimal but valid marker YAML and return its path."""
    cfg = root / "teleo.yaml"
    cfg.write_text(
        "\n".join(
            [
                "marker:",
                "  name: teleo",
                "  primers:",
                "    forward: ACACCGCCCGTCACTCT",
                "    reverse: CTTCCGGTACACTTACCATG",
                "paths:",
                f"  raw_data: {root}/raw",
                f"  output: {root}/out",
                f"  logs: {root}/logs",
                "taxonomy:",
                "  method: blast",
                "  databases:",
                "    blast:",
                f"      fasta: {root}/db.fasta",
                "",
            ]
        )
    )
    return cfg


def _build_orchestrator(tmp_path: Path) -> PipelineOrchestrator:
    cfg_path = _write_minimal_config(tmp_path)
    config = load_config(cfg_path)
    return PipelineOrchestrator(config)


def test_config_snapshot_written_to_output_tree(tmp_path) -> None:
    orch = _build_orchestrator(tmp_path)

    snapshot_path = orch.state.config_snapshot_path
    assert snapshot_path is not None
    snapshot_path = Path(snapshot_path)
    # Snapshot lives inside the configured output tree.
    assert snapshot_path.exists()
    assert orch.config.paths.output in snapshot_path.parents

    # It is parseable YAML and reflects the effective merged config (defaults included).
    loaded = yaml.safe_load(snapshot_path.read_text())
    assert loaded["marker"]["name"] == "teleo"
    # A defaulted key the YAML never set is present in the snapshot.
    assert "pipeline" in loaded and "steps" in loaded["pipeline"]


def test_config_snapshot_path_recorded_in_state_json(tmp_path) -> None:
    orch = _build_orchestrator(tmp_path)

    state_json = json.loads(Path(orch.state_file).read_text())
    recorded = state_json.get("config_snapshot_path")
    assert recorded is not None
    assert Path(recorded) == Path(orch.state.config_snapshot_path)
    assert Path(recorded).exists()


def test_state_json_carries_seednap_version(tmp_path) -> None:
    orch = _build_orchestrator(tmp_path)

    assert orch.state.seednap_version == SEEDNAP_VERSION

    state_json = json.loads(Path(orch.state_file).read_text())
    assert state_json.get("seednap_version") == SEEDNAP_VERSION
