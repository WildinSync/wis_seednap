"""The orchestrator's step dispatch must cover operational.VALID_STEPS.

The dispatch table (step name -> run_* method) used to be hardcoded in two places
and could drift from operational.VALID_STEPS. The orchestrator now builds the
dispatch once in __init__ and fails loudly if any valid stage lacks a handler. This
test pins that invariant and the resume version-mismatch [WARN].

The env's installed ``seednap`` may resolve to a different checkout than the one
under test, so we put this repo's ``src`` first on sys.path before importing.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from seednap.config.loader import load_config  # noqa: E402
from seednap.config.models.operational import VALID_STEPS  # noqa: E402
from seednap.pipeline import orchestrator as orchestrator_module  # noqa: E402
from seednap.pipeline.orchestrator import PipelineOrchestrator  # noqa: E402


class _ListHandler(logging.Handler):
    """Captures records on the orchestrator's own logger.

    setup_logging() clears the ROOT logger's handlers (which removes pytest's
    caplog handler), so we attach directly to the module logger instead.
    """

    def __init__(self) -> None:
        super().__init__()
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


def _build_orchestrator(root: Path) -> PipelineOrchestrator:
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
    return PipelineOrchestrator(load_config(cfg))


def test_dispatch_covers_every_valid_step(tmp_path) -> None:
    orch = _build_orchestrator(tmp_path)
    # Every valid stage has a dispatch handler (no missing keys).
    assert set(VALID_STEPS).issubset(set(orch._step_methods))
    # And the dispatch introduces no stage that is not a valid step.
    assert set(orch._step_methods) == set(VALID_STEPS)


def test_dispatch_handlers_are_callable(tmp_path) -> None:
    orch = _build_orchestrator(tmp_path)
    for name in VALID_STEPS:
        assert callable(orch._step_methods[name]), name


def test_resume_version_mismatch_emits_warning(tmp_path) -> None:
    # First run writes a state stamped with the current version.
    orch = _build_orchestrator(tmp_path)
    state_file = orch.state_file
    # Simulate a state written by a different seednap version.
    orch.state.seednap_version = "0.0.0-old"
    orch.state.save(state_file)

    handler = _ListHandler()
    orchestrator_module.logger.addHandler(handler)
    try:
        PipelineOrchestrator(orch.config, state_file=state_file, resume=True)
    finally:
        orchestrator_module.logger.removeHandler(handler)

    assert any(
        "[WARN] resume version mismatch" in m and "0.0.0-old" in m
        for m in handler.messages
    ), handler.messages
