"""Unit test for the export-staleness guard on --resume after a clean retry.

If run 1 leaves `clean` SKIPPED (transient error, no cleaned_table) and `export`
COMPLETED against the uncleaned table, a run-2 --resume can re-run clean and now
COMPLETE it (writing a fresh cleaned_table). But export is already COMPLETED, so
_should_run_step('export') returns False and export is NOT re-run -- the GBIF CSV
silently stays stale. The orchestrator must emit a [WARN] when clean completed
AFTER the already-completed export.

_warn_if_export_predates_clean only touches self.state, so we exercise it on a
lightweight stand-in carrying a real PipelineState rather than building a full
orchestrator (which needs a config + directories).
"""

from __future__ import annotations

from datetime import datetime, timedelta
from types import SimpleNamespace

from seednap.pipeline.orchestrator import PipelineOrchestrator
from seednap.pipeline.state import PipelineState, StepStatus


def _state_with_clean_after_export(*, cleaned_table: bool) -> PipelineState:
    state = PipelineState(marker="teleo")
    base = datetime(2026, 1, 1, 12, 0, 0)

    export = state.add_step("export")
    export.status = StepStatus.COMPLETED
    export.completed_at = base
    export.outputs = {"gbif_csv": "/tmp/teleo_blast_gbif.csv"}

    clean = state.add_step("clean")
    clean.status = StepStatus.COMPLETED
    clean.completed_at = base + timedelta(minutes=5)  # completed AFTER export
    clean.outputs = (
        {"cleaned_table": "/tmp/teleo_blast_cleaned.csv"} if cleaned_table else {}
    )
    return state


def _warn(stand_in) -> None:
    PipelineOrchestrator._warn_if_export_predates_clean(
        stand_in, stand_in.state.get_step("export")
    )


def test_warns_when_clean_completed_after_export(caplog) -> None:
    state = _state_with_clean_after_export(cleaned_table=True)
    stand_in = SimpleNamespace(state=state)
    with caplog.at_level("WARNING"):
        _warn(stand_in)
    msgs = [r.getMessage() for r in caplog.records]
    assert any("[WARN] export" in m and "stale" in m for m in msgs), msgs


def test_no_warn_when_clean_has_no_cleaned_table(caplog) -> None:
    # Clean re-ran but produced no cleaned_table (e.g. still skipped): nothing stale.
    state = _state_with_clean_after_export(cleaned_table=False)
    stand_in = SimpleNamespace(state=state)
    with caplog.at_level("WARNING"):
        _warn(stand_in)
    assert not any("[WARN] export" in r.getMessage() for r in caplog.records)


def test_no_warn_when_export_postdates_clean(caplog) -> None:
    # Normal ordering: export completed after clean -> the GBIF CSV is current.
    state = PipelineState(marker="teleo")
    base = datetime(2026, 1, 1, 12, 0, 0)
    clean = state.add_step("clean")
    clean.status = StepStatus.COMPLETED
    clean.completed_at = base
    clean.outputs = {"cleaned_table": "/tmp/teleo_blast_cleaned.csv"}
    export = state.add_step("export")
    export.status = StepStatus.COMPLETED
    export.completed_at = base + timedelta(minutes=5)
    stand_in = SimpleNamespace(state=state)
    with caplog.at_level("WARNING"):
        _warn(stand_in)
    assert not any("[WARN] export" in r.getMessage() for r in caplog.records)
