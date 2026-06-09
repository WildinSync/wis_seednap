"""Unit tests for pipeline.steps dependency validation.

pipeline.steps is the single source of truth for what runs and in what order. The
PipelineStepsConfig validator rejects unknown stages, duplicates, both feature paths at
once, and dependency-violating orders -- each with a self-explanatory message.
"""

from __future__ import annotations

import pytest

from seednap.config.models.operational import VALID_STEPS, PipelineStepsConfig


@pytest.mark.parametrize(
    "steps",
    [
        ["trim", "dada2", "taxonomy", "export", "report"],
        ["trim", "swarm", "taxonomy"],
        ["demultiplex", "trim", "dada2", "taxonomy", "clean", "export", "report"],
        ["trim"],
        ["trim", "swarm", "taxonomy", "clean", "report"],
    ],
)
def test_valid_step_orders_accepted(steps) -> None:
    assert PipelineStepsConfig(steps=steps).steps == steps


def test_default_includes_report_and_no_skip() -> None:
    cfg = PipelineStepsConfig()
    assert cfg.steps == ["trim", "dada2", "taxonomy", "export", "report"]
    assert not hasattr(cfg, "skip")  # the old skip list was removed


@pytest.mark.parametrize(
    "steps, needle",
    [
        (["trim", "dada2", "swarm", "taxonomy"], "mutually exclusive"),
        (["taxonomy", "trim"], "requires a feature step"),
        (["dada2", "trim", "taxonomy"], "must come after 'trim'"),
        (["trim", "dada2", "export"], "'export' requires 'taxonomy'"),
        (["trim", "dada2", "export", "taxonomy"], "must come after 'taxonomy'"),
        (["trim", "dada2", "taxonomy", "frobnicate"], "unknown stage"),
        (["trim", "trim", "dada2", "taxonomy"], "duplicate stage"),
        (["trim", "demultiplex", "dada2", "taxonomy"], "'demultiplex' must come before 'trim'"),
        (["trim", "clean"], "requires a feature step"),
    ],
)
def test_invalid_step_orders_rejected_with_message(steps, needle) -> None:
    with pytest.raises(Exception) as exc:
        PipelineStepsConfig(steps=steps)
    assert needle in str(exc.value)


def test_valid_steps_constant_is_complete() -> None:
    # Guard against drift: every stage the orchestrator dispatches must be a valid step.
    assert set(VALID_STEPS) == {
        "demultiplex", "trim", "dada2", "swarm", "taxonomy", "clean", "export", "report",
    }
