"""Unit test for `seednap monitor` (E4) and the E3 per-sample state metadata it consumes."""

from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from seednap.cli import main
from seednap.pipeline.state import PipelineState


def _state_with_counts(marker: str) -> PipelineState:
    state = PipelineState.from_config(marker=marker)
    state.add_step("trim")
    state.add_step("dada2")
    state.get_step("trim").complete({})
    dada = state.get_step("dada2")
    dada.complete({})
    # E3 substrate: a compact summary + the per-sample counts keyed on eventID.
    dada.metadata["read_tracking"] = {
        "n_samples": 2, "raw_reads_total": 100, "final_step": "nonchim",
        "final_reads_total": 80, "mean_retention_pct": 80.0, "n_warnings": 0,
    }
    dada.metadata["read_tracking_per_sample"] = {
        "DAR-1": {"raw": 50, "trimmed": 48, "nonchim": 40},
        "Blank-PCR-1": {"raw": 50, "trimmed": 49, "nonchim": 40},
    }
    return state


def test_monitor_writes_per_sample_summary(tmp_path):
    out = tmp_path
    state = _state_with_counts("demo")
    state.save(out / ".demo_state.json")

    result = CliRunner().invoke(main, ["monitor", "demo", "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert "demo" in result.output and "dada2" in result.output

    csv = out / "04_report" / "demo" / "monitoring_summary.csv"
    assert csv.exists()
    df = pd.read_csv(csv)
    assert set(df["eventID"]) == {"DAR-1", "Blank-PCR-1"}
    assert int(df.set_index("eventID").loc["DAR-1", "nonchim"]) == 40


def test_monitor_missing_state_errors(tmp_path):
    result = CliRunner().invoke(main, ["monitor", "nope", "-o", str(tmp_path)])
    assert result.exit_code != 0
    assert "State file not found" in result.output


def test_monitor_state_without_per_sample_counts_warns(tmp_path):
    """A state with no per-sample counts still summarises steps, but warns (no silent CSV)."""
    state = PipelineState.from_config(marker="bare")
    state.add_step("dada2")
    state.get_step("dada2").complete({})
    state.save(tmp_path / ".bare_state.json")
    result = CliRunner().invoke(main, ["monitor", "bare", "-o", str(tmp_path)])
    assert result.exit_code == 0
    assert not (tmp_path / "04_report" / "bare" / "monitoring_summary.csv").exists()
    assert "No per-sample" in result.output
