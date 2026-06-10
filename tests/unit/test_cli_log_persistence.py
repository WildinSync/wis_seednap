"""Unit tests that standalone CLI commands persist their logs (including [WARN]s).

Regression guard for the no-silent-fallback rule (CLAUDE.md section 4 and 7): main()
only configures a console handler, so before this fix the 16 standalone commands emitted
their safety [WARN]s to the console only -- nothing reached a log file, so a fallback that
quietly altered a biological dataset left no on-disk trace.

The fix gives each command that knows an output/marker a per-command log file at
<output-dir>/logs/<name>.log, and routes print_warning through the configured logger so a
warning is both visible on the console AND persisted. These tests assert both:

1. A standalone command writes a log file under the output tree.
2. A print_warning emitted by that command lands in the log with the [WARN] prefix.

`monitor` is used as the vehicle: a state JSON without per-sample read counts makes it
emit exactly one print_warning ("No per-sample ...") without needing any external tool.
"""

from pathlib import Path

from click.testing import CliRunner

from seednap.cli import main
from seednap.pipeline.state import PipelineState


def _bare_state(marker: str, out: Path) -> Path:
    """Write a minimal completed state with no per-sample read counts."""
    state = PipelineState.from_config(marker=marker)
    state.add_step("dada2")
    state.get_step("dada2").complete({})
    state_path = out / f".{marker}_state.json"
    state.save(state_path)
    return state_path


def test_standalone_command_persists_log_file(tmp_path: Path) -> None:
    """`monitor` writes <output-dir>/logs/<marker>.log so its logs are captured on disk."""
    _bare_state("demo", tmp_path)

    result = CliRunner().invoke(main, ["monitor", "demo", "-o", str(tmp_path)])
    assert result.exit_code == 0, result.output

    log_file = tmp_path / "logs" / "demo.log"
    assert log_file.exists(), f"expected a persisted log at {log_file}"
    assert log_file.read_text().strip(), "log file should not be empty"


def test_warning_reaches_the_log(tmp_path: Path) -> None:
    """A print_warning from a standalone command is persisted (with the [WARN] prefix)."""
    _bare_state("demo", tmp_path)

    result = CliRunner().invoke(main, ["monitor", "demo", "-o", str(tmp_path)])
    assert result.exit_code == 0, result.output

    # The warning must still be visible on the console for the user...
    assert "No per-sample" in result.output

    # ...and also captured in the persisted log with the standard [WARN] prefix.
    log_text = (tmp_path / "logs" / "demo.log").read_text()
    assert "[WARN]" in log_text
    assert "No per-sample" in log_text
