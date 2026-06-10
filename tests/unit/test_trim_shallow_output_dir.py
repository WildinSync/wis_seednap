"""Regression test for standalone trimming to a shallow output directory.

StandardTrimmer.trim_sample built the per-sample cutadapt log path from
``output_dir.parent.parent / "logs" / ...``. With a shallow standalone
``-o`` (e.g. ``-o /tmp/out``), ``parent.parent`` resolves to ``/`` and the
cutadapt runner's ``mkdir("/logs")`` raises PermissionError, killing the run
on the first sample (0 samples trimmed, orphan ``*_TEMPORARY.fastq`` left
behind). The fix derives the log dir from ``output_dir`` itself
(``output_dir / "logs"``), which is self-contained for both the standalone
CLI and the orchestrator (whose output_dir is ``<output>/01_trim/<marker>``).

These tests do not run cutadapt; CutadaptRunner.trim_primers is replaced with
a fake that mirrors the real behavior relevant here: it creates the requested
output FASTQs and writes the log file (which mkdir's the log dir's parent).
"""

from __future__ import annotations

from pathlib import Path

from seednap.steps.trimming.trimming_pipeline import StandardTrimmer

# Path to the production source, used only for the source-reading regression guard below.
_REPO_SRC = Path(__file__).resolve().parents[2] / "src"

# The regression guard reads the production source directly.
_REPO_TRIMMING = (
    _REPO_SRC / "seednap" / "steps" / "trimming" / "trimming_pipeline.py"
)


class _FakeCutadapt:
    """Stand-in for CutadaptRunner that mimics the file/log side effects.

    Real cutadapt creates the -o/-p outputs and the wrapper appends to the
    log file, creating the log dir's parent on the way. The fake reproduces
    exactly those filesystem touches so the log-path logic is exercised.
    """

    def trim_primers(self, *, r1_output, r2_output=None, log_file=None, **kwargs):
        if log_file is not None:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text("fake cutadapt log\n")
        Path(r1_output).write_text("")
        if r2_output is not None:
            Path(r2_output).write_text("")
        return ""


def _make_inputs(tmp_path: Path) -> tuple:
    raw = tmp_path / "raw"
    raw.mkdir()
    r1 = raw / "sampleA_R1.fastq"
    r2 = raw / "sampleA_R2.fastq"
    r1.write_text("@r\nACGT\n+\nIIII\n")
    r2.write_text("@r\nACGT\n+\nIIII\n")
    return r1, r2


def test_trim_to_shallow_output_dir_does_not_raise(tmp_path: Path) -> None:
    r1, r2 = _make_inputs(tmp_path)

    # A shallow output dir: output_dir.parent.parent would resolve toward the
    # filesystem root, which is the bug this guards against.
    output_dir = tmp_path / "out"

    trimmer = StandardTrimmer()
    trimmer.cutadapt = _FakeCutadapt()

    # Must not raise (previously crashed building/creating "/logs").
    r1_final, r2_final = trimmer.trim_sample(
        r1_input=r1,
        r2_input=r2,
        output_dir=output_dir,
        sample_name="sampleA",
        forward_primer="ACGT",
        reverse_primer="TGCA",
    )

    assert Path(r1_final).exists()
    assert Path(r2_final).exists()

    # Logs land under output_dir/logs, not via parent.parent.
    log_dir = output_dir / "logs"
    assert (log_dir / "sampleA_trim_pass1.txt").exists()
    assert (log_dir / "sampleA_trim_pass2.txt").exists()


def test_temp_files_cleaned_up_after_success(tmp_path: Path) -> None:
    r1, r2 = _make_inputs(tmp_path)
    output_dir = tmp_path / "out"

    trimmer = StandardTrimmer()
    trimmer.cutadapt = _FakeCutadapt()
    trimmer.trim_sample(
        r1_input=r1,
        r2_input=r2,
        output_dir=output_dir,
        sample_name="sampleA",
        forward_primer="ACGT",
        reverse_primer="TGCA",
    )

    # No orphan pass-1 temp files.
    assert not list(output_dir.glob("*_TEMPORARY.fastq"))


def test_temp_files_cleaned_up_on_failure(tmp_path: Path) -> None:
    """An aborted run must not leave orphan *_TEMPORARY.fastq files."""
    r1, r2 = _make_inputs(tmp_path)
    output_dir = tmp_path / "out"

    class _FailSecondPass(_FakeCutadapt):
        def __init__(self) -> None:
            self.calls = 0

        def trim_primers(self, **kwargs):
            self.calls += 1
            if self.calls == 1:
                # Pass 1 succeeds and writes the temp files.
                return super().trim_primers(**kwargs)
            raise RuntimeError("simulated cutadapt failure on pass 2")

    trimmer = StandardTrimmer()
    trimmer.cutadapt = _FailSecondPass()

    raised = False
    try:
        trimmer.trim_sample(
            r1_input=r1,
            r2_input=r2,
            output_dir=output_dir,
            sample_name="sampleA",
            forward_primer="ACGT",
            reverse_primer="TGCA",
        )
    except RuntimeError:
        raised = True

    assert raised, "expected the pass-2 failure to propagate"
    # The finally-block removed the orphan temp files despite the failure.
    assert not list(output_dir.glob("*_TEMPORARY.fastq"))


def test_production_does_not_use_parent_parent_for_logs() -> None:
    """Guard against a regression to the parent.parent log path."""
    src = _REPO_TRIMMING.read_text()
    assert "output_dir.parent.parent" not in src
