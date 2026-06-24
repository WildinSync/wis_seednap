"""Early diagnostic WARN when primer trimming discards most reads.

A run that feeds already-primer-trimmed FASTQs into the default
``trimming.discard_untrimmed=True`` path loses nearly every read at the trim step
(Cutadapt finds no primer to remove and discards the read). The final read-tracking
report already flags this per sample, but only after the long feature/taxonomy/export
steps have run. ``_warn_on_heavy_trim_loss`` surfaces it immediately after trimming,
naming the likely cause and the exact config fix.

These cover the run-level aggregation (``ReadTrackingBuilder.aggregate_trim_loss``)
and that the orchestrator emits the diagnostic only when the loss is catastrophic.
"""

from pathlib import Path

from seednap.config.loader import load_config
from seednap.pipeline.orchestrator import PipelineOrchestrator
from seednap.steps.report import ReadTrackingBuilder
from seednap.steps.trimming.trimming_pipeline import StandardTrimmer

_REPO = Path(__file__).resolve().parents[2]


def _write_trim_logs(logs_dir: Path, sample: str, raw: int, trimmed: int) -> None:
    """Write minimal two-pass Cutadapt logs in the real format (see read_tracking)."""
    logs_dir.mkdir(parents=True, exist_ok=True)
    (logs_dir / f"{sample}_trim_pass1.txt").write_text(
        f"Total read pairs processed:            {raw:,}\n"
        f"Pairs written (passing filters):       {raw:,} (100.0%)\n"
    )
    (logs_dir / f"{sample}_trim_pass2.txt").write_text(
        f"Total read pairs processed:            {raw:,}\n"
        f"Pairs written (passing filters):       {trimmed:,} (x%)\n"
    )


def test_aggregate_trim_loss_computes_run_level_loss(tmp_path):
    logs = tmp_path / "logs"
    _write_trim_logs(logs, "S1", raw=1000, trimmed=200)
    _write_trim_logs(logs, "S2", raw=1000, trimmed=600)

    raw_total, trimmed_total, loss_pct = ReadTrackingBuilder(
        "m", logs_dir=logs
    ).aggregate_trim_loss()

    assert (raw_total, trimmed_total) == (2000, 800)
    assert loss_pct == 60.0  # 1 - 800/2000


def test_aggregate_trim_loss_none_when_unmeasurable(tmp_path):
    empty = tmp_path / "logs"
    empty.mkdir()
    assert ReadTrackingBuilder("m", logs_dir=empty).aggregate_trim_loss() is None


def _orchestrator_with_one_sample(tmp_path, monkeypatch):
    """Load the teleo config, point it at a fresh raw/output tree with one sample."""
    monkeypatch.chdir(tmp_path)
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo_rhone.yaml"))
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "S1_R1.fastq.gz").write_bytes(b"")
    (raw / "S1_R2.fastq.gz").write_bytes(b"")
    cfg.paths.raw_data = raw
    cfg.paths.output = tmp_path / "out"
    cfg.paths.logs = tmp_path / "out" / "logs"
    cfg.report.warn_step_loss_pct = 70.0
    return cfg, PipelineOrchestrator(cfg)


def _fake_trimmer_writing_logs(raw: int, trimmed: int):
    """A StandardTrimmer.trim_sample stand-in that also writes the pass1/pass2 logs."""

    def fake_trim_sample(self, *, output_dir, sample_name, **_kw):
        output_dir.mkdir(parents=True, exist_ok=True)
        r1 = output_dir / f"{sample_name}.R1.fastq"
        r2 = output_dir / f"{sample_name}.R2.fastq"
        r1.write_text("x")
        r2.write_text("x")
        _write_trim_logs(output_dir / "logs", sample_name, raw=raw, trimmed=trimmed)
        return (r1, r2)

    return fake_trim_sample


def _run_log_text(cfg) -> str:
    """The orchestrator's run log (it clears caplog's handler, so we read its file)."""
    logs = list(Path(cfg.paths.logs).glob("*_pipeline_run.log"))
    assert logs, "the orchestrator must have written a run log"
    return logs[0].read_text()


def test_trim_warns_on_heavy_loss_naming_discard_untrimmed(tmp_path, monkeypatch):
    cfg, orch = _orchestrator_with_one_sample(tmp_path, monkeypatch)
    cfg.trimming.discard_untrimmed = True
    monkeypatch.setattr(
        StandardTrimmer, "trim_sample", _fake_trimmer_writing_logs(raw=100_000, trimmed=200)
    )

    orch.run_trim()

    log = _run_log_text(cfg)
    assert "[WARN] trim:" in log and "discarded" in log, "heavy trim loss must raise a [WARN]"
    assert "discard_untrimmed: false" in log  # names the exact fix
    assert "99.8%" in log  # 200/100000 retained -> 99.8% discarded


def test_trim_quiet_on_normal_loss(tmp_path, monkeypatch):
    cfg, orch = _orchestrator_with_one_sample(tmp_path, monkeypatch)
    cfg.trimming.discard_untrimmed = True
    monkeypatch.setattr(
        StandardTrimmer, "trim_sample", _fake_trimmer_writing_logs(raw=100_000, trimmed=80_000)
    )

    orch.run_trim()

    log = _run_log_text(cfg)
    assert "[WARN] trim:" not in log, "a normal trim loss (20%) must not raise the heavy-loss WARN"
