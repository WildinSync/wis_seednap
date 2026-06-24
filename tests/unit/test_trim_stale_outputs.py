"""Regression test for GitHub issue #2 (deeper cause): stale trimmed outputs reused.

A trim run wrote its per-sample FASTQs into <output>/01_trim/<marker>, and the
downstream feature step (SWARM/DADA2) discovered its inputs by scanning that
directory. So a re-run that produced a *different* sample set left the earlier
run's trimmed files behind, and the feature step silently processed those stale
samples -- exactly how a run configured for one dataset kept emitting results for
a previously-trimmed one. The trim step now clears stale outputs before writing.
"""

from pathlib import Path

from seednap.config.loader import load_config
from seednap.pipeline.orchestrator import PipelineOrchestrator
from seednap.steps.trimming.trimming_pipeline import StandardTrimmer

_REPO = Path(__file__).resolve().parents[2]


def test_trim_clears_previous_run_outputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo_rhone.yaml"))

    # New raw data: a single sample NEW (presence of the pair is enough for discovery).
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "NEW_R1.fastq.gz").write_bytes(b"")
    (raw / "NEW_R2.fastq.gz").write_bytes(b"")
    cfg.paths.raw_data = raw
    cfg.paths.output = tmp_path / "out"
    cfg.paths.logs = tmp_path / "out" / "logs"

    orch = PipelineOrchestrator(cfg)
    trim_dir = cfg.paths.output / "01_trim" / cfg.marker.name
    trim_dir.mkdir(parents=True, exist_ok=True)
    # Stale leftovers from a "previous run" with a different dataset.
    (trim_dir / "OLD.R1.fastq").write_text("stale")
    (trim_dir / "OLD.R2.fastq").write_text("stale")

    # Fake the trimmer so no external cutadapt is needed: it just writes the outputs.
    def fake_trim_sample(self, *, output_dir, sample_name, **_kw):
        output_dir.mkdir(parents=True, exist_ok=True)
        r1 = output_dir / f"{sample_name}.R1.fastq"
        r2 = output_dir / f"{sample_name}.R2.fastq"
        r1.write_text("x")
        r2.write_text("x")
        return (r1, r2)

    monkeypatch.setattr(StandardTrimmer, "trim_sample", fake_trim_sample)

    orch.run_trim()

    assert not (trim_dir / "OLD.R1.fastq").exists(), "stale previous-run output must be removed"
    assert not (trim_dir / "OLD.R2.fastq").exists(), "stale previous-run output must be removed"
    assert (trim_dir / "NEW.R1.fastq").exists(), "the current sample must be trimmed"
