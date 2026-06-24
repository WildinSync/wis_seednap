"""DADA2-by-library can derive the sample->library grouping from per-library subfolders.

When `dada2.per_library` is on but no metadata grouping is configured, and raw_data is
organized one folder per sequencing library/run (already-demultiplexed per-sample reads),
the library map is derived from the subfolder each sample lives in, so per_library works
without any lab metadata. Flat layouts, a single library, or `per_library` off all fall
back to the standard single-batch path (no map).
"""

import csv
from pathlib import Path

from seednap.config.loader import load_config
from seednap.pipeline.orchestrator import PipelineOrchestrator

_REPO = Path(__file__).resolve().parents[2]


def _orch(tmp_path, raw, per_library=True):
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo_rhone.yaml"))
    cfg.dada2.per_library = per_library
    cfg.report.sample_metadata = None
    cfg.demultiplex.metadata = None
    cfg.paths.raw_data = raw
    cfg.paths.output = tmp_path / "out"
    cfg.paths.logs = tmp_path / "out" / "logs"
    return PipelineOrchestrator(cfg)


def _pair(d: Path, sample: str) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{sample}_R1.fastq.gz").write_bytes(b"")
    (d / f"{sample}_R2.fastq.gz").write_bytes(b"")


def test_library_map_derived_from_subdirs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw / "LIB_A", "S1")
    _pair(raw / "LIB_A", "S2")
    _pair(raw / "LIB_B", "S3")
    out = _orch(tmp_path, raw)._build_library_map()
    assert out is not None and out.exists()
    mapping = {r["sample"]: r["library"] for r in csv.DictReader(open(out))}
    assert mapping == {"S1": "LIB_A", "S2": "LIB_A", "S3": "LIB_B"}


def test_flat_layout_has_no_library_map(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw, "S1")
    _pair(raw, "S2")
    assert _orch(tmp_path, raw)._build_library_map() is None


def test_single_subdir_is_single_batch(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw / "LIB_A", "S1")
    _pair(raw / "LIB_A", "S2")
    assert _orch(tmp_path, raw)._build_library_map() is None


def test_per_library_off_returns_none(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw / "LIB_A", "S1")
    _pair(raw / "LIB_B", "S2")
    assert _orch(tmp_path, raw, per_library=False)._build_library_map() is None
