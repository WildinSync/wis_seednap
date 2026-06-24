"""Sample discovery handles raw data organised into per-library subdirectories.

Already-demultiplexed runs are often kept one folder per sequencing library/run
rather than as one flat directory. Discovery searches the top level first, then
falls back to a recursive search of subdirectories, so these layouts work without
the user flattening them. Flat layouts are unchanged; a sample name that resolves
to several files across subfolders is ambiguous and raises rather than guessing.
"""

from pathlib import Path

import pytest

from seednap.config.loader import load_config
from seednap.pipeline.orchestrator import PipelineOrchestrator

_REPO = Path(__file__).resolve().parents[2]


def _orch(tmp_path, raw):
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo_rhone.yaml"))
    cfg.paths.raw_data = raw
    cfg.paths.output = tmp_path / "out"
    cfg.paths.logs = tmp_path / "out" / "logs"
    return PipelineOrchestrator(cfg)


def _pair(d: Path, sample: str) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{sample}_R1.fastq.gz").write_bytes(b"")
    (d / f"{sample}_R2.fastq.gz").write_bytes(b"")


def test_discovers_samples_in_per_library_subdirectories(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw / "LIB_A", "S1")
    _pair(raw / "LIB_A", "S2")
    _pair(raw / "LIB_B", "S3")
    orch = _orch(tmp_path, raw)
    assert orch._get_sample_list() == ["S1", "S2", "S3"]
    assert orch._find_read_file("S2", "R1") == raw / "LIB_A" / "S2_R1.fastq.gz"
    assert orch._find_read_file("S3", "R2") == raw / "LIB_B" / "S3_R2.fastq.gz"


def test_flat_layout_unchanged(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw, "S1")
    orch = _orch(tmp_path, raw)
    assert orch._get_sample_list() == ["S1"]
    assert orch._find_read_file("S1", "R1") == raw / "S1_R1.fastq.gz"


def test_ambiguous_sample_across_subdirs_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw / "LIB_A", "DUP")
    _pair(raw / "LIB_B", "DUP")
    orch = _orch(tmp_path, raw)
    with pytest.raises(FileNotFoundError, match="Ambiguous"):
        orch._find_read_file("DUP", "R1")


def test_no_fastqs_anywhere_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    (raw / "emptylib").mkdir(parents=True)
    orch = _orch(tmp_path, raw)
    with pytest.raises(FileNotFoundError, match="No forward-read FASTQ"):
        orch._get_sample_list()
