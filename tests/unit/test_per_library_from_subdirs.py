"""DADA2-by-library can derive the sample->library grouping from per-library subfolders.

When `dada2.per_library` is on and no metadata carries a library grouping (no metadata, or a
field CSV set for the report only), and raw_data is organized one folder per sequencing
library/run (already-demultiplexed per-sample reads), the library map is derived from the
subfolder each sample lives in. A library column in the metadata wins over the folders, and
a lab CSV shared by several markers is filtered on pcr_primer_forward. Flat layouts, a single
library, or `per_library` off all fall back to the standard single-batch path (no map).
"""

import csv
from pathlib import Path

from seednap.config.loader import load_config
from seednap.pipeline.orchestrator import PipelineOrchestrator

_REPO = Path(__file__).resolve().parents[2]


def _orch(tmp_path, raw, per_library=True):
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo.yaml"))
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


def _read_map(path):
    return {r["sample"]: r["library"] for r in csv.DictReader(open(path))}


def test_field_csv_without_library_still_uses_subdirs(tmp_path, monkeypatch):
    # A field CSV set for the report only (no library column) must not hide the folders.
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw / "LIB_A", "S1")
    _pair(raw / "LIB_B", "S2")
    field = tmp_path / "field.csv"
    field.write_text("eventID,eventDate\nS1,2024-05-01\nS2,2024-05-02\n")
    orch = _orch(tmp_path, raw)
    orch.config.report.sample_metadata = field
    assert _read_map(orch._build_library_map()) == {"S1": "LIB_A", "S2": "LIB_B"}


def test_field_csv_library_column_wins_over_subdirs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = tmp_path / "raw"
    _pair(raw / "LIB_A", "S1")
    _pair(raw / "LIB_B", "S2")
    field = tmp_path / "field.csv"
    field.write_text("eventID,library\nS1,run1\nS2,run2\n")
    orch = _orch(tmp_path, raw)
    orch.config.report.sample_metadata = field
    assert _read_map(orch._build_library_map()) == {"S1": "run1", "S2": "run2"}


def test_lab_csv_is_filtered_by_marker_primer(tmp_path, monkeypatch):
    # S1 is in libA for this marker (teleo) and in libB for another marker: the teleo map
    # must keep libA, whatever the row order.
    monkeypatch.chdir(tmp_path)
    lab = tmp_path / "lab.csv"
    lab.write_text(
        "eventID,library,tag_demultiplex,pcr_primer_forward\n"
        "S1,libA,ACGT,ACACCGCCCGTCACTCT\n"
        "S2,libC,TTGG,ACACCGCCCGTCACTCT\n"
        "S1,libB,ACGT,GTCGGTAAAACTCGTGCCAGC\n"
        "S3,libB,GGCC,GTCGGTAAAACTCGTGCCAGC\n"
    )
    orch = _orch(tmp_path, tmp_path / "raw")
    orch.config.demultiplex.metadata = lab
    assert _read_map(orch._build_library_map()) == {"S1": "libA", "S2": "libC"}


def test_lab_csv_filtered_when_field_csv_is_primary(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    lab = tmp_path / "lab.csv"
    lab.write_text(
        "eventID,library,tag_demultiplex,pcr_primer_forward\n"
        "S1,libA,ACGT,ACACCGCCCGTCACTCT\n"
        "S2,libC,TTGG,ACACCGCCCGTCACTCT\n"
        "S1,libB,ACGT,GTCGGTAAAACTCGTGCCAGC\n"
    )
    field = tmp_path / "field.csv"
    field.write_text("eventID,eventDate\nS1,2024-05-01\nS2,2024-05-02\n")
    orch = _orch(tmp_path, tmp_path / "raw")
    orch.config.demultiplex.metadata = lab
    orch.config.report.sample_metadata = field
    assert _read_map(orch._build_library_map()) == {"S1": "libA", "S2": "libC"}
