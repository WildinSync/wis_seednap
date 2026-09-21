"""Ligation demultiplex processes the metadata libraries and feeds trim its output.

Regression for two bugs: the demultiplex step globbed raw files by the marker name
(``teleo*_R1.fastq.gz``) instead of the metadata ``library`` values, and the trim
step kept reading ``paths.raw_data`` (the multiplexed library files) after
demultiplex had run.
"""

from pathlib import Path

import pytest

from seednap.config.loader import load_config
from seednap.pipeline.orchestrator import PipelineOrchestrator
from seednap.steps.trimming.trimming_pipeline import LigationTrimmer

_REPO = Path(__file__).resolve().parents[2]


def _orch(tmp_path, metadata_rows):
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo_rhone.yaml"))
    raw = tmp_path / "raw"
    raw.mkdir()
    meta = tmp_path / "meta.csv"
    meta.write_text(
        "eventID,tag_demultiplex,library,pcr_primer_forward\n"
        + "".join(f"{e},aacaagcc,{lib},{p}\n" for e, lib, p in metadata_rows)
    )
    cfg.paths.raw_data = raw
    cfg.paths.output = tmp_path / "out"
    cfg.paths.logs = tmp_path / "out" / "logs"
    cfg.demultiplex.protocol = "ligation"
    cfg.demultiplex.metadata = meta
    return PipelineOrchestrator(cfg), cfg


def _fake_process_library(calls):
    def fake(self, raw_reads_dir, library_name, metadata_csv, output_base_dir, **kw):
        calls.append(library_name)
        realigned = Path(output_base_dir) / "00_demultiplex_ligation" / "realigned"
        realigned.mkdir(parents=True)
        import csv

        with open(metadata_csv) as fh:
            for row in csv.DictReader(fh):
                if row["library"] == library_name:
                    for r in ("R1", "R2"):
                        (realigned / f"{row['eventID']}.{r}.fastq").write_text("")
        return realigned

    return fake


def test_demux_uses_metadata_libraries_and_trim_reads_demux_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fwd = load_config(
        str(_REPO / "config" / "markers" / "teleo_rhone.yaml")
    ).marker.primers.forward
    orch, cfg = _orch(
        tmp_path,
        [
            ("S1", "RUN_A_L001", fwd),
            ("S2", "RUN_A_L001", fwd),
            ("S3", "RUN_B_L002", fwd),
            ("X1", "OTHER_MARKER_LIB", "GGGGGGGGGG"),
        ],
    )
    calls = []
    monkeypatch.setattr(LigationTrimmer, "process_library", _fake_process_library(calls))

    outputs = orch.run_demultiplex()

    assert calls == ["RUN_A_L001", "RUN_B_L002"]  # library names, not marker.name
    samples_dir = Path(outputs["samples_dir"])
    assert sorted(p.name for p in samples_dir.iterdir()) == [
        "S1.R1.fastq", "S1.R2.fastq", "S2.R1.fastq", "S2.R2.fastq",
        "S3.R1.fastq", "S3.R2.fastq",
    ]
    # trim now discovers samples in the demux output, not in raw_data
    assert orch._trim_input_dir() == samples_dir
    assert orch._get_sample_list() == ["S1", "S2", "S3"]
    assert orch._find_read_file("S3", "R2") == samples_dir / "S3.R2.fastq"


def test_demux_rejects_sample_name_shared_by_two_libraries(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fwd = load_config(
        str(_REPO / "config" / "markers" / "teleo_rhone.yaml")
    ).marker.primers.forward
    orch, _ = _orch(tmp_path, [("S1", "LIB_A", fwd), ("S1", "LIB_B", fwd)])
    monkeypatch.setattr(LigationTrimmer, "process_library", _fake_process_library([]))
    with pytest.raises(ValueError, match="unique across libraries"):
        orch.run_demultiplex()


def test_trim_input_is_raw_data_without_demux(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    orch, cfg = _orch(tmp_path, [])
    assert orch._trim_input_dir() == cfg.paths.raw_data
