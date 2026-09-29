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
    cfg = load_config(str(_REPO / "config" / "markers" / "teleo.yaml"))
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
    def fake(self, raw_reads_dir, library_name, tag_file, output_dir, work_dir, log_dir, **kw):
        calls.append(library_name)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        samples = [
            line[1:].strip() for line in Path(tag_file).read_text().splitlines()
            if line.startswith(">")
        ]
        for sample in samples:
            for r in ("R1", "R2"):
                (Path(output_dir) / f"{sample}.{r}.fastq.gz").write_text("")
        return samples

    return fake


def test_demux_uses_metadata_libraries_and_trim_reads_demux_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fwd = load_config(
        str(_REPO / "config" / "markers" / "teleo.yaml")
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
        "S1.R1.fastq.gz", "S1.R2.fastq.gz", "S2.R1.fastq.gz", "S2.R2.fastq.gz",
        "S3.R1.fastq.gz", "S3.R2.fastq.gz",
    ]
    # tag files are written once, beside (not inside) the per-library work dirs,
    # and only for this marker's libraries
    demux_dir = Path(outputs["demux_dir"])
    assert sorted(p.name for p in (demux_dir / "cutadapt_tags").iterdir()) == [
        "RUN_A_L001.fasta", "RUN_B_L002.fasta",
    ]
    assert not (demux_dir / "work").exists()
    # trim now discovers samples in the demux output, not in raw_data
    assert orch._trim_input_dir() == samples_dir
    assert orch._get_sample_list() == ["S1", "S2", "S3"]
    assert orch._find_read_file("S3", "R2") == samples_dir / "S3.R2.fastq.gz"


def test_demux_samples_removed_only_after_trim_completes(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fwd = load_config(
        str(_REPO / "config" / "markers" / "teleo.yaml")
    ).marker.primers.forward
    orch, cfg = _orch(tmp_path, [("S1", "LIB_A", fwd)])
    monkeypatch.setattr(LigationTrimmer, "process_library", _fake_process_library([]))
    orch.state.add_step("demultiplex")
    orch.state.add_step("trim")
    samples_dir = Path(orch.run_demultiplex()["samples_dir"])

    orch._drop_demux_samples()  # trim not completed: resume still needs them
    assert samples_dir.is_dir()

    orch.state.start_step("trim")
    orch.state.complete_step("trim", {"trimmed_dir": str(tmp_path)})
    orch._drop_demux_samples()
    assert not samples_dir.exists()
    assert cfg.paths.raw_data.is_dir()  # raw data is never touched


def test_demux_rejects_sample_name_shared_by_two_libraries(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fwd = load_config(
        str(_REPO / "config" / "markers" / "teleo.yaml")
    ).marker.primers.forward
    orch, _ = _orch(tmp_path, [("S1", "LIB_A", fwd), ("S1", "LIB_B", fwd)])
    monkeypatch.setattr(LigationTrimmer, "process_library", _fake_process_library([]))
    with pytest.raises(ValueError, match="once per marker"):
        orch.run_demultiplex()


def test_demux_rejects_sample_listed_twice_in_one_library(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fwd = load_config(
        str(_REPO / "config" / "markers" / "teleo.yaml")
    ).marker.primers.forward
    orch, _ = _orch(tmp_path, [("S1", "LIB_A", fwd), ("S1", "LIB_A", fwd)])
    monkeypatch.setattr(LigationTrimmer, "process_library", _fake_process_library([]))
    with pytest.raises(ValueError, match="once per marker"):
        orch.run_demultiplex()


def test_demux_keeps_only_this_markers_rows(tmp_path, monkeypatch):
    """The same eventID may carry one row per marker, across or within libraries."""
    monkeypatch.chdir(tmp_path)
    fwd = load_config(
        str(_REPO / "config" / "markers" / "teleo.yaml")
    ).marker.primers.forward
    other = "GTCGGTAAAACTCGTGCCAGC"  # mifish
    orch, _ = _orch(
        tmp_path,
        [
            ("S1", "LIB_A", fwd),
            ("S1", "LIB_B", other),  # S1 again, other marker, other library
            ("S2", "LIB_B", fwd),
            ("S2", "LIB_B", other),  # S2 twice in one library, one row per marker
            ("S3", "LIB_B", other),  # other marker only
        ],
    )
    calls = []
    monkeypatch.setattr(LigationTrimmer, "process_library", _fake_process_library(calls))

    outputs = orch.run_demultiplex()

    assert calls == ["LIB_A", "LIB_B"]
    tags = Path(outputs["demux_dir"]) / "cutadapt_tags"
    assert (tags / "LIB_A.fasta").read_text().count(">") == 1
    assert (tags / "LIB_B.fasta").read_text().splitlines()[0] == ">S2"
    assert (tags / "LIB_B.fasta").read_text().count(">") == 1
    assert sorted(p.name for p in Path(outputs["samples_dir"]).iterdir()) == [
        "S1.R1.fastq.gz", "S1.R2.fastq.gz", "S2.R1.fastq.gz", "S2.R2.fastq.gz",
    ]


def test_trim_input_is_raw_data_without_demux(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    orch, cfg = _orch(tmp_path, [])
    assert orch._trim_input_dir() == cfg.paths.raw_data
