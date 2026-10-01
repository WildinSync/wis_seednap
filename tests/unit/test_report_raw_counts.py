"""Regression tests for GitHub issue #3 (missing raw read counts -> % retained NA).

The report recovers raw/trimmed read counts from the per-sample Cutadapt logs, which
the trim step writes under ``<output>/01_trim/<marker>/logs``. A logs-dir mismatch
(the report reading ``<output>/logs`` while trim wrote elsewhere) left raw/trimmed
empty and ``% retained`` as NA. These tests pin both the builder behaviour and the
wiring so the mismatch cannot silently return.
"""

import inspect

import pandas as pd

from seednap.steps.report.read_tracking import ReadTrackingBuilder


def test_builder_recovers_raw_trimmed_and_retention_from_trim_logs(tmp_path):
    """Given Cutadapt logs in <output>/01_trim/<marker>/logs, raw/trimmed and the
    computed % retained are populated (not NA)."""
    marker = "m"
    logs = tmp_path / "01_trim" / marker / "logs"
    logs.mkdir(parents=True)
    (logs / "S1_trim_pass1.txt").write_text("Total read pairs processed:            1,000\n")
    (logs / "S1_trim_pass2.txt").write_text("Pairs written (passing filters):           900 (90.0%)\n")

    otu = tmp_path / "02_swarm" / marker / "otu_table.csv"
    otu.parent.mkdir(parents=True)
    otu.write_text("OTU_ID,S1\nOTU1,800\n")

    df = ReadTrackingBuilder(marker=marker, logs_dir=logs, swarm_otu_table=otu).build()
    row = df[df["sample"] == "S1"].iloc[0]

    assert int(row["raw"]) == 1000
    assert int(row["trimmed"]) == 900
    assert not pd.isna(row["pct_retained"])  # the bug made this NA
    assert round(float(row["pct_retained"]), 1) == 80.0  # clustered 800 / raw 1000


def test_report_reads_logs_from_the_trim_output_dir_not_run_root():
    """The orchestrator and the standalone `report` command must read the Cutadapt
    logs from <output>/01_trim/<marker>/logs (where the trim step writes them), not
    from a bare <output>/logs. Guards against re-introducing the issue-#3 mismatch."""
    import seednap.cli as cli
    import seednap.pipeline.orchestrator as orch

    needle = '"01_trim" / marker / "logs"'
    assert needle in inspect.getsource(orch), "orchestrator report logs_dir must use 01_trim/<marker>/logs"
    assert needle in inspect.getsource(cli), "report command logs_dir must use 01_trim/<marker>/logs"
    # Every orchestrator report (the table AND the HTML) must build through the one
    # helper; the HTML report once read <output>/logs on its own and showed raw NA.
    src = inspect.getsource(orch)
    assert '"logs_dir": out / "logs"' not in src
    assert src.count("ReadTrackingBuilder(**self._read_tracking_kwargs(") == 2


def test_step_output_dirs_keep_marker_case(tmp_path, monkeypatch):
    """A mixed-case marker (e.g. ``euka_V9``) must land in the same subdirectory the
    report reads. The DADA2/SWARM/taxonomy steps used to lowercase it (``euka_v9``),
    so the report found no track_reads.csv and every DADA2 column came out NA."""
    import seednap.steps.dada2.processor as dada2_mod
    import seednap.steps.swarm.processor as swarm_mod
    from seednap.steps.taxonomic_assignment.assigner import TaxonomicAssigner

    monkeypatch.setattr(dada2_mod, "Dada2Runner", lambda **_: None)
    for name in ("VsearchRunner", "SwarmClusterer"):
        monkeypatch.setattr(swarm_mod, name, lambda **_: None)

    marker = "euka_V9"
    trimmed = tmp_path / "01_trim" / marker
    trimmed.mkdir(parents=True)

    dada2 = dada2_mod.Dada2Processor(marker, trimmed, tmp_path)
    swarm = swarm_mod.SwarmProcessor(marker, trimmed, tmp_path)
    taxo = TaxonomicAssigner("blast", marker, tmp_path)

    assert dada2.output_dir == tmp_path / "02_dada2" / marker
    assert swarm.output_dir == tmp_path / "02_swarm" / marker
    assert taxo.taxo_dir == tmp_path / "03_taxo" / marker


def test_dada2_r_script_keeps_marker_case():
    """dada2_process.R builds 02_dada2/<marker>/ itself; it must not lowercase the marker."""
    from importlib.resources import files

    src = files("seednap").joinpath("scripts/dada2_process.R").read_text()
    assert "marker <- args[1]" in src
    assert "tolower(args[1])" not in src


def test_raw_counted_from_fastq_when_trim_log_missing(tmp_path):
    """A sample without a pass-1 Cutadapt log gets its raw count from its raw R1 FASTQ,
    so raw and % retained are not NA."""
    import gzip

    logs = tmp_path / "logs"
    logs.mkdir()
    raw_dir = tmp_path / "raw" / "lib1"  # per-library subdirectory layout
    raw_dir.mkdir(parents=True)
    with gzip.open(raw_dir / "S1_R1.fastq.gz", "wt") as fh:
        fh.write("@r\nACGT\n+\nIIII\n" * 1000)
    (raw_dir / "S10_R1.fastq").write_text("@r\nACGT\n+\nIIII\n" * 7)
    d2 = tmp_path / "02_dada2"
    d2.mkdir()
    pd.DataFrame([{"sample": "S1", "input": 900, "filtered": 800, "denoised": 700,
                   "merged": 600, "nonchim": 500}]).to_csv(d2 / "track_reads.csv", index=False)

    b = ReadTrackingBuilder("m", logs_dir=logs, dada2_dir=d2, raw_dir=tmp_path / "raw")
    row = b.build().set_index("sample").loc["S1"]
    assert int(row["raw"]) == 1000  # S10's file must not be counted for S1
    assert float(row["pct_retained"]) == 50.0


def test_funnel_percentages_use_a_fully_measured_base():
    """Figure 1: with raw NA the percentages must not be computed against 1 (which
    printed read counts x 100 as percentages); they fall back to the first step
    measured for every sample."""
    from seednap.steps.report.html_report import HTMLReportBuilder

    df = pd.DataFrame({"sample": ["A", "B"], "raw": [pd.NA, 1000], "trimmed": [800, 900],
                       "clustered": [400, 450], "pct_retained": [pd.NA, 45.0]})
    rb = HTMLReportBuilder("m", df, steps=["raw", "trimmed", "clustered"])
    assert rb._funnel_base() == "trimmed"
    df["raw"] = [pd.NA, pd.NA]
    assert HTMLReportBuilder("m", df, steps=["raw", "trimmed", "clustered"])._funnel_base() == "trimmed"
    df["raw"] = [1000, 1000]
    assert HTMLReportBuilder("m", df, steps=["raw", "trimmed", "clustered"])._funnel_base() == "raw"
