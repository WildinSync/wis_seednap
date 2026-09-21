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
