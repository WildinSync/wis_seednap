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
