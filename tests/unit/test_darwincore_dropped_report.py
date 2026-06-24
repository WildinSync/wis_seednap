"""Increment 3: the DarwinCore deleted-entries report.

The 'darwincore' step records every occurrence the control + non-target filters remove, with
the reason, in a `<output>_dropped.csv` alongside the main output (for early-stage QA). The
report is written even when nothing is dropped (header only), so its presence confirms the
filters ran.
"""

import pandas as pd

from seednap.steps.formatting.darwincore_builder import DarwinCoreBuilder


def _builder(tmp_path):
    b = DarwinCoreBuilder.__new__(DarwinCoreBuilder)  # bypass __init__ (no real inputs needed)
    b.output_path = tmp_path / "teleo_blast_darwincore.csv"
    return b


def test_dropped_report_lists_removed_rows_with_reasons(tmp_path):
    b = _builder(tmp_path)
    snapshot = pd.DataFrame({
        "_dwc_row_id": [0, 1, 2],
        "eventID": ["Blank-PCR-1", "DAR-2025-0100", "DAR-2025-0101"],
        "taxon": ["Escherichia coli", "Salmo trutta", "Homo sapiens"],
        "nb_reads": [12, 200, 30],
    })
    dropped_reason = {0: "negative/positive control", 2: "non-target taxon (teleo)"}
    b._write_dropped_report(snapshot, dropped_reason)

    report = tmp_path / "teleo_blast_darwincore_dropped.csv"
    assert report.exists() and b.dropped_report_path == report
    out = pd.read_csv(report)
    assert len(out) == 2
    assert dict(zip(out["eventID"], out["drop_reason"])) == {
        "Blank-PCR-1": "negative/positive control",
        "DAR-2025-0101": "non-target taxon (teleo)",
    }
    # the kept row is absent
    assert "DAR-2025-0100" not in set(out["eventID"])


def test_dropped_report_written_even_when_nothing_dropped(tmp_path):
    b = _builder(tmp_path)
    snapshot = pd.DataFrame({
        "_dwc_row_id": [0, 1],
        "eventID": ["DAR-2025-0100", "DAR-2025-0101"],
        "taxon": ["Salmo trutta", "Perca fluviatilis"],
        "nb_reads": [200, 150],
    })
    b._write_dropped_report(snapshot, {})
    out = pd.read_csv(tmp_path / "teleo_blast_darwincore_dropped.csv")
    assert len(out) == 0
    assert "drop_reason" in out.columns
