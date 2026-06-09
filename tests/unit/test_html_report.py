"""HTMLReportBuilder correctness regressions.

Three confirmed bugs, each a behaviour that was wrong before the fix:

1. _section_taxonomy divided by len(tax) with no guard, so a taxonomy CSV with a
   header but zero feature rows (all-chimera / non-amplifying run) raised
   ZeroDivisionError and the whole report was lost. It must now render an
   explanatory empty-run note instead.
2. _richness's full-OTU-table fallback omitted the SWARM metadata columns
   cloud/amplicon/abundance from its exclusion set, so those per-OTU metadata
   columns were counted as phantom per-sample read-count columns. They must be
   excluded; only real sample columns survive.
3. Control detection was hardcoded to the literal 'blank' name prefix, silently
   counting legacy controls (CNEG/CEXT/water/etc.) as biological samples. It now
   uses config.manifest.classify_control.
"""

import pandas as pd
import pytest

from seednap.steps.report.html_report import HTMLReportBuilder, _is_negative_control


def _tracking_df(samples):
    """Minimal per-sample read-tracking frame for the report builder."""
    return pd.DataFrame(
        {
            "sample": samples,
            "raw": [100] * len(samples),
            "nonchim": [80] * len(samples),
            "pct_retained": [80.0] * len(samples),
        }
    )


def _builder(df, **kwargs):
    return HTMLReportBuilder(marker="teleo", tracking_df=df, **kwargs)


def test_section_taxonomy_zero_rows_does_not_divide_by_zero(tmp_path):
    """A header-only taxonomy CSV (0 feature rows) must render a note, not crash."""
    tax = tmp_path / "tax.csv"
    # Header present, zero data rows -> pd.read_csv returns a 0-row frame.
    pd.DataFrame(columns=["ASV_ID", "species", "genus", "sampleA"]).to_csv(tax, index=False)
    b = _builder(_tracking_df(["sampleA"]), taxonomy_csv=tax)
    html = b._section_taxonomy({})
    assert html is not None
    assert "No features survived" in html
    # And the full render must not raise.
    assert "<html" in b.render().lower()


def test_richness_excludes_swarm_metadata_columns(tmp_path):
    """cloud/amplicon/abundance are per-OTU metadata, never phantom samples."""
    otu = tmp_path / "otu_table_full.csv"
    pd.DataFrame(
        {
            "OTU": ["OTU_1", "OTU_2"],
            "total": [30, 10],
            "cloud": [5, 2],
            "amplicon": ["seedA", "seedB"],
            "length": [120, 121],
            "abundance": [25, 8],
            "chimera": ["N", "N"],
            "spread": [2, 1],
            "sequence": ["ACGT", "TTGC"],
            "sampleA": [20, 5],
            "sampleB": [10, 5],
        }
    ).to_csv(otu, index=False)
    # No taxonomy CSV, so _richness falls back to the full OTU table.
    b = _builder(_tracking_df(["sampleA", "sampleB"]), otu_table_full=otu)
    rich = b._richness()
    assert rich is not None
    # Only the two real sample columns appear; the metadata columns must not.
    assert set(rich.index) == {"sampleA", "sampleB"}
    for phantom in ("cloud", "amplicon", "abundance", "total", "spread"):
        assert phantom not in rich.index


def test_negative_control_detection_recognises_legacy_names():
    """Legacy control conventions are detected, not just the 'blank' prefix."""
    for ctrl in ("Blank-ext-1", "CNEG01", "CEXT_3", "EXT_NC", "PCR_NC", "water control"):
        assert _is_negative_control(ctrl) is True
    # Positive controls and real samples are not negative controls.
    assert _is_negative_control("CPOS_1") is False
    assert _is_negative_control("Site_A_rep1") is False


def test_n_controls_counts_legacy_controls():
    """_n_controls must count legacy-named controls a 'blank' prefix would miss."""
    df = _tracking_df(["Site_A", "Site_B", "CNEG01", "EXT_NC"])
    b = _builder(df)
    assert b._n_controls() == 2


def test_n_controls_warns_when_zero_controls(capsys, caplog):
    """A non-empty dataset with no detectable controls must emit a [WARN]."""
    df = _tracking_df(["Site_A", "Site_B"])
    b = _builder(df)
    import logging

    with caplog.at_level(logging.WARNING):
        assert b._n_controls() == 0
    assert any("0" in r.getMessage() and "control" in r.getMessage().lower()
               for r in caplog.records)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
