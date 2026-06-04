"""Unit tests for the read/sequence tracking builder + report config."""

import pandas as pd
import pytest

from seednap.config.models import ReportConfig
from seednap.steps.report.read_tracking import DADA2_STEPS, ReadTrackingBuilder


def _write_trim_logs(logs_dir, sample, raw, trimmed):
    """Write minimal Cutadapt pass1/pass2 logs in the real format."""
    (logs_dir / f"{sample}_trim_pass1.txt").write_text(
        f"Total read pairs processed:            {raw:,}\n"
        f"Pairs written (passing filters):       {raw:,} (100.0%)\n"
    )
    (logs_dir / f"{sample}_trim_pass2.txt").write_text(
        f"Total read pairs processed:            {raw:,}\n"
        f"Pairs written (passing filters):       {trimmed:,} (50.0%)\n"
    )


def _write_track(dada2_dir, rows):
    pd.DataFrame(rows).to_csv(dada2_dir / "track_reads.csv", index=False)


def test_cutadapt_parsing(tmp_path):
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=800)
    b = ReadTrackingBuilder("m", logs_dir=logs)
    df = b.build()
    row = df[df["sample"] == "S1"].iloc[0]
    assert row["raw"] == 1000
    assert row["trimmed"] == 800


def test_full_assembly_and_retention(tmp_path):
    logs = tmp_path / "logs"; logs.mkdir()
    d2 = tmp_path / "02_dada2"; d2.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    _write_track(d2, [dict(sample="S1", input=900, filtered=800,
                           denoised=700, merged=600, nonchim=550)])
    df = ReadTrackingBuilder("m", logs_dir=logs, dada2_dir=d2).build()
    r = df[df["sample"] == "S1"].iloc[0]
    assert [int(r[s]) for s in DADA2_STEPS] == [1000, 900, 800, 700, 600, 550]
    assert r["pct_retained"] == pytest.approx(55.0)


def test_absent_counts_are_NA_not_zero(tmp_path):
    """A sample with no DADA2 track entry must be NA (not a silent 0)."""
    logs = tmp_path / "logs"; logs.mkdir()
    d2 = tmp_path / "02_dada2"; d2.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    # S1 is in the logs but absent from the DADA2 track (only OTHER is present)
    _write_track(d2, [dict(sample="OTHER", input=10, filtered=9,
                           denoised=8, merged=7, nonchim=6)])
    b = ReadTrackingBuilder("m", logs_dir=logs, dada2_dir=d2)
    df = b.build()
    r = df[df["sample"] == "S1"].iloc[0]
    assert pd.isna(r["filtered"]) and pd.isna(r["nonchim"])
    msgs = b.warnings(df)
    assert any("absent (not measured)" in m for m in msgs)


def test_low_retention_and_step_loss_warnings(tmp_path):
    logs = tmp_path / "logs"; logs.mkdir()
    d2 = tmp_path / "02_dada2"; d2.mkdir()
    # 1000 raw -> only 5 nonchim = 0.5% retention; trimmed->filtered drops 99%
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    _write_track(d2, [dict(sample="S1", input=900, filtered=9,
                           denoised=7, merged=6, nonchim=5)])
    b = ReadTrackingBuilder("m", logs_dir=logs, dada2_dir=d2,
                            warn_below_retention_pct=30.0, warn_step_loss_pct=70.0)
    msgs = b.warnings(b.build())
    assert any("low overall retention" in m for m in msgs)
    assert any("trimmed->filtered dropped" in m for m in msgs)


def test_write_artifacts(tmp_path):
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    out = ReadTrackingBuilder("m", logs_dir=logs).write(tmp_path / "04_report")
    assert out["read_tracking_csv"].exists()
    assert out["read_tracking_txt"].exists()
    assert "Read tracking" in out["read_tracking_txt"].read_text()


def test_swarm_chain(tmp_path):
    """SWARM path reports raw -> trimmed -> clustered from otu_table column sums."""
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    _write_trim_logs(logs, "S2", raw=2000, trimmed=1800)
    otu = tmp_path / "otu_table.csv"
    pd.DataFrame({
        "sequence": ["AAA", "CCC"],
        "S1": [300, 100],   # clustered total = 400
        "S2": [800, 200],   # clustered total = 1000
    }).to_csv(otu, index=False)
    b = ReadTrackingBuilder("m", logs_dir=logs, swarm_otu_table=otu)
    assert b.steps == ["raw", "trimmed", "clustered"]
    df = b.build().set_index("sample")
    assert int(df.loc["S1", "clustered"]) == 400
    assert int(df.loc["S2", "clustered"]) == 1000
    assert df.loc["S2", "pct_retained"] == pytest.approx(50.0)  # 1000/2000


def test_html_report_is_self_contained(tmp_path):
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    d2 = tmp_path / "02_dada2"; d2.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    _write_track(d2, [dict(sample="S1", input=900, filtered=800,
                           denoised=700, merged=600, nonchim=550)])
    b = ReadTrackingBuilder("m", logs_dir=logs, dada2_dir=d2)
    df = b.build()
    out = HTMLReportBuilder("m", df, warnings=b.warnings(df),
                            summary={"n_features": 5}).write(tmp_path / "r.html")
    html = out.read_text()
    assert "<html" in html and "</html>" in html
    assert "http://" not in html and "https://" not in html  # self-contained
    assert "data:image/png;base64," in html  # embedded charts
    assert "S1" in html


def test_report_config_defaults_and_strictness():
    c = ReportConfig()
    assert c.read_tracking is True and c.html_report is False
    assert c.warn_below_retention_pct == 30.0 and c.warn_step_loss_pct == 70.0
    with pytest.raises(Exception):
        ReportConfig(unknown_field=1)  # extra="forbid"
