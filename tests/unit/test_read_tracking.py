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


def test_step_summary_swarm(tmp_path):
    """SWARM step summary: total reads per step + the OTU count at 'clustered'."""
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    _write_trim_logs(logs, "S2", raw=2000, trimmed=1800)
    otu = tmp_path / "otu_table.csv"
    pd.DataFrame({"sequence": ["AAA", "CCC", "GGG"],
                  "S1": [300, 100, 50], "S2": [800, 200, 100]}).to_csv(otu, index=False)
    ss = ReadTrackingBuilder("m", logs_dir=logs, swarm_otu_table=otu).step_summary().set_index("step")
    assert int(ss.loc["raw", "total_reads"]) == 3000
    assert int(ss.loc["trimmed", "total_reads"]) == 2700
    assert int(ss.loc["clustered", "total_reads"]) == 1550   # 450 + 1100
    assert int(ss.loc["clustered", "n_features"]) == 3       # 3 OTUs (rows)
    assert pd.isna(ss.loc["raw", "n_features"])              # no features at read-level steps
    assert pd.isna(ss.loc["trimmed", "n_features"])


def test_step_summary_dada2(tmp_path):
    """DADA2 step summary: ASV counts at merged/nonchim from feature_counts.csv, NA at the
    read-level steps, total reads = the per-sample column sums."""
    logs = tmp_path / "logs"; logs.mkdir()
    d2 = tmp_path / "02_dada2"; d2.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    _write_trim_logs(logs, "S2", raw=1000, trimmed=800)
    _write_track(d2, [
        {"sample": "S1", "input": 900, "filtered": 800, "denoised": 780, "merged": 700, "nonchim": 690},
        {"sample": "S2", "input": 800, "filtered": 700, "denoised": 690, "merged": 600, "nonchim": 590},
    ])
    pd.DataFrame({"step": ["merged", "nonchim"], "n_features": [120, 110]}).to_csv(
        d2 / "feature_counts.csv", index=False)
    ss = ReadTrackingBuilder("m", logs_dir=logs, dada2_dir=d2).step_summary().set_index("step")
    assert int(ss.loc["filtered", "total_reads"]) == 1500    # 800 + 700
    assert int(ss.loc["nonchim", "total_reads"]) == 1280     # 690 + 590
    assert int(ss.loc["merged", "n_features"]) == 120
    assert int(ss.loc["nonchim", "n_features"]) == 110
    for step in ("raw", "trimmed", "filtered", "denoised"):
        assert pd.isna(ss.loc[step, "n_features"])           # ASVs only exist from the merge stage


def test_step_summary_missing_feature_counts_is_na_not_guessed(tmp_path):
    """A DADA2 run without feature_counts.csv yields NA ASV counts (never a guessed value),
    while reads are still tracked. The builder also emits a [WARN] for the absent
    feature_counts.csv at runtime (see _feature_counts); this test asserts only the
    NA/reads outcome, not the log line."""
    logs = tmp_path / "logs"; logs.mkdir()
    d2 = tmp_path / "02_dada2"; d2.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    _write_track(d2, [{"sample": "S1", "input": 900, "filtered": 800,
                       "denoised": 780, "merged": 700, "nonchim": 690}])
    ss = ReadTrackingBuilder("m", logs_dir=logs, dada2_dir=d2).step_summary().set_index("step")
    assert pd.isna(ss.loc["merged", "n_features"]) and pd.isna(ss.loc["nonchim", "n_features"])
    assert int(ss.loc["nonchim", "total_reads"]) == 690      # reads are still reported


def test_step_summary_partial_na_sums_measured_samples(tmp_path):
    """A step measured for some samples but NA for others sums the measured ones (it is not
    blanked to NA, which would discard a usable run total over one dropped sample). The
    builder also emits a [WARN] naming the unmeasured samples at runtime (see step_summary);
    this test asserts only the summed total, not the log line."""
    logs = tmp_path / "logs"; logs.mkdir()
    d2 = tmp_path / "02_dada2"; d2.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    _write_trim_logs(logs, "S2", raw=1000, trimmed=800)   # in the logs ...
    # ... but only S1 reaches the DADA2 track, so S2 is NA from 'filtered' on.
    _write_track(d2, [{"sample": "S1", "input": 900, "filtered": 800,
                       "denoised": 780, "merged": 700, "nonchim": 690}])
    ss = ReadTrackingBuilder("m", logs_dir=logs, dada2_dir=d2).step_summary().set_index("step")
    assert int(ss.loc["raw", "total_reads"]) == 2000        # both samples measured
    assert int(ss.loc["trimmed", "total_reads"]) == 1700    # both measured
    assert int(ss.loc["nonchim", "total_reads"]) == 690     # only S1 measured -> summed, not NA
    assert pd.isna(ss.loc["nonchim", "n_features"])         # no feature_counts.csv here


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


def test_html_report_rich_sections(tmp_path):
    """With taxonomy + otu_full + state, the report gains taxonomy/contamination/timeline."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "DAR-1", raw=1000, trimmed=900)
    otu = tmp_path / "otu_table.csv"
    pd.DataFrame({"sequence": ["A"], "DAR-1": [800], "Blank-PCR-1": [2]}).to_csv(otu, index=False)
    b = ReadTrackingBuilder("m", logs_dir=logs, swarm_otu_table=otu)
    df = b.build()
    # taxonomy CSV (blast schema): ranks + is_contaminant_candidate + sample cols + Sequence
    tax = tmp_path / "m_blast.csv"
    pd.DataFrame({
        "ASV_ID": ["OTU_1", "OTU_2"], "pident": [100.0, ""],
        "kingdom": ["Metazoa", "Unassigned"], "phylum": ["Chordata", "Unassigned"],
        "class": ["Actinopteri", "Unassigned"], "order": ["Perciformes", "Unassigned"],
        "family": ["Percidae", "Unassigned"], "genus": ["Perca", "Unassigned"],
        "species": ["Perca_fluviatilis", "Unassigned"],
        "is_contaminant_candidate": [False, False],
        "DAR-1": [800, 50], "Blank-PCR-1": [2, 0], "Sequence": ["A", "C"],
    }).to_csv(tax, index=False)
    otu_full = tmp_path / "otu_table_full.csv"
    pd.DataFrame({"OTU": ["1", "2", "3"], "total": [800, 50, 5], "length": [82, 83, 200],
                  "chimera": ["N", "N", "Y"], "spread": [1, 1, 1]}).to_csv(otu_full, index=False)
    state = {"steps": {"trim": {"status": "completed", "duration_seconds": 10.0},
                       "swarm": {"status": "completed", "duration_seconds": 5.0}}}
    html = HTMLReportBuilder("m", df, steps=b.steps, state=state,
                             taxonomy_csv=tax, otu_table_full=otu_full,
                             warnings=b.warnings(df, log=False)).render()
    assert "Run provenance" in html and "completed" in html  # timeline table
    assert "Top species" in html and "Perca" in html
    assert "contamination" in html.lower()  # Blank-PCR-1 had reads
    assert "data:image/png;base64," in html
    # LaTeX-paper theme: sea-green accent, no neon dashboard colors
    assert "#2e8b57" in html and "#44f187" not in html and "#38bdf8" not in html
    assert "Figure 1." in html and "Table 1." in html  # numbered figs/tables


def test_html_report_degrades_without_extra_sources(tmp_path):
    """No taxonomy/state -> still renders (read-tracking only), no crash."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    b = ReadTrackingBuilder("m", logs_dir=logs)
    html = HTMLReportBuilder("m", b.build(), steps=b.steps).render()
    assert "<html" in html and "Read tracking" in html
    assert "Top species" not in html  # no taxonomy section


def test_html_report_scales_to_many_samples(tmp_path):
    """>50 samples must not produce a metres-tall chart; report renders fine."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    otu = {"sequence": ["A", "C"]}
    for i in range(120):
        s = f"S{i:03d}"
        _write_trim_logs(logs, s, raw=1000, trimmed=900)
        otu[s] = [500, 300]
    pd.DataFrame(otu).to_csv(tmp_path / "otu.csv", index=False)
    b = ReadTrackingBuilder("m", logs_dir=logs, swarm_otu_table=tmp_path / "otu.csv")
    df = b.build()
    assert len(df) == 120
    html = HTMLReportBuilder("m", df, steps=b.steps).render()
    # report stays small (histogram, not a 120-bar chart) and renders all rows
    assert len(html) < 1_500_000
    assert "Retention distribution" in html or "data:image/png;base64," in html


def test_dataset_provenance_section(tmp_path):
    """Dataset section surfaces location/marker/recorder from field+project metadata."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    field = tmp_path / "field.csv"
    pd.DataFrame({"eventID": ["S1", "S2", "Blank-ext-1"],
                  "decimalLatitude": ["46.30", "46.40", "NA"],
                  "decimalLongitude": ["7.17", "7.20", "NA"],
                  "eventDate": ["2025.08.18", "2025.08.19", "NA"],
                  "institution": ["ELE", "ELE", "ELE"], "body": ["river", "river", "NA"]}).to_csv(field, index=False)
    proj = tmp_path / "proj.csv"
    pd.DataFrame({"marker": ["teleo"], "recordedby": ["Ada L."], "seqmet": ["MiSeq"],
                  "otu_db": ["MIDORI2"]}).to_csv(proj, index=False)
    b = ReadTrackingBuilder("teleo_rhone", logs_dir=logs)
    html = HTMLReportBuilder("teleo_rhone", b.build(), steps=b.steps,
                             field_metadata_csv=field, project_metadata_csv=proj,
                             summary={"provenance": {"dataset_name": "teleo_rhone", "marker": "teleo_rhone"}}).render()
    assert "Dataset</h2>" in html
    assert "46.30" in html and "7.17" in html        # location (controls excluded)
    assert "Ada L." in html and "MiSeq" in html       # project provenance
    assert ">teleo<" in html                          # marker = project's "teleo", not "teleo_rhone"
    assert "2025.08.18" in html


def test_dataset_section_no_metadata_is_explicit(tmp_path):
    """No metadata -> explicit note, never a silent omission."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    b = ReadTrackingBuilder("m", logs_dir=logs)
    html = HTMLReportBuilder("m", b.build(), steps=b.steps).render()
    assert "Dataset</h2>" in html
    # absence is stated explicitly, never silently skipped
    assert "were not provided" in html


def _write_run_log(path, lines):
    """Write a file-logger-format run log (TIME | LEVEL | name:lineno | msg)."""
    rows = [f"2026-06-04 12:00:{i:02d} | {lvl:8s} | seednap.x:{i} | {msg}"
            for i, (lvl, msg) in enumerate(lines)]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_run_log_section_colorized(tmp_path):
    """The run-log section embeds the transcript with rich's exact level colors."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    run_log = tmp_path / "m_pipeline_run.log"
    _write_run_log(run_log, [
        ("INFO", "Starting seednap pipeline for marker: m"),
        ("WARNING", "Sample DAR-1: merged file is empty, skipping <x>"),
        ("ERROR", "Step 'taxonomy' failed"),
        ("INFO", "Pipeline complete"),
    ])
    b = ReadTrackingBuilder("m", logs_dir=logs)
    html = HTMLReportBuilder("m", b.build(), steps=b.steps, log_file=run_log).render()
    assert "Run log</h2>" in html        # the Run-log panel heading
    assert '<pre class="runlog">' in html
    assert '<div class="terminal">' in html and "term-dot" in html  # real terminal window
    # bright ANSI palette tuned for the dark terminal: info blue, warning amber, error red
    assert "#6cb6ff" in html and "#e3b341" in html and "#ff6b6b" in html
    assert "http://" not in html and "https://" not in html  # still self-contained
    assert "&lt;x&gt;" in html  # message HTML is escaped, never injected raw


def test_run_log_missing_is_explicit(tmp_path):
    """A missing or unprovided run log is stated, never silently omitted (section 4)."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    b = ReadTrackingBuilder("m", logs_dir=logs)
    df = b.build()
    # no log_file passed
    html_none = HTMLReportBuilder("m", df, steps=b.steps).render()
    assert "Run log</h2>" in html_none and "not embedded here" in html_none
    # log_file pointing at a nonexistent file
    html_missing = HTMLReportBuilder("m", df, steps=b.steps,
                                     log_file=tmp_path / "nope.log").render()
    assert "was not found" in html_missing


def test_run_log_truncation_keeps_all_events(tmp_path):
    """A long log is truncated but keeps every warning/error and marks omissions."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    lines = [("INFO", f"routine line {i}") for i in range(500)]
    lines[100] = ("WARNING", "important warning A")
    lines[300] = ("ERROR", "important error B")
    run_log = tmp_path / "m_pipeline_run.log"
    _write_run_log(run_log, lines)
    b = ReadTrackingBuilder("m", logs_dir=logs)
    builder = HTMLReportBuilder("m", b.build(), steps=b.steps, log_file=run_log, max_log_lines=120)
    items, truncated, total = builder._select_log_lines(run_log.read_text().splitlines())
    assert truncated and total == 500
    kept = "\n".join(t for _, t in items)
    assert "important warning A" in kept and "important error B" in kept
    assert any(i == -1 and "omitted" in t for i, t in items)  # explicit markers


def test_read_tracking_warnings_have_header(tmp_path):
    """Warnings in the report carry a header + explanation, not a bare dump."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "Blank-PCR-1", raw=10000, trimmed=50)  # heavy loss -> warns
    b = ReadTrackingBuilder("m", logs_dir=logs)
    df = b.build()
    warns = b.warnings(df, log=False)
    assert warns  # sanity: this sample triggers warnings
    html = HTMLReportBuilder("m", df, steps=b.steps, warnings=warns,
                             summary={"warn_below_retention_pct": 30.0,
                                      "warn_step_loss_pct": 70.0}).render()
    assert 'class="warn-head"' in html and "Read-tracking warnings" in html
    assert "not necessarily errors" in html  # the explanatory context
    # warnings render in the same colorized terminal style as the run log
    assert "read-tracking warnings</span>" in html      # terminal window title
    assert "term-body compact" in html                  # compact terminal body
    assert "#e3b341" in html                             # amber [WARN] colorization


def test_tables_wrap_long_text_and_scroll(tmp_path):
    """Tables wrap long values and give data tables natural, non-cramped widths."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    b = ReadTrackingBuilder("m", logs_dir=logs)
    html = HTMLReportBuilder("m", b.build(), steps=b.steps).render()
    assert "th, td{overflow-wrap:anywhere;}" in html          # long paths/names wrap
    assert ".scroll table{width:auto; min-width:100%;" in html  # data tables get full width
    assert ".scroll th, .scroll td{white-space:nowrap;}" in html  # one line + scroll, no cramping


def test_run_log_has_css_fullscreen_toggle(tmp_path):
    """The run-log terminal offers a pure-CSS Fullscreen toggle (no JS)."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    run_log = tmp_path / "m_pipeline_run.log"
    _write_run_log(run_log, [("INFO", "start"), ("INFO", "done")])
    b = ReadTrackingBuilder("m", logs_dir=logs)
    html = HTMLReportBuilder("m", b.build(), steps=b.steps, log_file=run_log).render()
    assert 'id="termmax"' in html and 'class="term-max-btn"' in html
    assert "Fullscreen" in html and "Exit fullscreen" in html
    assert "#termmax:checked ~ .terminal{position:fixed" in html  # CSS-only maximize
    assert "<script" not in html.lower()                          # still no JavaScript


def test_html_report_tabbed_panels(tmp_path):
    """Sections render as pure-CSS selectable tabs (no JS), one panel per section."""
    import re
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    run_log = tmp_path / "m_pipeline_run.log"
    _write_run_log(run_log, [("INFO", "go"), ("INFO", "done")])
    b = ReadTrackingBuilder("m", logs_dir=logs)
    html = HTMLReportBuilder("m", b.build(), steps=b.steps, log_file=run_log).render()
    radios = re.findall(r'class="tab-radio"[^>]*id="tab-(\d+)"', html)
    labels = re.findall(r'<label for="tab-(\d+)">([^<]+)</label>', html)
    panels = re.findall(r'<section class="panel" id="panel-(\d+)">', html)
    assert len(radios) == len(labels) == len(panels) >= 4   # several sections
    assert 'id="tab-0" checked' in html                     # first panel selected by default
    assert "#tab-0:checked ~ #panel-0" in html              # pure-CSS show rule
    assert "<script" not in html.lower()                    # no JavaScript
    assert "http://" not in html                            # self-contained
    # the run log is its own dedicated, selectable panel
    assert any(t == "Run log" for _, t in labels)
    assert "Notes & methods" in {t for _, t in labels}      # methods is a panel too


def test_summary_is_its_own_tab_not_repeated(tmp_path):
    """The summary is a single first tab, never duplicated above other panels."""
    import re
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "S1", raw=1000, trimmed=900)
    b = ReadTrackingBuilder("m", logs_dir=logs)
    html = HTMLReportBuilder("m", b.build(), steps=b.steps).render()
    labels = re.findall(r'<label for="tab-\d+">([^<]+)</label>', html)
    assert labels[0] == "Summary"                 # first tab
    assert html.count("Summary</h2>") == 1        # appears once, not on every panel
    assert html.count("Table 1.") == 1            # the run-summary table is not repeated
    # the old always-visible front matter is gone
    assert 'class="title-block"' not in html and 'class="abstract"' not in html


def test_report_has_no_em_or_curly_punctuation(tmp_path):
    """No em/en dashes or curly quotes in the rendered report (plain-ASCII house style)."""
    from seednap.steps.report import HTMLReportBuilder
    logs = tmp_path / "logs"; logs.mkdir()
    _write_trim_logs(logs, "Blank-PCR-1", raw=10000, trimmed=50)  # triggers warnings prose
    run_log = tmp_path / "m_pipeline_run.log"
    _write_run_log(run_log, [("INFO", "start"), ("WARNING", "low yield"), ("INFO", "done")])
    b = ReadTrackingBuilder("m", logs_dir=logs)
    df = b.build()
    html = HTMLReportBuilder("m", df, steps=b.steps, warnings=b.warnings(df, log=False),
                             log_file=run_log,
                             summary={"warn_below_retention_pct": 30.0,
                                      "warn_step_loss_pct": 70.0}).render()
    for bad in ("—", "–", "&mdash;", "&ndash;", "‘", "’", "&lsquo;", "&rsquo;", "&middot;"):
        assert bad not in html, f"found {bad!r} in report"


def test_report_config_defaults_and_strictness():
    c = ReportConfig()
    # The report step (listed in pipeline.steps) always writes the read-tracking table;
    # html_report toggles the HTML document and is on by default.
    assert c.html_report is True
    assert c.output_dir is None  # defaults to <paths.output>/04_report
    assert c.warn_below_retention_pct == 30.0 and c.warn_step_loss_pct == 70.0
    with pytest.raises(Exception):
        ReportConfig(read_tracking=True)  # removed: the report step is gated via pipeline.steps
    with pytest.raises(Exception):
        ReportConfig(unknown_field=1)  # extra="forbid"


def test_report_output_dir_expands_and_overrides(tmp_path):
    c = ReportConfig(output_dir=str(tmp_path / "reports"))
    assert c.output_dir == (tmp_path / "reports").resolve()
