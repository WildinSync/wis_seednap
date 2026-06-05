"""Unit tests for the FAIRe sample manifest: model, control classifier, date
normalisation, migrator, loader, and the cross-CSV eventID validator.

Fixtures are tiny synthetic CSVs that reproduce the real-world traps catalogued across the
lab datasets (BOM, per-marker dotted-date order, capitalised/unit-suffixed headers,
Blank-ext/Blank-PCR controls, the silent-ID-mismatch orphan column).
"""

import logging

import pandas as pd
import pytest

from seednap.config.manifest import (
    SampleManifestRow,
    classify_control,
    load_manifest,
    validate_against_abundance,
)
from seednap.config.manifest_migrate import migrate_to_manifest, normalise_event_dates


# --------------------------------------------------------------------------- #
# Model
# --------------------------------------------------------------------------- #
def test_valid_sample_row():
    r = SampleManifestRow(
        eventID="DAR-1", seq_run_id="run1", samp_category="sample",
        eventDate="2025-08-19", decimalLatitude="46.3", samp_size="15",
    )
    assert r.eventDate == "2025-08-19"
    assert r.decimalLatitude == pytest.approx(46.3)
    assert r.samp_size == pytest.approx(15.0)
    assert not r.is_control


def test_unknown_column_is_rejected():
    """extra='forbid': a typo'd column is a hard error (CLAUDE.md sec.5)."""
    with pytest.raises(Exception):
        SampleManifestRow(
            eventID="x", seq_run_id="r", samp_category="sample",
            eventDate="2025", typo_column="oops",
        )


def test_dotted_date_rejected_by_model():
    """Canonical manifests are ISO-only; the migrator normalises legacy dotted dates."""
    with pytest.raises(Exception):
        SampleManifestRow(eventID="x", seq_run_id="r", samp_category="sample", eventDate="2025.08.19")


@pytest.mark.parametrize("iso", ["2025", "2025-08", "2025-08-19"])
def test_iso_partial_dates_accepted(iso):
    r = SampleManifestRow(eventID="x", seq_run_id="r", samp_category="sample", eventDate=iso)
    assert r.eventDate == iso


def test_na_and_insdc_tokens_become_none():
    """'NA' and INSDC missing-value tokens normalise to None, not a literal/NaN."""
    r = SampleManifestRow(
        eventID="Blank-PCR-1", seq_run_id="r", samp_category="negative control",
        neg_cont_type="PCR negative", decimalLatitude="NA",
        extraction_ID="not applicable: control sample",
    )
    assert r.decimalLatitude is None
    assert r.extraction_ID is None


def test_invalid_samp_category_rejected():
    with pytest.raises(Exception):
        SampleManifestRow(eventID="x", seq_run_id="r", samp_category="blank", eventDate="2025")


def test_other_prefixed_categories_accepted():
    r = SampleManifestRow(eventID="x", seq_run_id="r", samp_category="other: mock", eventDate="2025")
    assert r.samp_category == "other: mock"


def test_negative_control_requires_neg_cont_type():
    with pytest.raises(Exception):
        SampleManifestRow(eventID="Blank-PCR-1", seq_run_id="r", samp_category="negative control")


def test_invalid_neg_cont_type_rejected():
    with pytest.raises(Exception):
        SampleManifestRow(eventID="b", seq_run_id="r", samp_category="negative control",
                          neg_cont_type="bogus")


def test_out_of_range_latitude_rejected():
    with pytest.raises(Exception):
        SampleManifestRow(eventID="x", seq_run_id="r", samp_category="sample",
                          eventDate="2025", decimalLatitude="200")


def test_sample_without_event_date_is_allowed_but_flagged(caplog):
    """eventDate is not a constructor hard-fail (demux-stage manifests have none) but
    check_completeness must WARN."""
    from seednap.config.manifest import SampleManifest

    row = SampleManifestRow(eventID="DAR-1", seq_run_id="r", samp_category="sample")
    m = SampleManifest(rows=[row])
    with caplog.at_level(logging.WARNING):
        m.check_completeness()
    assert any("eventDate" in rec.message and "DAR-1" in rec.message for rec in caplog.records)


# --------------------------------------------------------------------------- #
# Control classification (single source of truth, superset of the legacy regex)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "name,category,neg,warns",
    [
        ("Blank-ext-1", "negative control", "extraction negative", False),
        ("Blank-PCR-2", "negative control", "PCR negative", False),
        ("Blank_PCR-1", "negative control", "PCR negative", False),       # underscore form
        ("Blank-ext-2run2", "negative control", "extraction negative", False),  # run suffix
        ("CNEG01_03-MB1123A4", "negative control", "PCR negative", False),
        ("CEXT02", "negative control", "extraction negative", False),
        ("CMET01", "negative control", "process negative", True),          # inference -> warn
        ("CPCR03", "negative control", "PCR negative", False),             # MISSED by legacy regex
        ("EXT_NC", "negative control", "extraction negative", False),      # MISSED by legacy regex
        ("PCR_NC", "negative control", "PCR negative", False),             # MISSED by legacy regex
        ("water", "negative control", "other: water control", True),       # MISSED by legacy regex
        ("DAR-2025-1103", "sample", None, False),
        ("EVS1327", "sample", None, False),
    ],
)
def test_classify_control(name, category, neg, warns):
    c = classify_control(name)
    assert c.samp_category == category
    if category == "negative control":
        assert c.neg_cont_type == neg
    assert bool(c.warn_reason) == warns


def test_control_like_but_unclassified_warns():
    """A control-looking name that matches no rule is NOT silently treated as a sample."""
    c = classify_control("negctrl_weird")
    assert c.samp_category == "sample"
    assert c.warn_reason is not None
    assert c.rule == "unclassified-control-like"


def test_pcr_blank_flagged_for_null_extraction():
    assert classify_control("Blank-PCR-1").is_pcr_blank is True
    assert classify_control("Blank-ext-1").is_pcr_blank is False


# --------------------------------------------------------------------------- #
# Date normalisation (the silent-corruption guard)
# --------------------------------------------------------------------------- #
def test_year_first_dotted():
    m = normalise_event_dates(["2025.08.19", "2025.09.04"])
    assert m["2025.08.19"] == "2025-08-19" and m["2025.09.04"] == "2025-09-04"


def test_year_last_dotted_disambiguated_by_day_gt_12():
    m = normalise_event_dates(["19.08.2025", "04.09.2025"])
    assert m["19.08.2025"] == "2025-08-19"
    assert m["04.09.2025"] == "2025-09-04"  # 19>12 fixed the file as DD.MM.YYYY


def test_iso_passthrough_and_partials():
    m = normalise_event_dates(["2025-08-19", "2025.08", "2025"])
    assert m["2025-08-19"] == "2025-08-19"
    assert m["2025.08"] == "2025-08"
    assert m["2025"] == "2025"


def test_ambiguous_year_last_raises():
    with pytest.raises(ValueError, match="ambiguous"):
        normalise_event_dates(["04.09.2025", "05.06.2025"])


def test_mixed_year_first_and_last_raises():
    with pytest.raises(ValueError, match="mixes"):
        normalise_event_dates(["2025.08.19", "19.08.2025"])


def test_contradictory_order_raises():
    with pytest.raises(ValueError, match="contradictory"):
        normalise_event_dates(["19.08.2025", "04.13.2025"])


def test_unparseable_date_raises():
    with pytest.raises(ValueError, match="unrecognised"):
        normalise_event_dates(["not-a-date"])


def test_na_tokens_skipped_not_parsed():
    """A column mixing valid dates and 'NA' must parse the valid ones, not crash."""
    m = normalise_event_dates(["2024.04.30", "NA", ""])
    assert m["2024.04.30"] == "2024-04-30"
    assert "NA" not in m


def test_year_first_ddmm_order():
    """YYYY.DD.MM (year first, day before month) disambiguated by a day > 12."""
    m = normalise_event_dates(["2025.13.06", "2025.18.09"])
    assert m["2025.13.06"] == "2025-06-13"
    assert m["2025.18.09"] == "2025-09-18"


def test_slash_dates_with_time():
    """MM/DD/YYYY HH:MM:SS: slash separator + trailing clock time."""
    m = normalise_event_dates(["05/20/2021 00:00:00", "05/21/2021 00:00:00"])
    assert m["05/20/2021 00:00:00"] == "2021-05-20"  # 20>12 -> day -> MM.DD order


def test_two_digit_year_refused():
    with pytest.raises(ValueError, match="unrecognised|ambiguous"):
        normalise_event_dates(["24.08.22", "26.08.22"])


def test_ambiguous_dates_raise_without_order_but_resolve_with_order():
    """All-<=12 dotted dates are ambiguous: raise by default, parse with an explicit order."""
    ambiguous = ["2024.09.08", "2024.03.05"]
    with pytest.raises(ValueError, match="ambiguous"):
        normalise_event_dates(ambiguous)
    ymd = normalise_event_dates(ambiguous, order="YMD")
    assert ymd["2024.09.08"] == "2024-09-08"
    dmy_src = ["08.09.2024", "05.03.2024"]
    dmy = normalise_event_dates(dmy_src, order="DMY")
    assert dmy["08.09.2024"] == "2024-09-08"


# --------------------------------------------------------------------------- #
# Migrator
# --------------------------------------------------------------------------- #
def _write_field(path, rows, *, bom=False, header=None):
    """Write a synthetic field-metadata CSV (optionally with a BOM)."""
    df = pd.DataFrame(rows)
    if header:
        df = df.rename(columns=header)
    enc = "utf-8-sig" if bom else "utf-8"
    df.to_csv(path, index=False, encoding=enc)
    return path


def test_migrator_modern_field(tmp_path, caplog):
    field = _write_field(tmp_path / "metadata_field_demo.csv", [
        {"eventID": "DAR-1", "eventDate": "2025.08.19", "decimalLatitude": "46.3",
         "decimalLongitude": "6.9", "volume": "15", "depth": "0.25",
         "env_medium": "water", "size_frac": "0.45", "extraction_ID": "EXP1"},
        {"eventID": "DAR-2", "eventDate": "2025.08.20", "decimalLatitude": "46.4",
         "decimalLongitude": "6.8", "volume": "15", "depth": "0.25",
         "env_medium": "water", "size_frac": "0.45", "extraction_ID": "EXP1"},
        {"eventID": "Blank-ext-1", "eventDate": "", "decimalLatitude": "",
         "decimalLongitude": "", "volume": "", "depth": "",
         "env_medium": "", "size_frac": "", "extraction_ID": "EXP1"},
        {"eventID": "Blank-PCR-1", "eventDate": "", "decimalLatitude": "",
         "decimalLongitude": "", "volume": "", "depth": "",
         "env_medium": "", "size_frac": "", "extraction_ID": "NA"},
    ])
    with caplog.at_level(logging.WARNING):
        m = migrate_to_manifest(field, target_gene="teleo")

    assert len(m) == 4
    assert len(m.biological_samples()) == 2
    assert len(m.controls()) == 2
    # volume -> samp_size, depth -> maximumDepthInMeters, dotted date -> ISO
    s = {r.eventID: r for r in m.rows}
    assert s["DAR-1"].samp_size == pytest.approx(15.0)
    assert s["DAR-1"].maximumDepthInMeters == pytest.approx(0.25)
    assert s["DAR-1"].eventDate == "2025-08-19"
    assert s["DAR-1"].target_gene == "teleo"
    # control identity + extraction expectations
    assert s["Blank-ext-1"].neg_cont_type == "extraction negative"
    assert s["Blank-ext-1"].extraction_ID == "EXP1"
    assert s["Blank-PCR-1"].neg_cont_type == "PCR negative"
    assert s["Blank-PCR-1"].extraction_ID is None  # null is EXPECTED for PCR blanks
    # seq_run_id was synthesised (no library column) and warned
    assert m.seq_run_ids() == ["demo_teleo"]
    assert any("synthesised seq_run_id" in r.message for r in caplog.records)


def test_migrator_strips_bom_and_normalises_headers(tmp_path):
    """A BOM glued to eventID and capitalised/unit-suffixed headers must resolve."""
    field = _write_field(
        tmp_path / "metadata_field_bom.csv",
        [{"eventID": "DAR-1", "eventDate": "2025-05-15", "Site_names": "Foo",
          "Conductivity [µS]": "200", "volume": "15"}],
        bom=True,
    )
    m = migrate_to_manifest(field, target_gene="mam07")
    assert m.event_ids() == ["DAR-1"]
    assert m.rows[0].samp_size == pytest.approx(15.0)  # 'volume' still mapped despite BOM/extra cols


def test_migrator_ddmmyyyy_per_file(tmp_path):
    field = _write_field(tmp_path / "metadata_field_eu.csv", [
        {"eventID": "DAR-1", "eventDate": "19.08.2025", "extraction_ID": "EXP1"},
        {"eventID": "DAR-2", "eventDate": "04.09.2025", "extraction_ID": "EXP1"},
    ])
    m = migrate_to_manifest(field, target_gene="teleo")
    dates = {r.eventID: r.eventDate for r in m.rows}
    assert dates == {"DAR-1": "2025-08-19", "DAR-2": "2025-09-04"}


def test_migrator_zero_controls_warns(tmp_path, caplog):
    field = _write_field(tmp_path / "metadata_field_noctrl.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01"},
        {"eventID": "DAR-2", "eventDate": "2025-01-02"},
    ])
    with caplog.at_level(logging.WARNING):
        m = migrate_to_manifest(field, target_gene="teleo")
    assert len(m.controls()) == 0
    assert any("no" in r.message.lower() and "control" in r.message.lower() for r in caplog.records)


def test_migrator_unexpected_column_warns(tmp_path, caplog):
    """A stray/unknown column (possible shifted file) warns rather than dropping silently."""
    field = _write_field(tmp_path / "metadata_field_stray.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01", "Test": "junk"},
    ])
    with caplog.at_level(logging.WARNING):
        migrate_to_manifest(field, target_gene="teleo")
    assert any("'Test'" in r.message or "Test" in r.message for r in caplog.records
               if "unexpected" in r.message.lower())


def test_migrator_missing_event_id_column_raises(tmp_path):
    field = _write_field(tmp_path / "metadata_field_nokey.csv", [
        {"name": "DAR-1", "eventDate": "2025-01-01"},
    ])
    with pytest.raises(ValueError, match="eventID"):
        migrate_to_manifest(field, target_gene="teleo")


def test_migrator_strips_unit_from_numeric_value(tmp_path, caplog):
    """A numeric cell with a unit ('17 L') is cleaned with a WARN, not a crash."""
    field = _write_field(tmp_path / "metadata_field_unit.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01", "volume": "17 L"},
    ])
    with caplog.at_level(logging.WARNING):
        m = migrate_to_manifest(field, target_gene="teleo")
    assert m.rows[0].samp_size == pytest.approx(17.0)
    assert any("stripped non-numeric" in r.message for r in caplog.records)


def test_migrator_out_of_range_coord_nulled(tmp_path, caplog):
    """A typo'd out-of-range coordinate is nulled with a WARN; the row survives."""
    field = _write_field(tmp_path / "metadata_field_badcoord.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01", "decimalLongitude": "-56239"},
    ])
    with caplog.at_level(logging.WARNING):
        m = migrate_to_manifest(field, target_gene="teleo")
    assert len(m) == 1
    assert m.rows[0].decimalLongitude is None
    assert any("out of range" in r.message for r in caplog.records)


def test_migrator_negative_marine_depth_allowed(tmp_path):
    """Marine datasets encode depth as negative; it must be accepted, not rejected."""
    field = _write_field(tmp_path / "metadata_field_marine.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01", "depth": "-840"},
    ])
    m = migrate_to_manifest(field, target_gene="teleo")
    assert m.rows[0].maximumDepthInMeters == pytest.approx(-840.0)


def test_migrator_extraction_blank_na_token_fires_intended_warn(tmp_path, caplog):
    """Adversarial-audit regression: an extraction blank with extraction_ID='NA' must read
    as null and fire the 'got=none' WARN, not register a fabricated 'NA' batch."""
    field = _write_field(tmp_path / "metadata_field_na.csv", [
        {"eventID": "S1", "eventDate": "2024-01-27", "extraction_ID": "EXP1"},
        {"eventID": "Blank-ext-1", "eventDate": "2024-01-27", "extraction_ID": "NA"},
        {"eventID": "Blank-PCR-1", "eventDate": "2024-01-26", "extraction_ID": "NA"},
    ])
    with caplog.at_level(logging.WARNING):
        m = migrate_to_manifest(field, target_gene="mam07", seq_run_id="R1")
    rows = {r.eventID: r for r in m.rows}
    assert rows["Blank-ext-1"].extraction_ID is None
    msgs = " ".join(r.message for r in caplog.records)
    assert "extraction blank 'Blank-ext-1'" in msgs and "got=none" in msgs
    assert "orphan blank batch" not in msgs  # the misleading warn must NOT fire


def test_migrator_per_row_seq_run_fallback_warns(tmp_path, caplog):
    """Adversarial-audit regression: when a library column exists but a row's cell is empty,
    that row's fallback to the default batch must emit a per-row WARN naming the eventID."""
    field = _write_field(tmp_path / "metadata_field_partlib.csv", [
        {"eventID": "S1", "eventDate": "2024-01-01", "library": "LIB_A"},
        {"eventID": "S2", "eventDate": "2024-01-02", "library": ""},
    ])
    with caplog.at_level(logging.WARNING):
        m = migrate_to_manifest(field, target_gene="teleo")
    rows = {r.eventID: r for r in m.rows}
    assert rows["S1"].seq_run_id == "LIB_A"
    assert any("library/seq_run_id grouping for 'S2'" in r.message for r in caplog.records)


def test_migrator_header_alias_collision_warns(tmp_path, caplog):
    """Adversarial-audit regression: two columns mapping to the same field must WARN, not
    silently drop the second."""
    # volume and samp_size both -> samp_size
    path = tmp_path / "metadata_field_collide.csv"
    pd.DataFrame([{"eventID": "S1", "eventDate": "2024-01-01", "volume": "15", "samp_size": "20"}]).to_csv(path, index=False)
    with caplog.at_level(logging.WARNING):
        migrate_to_manifest(path, target_gene="teleo")
    assert any("map to 'samp_size'" in r.message for r in caplog.records)


def test_migrator_inf_nan_depth_nulled(tmp_path, caplog):
    """Adversarial-audit regression: a non-finite depth must be nulled with a WARN."""
    field = _write_field(tmp_path / "metadata_field_inf.csv", [
        {"eventID": "S1", "eventDate": "2024-01-01", "depth": "inf"},
    ])
    with caplog.at_level(logging.WARNING):
        m = migrate_to_manifest(field, target_gene="teleo")
    assert m.rows[0].maximumDepthInMeters is None
    assert any("inf/nan" in r.message for r in caplog.records)


def test_model_rejects_non_finite_float():
    """The strict model rejects inf/nan directly (allow_inf_nan=False)."""
    with pytest.raises(Exception):
        SampleManifestRow(eventID="x", seq_run_id="r", samp_category="sample",
                          eventDate="2024", maximumDepthInMeters="inf")


def test_positive_control_requires_pos_cont_type():
    """FAIRe Mandatory-if symmetry: a positive control needs pos_cont_type."""
    with pytest.raises(Exception):
        SampleManifestRow(eventID="POS1", seq_run_id="r", samp_category="positive control")
    ok = SampleManifestRow(eventID="POS1", seq_run_id="r", samp_category="positive control",
                           pos_cont_type="mock community")
    assert ok.pos_cont_type == "mock community"


def test_validator_keeps_sample_named_like_meta_token(tmp_path, caplog):
    """Adversarial-audit regression: a real sample named 'total' (an OTU-table meta token)
    must not be silently dropped from the cross-check."""
    field = _write_field(tmp_path / "metadata_field_meta.csv", [
        {"eventID": "total", "eventDate": "2024-01-01"},
        {"eventID": "S1", "eventDate": "2024-01-02"},
    ])
    m = migrate_to_manifest(field, target_gene="teleo")
    path = tmp_path / "otu_meta.csv"
    pd.DataFrame([["ACGT", "total", "S1"]], columns=["sequence", "total", "S1"]).to_csv(path, index=False)
    with caplog.at_level(logging.WARNING):
        res = validate_against_abundance(m, path)
    assert "total" in res.abundance_samples  # kept, not silently dropped
    assert res.ok


def test_migrator_legacy_lab_library_and_tags(tmp_path):
    """A demux lab CSV: library -> seq_run_id, tag_demultiplex -> mid_forward."""
    lab = tmp_path / "metadata_lab_demo.csv"
    pd.DataFrame([
        {"eventID": "SPY1", "tag_demultiplex": "acgt", "library": "LIB_A"},
        {"eventID": "CNEG01", "tag_demultiplex": "tgca", "library": "LIB_A"},
    ]).to_csv(lab, index=False)
    m = migrate_to_manifest(lab)
    rows = {r.eventID: r for r in m.rows}
    assert rows["SPY1"].seq_run_id == "LIB_A"
    assert rows["SPY1"].mid_forward == "acgt"
    assert rows["CNEG01"].samp_category == "negative control"


# --------------------------------------------------------------------------- #
# load_manifest (canonical CSV)
# --------------------------------------------------------------------------- #
def test_load_manifest_round_trip(tmp_path):
    field = _write_field(tmp_path / "metadata_field_rt.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01", "extraction_ID": "EXP1"},
        {"eventID": "Blank-ext-1", "eventDate": "", "extraction_ID": "EXP1"},
    ])
    m = migrate_to_manifest(field, target_gene="teleo")
    out = tmp_path / "manifest.csv"
    m.to_csv(out)
    reloaded = load_manifest(out)
    assert sorted(reloaded.event_ids()) == ["Blank-ext-1", "DAR-1"]
    assert len(reloaded.controls()) == 1


def test_load_manifest_duplicate_key_raises(tmp_path):
    path = tmp_path / "dup.csv"
    pd.DataFrame([
        {"eventID": "DAR-1", "seq_run_id": "r", "samp_category": "sample", "eventDate": "2025"},
        {"eventID": "DAR-1", "seq_run_id": "r", "samp_category": "sample", "eventDate": "2025"},
    ]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="duplicate"):
        load_manifest(path)


def test_load_manifest_unknown_column_raises(tmp_path):
    path = tmp_path / "bad.csv"
    pd.DataFrame([
        {"eventID": "DAR-1", "seq_run_id": "r", "samp_category": "sample",
         "eventDate": "2025", "bogus_col": "x"},
    ]).to_csv(path, index=False)
    with pytest.raises(ValueError):
        load_manifest(path)


def test_load_manifest_empty_raises(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("eventID,seq_run_id,samp_category\n")
    with pytest.raises(ValueError, match="empty"):
        load_manifest(path)


# --------------------------------------------------------------------------- #
# Cross-CSV eventID validator (E2)
# --------------------------------------------------------------------------- #
def _write_abundance(path, sample_cols, *, id_col="sequence"):
    cols = [id_col] + sample_cols
    pd.DataFrame([["ACGT"] + [1] * len(sample_cols)], columns=cols).to_csv(path, index=False)
    return path


def test_validator_detects_orphan_abundance_column(tmp_path, caplog):
    """The documented silent-ID-mismatch case: an abundance column with no manifest row."""
    field = _write_field(tmp_path / "metadata_field_v.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01"},
        {"eventID": "Blank-PCR-1", "eventDate": ""},
    ])
    m = migrate_to_manifest(field, target_gene="teleo")
    ab = _write_abundance(tmp_path / "otu.csv", ["DAR-1", "Blank-PCR-3"])  # 3, not 1
    with caplog.at_level(logging.WARNING):
        res = validate_against_abundance(m, ab)
    assert res.orphan_abundance_columns == ["Blank-PCR-3"]
    assert res.manifest_extra_rows == ["Blank-PCR-1"]
    assert not res.ok
    assert any("orphan" in r.message.lower() for r in caplog.records)


def test_validator_raise_on_orphan(tmp_path):
    field = _write_field(tmp_path / "metadata_field_v2.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01"},
    ])
    m = migrate_to_manifest(field, target_gene="teleo")
    ab = _write_abundance(tmp_path / "otu2.csv", ["DAR-1", "ORPHAN"])
    with pytest.raises(ValueError, match="absent from the manifest"):
        validate_against_abundance(m, ab, raise_on_orphan=True)


def test_validator_excludes_metadata_columns(tmp_path):
    """OTU-table metadata columns (sequence, total, chimera...) are not sample columns."""
    field = _write_field(tmp_path / "metadata_field_v3.csv", [
        {"eventID": "DAR-1", "eventDate": "2025-01-01"},
    ])
    m = migrate_to_manifest(field, target_gene="teleo")
    path = tmp_path / "otu3.csv"
    pd.DataFrame([["ACGT", 5, "N", "DAR-1"]],
                 columns=["sequence", "total", "chimera", "DAR-1"]).to_csv(path, index=False)
    res = validate_against_abundance(m, path)
    assert res.abundance_samples == ["DAR-1"]
    assert res.ok
