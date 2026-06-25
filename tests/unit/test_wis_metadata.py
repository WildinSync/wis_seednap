"""WIS database -> GBIF metadata bridge (`seednap wis-metadata`).

The bridge reads per-sample field metadata from the WIS PostgreSQL/PostGIS database and writes
the two CSVs the DarwinCore export consumes. SQLAlchemy + a Postgres driver are an optional
dependency, and the live PostGIS database is not reachable from CI, so these tests target the
pure transform layer (env_medium mapping, date formatting, column assembly) with synthetic dict
rows, and cross-check that the bridge's output satisfies DarwinCoreBuilder's own validation. The
thin Postgres/PostGIS query layer (ST_X/ST_Y over gis_point.geom) can only be validated against a
real WIS instance; here we cover its dependency-guard and zero-row behaviour with the SQL itself
stubbed.
"""

import logging
import sys
from datetime import date

import pandas as pd
import pytest

from seednap.steps.formatting.darwincore_builder import _ENVO_TERMS, DarwinCoreBuilder
from seednap.steps.formatting.wis_metadata import (
    WIS_SAMPLE_TYPE_TO_ENV_MEDIUM,
    WisMetadataExporter,
    _env_medium_from_sample_type,
    _format_event_date,
    build_project_metadata_df,
    build_sample_metadata_df,
)


def _sample_rows():
    """Two synthetic WIS sample rows in the shape fetch_sample_rows returns."""
    return [
        {
            "sample_id": "DAR2025005",
            "material_sample_id": "fw_ch_2025_dar2025005",
            "event_date": date(2025, 3, 15),
            "sample_depth": 1.5,
            "sample_size_frac": 0.45,
            "sample_size": 2000.0,
            "sample_type_code": "FW",
            "lat": 46.51,
            "lon": 6.63,
        },
        {
            "sample_id": "DAR2025006",
            "material_sample_id": "se_ch_2025_dar2025006",
            "event_date": date(2025, 3, 16),
            "sample_depth": None,
            "sample_size_frac": None,
            "sample_size": None,
            "sample_type_code": "SE",
            "lat": None,
            "lon": None,
        },
    ]


# --------------------------------------------------------------------------- env_medium mapping


@pytest.mark.parametrize(
    "code,expected",
    [
        ("FW", "water"),
        ("MA", "marine"),
        ("SE", "sediment"),
        ("SO", "soil"),
        ("SU", "water"),
        ("fw", "water"),
    ],  # case-insensitive
)
def test_env_medium_known_codes(code, expected):
    assert _env_medium_from_sample_type(code) == expected


def test_env_medium_missing_is_blank():
    assert _env_medium_from_sample_type(None) == ""
    assert _env_medium_from_sample_type(float("nan")) == ""


def test_env_medium_unmapped_passes_raw_label_and_warns(caplog):
    with caplog.at_level(logging.WARNING):
        out = _env_medium_from_sample_type("AI")  # Air: no aquatic/soil analogue
    assert out == "Air"  # raw label passed through, not silently mapped to water
    assert any("[WARN] wis-metadata" in r.message and "AI" in r.message for r in caplog.records)


def test_every_mapped_value_is_a_builder_env_medium_term():
    # Cross-check: the bridge must only ever emit env_medium labels the builder recognises,
    # so a mapped sample never trips the builder's G3 'unknown env_medium' failure.
    assert set(WIS_SAMPLE_TYPE_TO_ENV_MEDIUM.values()) <= set(_ENVO_TERMS)


# --------------------------------------------------------------------------- date formatting


@pytest.mark.parametrize(
    "value,expected",
    [
        (date(2024, 3, 15), "2024.03.15"),
        ("2024-03-15", "2024.03.15"),
        ("2024/03/15", "2024.03.15"),
        (None, ""),
    ],
)
def test_format_event_date(value, expected):
    assert _format_event_date(value) == expected


def test_format_event_date_unparseable_blank_and_warns(caplog):
    with caplog.at_level(logging.WARNING):
        assert _format_event_date("not-a-date") == ""
    assert any("[WARN] wis-metadata" in r.message for r in caplog.records)


# --------------------------------------------------------------------------- DataFrame assembly


def test_build_sample_metadata_df_shape_and_values():
    df = build_sample_metadata_df(_sample_rows())
    assert list(df.columns) == [
        "eventID",
        "eventDate",
        "env_medium",
        "decimalLatitude",
        "decimalLongitude",
        "depth",
        "size_frac",
        "samp_size",
    ]
    assert df.loc[0, "eventID"] == "DAR2025005"  # defaults to sample_id
    assert df.loc[0, "eventDate"] == "2025.03.15"
    assert df.loc[0, "env_medium"] == "water"
    assert df.loc[0, "decimalLatitude"] == 46.51
    assert df.loc[1, "env_medium"] == "sediment"
    assert pd.isna(df.loc[1, "decimalLatitude"])  # no COORDINATE point -> blank, not 0


def test_build_sample_metadata_df_event_id_field_selectable():
    df = build_sample_metadata_df(_sample_rows(), event_id_field="material_sample_id")
    assert df.loc[0, "eventID"] == "fw_ch_2025_dar2025005"


def test_build_sample_metadata_df_rejects_unknown_event_id_field():
    with pytest.raises(ValueError, match="event_id_field"):
        build_sample_metadata_df(_sample_rows(), event_id_field="nope")


def test_build_project_metadata_df_single_row():
    df = build_project_metadata_df(
        marker="teleo",
        recorded_by="ELE Lab",
        identification_remarks="BLAST LCA",
        identification_references="MitoFish 2025",
    )
    assert len(df) == 1
    assert df.loc[0, "marker"] == "teleo"
    assert df.loc[0, "recordedby"] == "ELE Lab"
    # otu_db / chimera_check are intentionally absent (filled by the darwincore step from config)
    assert "otu_db" not in df.columns and "chimera_check" not in df.columns


# ----------------------------------------------- cross-check: output satisfies builder contract


def test_bridge_output_passes_builder_validation(tmp_path):
    # The strongest check available without a DB: the sample metadata the bridge produces must
    # pass DarwinCoreBuilder's own _validate_sample_metadata (required columns, coordinate ranges,
    # env_medium vocabulary) and _validate_dates, so a DB-sourced export cannot fail validation
    # for a shape the bridge controls.
    df = build_sample_metadata_df(_sample_rows())
    DarwinCoreBuilder._validate_sample_metadata(df, tmp_path / "sample.csv")
    DarwinCoreBuilder._validate_dates(df["eventDate"])

    proj = build_project_metadata_df(
        marker="teleo",
        recorded_by="ELE Lab",
        identification_remarks="BLAST LCA",
        identification_references="MitoFish 2025",
    )
    DarwinCoreBuilder._validate_project_metadata(proj, tmp_path / "project.csv")


# --------------------------------------------------------------------------- query layer guards


def test_missing_sqlalchemy_raises_actionable_error(monkeypatch):
    # Simulate the optional extra not being installed; the bridge must say how to fix it.
    monkeypatch.setitem(sys.modules, "sqlalchemy", None)
    with pytest.raises(RuntimeError, match=r"seednap\[wis\]"):
        WisMetadataExporter._require_sqlalchemy()


def test_export_raises_on_zero_samples(tmp_path, monkeypatch):
    exporter = WisMetadataExporter("postgresql://x")
    monkeypatch.setattr(exporter, "fetch_sample_rows", lambda **_kw: [])
    with pytest.raises(ValueError, match="No samples found"):
        exporter.export(
            output_dir=tmp_path,
            marker="teleo",
            recorded_by="ELE",
            identification_remarks="r",
            identification_references="ref",
            monitoring="site-1",
        )


def test_export_round_trip_writes_both_csvs(tmp_path, monkeypatch):
    # Exercises export() end to end with the (Postgres-only) SQL stubbed: the two CSVs are written
    # with builder-valid content, and the sample CSV round-trips through the builder validation.
    exporter = WisMetadataExporter("postgresql://x")
    monkeypatch.setattr(exporter, "fetch_sample_rows", lambda **_kw: _sample_rows())
    sample_csv, project_csv = exporter.export(
        output_dir=tmp_path,
        marker="teleo",
        recorded_by="ELE Lab",
        identification_remarks="BLAST LCA",
        identification_references="MitoFish 2025",
        monitoring="site-1",
    )
    assert sample_csv.name == "teleo_sample_metadata.csv"
    assert project_csv.name == "teleo_project_metadata.csv"
    sm = pd.read_csv(sample_csv)
    assert list(sm["eventID"]) == ["DAR2025005", "DAR2025006"]
    DarwinCoreBuilder._validate_sample_metadata(sm, sample_csv)
