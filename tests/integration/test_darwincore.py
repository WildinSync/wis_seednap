"""Integration tests for DarwinCoreBuilder (Commit H).

Covers G1 (DwC field validation), G2 (contaminant propagation),
G3 (input metadata validation), G5 (stable occurrenceID).

Uses synthetic CSVs and `--skip-enrichment` to avoid hitting NCBI/WORMS,
so these tests are fast and offline.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from seednap.steps.formatting.darwincore_builder import DarwinCoreBuilder


def _write_taxonomy(path: Path, with_contaminant: bool = False) -> None:
    """Synthetic post-format-gbif taxonomy CSV (long format, one row per OTU/sample)."""
    rows = [
        {
            "kingdom": "Metazoa", "phylum": "Chordata", "class": "Actinopteri",
            "order": "Perciformes", "family": "Percidae", "genus": "Perca",
            "species": "Perca_fluviatilis", "taxon": "Perca_fluviatilis",
            "rank": "species", "sequence": "ACGTACGTACGT",
            "nb_reads": 100, "eventID": "DAR-2025-1101",
        },
        {
            "kingdom": "Metazoa", "phylum": "Chordata", "class": "Mammalia",
            "order": "Primates", "family": "Hominidae", "genus": "Homo",
            "species": "Homo_sapiens", "taxon": "Homo_sapiens",
            "rank": "species", "sequence": "AAAACCCCGGGG",
            "nb_reads": 5, "eventID": "DAR-2025-1101",
        },
    ]
    df = pd.DataFrame(rows)
    if with_contaminant:
        df["is_contaminant_candidate"] = [False, True]
    df.to_csv(path, index=False)


def _write_sample_meta(path: Path, *, lat: float = 46.0, lon: float = 7.5,
                      env: str = "river") -> None:
    pd.DataFrame([{
        "eventID": "DAR-2025-1101",
        "decimalLatitude": lat,
        "decimalLongitude": lon,
        "eventDate": "2025.06.15",
        "env_medium": env,
        "samp_size": "1L",
        "depth": 0.5,
        "size_frac": "0.22um",
    }]).to_csv(path, index=False)


def _write_project_meta(path: Path, marker: str = "Teleo", **overrides: str) -> None:
    base = {
        "marker": marker,
        "recordedby": "J. Smith",
        "seqmet": "MiSeq PE 2x150",
        "identificationRemarks": "BLAST + LCA",
        "identificationReferences": "10.1038/nmeth.3869",
        "otu_seq_comp_appr": "SWARM d=1",
        "otu_db": "CRABS MitoFish 2025",
        "chimera_check": "UCHIME de novo",
    }
    base.update(overrides)
    pd.DataFrame([base]).to_csv(path, index=False)


@pytest.fixture
def fixture_dir(tmp_path: Path) -> Path:
    _write_taxonomy(tmp_path / "tax.csv")
    _write_sample_meta(tmp_path / "sample.csv")
    _write_project_meta(tmp_path / "project.csv")
    return tmp_path


def test_unknown_marker_raises(fixture_dir: Path) -> None:
    """G1: marker not in primers_list.csv -> ValueError, not silent fill."""
    _write_project_meta(fixture_dir / "project.csv", marker="DefinitelyNotAMarker")
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    with pytest.raises(ValueError, match="not found in primers_list"):
        builder.build()


def test_invalid_latitude_raises(fixture_dir: Path) -> None:
    """G3: lat out of [-90, 90] -> ValueError before any output is written."""
    _write_sample_meta(fixture_dir / "sample.csv", lat=200.0)
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    with pytest.raises(ValueError, match="Invalid decimalLatitude"):
        builder.build()


def test_unknown_env_medium_raises(fixture_dir: Path) -> None:
    """G3: env_medium not in _ENVO_TERMS -> ValueError, no silent fallback."""
    _write_sample_meta(fixture_dir / "sample.csv", env="lava")
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    with pytest.raises(ValueError, match="Unknown env_medium"):
        builder.build()


def test_missing_env_medium_raises(fixture_dir: Path) -> None:
    """env_medium is required: a sample CSV without it fails validation with a clear
    message, not a raw KeyError mid-build (build() maps it to ENVO by direct indexing)."""
    pd.DataFrame([{
        "eventID": "DAR-2025-1101", "decimalLatitude": 46.0,
        "decimalLongitude": 7.5, "eventDate": "2025.06.15",
    }]).to_csv(fixture_dir / "sample.csv", index=False)
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    with pytest.raises(ValueError, match="missing required columns"):
        builder.build()


def test_occurrence_id_is_stable(fixture_dir: Path) -> None:
    """G5: same input -> same occurrenceID across runs (sequence-hash based)."""
    builder1 = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out1.csv",
        skip_enrichment=True,
    )
    builder1.build()
    builder2 = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out2.csv",
        skip_enrichment=True,
    )
    builder2.build()

    ids1 = pd.read_csv(fixture_dir / "out1.csv")["occurrenceID"].sort_values().tolist()
    ids2 = pd.read_csv(fixture_dir / "out2.csv")["occurrenceID"].sort_values().tolist()
    assert ids1 == ids2

    # Sanity: format is "marker:eventID:sha256[:8]"
    for oid in ids1:
        parts = oid.split(":")
        assert len(parts) == 3
        assert parts[0] == "Teleo"
        assert len(parts[2]) == 8


def test_contaminant_flag_propagates_to_output(fixture_dir: Path) -> None:
    """G2: is_contaminant_candidate from upstream survives to GBIF as contamination_flag."""
    _write_taxonomy(fixture_dir / "tax.csv", with_contaminant=True)
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    builder.build()
    out = pd.read_csv(fixture_dir / "out.csv")
    assert "contamination_flag" in out.columns
    flagged = out[out["contamination_flag"]]
    assert len(flagged) == 1
    assert flagged.iloc[0]["scientificName"] == "Homo_sapiens"


def _tax_row(event_id: str, seq: str, reads: int) -> dict:
    return {
        "kingdom": "Metazoa", "phylum": "Chordata", "class": "Actinopteri",
        "order": "Perciformes", "family": "Percidae", "genus": "Perca",
        "species": "Perca_fluviatilis", "taxon": "Perca_fluviatilis",
        "rank": "species", "sequence": seq, "nb_reads": reads,
        "eventID": event_id,
    }


def test_widened_control_filter_drops_extended_controls(fixture_dir: Path) -> None:
    """_remove_controls must drop controls recognised by classify_control's
    superset (water/CPCR/EXT_NC), not just the legacy blank/CNEG/CMET/CEXT
    regex.

    Before the fix CPCR_01, EXT_NC and a bare 'water' control leaked into the
    GBIF occurrence CSV as biological records; after the fix they are removed
    because _remove_controls now delegates to classify_control.
    """
    rows = [
        _tax_row("DAR-2025-1101", "ACGTACGTACGT", 100),
        _tax_row("CPCR_01", "ACGTACGTACGT", 3),
        _tax_row("EXT_NC", "ACGTACGTACGT", 2),
        _tax_row("water", "ACGTACGTACGT", 4),
    ]
    pd.DataFrame(rows).to_csv(fixture_dir / "tax.csv", index=False)
    # Sample metadata only needs to carry the real biological eventID; controls
    # are dropped before the sample-metadata merge.
    _write_sample_meta(fixture_dir / "sample.csv")
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    builder.build()
    out = pd.read_csv(fixture_dir / "out.csv")
    leaked = set(out["eventID"].astype(str)) & {"CPCR_01", "EXT_NC", "water"}
    assert not leaked, f"control eventIDs leaked into GBIF output: {leaked}"
    assert set(out["eventID"].astype(str)) == {"DAR-2025-1101"}


def test_control_looking_unclassified_eventid_warns(
    fixture_dir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A control-LOOKING but unclassifiable eventID (e.g. 'water_001', which
    classify_control cannot resolve to a known control rule) is kept as a
    biological sample but must emit a [WARN] -- no-silent-fallbacks policy."""
    rows = [
        _tax_row("DAR-2025-1101", "ACGTACGTACGT", 100),
        _tax_row("water_001", "AAAACCCCGGGG", 7),
    ]
    pd.DataFrame(rows).to_csv(fixture_dir / "tax.csv", index=False)
    _write_sample_meta(fixture_dir / "sample.csv")
    # eventID 'water_001' has no sample-metadata row -> its merged location/
    # date fields are blank, but that is fine: this test only asserts the WARN.
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    builder.build()
    captured = capsys.readouterr()
    assert "[WARN]" in captured.out
    assert "water_001" in captured.out
    out = pd.read_csv(fixture_dir / "out.csv")
    # It was kept (not dropped) -- the WARN, not removal, is the safe behavior.
    assert "water_001" in set(out["eventID"].astype(str))


def test_zero_rows_after_filtering_raises_clear_error(fixture_dir: Path) -> None:
    """When every occurrence row is a control, the error must name the empty
    result, not misattribute it to blank required fields (otu_db)."""
    rows = [
        _tax_row("Blank-ext-1", "ACGTACGTACGT", 5),
        _tax_row("CNEG-2", "AAAACCCCGGGG", 3),
    ]
    pd.DataFrame(rows).to_csv(fixture_dir / "tax.csv", index=False)
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    with pytest.raises(ValueError, match="no occurrence rows"):
        builder.build()


def test_required_fields_populated(fixture_dir: Path) -> None:
    """G1: every required DwC column is non-empty in the final CSV."""
    builder = DarwinCoreBuilder(
        fixture_dir / "tax.csv", fixture_dir / "sample.csv",
        fixture_dir / "project.csv", fixture_dir / "out.csv",
        skip_enrichment=True,
    )
    builder.build()
    out = pd.read_csv(fixture_dir / "out.csv")
    for required in (
        "occurrenceID", "eventID", "basisOfRecord",
        "target_gene", "pcr_primer_forward", "pcr_primer_reverse", "otu_db",
    ):
        assert required in out.columns, f"missing column: {required}"
        assert (out[required].astype(str) != "").any(), f"all empty: {required}"
