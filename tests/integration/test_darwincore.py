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
