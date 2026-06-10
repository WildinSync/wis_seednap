"""eventID normalization on the sample-metadata join (cli-fixes).

Regression: when a legacy taxonomy table carries R make.names()-dotted eventIDs
(DAR.2023.0025) while the sample-metadata sheet uses the canonical dashed form
(DAR-2023-0025), the raw `on="eventID"` left join matched ZERO rows. That left
eventDate / decimalLatitude / decimalLongitude / env_medium blank in 100% of the
GBIF output while the command still exited 0 -- a silent loss of all spatial /
temporal data.

These tests assert that:
  - dotted-vs-dashed eventIDs now join, coords/date/env populate, and the output
    eventID is written in the canonical dashed form;
  - a genuinely non-matching metadata set RAISES instead of silently emitting
    empty location/date fields;
  - a partial mismatch emits a [WARN] (no-silent-fallbacks) rather than passing
    quietly.

Uses synthetic CSVs and skip_enrichment=True to stay fast and offline.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from seednap.steps.formatting.darwincore_builder import DarwinCoreBuilder


def _tax_row(event_id: str, seq: str, reads: int) -> dict:
    return {
        "kingdom": "Metazoa", "phylum": "Chordata", "class": "Actinopteri",
        "order": "Perciformes", "family": "Percidae", "genus": "Perca",
        "species": "Perca_fluviatilis", "taxon": "Perca_fluviatilis",
        "rank": "species", "sequence": seq, "nb_reads": reads,
        "eventID": event_id,
    }


def _write_taxonomy(path: Path, event_ids: list[str]) -> None:
    """One taxonomy row per eventID, each with a distinct sequence."""
    rows = [
        _tax_row(eid, f"ACGTACGT{i:04d}".replace(" ", ""), 100 + i)
        for i, eid in enumerate(event_ids)
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def _write_sample_meta(path: Path, event_ids: list[str]) -> None:
    """Dashed (canonical) sample-metadata rows, distinct coords per eventID."""
    rows = [
        {
            "eventID": eid,
            "decimalLatitude": 46.0 + i * 0.01,
            "decimalLongitude": 7.5 + i * 0.01,
            "eventDate": "2025.06.15",
            "env_medium": "river",
            "samp_size": "1L",
            "depth": 0.5,
            "size_frac": "0.22um",
        }
        for i, eid in enumerate(event_ids)
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def _write_project_meta(path: Path, marker: str = "Teleo") -> None:
    pd.DataFrame([{
        "marker": marker,
        "recordedby": "J. Smith",
        "seqmet": "MiSeq PE 2x150",
        "identificationRemarks": "BLAST + LCA",
        "identificationReferences": "10.1038/nmeth.3869",
        "otu_seq_comp_appr": "SWARM d=1",
        "otu_db": "CRABS MitoFish 2025",
        "chimera_check": "UCHIME de novo",
    }]).to_csv(path, index=False)


def test_dotted_taxonomy_joins_dashed_metadata(tmp_path: Path) -> None:
    """Make.names-dotted taxonomy eventIDs join the dashed metadata: coords/date
    populate and the output eventID is written in the canonical dashed form."""
    _write_taxonomy(tmp_path / "tax.csv", ["DAR.2023.0025", "DAR.2023.0026"])
    _write_sample_meta(tmp_path / "sample.csv", ["DAR-2023-0025", "DAR-2023-0026"])
    _write_project_meta(tmp_path / "project.csv")

    builder = DarwinCoreBuilder(
        tmp_path / "tax.csv", tmp_path / "sample.csv",
        tmp_path / "project.csv", tmp_path / "out.csv",
        skip_enrichment=True,
    )
    builder.build()
    out = pd.read_csv(tmp_path / "out.csv")

    # Location / date / env populated for every row (not blank in 100% of rows).
    assert out["eventDate"].notna().all()
    assert out["decimalLatitude"].notna().all()
    assert out["decimalLongitude"].notna().all()
    assert (out["env_medium"].astype(str) != "").all()

    # Output eventID is the canonical dashed form, not the dotted source form.
    event_ids = set(out["eventID"].astype(str))
    assert event_ids == {"DAR-2023-0025", "DAR-2023-0026"}
    assert not any("." in e for e in event_ids)

    # occurrenceID carries the canonical eventID too (marker:eventID:hash).
    for oid in out["occurrenceID"].astype(str):
        assert "DAR-2023-002" in oid
        assert "DAR.2023" not in oid


def test_zero_match_metadata_raises(tmp_path: Path) -> None:
    """A genuinely non-matching metadata set must RAISE, not silently emit empty
    location/date fields while exiting 0."""
    _write_taxonomy(tmp_path / "tax.csv", ["DAR.2023.0025"])
    # Completely unrelated eventID -- no normalization can reconcile it.
    _write_sample_meta(tmp_path / "sample.csv", ["GREINA-2024-9999"])
    _write_project_meta(tmp_path / "project.csv")

    builder = DarwinCoreBuilder(
        tmp_path / "tax.csv", tmp_path / "sample.csv",
        tmp_path / "project.csv", tmp_path / "out.csv",
        skip_enrichment=True,
    )
    with pytest.raises(ValueError, match="matched ZERO"):
        builder.build()
    # Nothing must have been written when the join fails.
    assert not (tmp_path / "out.csv").exists()


def test_partial_match_warns(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A partial mismatch (one eventID has no metadata row) emits a [WARN] naming
    the unmatched eventID rather than passing silently."""
    _write_taxonomy(tmp_path / "tax.csv", ["DAR.2023.0025", "DAR.2023.0099"])
    # Only the first sample has metadata; the second is unmatched.
    _write_sample_meta(tmp_path / "sample.csv", ["DAR-2023-0025"])
    _write_project_meta(tmp_path / "project.csv")

    builder = DarwinCoreBuilder(
        tmp_path / "tax.csv", tmp_path / "sample.csv",
        tmp_path / "project.csv", tmp_path / "out.csv",
        skip_enrichment=True,
    )
    builder.build()
    captured = capsys.readouterr()
    assert "[WARN]" in captured.out
    assert "_merge_sample_metadata" in captured.out
    # The unmatched eventID is named so the operator can see what dropped.
    assert "DAR-2023-0099" in captured.out or "DAR.2023.0099" in captured.out

    out = pd.read_csv(tmp_path / "out.csv")
    # The matched row has coordinates; the unmatched one is blank (but warned).
    matched = out[out["eventID"] == "DAR-2023-0025"]
    assert matched["decimalLatitude"].notna().all()
