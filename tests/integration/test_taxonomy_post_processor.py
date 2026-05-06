"""Integration test for the canonical taxonomy post-processor.

`link_taxonomy_with_abundance` is the shared post-merge step used by ecotag,
DECIPHER, and DADA2 RDP. It is responsible for keeping their outputs honest
in exactly the same way the BLAST + LCA path is honest:

- LEFT-merge from abundance side -> every OTU survives.
- Cascade-null taxonomic ranks -> no orphan-rank rows.
- Empty taxonomy input -> all-Unassigned output, no crash.
- Contaminant flagging -> rows annotated, never deleted.
- Stable BLAST-compatible column order.

These tests run against synthetic CSV fixtures with no external tools
(no DECIPHER, no ecotag, no R), so they are fast and deterministic.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from seednap.utils.taxonomy import (
    CONTAMINANT_FLAG_COL,
    DEFAULT_RANK_COLUMNS,
    UNASSIGNED_LABEL,
    link_taxonomy_with_abundance,
)


@pytest.fixture
def abundance_csv(tmp_path: Path) -> Path:
    """4-OTU abundance table, sequences as index, 2 samples."""
    p = tmp_path / "abundance.csv"
    df = pd.DataFrame(
        {"S1": [10, 20, 30, 40], "S2": [5, 15, 25, 35]},
        index=["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"],
    )
    df.index.name = "sequence"
    df.to_csv(p)
    return p


@pytest.fixture
def taxonomy_csv_partial(tmp_path: Path) -> Path:
    """Taxonomy hits for 2 of the 4 OTUs only (the others should become Unassigned)."""
    p = tmp_path / "tax.csv"
    df = pd.DataFrame(
        {
            "sequence": ["AAAAAAAA", "CCCCCCCC"],
            "kingdom": ["Metazoa", "Metazoa"],
            "phylum": ["Chordata", "Chordata"],
            "class": ["Actinopteri", "Mammalia"],
            "order": ["Perciformes", "Primates"],
            "family": ["Percidae", "Hominidae"],
            "genus": ["Perca", "Homo"],
            "species": ["Perca_fluviatilis", "Homo_sapiens"],
        }
    )
    df.to_csv(p, index=False)
    return p


def test_left_merge_keeps_all_otus(abundance_csv: Path, taxonomy_csv_partial: Path, tmp_path: Path) -> None:
    """B1: OTUs without taxonomy hits survive the merge as Unassigned."""
    out = tmp_path / "out.csv"
    link_taxonomy_with_abundance(taxonomy_csv_partial, abundance_csv, out)
    result = pd.read_csv(out)
    assert len(result) == 4  # all 4 OTUs survive

    # OTU_3 (GGGGGGGG) and OTU_4 (TTTTTTTT) should be all-Unassigned
    unassigned_otus = result[result["kingdom"] == UNASSIGNED_LABEL]
    assert len(unassigned_otus) == 2
    for rank in DEFAULT_RANK_COLUMNS:
        assert (unassigned_otus[rank] == UNASSIGNED_LABEL).all()


def test_cascade_null_no_orphan_ranks(tmp_path: Path) -> None:
    """B3: a coarse Unassigned rank forces every finer rank to Unassigned."""
    abundance = tmp_path / "abd.csv"
    pd.DataFrame({"S1": [10]}, index=pd.Index(["AAAA"], name="sequence")).to_csv(abundance)

    # Taxonomy where phylum is missing but class through species are populated (orphan case)
    tax = tmp_path / "tax.csv"
    pd.DataFrame({
        "sequence": ["AAAA"],
        "kingdom": ["Metazoa"],
        "phylum": [None],            # missing -> should cascade
        "class": ["Mammalia"],       # would be orphan without cascade
        "order": ["Primates"],
        "family": ["Hominidae"],
        "genus": ["Homo"],
        "species": ["Homo_sapiens"],
    }).to_csv(tax, index=False)

    out = tmp_path / "out.csv"
    link_taxonomy_with_abundance(tax, abundance, out)
    row = pd.read_csv(out).iloc[0]
    assert row["kingdom"] == "Metazoa"
    assert row["phylum"] == UNASSIGNED_LABEL
    # Cascade: every rank below phylum must also be Unassigned
    for finer in ("class", "order", "family", "genus", "species"):
        assert row[finer] == UNASSIGNED_LABEL, (
            f"orphan-rank: {finer}={row[finer]} despite phylum=Unassigned"
        )


def test_empty_taxonomy_file_all_unassigned(abundance_csv: Path, tmp_path: Path) -> None:
    """B4: empty taxonomy file -> all OTUs Unassigned, no crash."""
    empty = tmp_path / "empty.csv"
    empty.write_text("")  # zero bytes
    out = tmp_path / "out.csv"
    link_taxonomy_with_abundance(empty, abundance_csv, out)
    result = pd.read_csv(out)
    assert len(result) == 4
    for rank in DEFAULT_RANK_COLUMNS:
        assert (result[rank] == UNASSIGNED_LABEL).all()


def test_contaminant_flagged_not_deleted(
    abundance_csv: Path, taxonomy_csv_partial: Path, tmp_path: Path
) -> None:
    """B5: contaminant species are flagged in is_contaminant_candidate, never removed."""
    out = tmp_path / "out.csv"
    link_taxonomy_with_abundance(
        taxonomy_csv_partial, abundance_csv, out,
        contaminants=["Homo_sapiens"],
    )
    result = pd.read_csv(out)
    assert len(result) == 4  # nothing deleted

    # Exactly one row flagged (the Homo_sapiens OTU)
    flagged = result[result[CONTAMINANT_FLAG_COL]]
    assert len(flagged) == 1
    assert flagged.iloc[0]["species"] == "Homo_sapiens"

    # Non-contaminants are False, not NaN
    non_flagged = result[~result[CONTAMINANT_FLAG_COL]]
    assert len(non_flagged) == 3


def test_blast_compatible_schema(
    abundance_csv: Path, taxonomy_csv_partial: Path, tmp_path: Path
) -> None:
    """B6: output column order is stable and matches the BLAST output schema."""
    out = tmp_path / "out.csv"
    link_taxonomy_with_abundance(
        taxonomy_csv_partial, abundance_csv, out,
        contaminants=["Homo_sapiens"],
    )
    result = pd.read_csv(out)
    expected_prefix = ["ASV_ID", "pident", *DEFAULT_RANK_COLUMNS, CONTAMINANT_FLAG_COL]
    assert list(result.columns)[: len(expected_prefix)] == expected_prefix
    # Sample columns come next, Sequence last
    assert "S1" in result.columns and "S2" in result.columns
    assert list(result.columns)[-1] == "Sequence"


def test_otu_id_format_and_sort(
    abundance_csv: Path, taxonomy_csv_partial: Path, tmp_path: Path
) -> None:
    """ASV_IDs are OTU_<n> generated from row order, output sorted by ASV number."""
    out = tmp_path / "out.csv"
    link_taxonomy_with_abundance(taxonomy_csv_partial, abundance_csv, out)
    result = pd.read_csv(out)
    assert list(result["ASV_ID"]) == ["OTU_1", "OTU_2", "OTU_3", "OTU_4"]


def test_pident_mapped_from_bootstrap_min(abundance_csv: Path, tmp_path: Path) -> None:
    """G/A2: DADA2 RDP path passes bootstrap_min as the pident column.

    Replicates the format the patched taxo_dada2_marker.R produces and
    verifies the post-processor surfaces it under `pident`.
    """
    tax = tmp_path / "tax.csv"
    pd.DataFrame({
        "sequence": ["AAAAAAAA", "CCCCCCCC"],
        "kingdom": ["Metazoa", "Metazoa"],
        "phylum": ["Chordata", "Chordata"],
        "class": ["Actinopteri", "Mammalia"],
        "order": ["Perciformes", "Primates"],
        "family": ["Percidae", "Hominidae"],
        "genus": ["Perca", "Homo"],
        "species": ["Perca_fluviatilis", "Homo_sapiens"],
        "bootstrap_min": [95.0, 88.0],
    }).to_csv(tax, index=False)

    out = tmp_path / "out.csv"
    link_taxonomy_with_abundance(
        tax, abundance_csv, out, pident_col="bootstrap_min"
    )
    result = pd.read_csv(out)

    # The two assigned OTUs carry the bootstrap as pident; the two unassigned
    # OTUs carry NaN.
    assigned = result[result["kingdom"] != UNASSIGNED_LABEL].sort_values("ASV_ID")
    assert assigned["pident"].tolist() == [95.0, 88.0]
    unassigned = result[result["kingdom"] == UNASSIGNED_LABEL]
    assert unassigned["pident"].isna().all()
