"""Unit tests for the obitab -> shared-schema mapping in EcotagRunner.

OBITools obitab output has NO kingdom/phylum/class columns, and its
order/family/genus/species columns hold numeric NCBI taxids -- the scientific
names live in order_name/family_name/genus_name/species_name. Before the fix,
EcotagRunner.link_with_abundance_table handed that raw schema straight to the
shared post-processor, which found no usable rank columns and silently marked
every ecotag OTU 'Unassigned'.

These tests synthesise the obitab schema (no real OBITools needed) and confirm:
  - the *_name columns are mapped to the shared rank names, so real taxonomy
    survives into the linked CSV (fails before the fix, passes after);
  - a schema with no rank-name column at all raises rather than silently
    zeroing out taxonomy.
"""

from __future__ import annotations

import stat
from pathlib import Path
from typing import List

import pandas as pd
import pytest

from seednap.steps.taxonomic_assignment.ecotag_runner import (
    EcotagError,
    EcotagRunner,
    _REQUIRED_OBITOOLS,
)


def _make_runner(tmp_path: Path) -> EcotagRunner:
    """Build an EcotagRunner with stub binaries (link_with_abundance_table
    does not invoke the binaries, so stubs are enough)."""
    bin_dir = tmp_path / "obitools_bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    for name in _REQUIRED_OBITOOLS:
        f = bin_dir / name
        f.write_text("#!/bin/sh\necho stub\n")
        f.chmod(f.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return EcotagRunner(bin_dir=bin_dir)


def _write_obitab_tsv(path: Path, sequences: List[str]) -> None:
    """Write a TSV matching the real obitab schema: lowercase `sequence`,
    numeric taxid columns family/genus/order/species, and the *_name columns
    holding the scientific names. NO kingdom/phylum/class columns."""
    n = len(sequences)
    df = pd.DataFrame({
        "id": [f"seq_{i}" for i in range(n)],
        "sequence": sequences,
        "order": [123, 456][:n],
        "order_name": ["Perciformes", "Primates"][:n],
        "family": [11, 22][:n],
        "family_name": ["Percidae", "Hominidae"][:n],
        "genus": [1, 2][:n],
        "genus_name": ["Perca", "Homo"][:n],
        "species": [9, 8][:n],
        "species_name": ["Perca_fluviatilis", "Homo_sapiens"][:n],
        "rank": ["species", "species"][:n],
        "scientific_name": ["Perca fluviatilis", "Homo sapiens"][:n],
    })
    df.to_csv(path, sep="\t", index=False)


def _write_abundance(path: Path, sequences: List[str]) -> None:
    df = pd.DataFrame(
        {"S1": list(range(10, 10 + len(sequences)))},
        index=pd.Index(sequences, name="sequence"),
    )
    df.to_csv(path)


def test_obitab_name_columns_mapped_to_ranks(tmp_path: Path) -> None:
    """The *_name scientific names survive into the linked CSV (not Unassigned)."""
    seqs = ["AAAAAAAA", "CCCCCCCC"]
    tsv = tmp_path / "query_ecotag.tsv"
    abd = tmp_path / "abundance.csv"
    out = tmp_path / "out.csv"
    _write_obitab_tsv(tsv, seqs)
    _write_abundance(abd, seqs)

    runner = _make_runner(tmp_path)
    runner.link_with_abundance_table(tsv, abd, out)

    result = pd.read_csv(out)
    assert len(result) == 2
    by_seq = result.set_index("Sequence")

    # Scientific names from *_name columns, not the numeric taxids.
    assert by_seq.loc["AAAAAAAA", "species"] == "Perca_fluviatilis"
    assert by_seq.loc["AAAAAAAA", "genus"] == "Perca"
    assert by_seq.loc["AAAAAAAA", "family"] == "Percidae"
    assert by_seq.loc["AAAAAAAA", "order"] == "Perciformes"
    assert by_seq.loc["CCCCCCCC", "species"] == "Homo_sapiens"

    # No row should be fully Unassigned (the pre-fix failure mode).
    from seednap.utils.taxonomy import UNASSIGNED_LABEL

    assert not (result["species"] == UNASSIGNED_LABEL).all()


def test_obitab_contaminant_flag_uses_name_column(tmp_path: Path) -> None:
    """Contaminant flagging keys off the mapped species name, not the taxid."""
    seqs = ["AAAAAAAA", "CCCCCCCC"]
    tsv = tmp_path / "query_ecotag.tsv"
    abd = tmp_path / "abundance.csv"
    out = tmp_path / "out.csv"
    _write_obitab_tsv(tsv, seqs)
    _write_abundance(abd, seqs)

    runner = _make_runner(tmp_path)
    runner.link_with_abundance_table(
        tsv, abd, out, contaminants=["Homo_sapiens"]
    )
    result = pd.read_csv(out)
    flagged = result[result["is_contaminant_candidate"]]
    assert len(flagged) == 1
    assert flagged.iloc[0]["species"] == "Homo_sapiens"


def test_no_rank_name_column_raises(tmp_path: Path) -> None:
    """A taxonomy table with no rank-name column raises instead of silently
    marking every OTU Unassigned (guards against obitab schema drift)."""
    seqs = ["AAAAAAAA"]
    tsv = tmp_path / "query_ecotag.tsv"
    abd = tmp_path / "abundance.csv"
    out = tmp_path / "out.csv"
    # Only numeric taxid columns + sequence, NO *_name columns.
    pd.DataFrame({
        "id": ["seq_0"],
        "sequence": seqs,
        "order": [123],
        "family": [11],
        "genus": [1],
        "species": [9],
    }).to_csv(tsv, sep="\t", index=False)
    _write_abundance(abd, seqs)

    runner = _make_runner(tmp_path)
    with pytest.raises(EcotagError, match="no rank-name column"):
        runner.link_with_abundance_table(tsv, abd, out)
