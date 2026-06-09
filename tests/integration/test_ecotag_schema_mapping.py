"""Integration test for EcotagRunner.link_with_abundance_table schema mapping.

OBITools obitab output has NO kingdom/phylum/class columns, and its
order/family/genus/species columns hold numeric NCBI taxids -- the scientific
names live in order_name/family_name/genus_name/species_name. The shared
taxonomy post-processor (link_taxonomy_with_abundance) expects rank columns
named kingdom..species. Without mapping the obitab schema first, every ecotag
OTU would be silently emitted as fully Unassigned.

These tests synthesise obitab-style TSVs + a tiny abundance CSV; no OBITools
binaries are executed (link_with_abundance_table only reads files), so a stub
bin dir satisfies EcotagRunner construction.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pandas as pd
import pytest

from seednap.steps.taxonomic_assignment.ecotag_runner import (
    EcotagError,
    EcotagRunner,
)


def _stub_runner(tmp_path: Path) -> EcotagRunner:
    """Build an EcotagRunner whose bin_dir holds no-op stubs (never executed here)."""
    bin_dir = tmp_path / "obibin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    for n in ("ecotag", "obiannotate", "obitab"):
        f = bin_dir / n
        f.write_text("#!/bin/sh\necho stub\n")
        f.chmod(f.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return EcotagRunner(bin_dir=bin_dir)


@pytest.fixture
def abundance_csv(tmp_path: Path) -> Path:
    """2-OTU abundance table, sequences as index, 2 samples."""
    p = tmp_path / "abundance.csv"
    df = pd.DataFrame(
        {"S1": [10, 20], "S2": [5, 15]},
        index=["AAAAAAAA", "CCCCCCCC"],
    )
    df.index.name = "sequence"
    df.to_csv(p)
    return p


def test_obitab_schema_mapped_not_all_unassigned(
    tmp_path: Path, abundance_csv: Path
) -> None:
    """obitab *_name columns must reach the output as kingdom..species ranks."""
    # obitab TSV: lowercase `sequence`, numeric-taxid rank columns, *_name names.
    tax = tmp_path / "query_ecotag.tsv"
    pd.DataFrame(
        {
            "sequence": ["AAAAAAAA", "CCCCCCCC"],
            "order": [123, 456],          # numeric taxids (must be dropped)
            "order_name": ["Perciformes", "Primates"],
            "family": [789, 321],
            "family_name": ["Percidae", "Hominidae"],
            "genus": [11, 22],
            "genus_name": ["Perca", "Homo"],
            "species": [111, 222],
            "species_name": ["Perca_fluviatilis", "Homo_sapiens"],
            "rank": ["species", "species"],
            "scientific_name": ["Perca fluviatilis", "Homo sapiens"],
        }
    ).to_csv(tax, sep="\t", index=False)

    runner = _stub_runner(tmp_path)
    out = tmp_path / "ecotag_linked.csv"
    runner.link_with_abundance_table(
        taxonomy_tsv=tax, abundance_csv=abundance_csv, output_csv=out
    )

    result = pd.read_csv(out)
    assert len(result) == 2

    # Pre-fix bug: every rank Unassigned. Post-fix: the *_name values survive.
    species = set(result["species"])
    assert "Perca_fluviatilis" in species
    assert "Homo_sapiens" in species
    # Names, not numeric taxids, in the rank columns.
    assert set(result["genus"]) == {"Perca", "Homo"}
    # No row collapsed to fully Unassigned (order is the coarsest mapped rank).
    assert (result["order"] != "Unassigned").all()
    # The full BLAST-compatible 7-rank schema is kept (GBIF formatting needs it),
    # but obitab does not resolve kingdom/phylum/class, so they are Unassigned
    # placeholders. The key correctness property: the absent coarse ranks must
    # NOT cascade-null the resolved order..species ranks (the pre-fix bug).
    for placeholder in ("kingdom", "phylum", "class"):
        assert placeholder in result.columns
        assert (result[placeholder] == "Unassigned").all()


def test_no_rank_name_column_raises(tmp_path: Path, abundance_csv: Path) -> None:
    """A drifted obitab schema with no rank-name column must fail loudly."""
    tax = tmp_path / "query_ecotag.tsv"
    # Only numeric-taxid rank columns + metadata, no *_name columns at all.
    pd.DataFrame(
        {
            "sequence": ["AAAAAAAA", "CCCCCCCC"],
            "order": [123, 456],
            "family": [789, 321],
            "rank": ["species", "species"],
        }
    ).to_csv(tax, sep="\t", index=False)

    runner = _stub_runner(tmp_path)
    out = tmp_path / "ecotag_linked.csv"
    with pytest.raises(EcotagError, match=r"no rank-name column"):
        runner.link_with_abundance_table(
            taxonomy_tsv=tax, abundance_csv=abundance_csv, output_csv=out
        )
