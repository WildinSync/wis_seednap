"""Integration test for BlastTaxonomicAssigner correctness.

Builds a tiny synthetic BLAST TSV + reference DB + ASV count table, then
verifies the four behaviors that Commits A, B, and C had to get right:

1. (I-1) OTUs with no BLAST hits survive the merge as 'Unassigned' rows.
2. (I-2) The top-bitscore band collapses near-best disagreeing hits via LCA
   (a single near-best hit alone does not lock in a confident species call
   when the runner-up disagrees).
3. (I-3) Cascade nulling: hits below a rank threshold null that rank and
   every finer rank.
4. (I-7) Candidate contaminants get is_contaminant_candidate=True without
   being deleted.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from seednap.steps.taxonomic_assignment.blast_runner import BlastTaxonomicAssigner


# Reference DB headers in CRABS format: >ACC<TAB>kingdom;phylum;class;order;family;genus;species
REF_FASTA = """\
>REF1\tMetazoa;Chordata;Actinopteri;Perciformes;Percidae;Perca;Perca_fluviatilis
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>REF2\tMetazoa;Chordata;Actinopteri;Perciformes;Percidae;Perca;Perca_flavescens
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>REF3\tMetazoa;Chordata;Mammalia;Primates;Hominidae;Homo;Homo_sapiens
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>REF4\tMetazoa;Chordata;Actinopteri;Cypriniformes;Cyprinidae;Cyprinus;Cyprinus_carpio
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
"""

# OTU sequences (each unique to avoid asv_sequences/abundance merge fan-out).
QUERY_FASTA = """\
>OTU_1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>OTU_2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
>OTU_3
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>OTU_4
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
"""


def _abundance_csv(query_fasta: str, samples=("S1", "S2")) -> str:
    """Make an abundance CSV (sequences as rows, samples as columns)."""
    rows = []
    for line in query_fasta.strip().splitlines():
        if line.startswith(">"):
            continue
        rows.append(line.strip())
    header = "sequence," + ",".join(samples)
    body_lines = []
    for i, seq in enumerate(rows, start=1):
        body_lines.append(f"{seq},{i*10},{i*5}")
    return header + "\n" + "\n".join(body_lines) + "\n"


@pytest.fixture
def fixture_dir(tmp_path: Path) -> Path:
    """Write the synthetic ref DB, query FASTA, and ASV count CSV into tmp_path."""
    ref = tmp_path / "ref.fasta"
    ref.write_text(REF_FASTA)
    query = tmp_path / "query.fasta"
    query.write_text(QUERY_FASTA)
    counts = tmp_path / "abundance.csv"
    counts.write_text(_abundance_csv(QUERY_FASTA))
    return tmp_path


def _write_blast_tsv(path: Path, rows: list[tuple]) -> None:
    """Write a BLAST -outfmt 6 TSV with 12 columns (no qseq/sseq)."""
    cols_order = (
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    )
    lines = []
    for r in rows:
        assert len(r) == len(cols_order), f"need {len(cols_order)} cols, got {len(r)}"
        lines.append("\t".join(str(x) for x in r))
    path.write_text("\n".join(lines) + "\n")


def test_unassigned_otus_survive_merge(fixture_dir: Path) -> None:
    """I-1: OTUs with no BLAST hits appear in the output as Unassigned."""
    blast_tsv = fixture_dir / "blast.tsv"
    # Only OTU_1 gets a hit. OTU_2, OTU_3, OTU_4 should appear as Unassigned.
    _write_blast_tsv(
        blast_tsv,
        [("OTU_1", "REF1", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119)],
    )
    ass = BlastTaxonomicAssigner(reference_fasta=fixture_dir / "ref.fasta")
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv,
        asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    # All 4 OTUs survive (none silently dropped)
    assert len(df) == 4
    # OTU_1 is assigned, the others are Unassigned at every rank
    assigned = df[df["ASV_ID"] == "OTU_1"].iloc[0]
    assert assigned["species"] == "Perca_fluviatilis"
    for unassigned_id in ("OTU_2", "OTU_3", "OTU_4"):
        row = df[df["ASV_ID"] == unassigned_id].iloc[0]
        for rank in ("kingdom", "phylum", "class", "order", "family", "genus", "species"):
            assert row[rank] == "Unassigned", (
                f"{unassigned_id}/{rank} = {row[rank]!r}, expected 'Unassigned'"
            )


def test_top_bitscore_band_collapses_near_best_disagreement(fixture_dir: Path) -> None:
    """I-2: When equally-good (same-identity) hits disagree on family, LCA collapses below class.

    OTU_1 has two 100%-identity hits, Perca_fluviatilis (bitscore 119) and
    Cyprinus_carpio (bitscore 113, within MEGAN's 10% band). They agree on class
    (Actinopteri) but disagree on order/family/genus/species. Both are at the top
    identity, so the lca_pident_delta floor keeps both and the combined LCA row keeps
    class but nulls order onwards.
    """
    blast_tsv = fixture_dir / "blast.tsv"
    _write_blast_tsv(
        blast_tsv,
        [
            ("OTU_1", "REF1", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
            ("OTU_1", "REF4", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-28, 113),
        ],
    )
    ass = BlastTaxonomicAssigner(
        reference_fasta=fixture_dir / "ref.fasta",
        threshold_species=0.0,
        threshold_genus=0.0,
        threshold_family=0.0,
        threshold_order=0.0,
        threshold_class=0.0,
        top_bitscore_pct=10.0,  # 119 * 0.9 = 107.1, so 113 is in band
    )
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv,
        asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    otu1 = df[df["ASV_ID"] == "OTU_1"].iloc[0]
    # Class agreed (both Actinopteri) so it must be preserved
    assert otu1["class"] == "Actinopteri"
    # Order, family, genus, species all disagreed -> Unassigned via LCA combined row
    assert otu1["order"] == "Unassigned"
    assert otu1["family"] == "Unassigned"
    assert otu1["genus"] == "Unassigned"
    assert otu1["species"] == "Unassigned"


def test_lca_pident_floor_excludes_lower_identity_offtarget(fixture_dir: Path) -> None:
    """D1 regression: a hit inside the bitscore band but >lca_pident_delta below the best
    identity must NOT collapse the LCA. Mirrors the real greina case where a 98.6% peanut
    worm sat in the bitscore band of seven 100% Bos hits and nulled them to kingdom.

    Best hit Perca_fluviatilis at 100% (bitscore 119); a disagreeing Cyprinus_carpio at
    98.4% (bitscore 113, inside the 10% band). With the default delta=1.0 the 98.4% hit is
    excluded, so the confident 100% Perca call survives. With delta=0 (floor disabled) the
    old over-collapse to class returns.
    """
    blast_tsv = fixture_dir / "blast.tsv"
    _write_blast_tsv(
        blast_tsv,
        [
            ("OTU_1", "REF1", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
            ("OTU_1", "REF4", 98.4, 64, 1, 0, 1, 64, 1, 64, 1e-28, 113),
        ],
    )
    kw = dict(
        reference_fasta=fixture_dir / "ref.fasta",
        threshold_species=0.0, threshold_genus=0.0, threshold_family=0.0,
        threshold_order=0.0, threshold_class=0.0, top_bitscore_pct=10.0,
    )
    # Default floor (1.0): the 98.4% off-target is excluded -> Perca survives.
    kept = BlastTaxonomicAssigner(lca_pident_delta=1.0, **kw).assign_taxonomy(
        blast_tsv=blast_tsv, asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    row = kept[kept["ASV_ID"] == "OTU_1"].iloc[0]
    assert row["genus"] == "Perca"
    assert row["species"] == "Perca_fluviatilis"

    # Floor disabled (0.0): the 98.4% hit re-enters the band and collapses to class.
    collapsed = BlastTaxonomicAssigner(lca_pident_delta=0.0, **kw).assign_taxonomy(
        blast_tsv=blast_tsv, asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    crow = collapsed[collapsed["ASV_ID"] == "OTU_1"].iloc[0]
    assert crow["class"] == "Actinopteri"
    assert crow["genus"] == "Unassigned"


def test_lca_algorithm_factory():
    """B1 seam: 'cascade' yields the current resolver; taxid methods raise loudly until a
    taxid-mapped DB is provisioned; an unknown value errors."""
    from seednap.steps.taxonomic_assignment.blast_runner import (
        BlastLCAResolver,
        BlastTaxonomicAssigner,
    )

    r = BlastTaxonomicAssigner._make_lca_resolver(
        "cascade", top_bitscore_pct=10.0, lca_pident_delta=1.0
    )
    assert isinstance(r, BlastLCAResolver)
    for algo in ("collapsed_taxonomy", "fishbase_tiered"):
        with pytest.raises(NotImplementedError, match="taxid"):
            BlastTaxonomicAssigner._make_lca_resolver(
                algo, top_bitscore_pct=10.0, lca_pident_delta=1.0
            )
    with pytest.raises(ValueError, match="unknown"):
        BlastTaxonomicAssigner._make_lca_resolver(
            "bogus", top_bitscore_pct=10.0, lca_pident_delta=1.0
        )


def test_top_bitscore_band_excludes_far_hit(fixture_dir: Path) -> None:
    """I-2 (negative): Hits outside the band don't trigger LCA collapse.

    Same as above but the runner-up bitscore (90) is outside the 10% band
    (threshold = 119 * 0.9 = 107.1). LCA should ignore it and keep the best
    hit's species.
    """
    blast_tsv = fixture_dir / "blast.tsv"
    _write_blast_tsv(
        blast_tsv,
        [
            ("OTU_1", "REF1", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
            ("OTU_1", "REF4", 92.0, 64, 5, 0, 1, 64, 1, 64, 1e-22, 90),
        ],
    )
    ass = BlastTaxonomicAssigner(
        reference_fasta=fixture_dir / "ref.fasta",
        threshold_species=0.0,
        top_bitscore_pct=10.0,
    )
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv,
        asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    otu1 = df[df["ASV_ID"] == "OTU_1"].iloc[0]
    assert otu1["species"] == "Perca_fluviatilis"


def test_cascade_null_below_threshold(fixture_dir: Path) -> None:
    """I-3: Below threshold for rank R, R and every finer rank are nulled."""
    blast_tsv = fixture_dir / "blast.tsv"
    # pident = 87 -> below family threshold (90), above order (80).
    # Expect family/genus/species nulled, order/class kept.
    _write_blast_tsv(
        blast_tsv,
        [("OTU_1", "REF1", 87.0, 64, 8, 0, 1, 64, 1, 64, 1e-25, 100)],
    )
    ass = BlastTaxonomicAssigner(
        reference_fasta=fixture_dir / "ref.fasta",
        threshold_species=99.0,
        threshold_genus=96.0,
        threshold_family=90.0,
        threshold_order=80.0,
        threshold_class=70.0,
    )
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv,
        asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    otu1 = df[df["ASV_ID"] == "OTU_1"].iloc[0]
    assert otu1["kingdom"] == "Metazoa"
    assert otu1["phylum"] == "Chordata"
    assert otu1["class"] == "Actinopteri"
    assert otu1["order"] == "Perciformes"
    # Below family threshold -> family + genus + species nulled
    assert otu1["family"] == "Unassigned"
    assert otu1["genus"] == "Unassigned"
    assert otu1["species"] == "Unassigned"


def test_contaminant_flagged_not_deleted(fixture_dir: Path) -> None:
    """I-7: Contaminant species are flagged in is_contaminant_candidate, never removed."""
    blast_tsv = fixture_dir / "blast.tsv"
    _write_blast_tsv(
        blast_tsv,
        [
            ("OTU_1", "REF1", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
            ("OTU_2", "REF3", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
        ],
    )
    ass = BlastTaxonomicAssigner(
        reference_fasta=fixture_dir / "ref.fasta",
        threshold_species=0.0,
        contaminants=["Homo_sapiens"],
    )
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv,
        asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    # Both OTUs survive
    otu1 = df[df["ASV_ID"] == "OTU_1"].iloc[0]
    otu2 = df[df["ASV_ID"] == "OTU_2"].iloc[0]
    assert otu1["species"] == "Perca_fluviatilis"
    assert otu1["is_contaminant_candidate"] is False or otu1["is_contaminant_candidate"] == 0
    assert otu2["species"] == "Homo_sapiens"
    assert otu2["is_contaminant_candidate"] is True or otu2["is_contaminant_candidate"] == 1


def test_empty_blast_output_all_unassigned(fixture_dir: Path) -> None:
    """Edge case: empty BLAST TSV -> all OTUs Unassigned, no crash."""
    blast_tsv = fixture_dir / "blast.tsv"
    blast_tsv.write_text("")
    ass = BlastTaxonomicAssigner(reference_fasta=fixture_dir / "ref.fasta")
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv,
        asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    assert len(df) == 4  # all four OTUs
    assert (df["species"] == "Unassigned").all()
    assert df["pident"].isna().all()
