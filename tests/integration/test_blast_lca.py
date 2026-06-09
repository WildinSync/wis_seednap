"""Integration tests for BlastTaxonomicAssigner correctness.

Builds a tiny synthetic BLAST TSV + reference DB + ASV count table, then
exercises both LCA strategies (the default 'cascade' resolver and the
'collapsed_taxonomy' resolver) plus the surrounding merge/contaminant logic.
The behaviors covered:

- Unassigned survival: OTUs with no BLAST hits survive the merge as
  'Unassigned' rows rather than being silently dropped (also the empty-TSV
  edge case).
- Top-bitscore band collapse: equally-good (near-best) hits inside the
  bitscore band that disagree on a rank collapse to their LCA, while hits
  outside the band are ignored and the best hit's call survives.
- lca_pident_delta floor: a hit inside the bitscore band but more than
  lca_pident_delta below the best identity is excluded so it cannot
  over-collapse a confident high-identity call.
- Cascade nulling: a hit below the threshold for rank R nulls R and every
  finer rank, keeping coarser ranks.
- collapsed_taxonomy resolver: its identity window excludes lower-%id
  off-targets, equally-good disagreeing hits collapse to their LCA, and a
  hit below lca_pid yields Unassigned.
- LCA-algorithm factory: 'cascade' and 'collapsed_taxonomy' build their
  resolvers, 'fishbase_tiered' raises NotImplementedError, unknown values
  raise ValueError.
- CRABS 'NA'-sentinel handling for both resolvers: a missing-rank 'NA' is
  treated as missing, never as a taxon and never leaked as the literal
  string 'NA'.
- Contaminant flagging: contaminant species get is_contaminant_candidate
  set without being removed from the table.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from seednap.steps.taxonomic_assignment.blast_runner import (
    BlastTaxonomicAssigner,
    CollapsedTaxonomyLCAResolver,
)


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
>REF5\tMetazoa;Chordata;Actinopteri;NA;Percidae;Perca;Perca_fluviatilis
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
    """OTUs with no BLAST hits appear in the output as Unassigned."""
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
    """When equally-good (same-identity) hits disagree on family, LCA collapses below class.

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
    """Regression: a hit inside the bitscore band but >lca_pident_delta below the best
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
    """Factory: 'cascade' and 'collapsed_taxonomy' build their resolvers; 'fishbase_tiered'
    raises (not implemented); an unknown value errors."""
    from seednap.steps.taxonomic_assignment.blast_runner import (
        BlastLCAResolver,
        BlastTaxonomicAssigner,
        CollapsedTaxonomyLCAResolver,
    )

    assert isinstance(
        BlastTaxonomicAssigner._make_lca_resolver("cascade", top_bitscore_pct=10.0, lca_pident_delta=1.0),
        BlastLCAResolver,
    )
    assert isinstance(
        BlastTaxonomicAssigner._make_lca_resolver(
            "collapsed_taxonomy", top_bitscore_pct=10.0, lca_pident_delta=1.0, lca_pid=90.0, lca_diff=1.0
        ),
        CollapsedTaxonomyLCAResolver,
    )
    with pytest.raises(NotImplementedError, match="fishbase"):
        BlastTaxonomicAssigner._make_lca_resolver(
            "fishbase_tiered", top_bitscore_pct=10.0, lca_pident_delta=1.0
        )
    with pytest.raises(ValueError, match="unknown"):
        BlastTaxonomicAssigner._make_lca_resolver(
            "bogus", top_bitscore_pct=10.0, lca_pident_delta=1.0
        )


def test_collapsed_taxonomy_window_and_collapse(fixture_dir: Path) -> None:
    """collapsed_taxonomy: the identity window excludes a lower-%id off-target (keeps the
    confident call), and equally-good disagreeing hits collapse to their LCA."""
    blast_tsv = fixture_dir / "blast.tsv"
    # OTU_1: Perca 100% + Cyprinus 98.4% -> window [99,100] drops Cyprinus -> Perca survives.
    # OTU_2: Perca 100% + Cyprinus 100% -> both in window, disagree at order -> collapse to class.
    _write_blast_tsv(
        blast_tsv,
        [
            ("OTU_1", "REF1", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
            ("OTU_1", "REF4", 98.4, 64, 1, 0, 1, 64, 1, 64, 1e-28, 113),
            ("OTU_2", "REF1", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
            ("OTU_2", "REF4", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 113),
        ],
    )
    ass = BlastTaxonomicAssigner(
        reference_fasta=fixture_dir / "ref.fasta",
        lca_algorithm="collapsed_taxonomy", lca_pid=90.0, lca_diff=1.0,
    )
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv, asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    o1 = df[df["ASV_ID"] == "OTU_1"].iloc[0]
    assert o1["genus"] == "Perca" and o1["species"] == "Perca_fluviatilis"  # off-target excluded
    o2 = df[df["ASV_ID"] == "OTU_2"].iloc[0]
    assert o2["class"] == "Actinopteri" and o2["genus"] == "Unassigned"  # tie -> LCA collapse


def _lineage_row(otu: str, ref: str, pident: float, **ranks: object) -> dict:
    base = dict(qseqid=otu, sseqid=ref, pident=pident, bitscore=120, blast_rank=1)
    base.update(ranks)
    return base


def test_collapsed_taxonomy_na_sentinel_is_not_a_taxon() -> None:
    """The CRABS missing-rank sentinel 'NA' must be treated as missing, not as a taxon:
    it must neither over-collapse an otherwise-agreeing call nor leak a literal 'NA'."""
    resolver = CollapsedTaxonomyLCAResolver(lca_pid=90.0, lca_diff=1.0)
    full = dict(kingdom="Metazoa", phylum="Chordata", **{"class": "Actinopteri"},
                family="Centropomidae", genus="Lates", species="Lates_calcarifer")
    # Two equally-good hits for the SAME species; one ref carries order='NA'. Correct LCA is
    # the full species -- 'NA' is missing, not a disagreeing order.
    g = pd.DataFrame([
        _lineage_row("OTU_X", "A", 100.0, order="Perciformes", **full),
        _lineage_row("OTU_X", "B", 100.0, order="NA", **full),
    ])
    kept = resolver.resolve_ambiguous_hits(g.copy())
    row = kept[kept["keep_for_analysis"]].iloc[0]
    assert row["species"] == "Lates_calcarifer"  # not over-collapsed
    assert row["genus"] == "Lates" and row["family"] == "Centropomidae"
    assert row["order"] == "Perciformes"  # the one real order wins; 'NA' ignored

    # A single hit whose order is 'NA' must yield a NULL order, never the literal string 'NA'.
    g1 = pd.DataFrame([_lineage_row("OTU_Y", "B", 100.0, order="NA", **full)])
    row1 = resolver.resolve_ambiguous_hits(g1.copy())
    keep1 = row1[row1["keep_for_analysis"]].iloc[0]
    assert keep1["order"] is None
    assert keep1["species"] == "Lates_calcarifer"


def test_formatter_normalizes_na_for_cascade(fixture_dir: Path) -> None:
    """Option A: the formatter maps the CRABS 'NA' sentinel to missing, so the CASCADE path
    (the shipping default) neither leaks a literal 'NA' on a lone hit nor over-collapses when
    one in-band hit carries 'NA' at a rank where the others agree. REF5 is Perca_fluviatilis
    with order='NA'."""
    blast_tsv = fixture_dir / "blast.tsv"
    _write_blast_tsv(
        blast_tsv,
        [
            # OTU_1: lone hit with order='NA' -> order must become Unassigned, never literal 'NA'.
            ("OTU_1", "REF5", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
            # OTU_2: REF1 (order=Perciformes) + REF5 (order='NA'), same species, both best ->
            # 'NA' must not count as a disagreeing order; cascade keeps the full Perciformes call.
            ("OTU_2", "REF1", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
            ("OTU_2", "REF5", 100.0, 64, 0, 0, 1, 64, 1, 64, 1e-30, 119),
        ],
    )
    ass = BlastTaxonomicAssigner(
        reference_fasta=fixture_dir / "ref.fasta",
        threshold_species=99, threshold_genus=96, threshold_family=90,
        threshold_order=80, threshold_class=70,
    )
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv, asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    assert int((df[ranks].astype(str) == "NA").sum().sum()) == 0  # no literal 'NA' anywhere
    o1 = df[df["ASV_ID"] == "OTU_1"].iloc[0]
    assert o1["order"] == "Unassigned"  # the 'NA' sentinel became missing, not a taxon
    assert o1["genus"] == "Perca" and o1["species"] == "Perca_fluviatilis"
    o2 = df[df["ASV_ID"] == "OTU_2"].iloc[0]
    assert o2["order"] == "Perciformes"  # not over-collapsed by the 'NA'-bearing hit
    assert o2["species"] == "Perca_fluviatilis"


def test_collapsed_taxonomy_below_floor_is_unassigned(fixture_dir: Path) -> None:
    """A hit below lca_pid yields Unassigned (no rank assigned)."""
    blast_tsv = fixture_dir / "blast.tsv"
    _write_blast_tsv(blast_tsv, [("OTU_1", "REF1", 85.0, 64, 9, 0, 1, 64, 1, 64, 1e-20, 100)])
    ass = BlastTaxonomicAssigner(
        reference_fasta=fixture_dir / "ref.fasta", lca_algorithm="collapsed_taxonomy",
        lca_pid=90.0, lca_diff=1.0,
    )
    df = ass.assign_taxonomy(
        blast_tsv=blast_tsv, asv_count_csv=fixture_dir / "abundance.csv",
        asv_fasta=fixture_dir / "query.fasta",
    )
    o1 = df[df["ASV_ID"] == "OTU_1"].iloc[0]
    assert o1["kingdom"] == "Unassigned" and o1["species"] == "Unassigned"


def test_top_bitscore_band_excludes_far_hit(fixture_dir: Path) -> None:
    """Negative case: hits outside the band don't trigger LCA collapse.

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
    """Below threshold for rank R, R and every finer rank are nulled."""
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
    """Contaminant species are flagged in is_contaminant_candidate, never removed."""
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
