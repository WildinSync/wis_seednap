"""Canonical taxonomy schema is the single source of truth.

`seednap.utils.taxonomy` defines the one 7-rank list (TAXONOMIC_RANKS) and the
one missing-taxonomy sentinel set (MISSING_TAXON_VALUES). Every module that used
to re-declare its own rank list or placeholder set must now point at these, so
the schema and missing-value matching can never drift between methods.

These tests assert:

1. The canonical rank list is the expected 7 ranks in order.
2. Every module that consumes the rank list resolves to the same canonical list.
3. The canonical sentinel set is the COMPLETE, case-insensitive union of every
   variant that was previously scattered across blast_runner, the taxonomy
   post-processor, and the HTML report (so no call site newly mis-classifies a
   value), and no real taxon name normalises to a sentinel.
4. The blast_runner / html_report aliases ARE the canonical objects (identity,
   not a copy), so they cannot silently fall out of sync.

No external tools are touched, so these run fast and deterministically.
"""

from __future__ import annotations

import pytest

from seednap.utils.taxonomy import (
    DEFAULT_RANK_COLUMNS,
    MISSING_TAXON_VALUES,
    TAXONOMIC_RANKS,
    is_missing_taxon,
)

# The exact 7-rank lineage, coarse-to-fine.
EXPECTED_RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")

# Every missing-taxonomy literal that used to live in a separate module, kept
# here verbatim so this test fails if the canonical union ever stops covering
# one of them.
#   blast_runner.MISSING_RANK_SENTINELS:        ("", "na", "nan")
#   taxonomy._normalize_unassigned placeholders: {"", "NA", "None", "nan", "NaN", "Na"}
#   html_report._UNASSIGNED:                     {"Unassigned", "unassigned", "", "NA", "nan", "None"}
PREVIOUSLY_SCATTERED_SENTINELS = [
    "", "na", "nan",
    "", "NA", "None", "nan", "NaN", "Na",
    "Unassigned", "unassigned", "", "NA", "nan", "None",
]


def test_canonical_rank_list_is_the_seven_ranks() -> None:
    """TAXONOMIC_RANKS is exactly the 7 ranks in coarse-to-fine order."""
    assert tuple(TAXONOMIC_RANKS) == EXPECTED_RANKS
    # The post-processor alias must be the very same object.
    assert DEFAULT_RANK_COLUMNS is TAXONOMIC_RANKS


def test_every_module_uses_the_canonical_rank_list() -> None:
    """blast_runner, gbif_formatter, and html_report all resolve to the same ranks."""
    from seednap.steps.report import html_report
    from seednap.steps.formatting.gbif_formatter import GBIFFormatter
    from seednap.steps.taxonomic_assignment.blast_runner import (
        BlastLCAResolver,
        BlastOutputFormatter,
        CollapsedTaxonomyLCAResolver,
    )

    canonical = list(TAXONOMIC_RANKS)
    assert BlastOutputFormatter.TAXONOMIC_RANKS == canonical
    assert BlastLCAResolver.TAXONOMIC_RANKS == canonical
    assert CollapsedTaxonomyLCAResolver.TAXONOMIC_RANKS == canonical
    assert GBIFFormatter().taxonomic_ranks == canonical
    assert html_report._RANKS == canonical


def test_canonical_sentinel_set_is_the_complete_union() -> None:
    """MISSING_TAXON_VALUES is the full case-insensitive union, no more, no less."""
    expected = {"", "na", "nan", "none", "unassigned"}
    assert set(MISSING_TAXON_VALUES) == expected
    # Stored lowercased so matching is unambiguous.
    assert all(v == v.lower() for v in MISSING_TAXON_VALUES)


def test_sentinel_set_covers_every_previously_scattered_value() -> None:
    """No call site can newly mis-classify a value it used to treat as missing.

    Each previously-scattered literal must still be classified as missing,
    matched case-insensitively via the canonical predicate.
    """
    for value in PREVIOUSLY_SCATTERED_SENTINELS:
        assert is_missing_taxon(value), f"sentinel {value!r} no longer matches"
        assert value.strip().lower() in MISSING_TAXON_VALUES


def test_is_missing_taxon_handles_case_whitespace_and_nan() -> None:
    """The canonical predicate is case-insensitive, strips whitespace, and
    treats None / NaN as missing."""
    import math

    for value in ("NONE", " na ", "  ", "NaN", "UNASSIGNED", None, float("nan"), math.nan):
        assert is_missing_taxon(value), f"{value!r} should be missing"


def test_real_taxa_are_not_classified_missing() -> None:
    """A genuine taxon name must never be read as a missing-taxonomy sentinel,
    including near-misses like 'Naja' that merely start with a sentinel."""
    for taxon in ("Homo_sapiens", "Perca", "Actinopteri", "Metazoa", "Naja", "Nanger"):
        assert not is_missing_taxon(taxon), f"{taxon!r} wrongly read as missing"


def test_blast_and_html_aliases_are_the_canonical_objects() -> None:
    """The historical names re-export the canonical set by identity, so they
    cannot drift from the single source."""
    from seednap.steps.report import html_report
    from seednap.steps.taxonomic_assignment.blast_runner import MISSING_RANK_SENTINELS

    assert MISSING_RANK_SENTINELS is MISSING_TAXON_VALUES
    assert html_report._UNASSIGNED is MISSING_TAXON_VALUES


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
