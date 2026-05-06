"""Shared taxonomy utilities for linking taxonomy results with abundance tables.

This module provides `link_taxonomy_with_abundance`, the canonical post-processor
for any taxonomy method (DECIPHER, ecotag, DADA2 RDP) that produces a
(sequence, rank_columns) table. It guarantees the same output schema and the
same correctness properties as the BLAST + LCA path:

- LEFT merge from abundance side -> every OTU survives, no silent dropout.
- Missing taxonomy filled with `Unassigned`.
- Cascade-null taxonomic ranks: if a coarse rank is `Unassigned`, every finer
  rank is forced to `Unassigned` too (no orphan-rank rows).
- Empty taxonomy input -> all OTUs marked `Unassigned`, no crash.
- `is_contaminant_candidate` boolean column for downstream filtering.
- Stable BLAST-compatible column order: ASV_ID, pident, kingdom..species,
  is_contaminant_candidate, sample_cols..., Sequence.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Sequence, Union

import pandas as pd

logger = logging.getLogger(__name__)


# Match BLAST's TAXONOMIC_RANKS to keep schemas aligned across methods.
DEFAULT_RANK_COLUMNS: tuple = (
    "kingdom", "phylum", "class", "order", "family", "genus", "species",
)
UNASSIGNED_LABEL = "Unassigned"
CONTAMINANT_FLAG_COL = "is_contaminant_candidate"


def _normalize_unassigned(series: pd.Series, unassigned_label: str = UNASSIGNED_LABEL) -> pd.Series:
    """Treat NaN, empty string, 'NA', 'None', 'nan' as unassigned.

    Reference DBs and R scripts use a mix of these for missing taxonomy; we
    collapse them to a single label so cascade nulling is unambiguous.
    """
    s = series.astype(object).where(series.notna(), unassigned_label)
    s = s.astype(str)
    placeholders = {"", "NA", "None", "nan", "NaN", "Na"}
    return s.where(~s.isin(placeholders), unassigned_label)


def link_taxonomy_with_abundance(
    taxonomy_path: Union[str, Path],
    abundance_path: Union[str, Path],
    output_path: Union[str, Path],
    sequence_col: str = "sequence",
    taxonomy_sep: str = ",",
    *,
    rank_columns: Sequence[str] = DEFAULT_RANK_COLUMNS,
    contaminants: Optional[List[str]] = None,
    pident_col: Optional[str] = None,
    unassigned_label: str = UNASSIGNED_LABEL,
) -> Path:
    """Merge a taxonomy table with a DADA2/SWARM abundance table on sequence.

    Used by ecotag (TSV taxonomy), DECIPHER (CSV taxonomy), and DADA2 RDP
    (CSV taxonomy). Produces a CSV with the same schema as BLAST output so
    downstream consumers (GBIF export, plotting, reporting) see one shape.

    Args:
        taxonomy_path: Path to taxonomy file (CSV or TSV).
        abundance_path: Path to abundance table (sequences as index, samples as columns).
        output_path: Path to output merged CSV.
        sequence_col: Name of sequence column for the join (default: 'sequence').
        taxonomy_sep: Delimiter for the taxonomy file (default: ',').
        rank_columns: Taxonomic ranks ordered coarse-to-fine. Used for cascade
            nulling and for the output column order.
        contaminants: Optional list of species names (CRABS underscore format)
            to flag in `is_contaminant_candidate`. Rows are flagged, never
            deleted -- downstream decides.
        pident_col: Name of an identity / confidence column in the taxonomy file
            to copy to the output's `pident` column. If not provided (or
            absent), the output's pident is NaN. For DECIPHER you can pass
            `confidence_species`; for ecotag and methods without a single
            confidence number, leave as None.
        unassigned_label: Label for missing taxonomy entries.

    Returns:
        Path to the output CSV.

    Raises:
        FileNotFoundError: If abundance_path doesn't exist. (Empty taxonomy
            file is allowed and produces an all-Unassigned output.)
    """
    taxonomy_path = Path(taxonomy_path)
    abundance_path = Path(abundance_path)
    output_path = Path(output_path)

    if not abundance_path.exists():
        raise FileNotFoundError(f"Abundance file not found: {abundance_path}")

    # Load abundance: sequences as index, samples as columns
    abundance_df = pd.read_csv(abundance_path, index_col=0)
    sample_cols = list(abundance_df.columns)
    abundance_df = abundance_df.reset_index().rename(columns={"index": sequence_col})

    # Generate ASV_IDs from row order (matches BLAST/SWARM convention)
    abundance_df["ASV_ID"] = [f"OTU_{i + 1}" for i in range(len(abundance_df))]

    # Load taxonomy. Empty taxonomy file -> empty df with sequence column;
    # the LEFT merge below will leave every OTU at NaN -> Unassigned.
    if not taxonomy_path.exists() or taxonomy_path.stat().st_size == 0:
        logger.warning(
            f"Taxonomy file {taxonomy_path} is empty or missing; all OTUs will "
            f"be marked '{unassigned_label}'."
        )
        taxo_df = pd.DataFrame(columns=[sequence_col, *rank_columns])
    else:
        try:
            taxo_df = pd.read_csv(taxonomy_path, sep=taxonomy_sep)
        except pd.errors.EmptyDataError:
            logger.warning(
                f"Taxonomy file {taxonomy_path} has no rows; all OTUs will be "
                f"marked '{unassigned_label}'."
            )
            taxo_df = pd.DataFrame(columns=[sequence_col, *rank_columns])

    if sequence_col not in taxo_df.columns:
        # Method's taxonomy file uses a different column name for the sequence.
        # Best effort: rename the first string-y column. Fail loudly if we can't.
        candidates = [c for c in taxo_df.columns if c.lower() in ("sequence", "seq", "asv", "asv_seq")]
        if candidates:
            taxo_df = taxo_df.rename(columns={candidates[0]: sequence_col})
        elif len(taxo_df) == 0:
            taxo_df[sequence_col] = pd.Series(dtype=str)
        else:
            raise ValueError(
                f"Taxonomy file {taxonomy_path} has no '{sequence_col}' column; "
                f"available columns: {list(taxo_df.columns)}"
            )

    # LEFT merge from abundance side -> every OTU survives (fix for I-1 / B1)
    n_with_taxo = int(taxo_df[sequence_col].isin(abundance_df[sequence_col]).sum())
    n_total = len(abundance_df)
    if n_with_taxo < n_total:
        logger.warning(
            f"{n_total - n_with_taxo} of {n_total} OTUs had no taxonomy hit and "
            f"will be marked '{unassigned_label}' in the output."
        )

    result = pd.merge(abundance_df, taxo_df, on=sequence_col, how="left")

    # Normalize and fill ranks
    for rank in rank_columns:
        if rank not in result.columns:
            result[rank] = unassigned_label
        else:
            result[rank] = _normalize_unassigned(result[rank], unassigned_label)

    # Cascade null: coarse Unassigned -> all finer Unassigned (fix for I-3 / B3)
    for i, rank in enumerate(rank_columns):
        is_unassigned = result[rank] == unassigned_label
        for finer in list(rank_columns)[i + 1:]:
            result.loc[is_unassigned, finer] = unassigned_label

    # Optional pident column (copied from the method's confidence column, if any)
    if pident_col and pident_col in result.columns:
        result["pident"] = pd.to_numeric(result[pident_col], errors="coerce")
    else:
        result["pident"] = pd.NA

    # Contaminant flag (fix for I-7 / B5)
    if contaminants:
        contam_set = set(contaminants)
        is_contam = result["species"].astype(str).isin(contam_set)
        result[CONTAMINANT_FLAG_COL] = is_contam
        n_flagged = int(is_contam.sum())
        if n_flagged > 0:
            breakdown = result.loc[is_contam, "species"].value_counts().to_dict()
            logger.warning(
                f"Flagged {n_flagged} OTUs as candidate contaminants "
                f"(by species match): {breakdown}"
            )
    else:
        result[CONTAMINANT_FLAG_COL] = False

    # Stable column order: ASV_ID, pident, ranks..., contam flag, samples..., Sequence
    if "Sequence" not in result.columns:
        result["Sequence"] = result[sequence_col]
    ordered = (
        ["ASV_ID", "pident"]
        + list(rank_columns)
        + [CONTAMINANT_FLAG_COL]
        + sample_cols
        + ["Sequence"]
    )
    result = result[[c for c in ordered if c in result.columns]]

    # Sort by ASV number for deterministic output
    result["asv_num"] = (
        result["ASV_ID"].astype(str).str.extract(r"(\d+)").astype("Int64")
    )
    result = result.sort_values("asv_num").drop(columns="asv_num").reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)

    logger.info(
        f"Linked taxonomy with abundances: {output_path} "
        f"({len(result)} OTUs, {n_total - n_with_taxo} unassigned)"
    )
    return output_path
