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
from typing import List, Optional, Sequence, Tuple, Union

import pandas as pd

logger = logging.getLogger(__name__)


# Match BLAST's TAXONOMIC_RANKS to keep schemas aligned across methods.
DEFAULT_RANK_COLUMNS: Tuple[str, ...] = (
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
    cascade_rank_columns: Optional[Sequence[str]] = None,
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
        rank_columns: Taxonomic ranks ordered coarse-to-fine. Used for the
            output column order (and, by default, for cascade nulling).
        cascade_rank_columns: Subset of `rank_columns` (coarse-to-fine) over
            which cascade nulling applies. Defaults to `rank_columns`. The
            ecotag path passes the obitab-resolvable ranks (order..species)
            here while keeping kingdom/phylum/class in `rank_columns` for
            schema parity: obitab never resolves kingdom/phylum/class, so they
            arrive as placeholders and must NOT trigger the cascade (which,
            keyed on the coarsest rank, would otherwise force every finer rank
            to Unassigned and silently zero out all ecotag taxonomy). Those
            coarse ranks are enriched downstream (DarwinCore NCBI/WORMS lookup).
        contaminants: Optional list of species names (CRABS underscore format)
            to flag in `is_contaminant_candidate`. Rows are flagged, never
            deleted -- downstream decides.
        pident_col: Name of an identity / confidence column in the taxonomy file
            to copy to the output's `pident` column (for schema parity with the
            BLAST path). If not provided (or absent), the output's pident is NaN.
            The DADA2 RDP path wires this as `bootstrap_min` (see
            _assign_dada2 in assigner.py), exposing the RDP bootstrap as pident;
            the DECIPHER and ecotag paths do not pass it, so their pident is NaN.
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
        raise FileNotFoundError(
            f"Abundance table not found: {abundance_path}. Taxonomy linking needs "
            f"the feature-count table written by the clustering step (DADA2 writes "
            f"02_dada2/<marker>/seqtab_clean_t.csv; SWARM writes "
            f"02_swarm/<marker>/otu_table.csv). The clustering step is marked "
            f"complete in the run state but its output file is missing from disk "
            f"(deleted or moved). Re-run the dada2 or swarm step to regenerate it "
            f"before re-running taxonomy; a plain --resume will not help because it "
            f"skips steps already marked complete in the state JSON without checking "
            f"that their outputs still exist."
        )

    # Load abundance: sequences as index, samples as columns
    try:
        abundance_df = pd.read_csv(abundance_path, index_col=0)
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
        raise ValueError(
            f"Could not read the abundance/count table: {abundance_path} ({exc}). "
            f"The file exists but is empty or not valid CSV. If you passed it to "
            f"`seednap assign-taxonomy` by hand, check you pointed at the right "
            f"seqtab_clean_t.csv / otu_table.csv and that it is not a stale or "
            f"truncated file; otherwise the dada2/swarm step that wrote it may have "
            f"been interrupted -- re-run that step (e.g. `--resume`) or delete the "
            f"file and regenerate it. Expected layout: sequences in the first "
            f"(index) column, one column per sample."
        ) from exc
    sample_cols = list(abundance_df.columns)
    abundance_df = abundance_df.reset_index().rename(columns={"index": sequence_col})

    # Generate feature IDs from row order. The values are OTU_-prefixed, but the
    # column is deliberately named ASV_ID: BLAST output uses an ASV_ID column for
    # both ASV_- and OTU_-prefixed IDs, and downstream consumers (GBIF export,
    # plotting, reporting) key on that one column name. Do not rename the column
    # to "OTU_ID" to match the prefix; that would break the shared cross-method
    # schema even though it looks more consistent.
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
    # Only ranks that the method actually resolves should drive the cascade;
    # placeholder ranks (e.g. ecotag's absent kingdom/phylum/class) are kept in
    # the output for schema parity but excluded here via cascade_rank_columns.
    cascade_cols = list(cascade_rank_columns) if cascade_rank_columns is not None else list(rank_columns)
    for i, rank in enumerate(cascade_cols):
        if rank not in result.columns:
            continue
        is_unassigned = result[rank] == unassigned_label
        for finer in cascade_cols[i + 1:]:
            if finer in result.columns:
                result.loc[is_unassigned, finer] = unassigned_label

    # Optional pident column (copied from the method's confidence column, if any)
    if pident_col and pident_col in result.columns:
        result["pident"] = pd.to_numeric(result[pident_col], errors="coerce")
    else:
        result["pident"] = pd.NA

    # Contaminant flag (fix for I-7 / B5)
    if contaminants:
        contam_set = set(contaminants)

        # DADA2 addSpecies runs with allowMultiple=TRUE, so a multi-hit cell
        # holds a '/'-joined value (e.g. 'Salmo_trutta/Salmo_salar'). Match if
        # ANY component is a configured contaminant; an exact equality test
        # would miss a contaminant hidden inside an ambiguous multi-match.
        def _row_is_contam(value: object) -> bool:
            return any(part in contam_set for part in str(value).split("/"))

        species_str = result["species"].astype(str)
        is_contam = species_str.map(_row_is_contam)
        result[CONTAMINANT_FLAG_COL] = is_contam
        n_flagged = int(is_contam.sum())
        if n_flagged > 0:
            # Break down by the matched components, not the raw (possibly
            # multi-hit) cell, so the [WARN] names the actual contaminants.
            matched_components: List[str] = []
            for value in species_str[is_contam]:
                matched_components.extend(
                    part for part in str(value).split("/") if part in contam_set
                )
            breakdown = pd.Series(matched_components).value_counts().to_dict()
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
