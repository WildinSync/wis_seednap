"""Shared taxonomy utilities for linking taxonomy results with abundance tables."""

import logging
from pathlib import Path
from typing import Union

import pandas as pd

logger = logging.getLogger(__name__)


def link_taxonomy_with_abundance(
    taxonomy_path: Union[str, Path],
    abundance_path: Union[str, Path],
    output_path: Union[str, Path],
    sequence_col: str = "sequence",
    taxonomy_sep: str = ",",
) -> Path:
    """
    Merge a taxonomy table with a DADA2 abundance table on a shared sequence column.

    Used by both ecotag (TSV taxonomy) and DECIPHER (CSV taxonomy) runners.

    Args:
        taxonomy_path: Path to taxonomy file (CSV or TSV).
        abundance_path: Path to DADA2 abundance table (seqtab_clean_t.csv).
        output_path: Path to output merged CSV.
        sequence_col: Name of sequence column for the join (default: 'sequence').
        taxonomy_sep: Delimiter for the taxonomy file (default: ',').

    Returns:
        Path to the output CSV.

    Raises:
        FileNotFoundError: If input files don't exist.
    """
    taxonomy_path = Path(taxonomy_path)
    abundance_path = Path(abundance_path)
    output_path = Path(output_path)

    if not taxonomy_path.exists():
        raise FileNotFoundError(f"Taxonomy file not found: {taxonomy_path}")
    if not abundance_path.exists():
        raise FileNotFoundError(f"Abundance file not found: {abundance_path}")

    taxo_df = pd.read_csv(taxonomy_path, sep=taxonomy_sep)

    abundance_df = pd.read_csv(abundance_path, index_col=0)
    abundance_df = abundance_df.reset_index().rename(columns={"index": sequence_col})

    result = pd.merge(taxo_df, abundance_df, on=sequence_col, how="left")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)

    logger.info(f"Linked taxonomy with abundances: {output_path}")
    return output_path
