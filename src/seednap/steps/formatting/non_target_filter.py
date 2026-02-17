"""Filter non-target taxa from eDNA metabarcoding results."""

import logging
from typing import Dict, List

import pandas as pd

logger = logging.getLogger(__name__)

# Per-marker non-target taxa definitions.
# Each key maps a taxonomic rank to a list of taxa to remove.
NON_TARGET_TAXA: Dict[str, Dict[str, List[str]]] = {
    "teleo": {
        "class": ["Gastropoda", "Holothuroidea", "Deinococci", "Branchiopoda"],
        "order": ["Primates", "Galliformes"],
        "family": ["Hominidae"],
        "genus": ["Gorilla", "Homo", "Pan", "Macaca"],
    },
}


class NonTargetFilter:
    """Remove non-target taxa from taxonomy results based on marker type."""

    def filter(self, df: pd.DataFrame, marker: str) -> pd.DataFrame:
        """
        Remove rows matching known non-target taxa for the given marker.

        Args:
            df: DataFrame with taxonomic columns (class, order, family, genus).
            marker: Marker name (e.g. 'teleo').

        Returns:
            Filtered DataFrame with non-target rows removed.
        """
        rules = NON_TARGET_TAXA.get(marker)
        if rules is None:
            logger.debug(f"No non-target filter defined for marker '{marker}'")
            return df

        mask = pd.Series(True, index=df.index)
        for rank, taxa in rules.items():
            if rank in df.columns:
                rank_mask = ~df[rank].isin(taxa)
                mask = mask & rank_mask

        n_removed = (~mask).sum()
        if n_removed > 0:
            logger.info(
                f"Removed {n_removed} non-target occurrence(s) for marker '{marker}'"
            )

        return df.loc[mask].reset_index(drop=True)
