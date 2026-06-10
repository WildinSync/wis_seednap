"""Filter non-target taxa from eDNA metabarcoding results.

Sits in the formatting stage of the pipeline, between taxonomic assignment and
the GBIF/DarwinCore export. A metabarcoding primer (the short PCR primer pair
that defines a marker, e.g. ``teleo`` for fish) is never perfectly specific: it
co-amplifies and the reference DB then assigns reads to organisms outside the
intended target group. Common contaminants are human and other primate DNA from
handling, and incidental hits to unrelated classes. This module drops those
known off-target assignments per marker so they do not reach the published
occurrence records.
"""

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
    """Remove non-target taxa from taxonomy results based on marker type.

    Looks up the per-marker rules in ``NON_TARGET_TAXA`` and removes any
    occurrence row whose assigned taxonomy matches a banned taxon. A marker with
    no rules (not present in the table) is passed through unchanged.
    """

    def filter(self, df: pd.DataFrame, marker: str) -> pd.DataFrame:
        """
        Remove rows matching known non-target taxa for the given marker.

        Drops occurrence rows whose ``class``/``order``/``family``/``genus``
        assignment is on the marker's non-target list (e.g. human and other
        primate hits for the fish marker ``teleo``). A rank is only checked when
        its column is present in ``df``; unknown markers are returned unchanged
        with a debug log.

        Args:
            df: Taxonomy results, one row per occurrence, carrying taxonomic
                columns (any of ``class``, ``order``, ``family``, ``genus``).
            marker: Marker name used to select the rule set (e.g. ``'teleo'``).

        Returns:
            A copy of ``df`` with non-target rows removed and the index reset.
            When the marker has no rules, the input ``df`` is returned as-is.
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
