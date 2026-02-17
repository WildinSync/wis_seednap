"""Enrich taxonomy DataFrames with missing kingdom/phylum via NCBI and WORMS."""

import json
import logging
import os
import time
import urllib.error
import urllib.request
from typing import Dict, Optional, Tuple

import pandas as pd
from dotenv import load_dotenv

logger = logging.getLogger(__name__)

# Load .env from the project root (or wherever the user placed it).
load_dotenv()

# Rate-limit: NCBI allows 10 req/s with API key, 3 req/s without.
_NCBI_DELAY_WITH_KEY = 0.1  # seconds between requests
_NCBI_DELAY_NO_KEY = 0.34


class TaxonomyEnricher:
    """
    Fill missing kingdom and phylum values by querying NCBI Entrez and WORMS.

    Reads ``NCBI_API_KEY`` from the environment (or a ``.env`` file via
    python-dotenv).  If the key is absent the enrichment step is skipped.
    """

    def __init__(self, email: str = "seednap@ethz.ch") -> None:
        self._cache: Dict[str, Tuple[Optional[str], Optional[str], Optional[str]]] = {}
        self._api_key: Optional[str] = os.environ.get("NCBI_API_KEY")
        self._email = email
        self._delay = _NCBI_DELAY_WITH_KEY if self._api_key else _NCBI_DELAY_NO_KEY

    # ------------------------------------------------------------------
    # Public
    # ------------------------------------------------------------------

    def enrich(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Fill missing *kingdom*, *phylum*, and *class* columns in *df*.

        Strategy (mirrors the original R ``add_kingdom_phylum``):
        1. Collect unique non-NA ``class`` values → query for phylum/kingdom.
        2. Collect unique ``scientificName`` values where class is NA →
           query for class/phylum/kingdom.
        3. Left-join results back into the DataFrame.

        Args:
            df: Must contain a ``scientificName`` column.  ``class``,
                ``phylum``, and ``kingdom`` columns are created if absent.

        Returns:
            DataFrame with filled taxonomy columns.
        """
        if self._api_key is None:
            logger.warning(
                "NCBI_API_KEY not set — skipping taxonomy enrichment. "
                "Set the variable in your .env file or environment to enable "
                "kingdom/phylum lookup."
            )
            return df

        # Lazy import so Biopython is not required at module load time.
        from Bio import Entrez

        Entrez.email = self._email
        Entrez.api_key = self._api_key

        df = df.copy()
        for col in ("class", "phylum", "kingdom"):
            if col not in df.columns:
                df[col] = pd.NA

        # 1. Queries derived from existing class values
        class_queries = df["class"].dropna().unique().tolist()
        class_queries = [q for q in class_queries if q != ""]

        # 2. Queries from scientificName where class is missing
        sci_queries = (
            df.loc[df["class"].isna() | (df["class"] == ""), "scientificName"]
            .dropna()
            .unique()
            .tolist()
        )

        all_queries = list(set(class_queries + sci_queries))
        if not all_queries:
            return df

        logger.info(f"Enriching taxonomy for {len(all_queries)} unique name(s)…")
        for query in all_queries:
            self._fetch(query)

        # Build lookup tables
        lookup_by_class: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
        for q in class_queries:
            if q in self._cache:
                _, phylum, kingdom = self._cache[q]
                lookup_by_class[q] = (phylum, kingdom)

        lookup_by_sci: Dict[str, Tuple[Optional[str], Optional[str], Optional[str]]] = {}
        for q in sci_queries:
            if q in self._cache:
                lookup_by_sci[q] = self._cache[q]

        # 3. Apply lookups
        # By class → fill phylum/kingdom
        def _fill_from_class(row: pd.Series) -> pd.Series:
            cls = row.get("class")
            if pd.notna(cls) and cls != "" and cls in lookup_by_class:
                phylum, kingdom = lookup_by_class[cls]
                if pd.isna(row.get("phylum")) or row["phylum"] == "":
                    row["phylum"] = phylum
                if pd.isna(row.get("kingdom")) or row["kingdom"] == "":
                    row["kingdom"] = kingdom
            return row

        df = df.apply(_fill_from_class, axis=1)

        # By scientificName → fill class/phylum/kingdom where class was missing
        def _fill_from_sci(row: pd.Series) -> pd.Series:
            if pd.notna(row.get("class")) and row["class"] != "":
                return row
            sci = row.get("scientificName")
            if pd.notna(sci) and sci in lookup_by_sci:
                cls, phylum, kingdom = lookup_by_sci[sci]
                if cls and (pd.isna(row.get("class")) or row["class"] == ""):
                    row["class"] = cls
                if phylum and (pd.isna(row.get("phylum")) or row["phylum"] == ""):
                    row["phylum"] = phylum
                if kingdom and (pd.isna(row.get("kingdom")) or row["kingdom"] == ""):
                    row["kingdom"] = kingdom
            return row

        df = df.apply(_fill_from_sci, axis=1)

        return df

    # ------------------------------------------------------------------
    # Private: fetch taxonomy for a single name
    # ------------------------------------------------------------------

    def _fetch(self, name: str) -> None:
        """Query NCBI then WORMS for *name*; cache the result."""
        if name in self._cache:
            return

        result = self._query_ncbi(name)
        if result is None:
            result = self._query_worms(name)
        if result is None:
            result = (None, None, None)

        self._cache[name] = result

    def _query_ncbi(
        self, name: str
    ) -> Optional[Tuple[Optional[str], Optional[str], Optional[str]]]:
        """Return (class, phylum, kingdom) from NCBI Taxonomy, or None."""
        from Bio import Entrez

        try:
            time.sleep(self._delay)
            handle = Entrez.esearch(db="taxonomy", term=name, retmax=1)
            search = Entrez.read(handle)
            handle.close()

            if not search.get("IdList"):
                return None

            tax_id = search["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            if not records:
                return None

            lineage = records[0].get("LineageEx", [])
            cls = phylum = kingdom = None
            for entry in lineage:
                rank = entry.get("Rank", "")
                sci_name = entry.get("ScientificName", "")
                if rank == "class":
                    cls = sci_name
                elif rank == "phylum":
                    phylum = sci_name
                elif rank == "kingdom":
                    kingdom = sci_name

            return (cls, phylum, kingdom)

        except Exception as exc:
            logger.debug(f"NCBI query failed for '{name}': {exc}")
            return None

    def _query_worms(
        self, name: str
    ) -> Optional[Tuple[Optional[str], Optional[str], Optional[str]]]:
        """Return (class, phylum, kingdom) from WORMS REST API, or None."""
        url = (
            f"https://www.marinespecies.org/rest/AphiaRecordsByName/"
            f"{urllib.request.quote(name)}?like=false&marine_only=false"
        )
        try:
            time.sleep(0.5)  # be polite to WORMS
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())

            if not data or not isinstance(data, list):
                return None

            record = data[0]
            return (
                record.get("class"),
                record.get("phylum"),
                record.get("kingdom"),
            )

        except Exception as exc:
            logger.debug(f"WORMS query failed for '{name}': {exc}")
            return None
