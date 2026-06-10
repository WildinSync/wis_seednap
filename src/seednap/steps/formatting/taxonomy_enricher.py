"""Enrich taxonomy DataFrames with missing kingdom/phylum via NCBI and WORMS.

Runs inside the formatting stage, called by ``DarwinCoreBuilder`` while building
the GBIF occurrence table. Taxonomic assignment (DADA2 RDP, BLAST, DECIPHER,
ecotag) usually resolves a read to a low rank (species/genus) but leaves the
higher ranks (kingdom, phylum, sometimes class) blank, because the marker's
reference database does not store the full lineage. GBIF expects those higher
ranks populated. This module looks each name up in two public taxonomic
authorities -- NCBI Taxonomy (via the Entrez API) and, as a fallback, the World
Register of Marine Species (WORMS) -- and fills only the empty cells, so the
denser assignment already produced upstream is preserved.
"""

import json
import logging
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Dict, Optional, Tuple

import pandas as pd
from dotenv import load_dotenv

logger = logging.getLogger(__name__)

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
        """Configure the enricher and read the NCBI API key from the environment.

        Loads a ``.env`` file (via python-dotenv) so the ``NCBI_API_KEY`` can be
        supplied either as a real environment variable or from that file. The
        presence of a key selects the faster NCBI rate limit (10 req/s with key,
        3 req/s without); its absence disables enrichment entirely in ``enrich``.

        Args:
            email: Contact address sent to NCBI Entrez as required by their
                usage policy. Defaults to the lab address.
        """
        load_dotenv()
        self._cache: Dict[str, Tuple[Optional[str], Optional[str], Optional[str]]] = {}
        self._last_worms_error: Optional[str] = None
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
        3. Fill blanks per row from the query results (two row-wise passes,
           not a merge). Each assignment only writes when the target cell is
           empty (NA or ""); existing values are never overwritten, so
           DADA2/BLAST-supplied taxonomy is preserved.

        Args:
            df: Must contain a ``scientificName`` column.  ``class``,
                ``phylum``, and ``kingdom`` columns are created if absent.

        Returns:
            A copy of ``df`` with ``class``, ``phylum``, and ``kingdom`` columns
            present and any previously empty cells filled where a lookup
            succeeded. When ``NCBI_API_KEY`` is unset the input ``df`` is
            returned unchanged (a ``[WARN]`` is logged); names that neither NCBI
            nor WORMS resolves keep their empty higher ranks (also warned).

        Raises:
            ImportError: If Biopython (the ``Bio.Entrez`` module) is not
                installed; it is imported lazily only when a key is present.
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

        # Biopython ships no type stubs; mypy infers these module attributes as
        # NoneType from their defaults. Assigning the configured email/key is the
        # documented Entrez usage, so the assignment is correct.
        Entrez.email = self._email  # type: ignore[assignment]
        Entrez.api_key = self._api_key  # type: ignore[assignment]

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

        # kingdom/phylum are not in darwincore_builder._DWC_REQUIRED_FIELDS, so
        # the pre-write guard cannot catch a wholly-unenriched name. Surface the
        # count here so missing taxonomy is visible before GBIF submission.
        n_unenriched = sum(
            1 for q in all_queries if self._cache.get(q) == (None, None, None)
        )
        if n_unenriched:
            logger.warning(
                f"[WARN] taxonomy enrichment: expected=kingdom/phylum for "
                f"{len(all_queries)} name(s), got={n_unenriched} unenriched, "
                f"fallback=those rows ship to GBIF with empty kingdom/phylum"
            )

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
            """Fill empty phylum/kingdom from the row's class lookup.

            Uses the row's existing ``class`` value as the lookup key. Only
            writes ``phylum``/``kingdom`` cells that are currently empty (NA or
            ""), so an upstream assignment is never overwritten.

            Args:
                row: One occurrence row; read keys ``class``, ``phylum``,
                    ``kingdom``.

            Returns:
                The same row with empty phylum/kingdom filled where the class
                lookup provided a value.
            """
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
            """Fill empty class/phylum/kingdom from the scientificName lookup.

            Applies only to rows whose ``class`` is still missing; uses the
            row's ``scientificName`` as the lookup key. Each of class, phylum,
            and kingdom is written only when the lookup returned a value and the
            target cell is empty, so existing assignments are preserved.

            Args:
                row: One occurrence row; read keys ``class``, ``scientificName``,
                    ``phylum``, ``kingdom``.

            Returns:
                The same row, unchanged if ``class`` was already set, otherwise
                with empty class/phylum/kingdom filled from the name lookup.
            """
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
        """Query NCBI then WORMS for *name*; cache the result.

        Tries NCBI Taxonomy first and, only if that returns nothing, falls back
        to WORMS. Whatever is found (or ``(None, None, None)`` if both fail) is
        stored in the in-memory cache keyed by ``name`` so the same name is never
        queried twice in one run. A total miss is reported with a ``[WARN]``.

        Args:
            name: A taxon name to resolve (either a ``class`` value or a
                ``scientificName``).

        Returns:
            None. The (class, phylum, kingdom) tuple is written to ``self._cache``
            as a side effect.
        """
        if name in self._cache:
            return

        result = self._query_ncbi(name)
        if result is None:
            result = self._query_worms(name)
        if result is None:
            worms_err = self._last_worms_error or "no record found"
            logger.warning(
                f"[WARN] taxonomy enrichment: expected=kingdom/phylum for "
                f"'{name}' from NCBI or WORMS, got=both lookups returned no "
                f"result (last WORMS error: {worms_err}), fallback=DarwinCore "
                f"row exported with blank kingdom/phylum"
            )
            result = (None, None, None)

        self._cache[name] = result

    def _query_ncbi(
        self, name: str
    ) -> Optional[Tuple[Optional[str], Optional[str], Optional[str]]]:
        """Return (class, phylum, kingdom) from NCBI Taxonomy, or None.

        Runs an Entrez ``esearch`` to map the name to a taxonomy ID, then an
        ``efetch`` to pull the lineage and reads the class/phylum/kingdom ranks
        out of it. Sleeps ``self._delay`` seconds first to respect the NCBI rate
        limit. Any network or parsing error is caught and logged as a ``[WARN]``,
        and the method returns None so the caller can fall back to WORMS.

        Args:
            name: Taxon name to search in the NCBI Taxonomy database.

        Returns:
            A ``(class, phylum, kingdom)`` tuple, with each element None if that
            rank is absent from the lineage. Returns None when the name has no
            NCBI match, the fetch is empty, or any error occurs.
        """
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
            logger.warning(
                f"[WARN] taxonomy enrichment: expected=NCBI Taxonomy lineage "
                f"for '{name}', got=error ({exc}), fallback=trying WORMS next"
            )
            return None

    def _query_worms(
        self, name: str
    ) -> Optional[Tuple[Optional[str], Optional[str], Optional[str]]]:
        """Return (class, phylum, kingdom) from WORMS REST API, or None.

        Calls the WORMS ``AphiaRecordsByName`` endpoint (exact, non-marine-only
        match) and reads the ranks from the first returned record. Sleeps 0.5 s
        first to be polite to the service. On any error the message is stored in
        ``self._last_worms_error`` (for the caller's warning), logged at debug
        level, and None is returned.

        Args:
            name: Taxon name to look up in the World Register of Marine Species.

        Returns:
            A ``(class, phylum, kingdom)`` tuple, with each element None if WORMS
            omits that rank. Returns None when there is no match or any error
            occurs.
        """
        url = (
            f"https://www.marinespecies.org/rest/AphiaRecordsByName/"
            f"{urllib.parse.quote(name)}?like=false&marine_only=false"
        )
        self._last_worms_error = None
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
            self._last_worms_error = str(exc)
            logger.debug(f"WORMS query failed for '{name}': {exc}")
            return None
