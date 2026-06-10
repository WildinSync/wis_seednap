"""DADA2 R package wrapper for taxonomic assignment.

This module provides a Python wrapper around DADA2's naive Bayesian classifier
(assignTaxonomy + addSpecies) for taxonomic assignment of eDNA sequences.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Union

from seednap.utils.r_runner import RScriptRunner, r_script_path

logger = logging.getLogger(__name__)


class Dada2TaxonomyError(Exception):
    """Exception raised for DADA2 taxonomy assignment errors."""

    pass


class Dada2TaxonomyRunner(RScriptRunner):
    """
    Run DADA2 taxonomic assignment via Rscript.

    DADA2 uses a naive Bayesian classifier with RDP training set
    for genus-level assignment and exact matching for species-level.
    """

    _error_class = Dada2TaxonomyError

    def __init__(self, timeout: int = 7200) -> None:
        """
        Initialize DADA2 taxonomy runner.

        Args:
            timeout: Command timeout in seconds (default: 7200 = 2 hours)
        """
        super().__init__(timeout=timeout)

    def run_dada2_taxonomy(
        self,
        marker: str,
        output_dir: Union[str, Path],
        rdp_db_path: Union[str, Path],
        species_db_path: Union[str, Path],
        query_fasta: Union[str, Path],
        multithread: bool = True,
        bootstrap_threshold: int = 80,
        script_path: Optional[Union[str, Path]] = None,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Dict[str, Path]:
        """
        Run DADA2 taxonomic assignment.

        Reads sequences from a query FASTA (produced by either DADA2 ASV or
        SWARM OTU clustering), applies the RDP naive Bayesian classifier
        with a Wang 2007 bootstrap threshold (default 80%), cascade-nulls
        finer ranks below the threshold, and writes a per-sequence taxonomy
        CSV. The merge with the abundance table is done by the Python
        caller via `seednap.utils.taxonomy.link_taxonomy_with_abundance`.

        The script is cluster-method agnostic by design -- it does not read
        seqtab_clean.rds and does not assume a particular output directory
        layout. Pass the query FASTA explicitly.

        Args:
            marker: Marker name (used for log messages and output filename)
            output_dir: Base output directory; the per-sequence taxonomy CSV
                is written under output_dir/02_dada2/{marker}/ for backward
                compatibility, and the merged final_table goes to
                output_dir/{marker}_dada2RDP.csv.
            rdp_db_path: Path to RDP-formatted database (kingdom..genus)
            species_db_path: Path to species-level database (exact match)
            query_fasta: Path to query.fasta (sequences to assign taxonomy to)
            multithread: Use multithreading (default: True)
            bootstrap_threshold: Min bootstrap (%) for a rank to be retained
                (default 80, per Wang 2007 RDP standard for short rRNA reads)
            script_path: Path to R script (default: seednap/scripts/taxo_dada2_marker.R)
            log_file: Path to log file

        Returns:
            Dictionary with paths to output files:
            - taxonomy: Taxonomy table CSV (per-sequence, no abundances)
            - final_table: Where the merged output WILL go (Python writes it)

        Raises:
            Dada2TaxonomyError: If taxonomy assignment fails
            FileNotFoundError: If any input file is missing
        """
        if script_path is None:
            script_path = r_script_path("taxo_dada2_marker.R")

        rdp_db_path = Path(rdp_db_path)
        species_db_path = Path(species_db_path)
        query_fasta = Path(query_fasta)
        output_dir = Path(output_dir)

        if not rdp_db_path.exists():
            raise FileNotFoundError(
                f"RDP training database not found: {rdp_db_path}. "
                f"DADA2's naive-Bayes classifier needs an RDP-formatted training "
                f"FASTA (kingdom..genus ranks) for this marker; this is the 'all' "
                f"database, separate from the species-level DB. For run-pipeline, "
                f"set taxonomy.databases.dada2.all in the marker YAML; for the "
                f"standalone `dada2 --assign-taxonomy` command, pass --rdp-db. "
                f"Confirm the file actually exists on this host -- DB paths differ "
                f"between the eDNA server and local checkouts, and `seednap validate` "
                f"flags a MISSING path but does not block the run."
            )
        if not species_db_path.exists():
            raise FileNotFoundError(
                f"Species-level database not found: {species_db_path}. "
                f"DADA2 addSpecies needs the exact-match species-assignment training "
                f"FASTA, configured under taxonomy.databases.dada2.species in the "
                f"marker YAML, but that path does not exist on this host (config "
                f"validation checks the path string, not its presence on disk). Point "
                f"taxonomy.databases.dada2.species at the species-assignment FASTA for "
                f"this marker and confirm the file exists on this server. The species "
                f"DB is required for the dada2 method (addSpecies is always run); it "
                f"cannot be skipped by omitting the key."
            )
        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")

        # The R script writes per-sequence taxonomy here; Python does the
        # merge into final_table.
        marker_dir = output_dir / "02_dada2" / marker
        marker_dir.mkdir(parents=True, exist_ok=True)
        taxonomy_csv = marker_dir / "taxonomy_dada2RDP.csv"

        logger.info(
            f"Running DADA2 taxonomic assignment for {marker} "
            f"(bootstrap_threshold={bootstrap_threshold}%, query={query_fasta})"
        )

        self._run_r_script(
            script_path=script_path,
            args=[
                marker,
                str(rdp_db_path),
                str(species_db_path),
                str(query_fasta),
                str(taxonomy_csv),
                str(multithread).upper(),
                str(bootstrap_threshold),
            ],
            log_file=log_file,
        )

        return {
            "taxonomy": taxonomy_csv,
            "final_table": output_dir / f"{marker}_dada2RDP.csv",
        }
