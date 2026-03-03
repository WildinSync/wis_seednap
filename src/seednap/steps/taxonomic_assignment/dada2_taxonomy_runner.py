"""DADA2 R package wrapper for taxonomic assignment.

This module provides a Python wrapper around DADA2's naive Bayesian classifier
(assignTaxonomy + addSpecies) for taxonomic assignment of eDNA sequences.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Union

from seednap.utils.r_runner import RScriptRunner

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

    def __init__(self, timeout: int = 7200):
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
        multithread: bool = True,
        script_path: Optional[Union[str, Path]] = None,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Dict[str, Path]:
        """
        Run DADA2 taxonomic assignment.

        Uses DADA2's naive Bayesian classifier with RDP training set.

        Args:
            marker: Marker name
            output_dir: Base output directory
            rdp_db_path: Path to RDP-formatted database (genus-level)
            species_db_path: Path to species-level database
            multithread: Use multithreading (default: True)
            script_path: Path to R script (default: scripts/taxo_dada2_marker.R)
            log_file: Path to log file

        Returns:
            Dictionary with paths to output files:
            - taxonomy: Taxonomy table CSV
            - final_table: Complete table with taxonomy and abundances

        Raises:
            Dada2TaxonomyError: If taxonomy assignment fails
        """
        if script_path is None:
            script_path = Path("scripts/taxo_dada2_marker.R")

        rdp_db_path = Path(rdp_db_path)
        species_db_path = Path(species_db_path)

        if not rdp_db_path.exists():
            raise FileNotFoundError(f"RDP database not found: {rdp_db_path}")
        if not species_db_path.exists():
            raise FileNotFoundError(f"Species database not found: {species_db_path}")

        # Check that sequence table exists
        output_dir = Path(output_dir)
        seqtab_rds = output_dir / "02_dada2" / marker / "seqtab_clean.rds"
        if not seqtab_rds.exists():
            raise FileNotFoundError(
                f"Sequence table not found: {seqtab_rds}. "
                "Run DADA2 processing first."
            )

        logger.info(f"Running DADA2 taxonomic assignment for {marker}")

        # Run taxonomy assignment
        self._run_r_script(
            script_path=script_path,
            args=[
                marker,
                str(rdp_db_path),
                str(species_db_path),
                str(output_dir),
                str(multithread).upper(),
            ],
            log_file=log_file,
        )

        # Construct output paths
        marker_dir = output_dir / "02_dada2" / marker

        return {
            "taxonomy": marker_dir / "taxonomy_dada2RDP.csv",
            "final_table": output_dir / f"{marker}_dada2RDP.csv",
        }
