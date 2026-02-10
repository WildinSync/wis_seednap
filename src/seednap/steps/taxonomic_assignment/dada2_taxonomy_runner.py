"""DADA2 R package wrapper for taxonomic assignment.

This module provides a Python wrapper around DADA2's naive Bayesian classifier
(assignTaxonomy + addSpecies) for taxonomic assignment of eDNA sequences.
"""

import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)


class Dada2TaxonomyError(Exception):
    """Exception raised for DADA2 taxonomy assignment errors."""

    pass


class Dada2TaxonomyRunner:
    """
    Run DADA2 taxonomic assignment via Rscript.

    DADA2 uses a naive Bayesian classifier with RDP training set
    for genus-level assignment and exact matching for species-level.
    """

    def __init__(self, timeout: int = 7200):
        """
        Initialize DADA2 taxonomy runner.

        Args:
            timeout: Command timeout in seconds (default: 7200 = 2 hours)
        """
        self.timeout = timeout
        self._check_r_availability()

    def _check_r_availability(self) -> None:
        """
        Check if R and DADA2 package are available.

        Raises:
            Dada2TaxonomyError: If R or DADA2 is not found
        """
        try:
            result = subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
            logger.debug(f"Found Rscript: {result.stderr.strip()}")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            raise Dada2TaxonomyError("Rscript not found. Is R installed?") from e

    def _run_r_script(
        self,
        script_path: Union[str, Path],
        args: List[str],
        log_file: Optional[Union[str, Path]] = None,
    ) -> str:
        """
        Execute R script via Rscript.

        Args:
            script_path: Path to R script file
            args: List of arguments to pass to script
            log_file: Optional path to log file for stdout/stderr

        Returns:
            stdout from Rscript

        Raises:
            Dada2TaxonomyError: If Rscript command fails
        """
        script_path = Path(script_path)
        if not script_path.exists():
            raise FileNotFoundError(f"R script not found: {script_path}")

        cmd = ["Rscript", str(script_path)] + [str(arg) for arg in args]
        logger.info(f"Running R script: {script_path.name} {' '.join([str(a) for a in args])}")

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True, timeout=self.timeout
            )

            # Write to log file if specified
            if log_file:
                log_path = Path(log_file)
                log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(log_path, "w") as f:
                    f.write(f"Command: {' '.join(cmd)}\n")
                    f.write(f"\n{'='*80}\n")
                    f.write("STDOUT:\n")
                    f.write(result.stdout)
                    f.write(f"\n{'='*80}\n")
                    f.write("STDERR:\n")
                    f.write(result.stderr)

            logger.debug(f"R script completed successfully: {script_path.name}")
            return result.stdout

        except subprocess.CalledProcessError as e:
            error_msg = f"R script failed: {script_path.name}\n{e.stderr}"
            logger.error(error_msg)
            raise Dada2TaxonomyError(error_msg) from e

        except subprocess.TimeoutExpired as e:
            error_msg = f"R script timed out after {self.timeout} seconds: {script_path.name}"
            logger.error(error_msg)
            raise Dada2TaxonomyError(error_msg) from e

    def run_dada2_taxonomy(
        self,
        marker: str,
        output_dir: Union[str, Path],
        rdp_db_path: Union[str, Path],
        species_db_path: Union[str, Path],
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
            script_path: Path to R script (default: scripts/taxo_dada2_marker.R)
            log_file: Path to log file

        Returns:
            Dictionary with paths to output files:
            - taxonomy: Taxonomy table CSV
            - complete: Complete table with taxonomy and abundances

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
            args=[marker, str(rdp_db_path), str(species_db_path)],
            log_file=log_file,
        )

        # Construct output paths
        marker_dir = output_dir / "02_dada2" / marker

        return {
            "taxonomy": marker_dir / "taxonomy_dada2RDP.csv",
            "complete": output_dir / f"{marker}_dada2RDP.csv",
        }
