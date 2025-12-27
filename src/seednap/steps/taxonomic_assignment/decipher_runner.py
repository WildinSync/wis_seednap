"""DECIPHER R package wrapper for taxonomic assignment.

This module provides a Python wrapper around the DECIPHER R package
for taxonomic assignment of eDNA sequences.
"""

import logging
import subprocess
from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


class DecipherError(Exception):
    """Exception raised for DECIPHER command errors."""

    pass


class DecipherRunner:
    """
    Run DECIPHER taxonomic assignment via Rscript.

    DECIPHER uses a trained classifier (created with DECIPHER::LearnTaxa)
    to assign taxonomy to sequences with confidence scores.
    """

    def __init__(self, timeout: int = 7200):
        """
        Initialize DECIPHER runner.

        Args:
            timeout: Command timeout in seconds (default: 7200 = 2 hours)
        """
        self.timeout = timeout
        self._check_r_availability()

    def _check_r_availability(self) -> None:
        """
        Check if R and DECIPHER package are available.

        Raises:
            DecipherError: If R or DECIPHER is not found
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
            raise DecipherError("Rscript not found. Is R installed?") from e

        # Check for DECIPHER package
        r_code = """
        if (!requireNamespace("DECIPHER", quietly = TRUE)) {
            stop("DECIPHER package not installed")
        }
        cat(as.character(packageVersion("DECIPHER")))
        """

        try:
            result = subprocess.run(
                ["Rscript", "-e", r_code],
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )
            version = result.stdout.strip()
            logger.info(f"Found DECIPHER version: {version}")
        except subprocess.CalledProcessError as e:
            raise DecipherError(
                "DECIPHER R package not installed. "
                "Install with: install.packages('DECIPHER')"
            ) from e

    def _run_r_script(
        self,
        script_path: Union[str, Path],
        args: list,
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
            DecipherError: If Rscript command fails
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
            raise DecipherError(error_msg) from e

        except subprocess.TimeoutExpired as e:
            error_msg = f"R script timed out after {self.timeout} seconds: {script_path.name}"
            logger.error(error_msg)
            raise DecipherError(error_msg) from e

    def run_decipher_assignment(
        self,
        marker: str,
        output_dir: Union[str, Path],
        trained_classifier_path: Union[str, Path],
        threshold: int = 60,
        processors: int = 8,
        script_path: Optional[Union[str, Path]] = None,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Dict[str, Path]:
        """
        Run DECIPHER taxonomic assignment.

        Uses the IdTaxa function from DECIPHER R package to assign taxonomy
        with confidence scores.

        Args:
            marker: Marker name
            output_dir: Base output directory
            trained_classifier_path: Path to trained DECIPHER classifier (.rds file)
            threshold: Minimum confidence threshold (0-100, default: 60)
            processors: Number of CPU cores to use (default: 8)
            script_path: Path to R script (default: scripts/taxo_decipher_marker.R)
            log_file: Optional path to log file

        Returns:
            Dictionary with paths to output files:
            - taxonomy: Taxonomy table CSV (with confidence scores)
            - complete: Complete table with taxonomy and abundances

        Raises:
            DecipherError: If DECIPHER assignment fails
            FileNotFoundError: If required inputs are missing
        """
        if script_path is None:
            script_path = Path("scripts/taxo_decipher_marker.R")

        trained_classifier_path = Path(trained_classifier_path)

        if not trained_classifier_path.exists():
            raise FileNotFoundError(
                f"Trained classifier not found: {trained_classifier_path}"
            )

        # Check that sequence table exists
        output_dir = Path(output_dir)
        seqtab_rds = output_dir / "02_dada2" / marker / "seqtab_clean.rds"
        if not seqtab_rds.exists():
            raise FileNotFoundError(
                f"Sequence table not found: {seqtab_rds}. "
                "Run DADA2 processing first."
            )

        logger.info(f"Running DECIPHER taxonomic assignment for {marker}")

        # Run R script
        output = self._run_r_script(
            script_path=script_path,
            args=[marker, str(trained_classifier_path)],
            log_file=log_file,
        )

        # Construct output paths
        marker_dir = output_dir / "02_dada2" / marker

        return {
            "taxonomy": marker_dir / "taxo_assigned_decipher.csv",
            "complete": output_dir / f"{marker}_decipher.csv",
        }

    def link_with_abundance_table(
        self,
        taxonomy_csv: Union[str, Path],
        abundance_csv: Union[str, Path],
        output_csv: Union[str, Path],
        sequence_col: str = "sequence",
    ) -> Path:
        """
        Link DECIPHER taxonomy with DADA2 abundance table.

        Args:
            taxonomy_csv: Path to DECIPHER taxonomy CSV
            abundance_csv: Path to DADA2 abundance table (seqtab_clean_t.csv)
            output_csv: Path to output CSV file
            sequence_col: Name of sequence column (default: 'sequence')

        Returns:
            Path to output CSV file with merged taxonomy and abundances
        """
        taxonomy_csv = Path(taxonomy_csv)
        abundance_csv = Path(abundance_csv)
        output_csv = Path(output_csv)

        if not taxonomy_csv.exists():
            raise FileNotFoundError(f"Taxonomy CSV not found: {taxonomy_csv}")
        if not abundance_csv.exists():
            raise FileNotFoundError(f"Abundance CSV not found: {abundance_csv}")

        # Read taxonomy
        taxo_df = pd.read_csv(taxonomy_csv)

        # Read abundance table
        abundance_df = pd.read_csv(abundance_csv, index_col=0)
        abundance_df = abundance_df.reset_index().rename(columns={"index": sequence_col})

        # Merge
        result = pd.merge(taxo_df, abundance_df, on=sequence_col, how="left")

        # Write output
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_csv, index=False)

        logger.info(f"Linked DECIPHER taxonomy with abundances: {output_csv}")
        return output_csv
