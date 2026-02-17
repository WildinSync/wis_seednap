"""DECIPHER R package wrapper for taxonomic assignment.

This module provides a Python wrapper around the DECIPHER R package
for taxonomic assignment of eDNA sequences.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Union

from seednap.utils.r_runner import RScriptRunner
from seednap.utils.subprocess import run_subprocess

logger = logging.getLogger(__name__)


class DecipherError(Exception):
    """Exception raised for DECIPHER command errors."""

    pass


class DecipherRunner(RScriptRunner):
    """
    Run DECIPHER taxonomic assignment via Rscript.

    DECIPHER uses a trained classifier (created with DECIPHER::LearnTaxa)
    to assign taxonomy to sequences with confidence scores.
    """

    _error_class = DecipherError

    def __init__(self, timeout: int = 7200):
        """
        Initialize DECIPHER runner.

        Args:
            timeout: Command timeout in seconds (default: 7200 = 2 hours)
        """
        super().__init__(timeout=timeout)
        self._check_decipher_package()

    def _check_decipher_package(self) -> None:
        """
        Check if DECIPHER R package is installed.

        Raises:
            DecipherError: If DECIPHER is not found
        """
        r_code = """
        if (!requireNamespace("DECIPHER", quietly = TRUE)) {
            stop("DECIPHER package not installed")
        }
        cat(as.character(packageVersion("DECIPHER")))
        """

        stdout = run_subprocess(
            ["Rscript", "-e", r_code],
            timeout=30,
            error_class=DecipherError,
        )
        version = stdout.strip()
        logger.info(f"Found DECIPHER version: {version}")

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
            - final_table: Complete table with taxonomy and abundances

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
        self._run_r_script(
            script_path=script_path,
            args=[marker, str(trained_classifier_path), str(threshold), str(processors), str(output_dir)],
            log_file=log_file,
        )

        # Construct output paths
        marker_dir = output_dir / "02_dada2" / marker

        return {
            "taxonomy": marker_dir / "taxo_assigned_decipher.csv",
            "final_table": output_dir / f"{marker}_decipher.csv",
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
        from seednap.utils.taxonomy import link_taxonomy_with_abundance

        return link_taxonomy_with_abundance(
            taxonomy_path=taxonomy_csv,
            abundance_path=abundance_csv,
            output_path=output_csv,
            sequence_col=sequence_col,
        )
