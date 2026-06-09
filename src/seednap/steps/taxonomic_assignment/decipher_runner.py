"""DECIPHER R package wrapper for taxonomic assignment.

This module provides a Python wrapper around the DECIPHER R package
for taxonomic assignment of eDNA sequences.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Union

from seednap.utils.r_runner import RScriptRunner, r_script_path
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
            stop("DECIPHER R package not found in the active environment. ",
                 "The 'decipher' taxonomy method needs the DECIPHER Bioconductor ",
                 "package under the same R interpreter the pipeline calls. ",
                 "Fix: activate the project env (conda activate metabarcoding), ",
                 "install it the way this repo pins it ",
                 "(conda install -c bioconda bioconductor-decipher; ",
                 "it is also in environment.yml), ",
                 "then verify with Rscript -e 'packageVersion(\\"DECIPHER\\")'.")
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
        query_fasta: Union[str, Path],
        threshold: int = 60,
        processors: int = 8,
        script_path: Optional[Union[str, Path]] = None,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Dict[str, Path]:
        """
        Run DECIPHER taxonomic assignment.

        Reads sequences from a query FASTA (produced by either DADA2 ASV or
        SWARM OTU clustering), runs the IdTaxa classifier with a confidence
        threshold, and writes a per-sequence taxonomy CSV. The merge with
        the abundance table is done by the Python caller via
        `seednap.utils.taxonomy.link_taxonomy_with_abundance`.

        The script is cluster-method agnostic by design -- it does not read
        seqtab_clean.rds.

        Args:
            marker: Marker name (used for log messages and output filenames)
            output_dir: Base output directory
            trained_classifier_path: Path to trained DECIPHER classifier (.rds)
            query_fasta: Path to query.fasta (sequences to assign taxonomy to)
            threshold: Minimum confidence threshold (0-100, default: 60)
            processors: Number of CPU cores to use (default: 8)
            script_path: Path to R script (default: scripts/taxo_decipher_marker.R)
            log_file: Optional path to log file

        Returns:
            Dictionary with paths to output files:
            - taxonomy: Per-sequence taxonomy CSV (with confidence scores)
            - final_table: Where the merged output WILL go (Python writes it)

        Raises:
            DecipherError: If DECIPHER assignment fails
            FileNotFoundError: If required inputs are missing
        """
        if script_path is None:
            script_path = r_script_path("taxo_decipher_marker.R")

        trained_classifier_path = Path(trained_classifier_path)
        query_fasta = Path(query_fasta)
        output_dir = Path(output_dir)

        if not trained_classifier_path.exists():
            raise FileNotFoundError(
                f"Trained DECIPHER classifier not found: {trained_classifier_path}. "
                f"DECIPHER IdTaxa needs a trained classifier .rds (built with "
                f"DECIPHER::LearnTaxa), set via taxonomy.databases.decipher.trained in "
                f"the marker YAML (or --trained-classifier for the standalone "
                f"assign-taxonomy command), but that path is missing or wrong for this "
                f"host. Config load and `seednap validate` resolve but do not require "
                f"this path to exist, so a typo or a config copied from another server "
                f"reaches this point. Fix: point taxonomy.databases.decipher.trained "
                f"(or --trained-classifier) at the trained .rds for this marker and "
                f"confirm the file exists on this machine."
            )
        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")

        marker_dir = output_dir / "02_dada2" / marker
        marker_dir.mkdir(parents=True, exist_ok=True)
        taxonomy_csv = marker_dir / "taxo_assigned_decipher.csv"

        logger.info(
            f"Running DECIPHER taxonomic assignment for {marker} "
            f"(threshold={threshold}, query={query_fasta})"
        )

        self._run_r_script(
            script_path=script_path,
            args=[
                marker,
                str(trained_classifier_path),
                str(query_fasta),
                str(taxonomy_csv),
                str(threshold),
                str(processors),
            ],
            log_file=log_file,
        )

        return {
            "taxonomy": taxonomy_csv,
            "final_table": output_dir / f"{marker}_decipher.csv",
        }

    def link_with_abundance_table(
        self,
        taxonomy_csv: Union[str, Path],
        abundance_csv: Union[str, Path],
        output_csv: Union[str, Path],
        sequence_col: str = "sequence",
        contaminants: Optional[list] = None,
    ) -> Path:
        """
        Link DECIPHER taxonomy with the DADA2/SWARM abundance table.

        Delegates to the shared taxonomy post-processor so DECIPHER, ecotag,
        DADA2 RDP, and BLAST all share the same output schema and the same
        correctness guarantees (left-merge from abundance side, cascade null,
        contaminant flagging, stable column order).

        Args:
            taxonomy_csv: Path to DECIPHER taxonomy CSV
            abundance_csv: Path to abundance table (seqtab_clean_t.csv)
            output_csv: Path to output CSV file
            sequence_col: Name of sequence column (default: 'sequence')
            contaminants: Optional list of species to flag as contaminants

        Returns:
            Path to output CSV file with merged taxonomy and abundances
        """
        from seednap.utils.taxonomy import link_taxonomy_with_abundance

        return link_taxonomy_with_abundance(
            taxonomy_path=taxonomy_csv,
            abundance_path=abundance_csv,
            output_path=output_csv,
            sequence_col=sequence_col,
            contaminants=contaminants,
        )
