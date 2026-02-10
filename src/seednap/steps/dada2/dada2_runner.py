"""DADA2 subprocess wrapper for R script execution.

This module provides Python wrappers around DADA2 R scripts for amplicon
sequence denoising, error correction, and ASV generation.
"""

import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)


class Dada2Error(Exception):
    """Exception raised for DADA2 R script errors."""

    pass


class Dada2Runner:
    """Run DADA2 R scripts via Rscript subprocess.

    This class wraps the DADA2 R scripts and provides methods for:
    - Quality control and filtering
    - Error learning and sample inference
    - Chimera removal
    - ASV table generation
    """

    def __init__(self, timeout: int = 14400):
        """
        Initialize DADA2 runner.

        Args:
            timeout: Command timeout in seconds (default: 14400 = 4 hours)
        """
        self.timeout = timeout
        self._check_r_availability()

    def _check_r_availability(self) -> None:
        """
        Check if R and required packages are available.

        Raises:
            Dada2Error: If R or Rscript is not found
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
            raise Dada2Error("Rscript not found. Is R installed?") from e

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
            Dada2Error: If Rscript command fails
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
            raise Dada2Error(error_msg) from e

        except subprocess.TimeoutExpired as e:
            error_msg = f"R script timed out after {self.timeout} seconds: {script_path.name}"
            logger.error(error_msg)
            raise Dada2Error(error_msg) from e

    def run_dada2_process(
        self,
        marker: str,
        trimmed_reads_dir: Union[str, Path],
        output_dir: Union[str, Path],
        max_ee: int = 2,
        trunc_q: int = 11,
        min_overlap: int = 20,
        script_path: Optional[Union[str, Path]] = None,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Dict[str, Path]:
        """
        Run complete DADA2 processing workflow.

        This executes the main DADA2 pipeline:
        1. Quality plots (pre-filtering)
        2. Filter and trim
        3. Quality plots (post-filtering)
        4. Learn error rates
        5. Sample inference
        6. Merge paired-end reads
        7. Make sequence table
        8. Remove chimeras

        Args:
            marker: Marker name (e.g., 'teleo', 'amph')
            trimmed_reads_dir: Directory with primer-trimmed FASTQ files
            output_dir: Base output directory
            max_ee: Maximum expected errors after filtering (default: 2)
            trunc_q: Truncate reads at first quality score below this (default: 11)
            min_overlap: Minimum overlap for merging paired reads (default: 20)
            script_path: Path to R script (default: scripts/dada2_process.R)
            log_file: Path to log file

        Returns:
            Dictionary with paths to output files:
            - seqtab: Sequence table CSV
            - seqtab_clean: Chimera-free sequence table
            - query_fasta: Query FASTA for taxonomic assignment
            - corresp_seq: ASV correspondence CSV
            - metrics_dir: Directory with QC/metrics plots

        Raises:
            Dada2Error: If processing fails
        """
        if script_path is None:
            # Default to the legacy R script
            script_path = Path("scripts/dada2_process.R")

        # For now, use the existing R script with marker argument
        # In the future, we'll create a parameterized version
        output = self._run_r_script(
            script_path=script_path, args=[marker], log_file=log_file
        )

        # Construct output paths
        output_dir = Path(output_dir)
        marker_dir = output_dir / "02_dada2" / marker

        return {
            "seqtab": marker_dir / "seqtab.rds",
            "seqtab_clean": marker_dir / "seqtab_clean.csv",
            "seqtab_clean_rds": marker_dir / "seqtab_clean.rds",
            "seqtab_clean_t": marker_dir / "seqtab_clean_t.csv",
            "query_fasta": marker_dir / "query.fasta",
            "corresp_seq": marker_dir / "corresp_seq.csv",
            "metrics_dir": marker_dir / "QC",  # Will rename to metrics
        }

    def check_dada2_packages(self) -> Dict[str, str]:
        """
        Check DADA2 and related R package versions.

        Returns:
            Dictionary mapping package names to versions

        Raises:
            Dada2Error: If packages are not installed
        """
        r_code = """
        packages <- c("dada2", "Biostrings", "DECIPHER", "dplyr", "ggplot2", "patchwork")
        versions <- sapply(packages, function(pkg) {
            if (requireNamespace(pkg, quietly = TRUE)) {
                as.character(packageVersion(pkg))
            } else {
                "NOT_INSTALLED"
            }
        })
        cat(paste(names(versions), versions, sep=":", collapse="\\n"))
        """

        try:
            result = subprocess.run(
                ["Rscript", "-e", r_code],
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )

            versions = {}
            for line in result.stdout.strip().split("\n"):
                if ":" in line:
                    pkg, ver = line.split(":", 1)
                    versions[pkg.strip()] = ver.strip()

            # Check for missing packages
            missing = [pkg for pkg, ver in versions.items() if ver == "NOT_INSTALLED"]
            if missing:
                raise Dada2Error(
                    f"Required R packages not installed: {', '.join(missing)}"
                )

            logger.info(f"Found R packages: {versions}")
            return versions

        except subprocess.CalledProcessError as e:
            raise Dada2Error(f"Failed to check R packages: {e.stderr}") from e
