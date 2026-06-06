"""DADA2 subprocess wrapper for R script execution.

This module provides Python wrappers around DADA2 R scripts for amplicon
sequence denoising, error correction, and ASV generation.
"""

import logging
import subprocess
from pathlib import Path
from typing import Dict, Optional, Union

from seednap.utils.r_runner import RScriptRunner

logger = logging.getLogger(__name__)


class Dada2Error(Exception):
    """Exception raised for DADA2 R script errors."""

    pass


class Dada2Runner(RScriptRunner):
    """Run DADA2 R scripts via Rscript subprocess.

    This class wraps the DADA2 R scripts and provides methods for:
    - Quality control and filtering
    - Error learning and sample inference
    - Chimera removal
    - ASV table generation
    """

    _error_class = Dada2Error

    def __init__(self, timeout: int = 14400):
        """
        Initialize DADA2 runner.

        Args:
            timeout: Command timeout in seconds (default: 14400 = 4 hours)
        """
        super().__init__(timeout=timeout)

    def run_dada2_process(
        self,
        marker: str,
        trimmed_reads_dir: Union[str, Path],
        output_dir: Union[str, Path],
        max_ee: float = 2.0,
        trunc_q: int = 11,
        min_overlap: int = 20,
        max_n: int = 0,
        rm_phix: bool = True,
        multithread: bool = True,
        chimera_method: str = "consensus",
        max_mismatch: int = 0,
        pool: bool = False,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
        library_map: Optional[Union[str, Path]] = None,
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
            max_ee: Maximum expected errors after filtering (default: 2.0)
            trunc_q: Truncate reads at first quality score below this (default: 11)
            min_overlap: Minimum overlap for merging paired reads (default: 20)
            max_n: Maximum number of N bases allowed (default: 0)
            rm_phix: Remove PhiX reads (default: True)
            multithread: Use multithreading (default: True)
            chimera_method: Chimera detection method (default: "consensus")
            max_mismatch: Maximum mismatches in overlap region (default: 0)
            pool: Pool samples for denoising (default: False)
            min_len: Minimum read length (None = no filter)
            max_len: Maximum read length (None = no filter)
            library_map: Optional path to a sample-to-library map used for
                per-library error learning (None = no map)
            script_path: Path to R script (default: scripts/dada2_process.R)
            log_file: Path to log file

        Returns:
            Dictionary with paths to output files:
            - seqtab: Sequence table RDS
            - seqtab_clean: Chimera-free sequence table CSV
            - seqtab_clean_rds: Chimera-free sequence table RDS
            - seqtab_clean_t: Transposed chimera-free sequence table CSV
            - query_fasta: Query FASTA for taxonomic assignment
            - corresp_seq: ASV correspondence CSV
            - metrics_dir: Directory with QC/metrics plots

        Raises:
            Dada2Error: If processing fails
        """
        if script_path is None:
            # Default to the legacy R script
            script_path = Path("scripts/dada2_process.R")

        self._run_r_script(
            script_path=script_path,
            args=[
                marker,
                str(trimmed_reads_dir),
                str(output_dir),
                str(max_ee),
                str(trunc_q),
                str(min_overlap),
                str(max_n),
                str(rm_phix).upper(),
                str(multithread).upper(),
                chimera_method,
                str(max_mismatch),
                str(pool).upper(),
                str(min_len if min_len is not None else 0),
                str(max_len if max_len is not None else 0),
                str(library_map) if library_map is not None else "",
            ],
            log_file=log_file,
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
            "metrics_dir": marker_dir / "QC",  # DADA2 QC plot directory
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
