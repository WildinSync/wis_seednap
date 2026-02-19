"""SWARM clustering algorithm wrapper.

Provides a Python interface around the SWARM binary for OTU clustering
of dereplicated amplicon sequences.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Union

from seednap.utils.subprocess import run_subprocess

logger = logging.getLogger(__name__)


class SwarmError(Exception):
    """Exception raised for SWARM clustering errors."""

    pass


class SwarmClusterer:
    """
    Wrapper around the SWARM clustering binary.

    SWARM uses a local linking algorithm to cluster amplicon sequences
    into OTUs based on a distance threshold.
    """

    def __init__(self, timeout: int = 3600):
        """
        Initialize SWARM clusterer.

        Args:
            timeout: Command timeout in seconds (default: 1 hour)

        Raises:
            SwarmError: If swarm is not installed
        """
        self.timeout = timeout
        self._check_availability()

    def _check_availability(self) -> None:
        """Check that swarm is installed."""
        run_subprocess(
            ["swarm", "--version"],
            timeout=10,
            error_class=SwarmError,
        )

    def cluster(
        self,
        input_fasta: Union[str, Path],
        output_dir: Union[str, Path],
        *,
        d: int = 1,
        fastidious: bool = True,
        boundary: int = 3,
        threads: int = 4,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Dict[str, Path]:
        """
        Run SWARM clustering on dereplicated sequences.

        Args:
            input_fasta: Input FASTA with ;size=N; abundance annotations
            output_dir: Directory for SWARM output files
            d: Clustering distance threshold (default: 1)
            fastidious: Enable fastidious mode to refine singletons (default: True)
            boundary: Min mass for large OTUs in fastidious mode (default: 3)
            threads: Number of threads (default: 4)
            log_file: Optional log file path

        Returns:
            Dictionary with paths to output files:
            - swarm_file: Cluster membership file
            - stats_file: Clustering statistics
            - representatives: Seed/representative sequences
            - struct_file: Internal cluster structure

        Raises:
            SwarmError: If clustering fails
        """
        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

        swarm_file = output_dir / "all.swarm"
        stats_file = output_dir / "all.stats"
        representatives = output_dir / "all.representatives"
        struct_file = output_dir / "all.struct"

        cmd = [
            "swarm",
            str(input_fasta),
            "-d", str(d),
            "-t", str(threads),
            "--usearch-abundance",
            "--internal-structure", str(struct_file),
            "-s", str(stats_file),
            "--seeds", str(representatives),
            "-o", str(swarm_file),
        ]

        if fastidious:
            cmd.extend(["--fastidious", "--boundary", str(boundary)])

        run_subprocess(cmd, timeout=self.timeout, log_file=log_file, error_class=SwarmError)

        logger.info(f"SWARM clustering completed (d={d}, fastidious={fastidious})")

        return {
            "swarm_file": swarm_file,
            "stats_file": stats_file,
            "representatives": representatives,
            "struct_file": struct_file,
        }
