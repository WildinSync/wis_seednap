"""SWARM clustering algorithm wrapper.

Provides a Python interface around the SWARM binary for OTU clustering
of dereplicated amplicon sequences.
"""

import logging
import shutil
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
        # Preflight with shutil.which so a missing binary is reported with
        # actionable guidance regardless of how a given swarm build handles a
        # `--version` exit code (some builds print the version then exit
        # non-zero, which would otherwise surface as a raw "Command failed"
        # message). Matches the shutil.which preflight pattern in ecotag_runner.
        if shutil.which("swarm") is None:
            raise SwarmError(
                "swarm not found on PATH. SWARM OTU clustering needs the 'swarm' "
                "binary, which is not on your PATH. This almost always means the "
                "seednap conda environment is not activated (it ships swarm alongside "
                "the pipeline).\n"
                "  Fix: activate the environment before running seednap, then re-run:\n"
                "    conda activate /home/shared/edna/envs/seednap   # ETH ELE eDNA server\n"
                "    conda activate seednap                          # local development (environment.yml)\n"
                "  If swarm is still missing after activation, install it with "
                "'conda install -c bioconda swarm'. Confirm with 'swarm --version'."
            )
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
            boundary: Mass threshold for fastidious mode (swarm --boundary,
                default: 3). OTUs whose total abundance (mass) exceeds this
                value are "large"/heavy; OTUs with mass at or below it are
                "small"/light and are the candidates grafted onto large OTUs.
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

        try:
            run_subprocess(cmd, timeout=self.timeout, log_file=log_file, error_class=SwarmError)
        except SwarmError as e:
            log_hint = f" Full swarm output is at {log_file}." if log_file else ""
            raise SwarmError(
                "swarm clustering failed (non-zero exit). The most common cause is "
                "a config mismatch: swarm's fastidious mode requires d=1. If you set "
                "swarm.clustering.d > 1 in the marker config, either set "
                "swarm.clustering.d=1 or set swarm.clustering.fastidious=false. Less "
                f"commonly, the dereplicated input ({input_fasta}) is missing the "
                ";size=N abundance annotations swarm needs.\n"
                f"  Fix: check the swarm.clustering.d and swarm.clustering.fastidious "
                f"keys in the marker config, then re-run with --resume.{log_hint}\n"
                f"  Underlying error:\n{e}"
            ) from e

        logger.info(f"SWARM clustering completed (d={d}, fastidious={fastidious})")

        return {
            "swarm_file": swarm_file,
            "stats_file": stats_file,
            "representatives": representatives,
            "struct_file": struct_file,
        }
