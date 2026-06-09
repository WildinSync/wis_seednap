"""DADA2 workflow orchestration for eDNA metabarcoding.

This module provides high-level orchestration for the complete DADA2 workflow,
integrating R script execution with metrics collection and quality reporting.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Union

from seednap.steps.dada2.dada2_runner import Dada2Runner
from seednap.steps.dada2.metrics import MetricsCollector

logger = logging.getLogger(__name__)


class Dada2Processor:
    """
    Orchestrate complete DADA2 workflow for amplicon sequencing.

    This class coordinates:
    - DADA2 processing (quality control, filtering, denoising, chimera removal)
    - Metrics collection and tracking
    - Optional taxonomic assignment
    - Report generation
    """

    def __init__(
        self,
        marker: str,
        trimmed_reads_dir: Union[str, Path],
        output_base_dir: Union[str, Path],
        timeout: int = 14400,
    ):
        """
        Initialize DADA2 processor.

        Args:
            marker: Marker name (e.g., 'teleo', 'amph')
            trimmed_reads_dir: Directory with primer-trimmed FASTQ files
            output_base_dir: Base output directory
            timeout: Timeout for R scripts in seconds (default: 4 hours)
        """
        self.marker = marker.lower()
        self.trimmed_reads_dir = Path(trimmed_reads_dir)
        self.output_base_dir = Path(output_base_dir)

        if not self.trimmed_reads_dir.exists():
            raise FileNotFoundError(
                f"Trimmed reads directory not found: {self.trimmed_reads_dir}. "
                f"DADA2 needs primer-trimmed paired-end FASTQ files as input. "
                f"If you ran `run-pipeline`, this points at the trim step's output "
                f"(when trim ran) or falls back to paths.raw_data (when trim was "
                f"skipped), so either include \"trim\" before \"dada2\" in "
                f"pipeline.steps, or fix paths.raw_data in the marker YAML to a "
                f"directory that exists. If you invoked the `dada2` subcommand "
                f"directly, pass the directory holding the already-trimmed FASTQs "
                f"as the TRIMMED_READS_DIR argument."
            )

        self.runner = Dada2Runner(timeout=timeout)
        self.metrics = MetricsCollector(marker=self.marker, output_dir=self.output_base_dir)

        self.output_dir = self.output_base_dir / "02_dada2" / self.marker
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Initialized DADA2 processor for marker: {self.marker}")

    def process(
        self,
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
        collect_metrics: bool = True,
    ) -> Dict[str, Path]:
        """
        Run complete DADA2 processing workflow.

        Args:
            max_ee: Maximum expected errors for filtering (default: 2.0)
            trunc_q: Truncate reads at first base with quality < this (default: 11)
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
            collect_metrics: Collect and export metrics (default: True)

        Returns:
            Dictionary with paths to output files

        Raises:
            FileNotFoundError: If required input files are missing
            Dada2Error: If DADA2 processing fails
        """
        logger.info(f"Starting DADA2 processing for {self.marker}")
        logger.info(
            f"Parameters: maxEE={max_ee}, truncQ={trunc_q}, minOverlap={min_overlap}, "
            f"maxN={max_n}, rmPhix={rm_phix}, chimera={chimera_method}, "
            f"maxMismatch={max_mismatch}, pool={pool}, multithread={multithread}"
        )

        # Verify trimmed reads directory has FASTQ files
        fastq_files = list(self.trimmed_reads_dir.glob("*R1*.fastq*"))
        if not fastq_files:
            fastq_files = list(self.trimmed_reads_dir.glob("*R1*.fq*"))
        if not fastq_files:
            raise FileNotFoundError(
                f"No FASTQ files found in trimmed reads directory: {self.trimmed_reads_dir}. "
                f"Check that the trim step produced output or that paths.raw_data is correct."
            )
        logger.info(f"Found {len(fastq_files)} R1 files in {self.trimmed_reads_dir}")

        # Check for required R packages
        logger.info("Checking R package dependencies...")
        packages = self.runner.check_dada2_packages()
        logger.info(f"Found DADA2 version: {packages.get('dada2', 'unknown')}")

        # Run DADA2 processing
        log_file = self.output_dir / "dada2_processing.log"
        outputs = self.runner.run_dada2_process(
            marker=self.marker,
            trimmed_reads_dir=self.trimmed_reads_dir,
            output_dir=self.output_base_dir,
            max_ee=max_ee,
            trunc_q=trunc_q,
            min_overlap=min_overlap,
            max_n=max_n,
            rm_phix=rm_phix,
            multithread=multithread,
            chimera_method=chimera_method,
            max_mismatch=max_mismatch,
            pool=pool,
            min_len=min_len,
            max_len=max_len,
            library_map=library_map,
            log_file=log_file,
        )

        logger.info("DADA2 processing completed successfully")
        logger.info(f"Outputs saved to: {self.output_dir}")

        # Collect metrics if requested
        if collect_metrics:
            logger.info("Collecting pipeline metrics...")
            self._collect_metrics(outputs)

            # Generate and save summary report
            summary = self.metrics.generate_summary_report()
            print(summary)

            summary_file = self.output_dir / "metrics" / "summary.txt"
            summary_file.parent.mkdir(parents=True, exist_ok=True)
            summary_file.write_text(summary)

            # Export metrics
            self.metrics.export_to_json()
            self.metrics.export_to_csv()

        return outputs

    def _collect_metrics(self, outputs: Dict[str, Path]) -> None:
        """
        Collect metrics from DADA2 outputs.

        Args:
            outputs: Dictionary of output paths from DADA2 processing
        """
        # Note: For full metrics collection, we'd need to parse the DADA2 log files
        # or run additional R commands to extract intermediate read counts.
        # For now, we collect what's available from the outputs.

        # Collect ASV metrics
        if outputs["seqtab_clean_t"].exists():
            self.metrics.collect_asv_metrics(
                seqtab_path=outputs["seqtab_clean_t"],
                corresp_seq_path=outputs.get("corresp_seq"),
            )

        logger.info("Metrics collection completed")

    def get_metrics_summary(self) -> str:
        """
        Get human-readable metrics summary.

        Returns:
            Formatted summary string
        """
        return self.metrics.generate_summary_report()
