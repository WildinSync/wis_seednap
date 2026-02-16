"""DADA2 workflow orchestration for eDNA metabarcoding.

This module provides high-level orchestration for the complete DADA2 workflow,
integrating R script execution with metrics collection and quality reporting.
"""

import logging
from pathlib import Path
from typing import Dict, Union

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
            raise FileNotFoundError(f"Trimmed reads directory not found: {self.trimmed_reads_dir}")

        self.runner = Dada2Runner(timeout=timeout)
        self.metrics = MetricsCollector(marker=self.marker, output_dir=self.output_base_dir)

        self.output_dir = self.output_base_dir / "02_dada2" / self.marker
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Initialized DADA2 processor for marker: {self.marker}")

    def process(
        self,
        max_ee: int = 2,
        trunc_q: int = 11,
        min_overlap: int = 20,
        collect_metrics: bool = True,
    ) -> Dict[str, Path]:
        """
        Run complete DADA2 processing workflow.

        Workflow steps:
        1. Quality control (pre-filtering plots)
        2. Filter and trim
        3. Quality control (post-filtering plots)
        4. Learn error rates
        5. Denoise (sample inference)
        6. Merge paired-end reads
        7. Make sequence table
        8. Remove chimeras
        9. Generate outputs (CSV, FASTA, RDS)
        10. Collect metrics (optional)

        Args:
            max_ee: Maximum expected errors for filtering (default: 2)
            trunc_q: Truncate reads at first base with quality < this (default: 11)
            min_overlap: Minimum overlap for merging paired reads (default: 20)
            collect_metrics: Collect and export metrics (default: True)

        Returns:
            Dictionary with paths to output files

        Raises:
            FileNotFoundError: If required input files are missing
            Dada2Error: If DADA2 processing fails
        """
        logger.info(f"Starting DADA2 processing for {self.marker}")
        logger.info(f"Parameters: maxEE={max_ee}, truncQ={trunc_q}, minOverlap={min_overlap}")

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
            log_file=log_file,
        )

        logger.info(f"DADA2 processing completed successfully")
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

    def assign_taxonomy(
        self,
        rdp_db_path: Union[str, Path],
        species_db_path: Union[str, Path],
    ) -> Dict[str, Path]:
        """
        Perform taxonomic assignment using DADA2's naive Bayesian classifier.

        This should be run after process() completes successfully.

        Args:
            rdp_db_path: Path to RDP-formatted taxonomy database (genus-level)
            species_db_path: Path to species-level database

        Returns:
            Dictionary with paths to taxonomy output files

        Raises:
            FileNotFoundError: If database files or required inputs are missing
            Dada2Error: If taxonomy assignment fails
        """
        logger.info(f"Starting DADA2 taxonomic assignment for {self.marker}")

        # Check that sequence table exists
        seqtab_rds = self.output_dir / "seqtab_clean.rds"
        if not seqtab_rds.exists():
            raise FileNotFoundError(
                f"Sequence table not found: {seqtab_rds}. "
                "Run process() before assign_taxonomy()."
            )

        from seednap.steps.taxonomic_assignment.dada2_taxonomy_runner import Dada2TaxonomyRunner

        log_file = self.output_dir / "dada2_taxonomy.log"
        taxonomy_runner = Dada2TaxonomyRunner()
        outputs = taxonomy_runner.run_dada2_taxonomy(
            marker=self.marker,
            output_dir=self.output_base_dir,
            rdp_db_path=rdp_db_path,
            species_db_path=species_db_path,
            log_file=log_file,
        )

        logger.info(f"Taxonomic assignment completed successfully")
        logger.info(f"Taxonomy table: {outputs['final_table']}")

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
