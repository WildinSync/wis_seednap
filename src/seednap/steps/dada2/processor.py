"""DADA2 workflow orchestration for eDNA metabarcoding.

This module provides high-level orchestration for the complete DADA2 workflow,
integrating R script execution with ASV metrics collection and a short metrics
summary. The QC plots are produced by the R script (dada2_process.R), not by
this Python layer.
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
      by invoking the dada2_process.R script via Dada2Runner
    - ASV metrics collection and a short metrics summary (summary.txt + JSON/CSV)

    Taxonomic assignment is a separate orchestrator step
    (steps/taxonomic_assignment/), and the pipeline run report (04_report) is
    produced elsewhere; neither happens in this module.
    """

    def __init__(
        self,
        marker: str,
        trimmed_reads_dir: Union[str, Path],
        output_base_dir: Union[str, Path],
        timeout: int = 14400,
    ) -> None:
        """
        Initialize DADA2 processor.

        Args:
            marker: Marker name (e.g., 'teleo', 'amph'); lowercased internally.
            trimmed_reads_dir: Directory with primer-trimmed paired-end FASTQ files
                (the DADA2 input).
            output_base_dir: Base output directory; DADA2 outputs are written under
                ``<output_base_dir>/02_dada2/<marker>/``.
            timeout: Timeout for the R scripts in seconds (default: 14400 = 4 hours).

        Raises:
            FileNotFoundError: If ``trimmed_reads_dir`` does not exist.
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
            pool: Pool samples for denoising (default: False). Ignored when
                library_map resolves to >=2 libraries: the per-library branch
                in dada2_process.R always denoises each sample independently
                regardless of this flag.
            min_len: Minimum read length (None = no filter)
            max_len: Maximum read length (None = no filter)
            library_map: Optional path to a sample-to-library map used for
                per-library error learning (None = no map). When it groups
                samples into >=2 libraries, errors are learned per library and
                samples are denoised per-sample within each library (pool is
                not honored on this path). With 0 or 1 library the standard
                single-batch path runs and pool applies normally.
            collect_metrics: Collect and export ASV metrics (summary.txt +
                JSON/CSV) after processing (default: True)

        Returns:
            Dictionary mapping output names to paths, as produced by
            Dada2Runner.run_dada2_process: "seqtab", "seqtab_clean",
            "seqtab_clean_rds", "seqtab_clean_t", "query_fasta", "corresp_seq",
            and "metrics_dir".

        Raises:
            FileNotFoundError: If no R1 FASTQ files are found in the trimmed
                reads directory.
            Dada2Error: If the R package check fails or DADA2 processing fails.
            SeednapError: If metrics collection runs but the sequence table is
                empty (0 ASVs).
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
        Collect ASV-level metrics from the final DADA2 sequence table.

        Reads the transposed chimera-free sequence table (one row per ASV) and,
        when present, the ASV-to-sequence correspondence file, and stores the
        derived statistics on the internal MetricsCollector. Per-step read
        counts are not recomputed here; the R script writes those to
        track_reads.csv / feature_counts.csv for the run report.

        Args:
            outputs: Dictionary of output paths from DADA2 processing; uses the
                "seqtab_clean_t" (transposed chimera-free table) and optional
                "corresp_seq" (ASV correspondence) entries.

        Returns:
            None. Results are stored on ``self.metrics``.

        Raises:
            SeednapError: Propagated from collect_asv_metrics if the sequence
                table exists but is empty (0 bytes or zero ASVs).
        """
        # This collects only ASV-level metrics from the final sequence table.
        # Per-step intermediate read counts are written by dada2_process.R to
        # track_reads.csv / feature_counts.csv and surfaced by the run report;
        # they are not recomputed here.

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
