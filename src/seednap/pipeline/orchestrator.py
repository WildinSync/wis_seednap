"""Pipeline orchestrator for end-to-end workflow execution.

This module orchestrates the complete seednap pipeline:
1. Demultiplexing (optional)
2. Primer trimming (cutadapt)
3. DADA2 processing
4. Taxonomic assignment
5. Export to GBIF format
"""

from pathlib import Path
from typing import Dict, List, Optional, Union

from seednap.config.loader import load_config
from seednap.config.models import PipelineConfig
from seednap.pipeline.state import PipelineState
from seednap.steps.dada2.processor import Dada2Processor
from seednap.steps.formatting.gbif_formatter import GBIFFormatter
from seednap.steps.taxonomic_assignment.assigner import TaxonomicAssigner
from seednap.steps.trimming.tag_generator import TagFileGenerator
from seednap.steps.trimming.trimming_pipeline import LigationTrimmer, StandardTrimmer
from seednap.utils.logging import get_logger, log_pipeline_step, setup_logging

logger = get_logger(__name__)


class PipelineOrchestrator:
    """
    Orchestrate complete SeeDNAP eDNA metabarcoding pipeline.

    This class coordinates all pipeline steps, handles state management,
    and enables resumability after failures.
    """

    def __init__(
        self,
        config: Union[str, Path, PipelineConfig],
        state_file: Optional[Union[str, Path]] = None,
        resume: bool = False,
    ):
        """
        Initialize pipeline orchestrator.

        Args:
            config: Pipeline configuration (YAML path or PipelineConfig object)
            state_file: Path to state file for tracking progress (default: auto-generated)
            resume: Whether to resume from previous run (default: False)

        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If resume=True but no state file exists
        """
        # Load configuration
        if isinstance(config, (str, Path)):
            self.config_path = Path(config)
            self.config = load_config(self.config_path)
        else:
            self.config_path = None
            self.config = config

        # Setup logging
        self._setup_logging()

        # Create output directory structure
        self._create_directories()

        # Setup state management
        if state_file is None:
            state_file = (
                self.config.paths.output / f".{self.config.marker.name}_state.json"
            )
        self.state_file = Path(state_file)

        # Load or create pipeline state
        if resume:
            if not self.state_file.exists():
                raise ValueError(
                    f"Cannot resume: state file not found at {self.state_file}"
                )
            self.state = PipelineState.load(self.state_file)
            logger.info(f"Resuming pipeline from {self.state_file}")
            logger.info(
                f"Completed steps: {', '.join(self.state.get_completed_steps())}"
            )
        else:
            self.state = PipelineState.from_config(
                marker=self.config.marker.name, config_path=self.config_path
            )
            # Initialize all steps as pending
            for step in self.config.pipeline.steps:
                self.state.add_step(step)
            logger.info(
                f"Starting new pipeline run for marker: {self.config.marker.name}"
            )

    def _setup_logging(self) -> None:
        """Configure logging using existing logging utilities."""
        log_config = self.config.logging

        # Determine log file path
        log_file = None
        if log_config.file:
            if hasattr(self, "state") and self.state.started_at:
                timestamp = self.state.started_at.strftime("%Y%m%d_%H%M%S")
            else:
                timestamp = "run"
            log_file = (
                self.config.paths.logs
                / f"{self.config.marker.name}_pipeline_{timestamp}.log"
            )

        # Setup logging using utility function
        setup_logging(
            level=log_config.level,
            log_file=log_file,
            format_style=log_config.format,
            console_output=log_config.console,
        )

        if log_file:
            logger.info(f"Logging to file: {log_file}")

    def _create_directories(self) -> None:
        """Create output directory structure."""
        marker = self.config.marker.name
        base = self.config.paths.output

        # Create main output directories
        (base / "01_trim" / marker).mkdir(parents=True, exist_ok=True)
        (base / "02_dada2" / marker).mkdir(parents=True, exist_ok=True)
        (base / "03_taxo" / marker).mkdir(parents=True, exist_ok=True)

        # Create logs directory
        self.config.paths.logs.mkdir(parents=True, exist_ok=True)

        logger.debug(f"Created output directory structure in {base}")

    def _save_state(self) -> None:
        """Save current pipeline state."""
        self.state.save(self.state_file)

    def _should_run_step(self, step_name: str) -> bool:
        """
        Determine if a step should be run.

        Args:
            step_name: Name of the step

        Returns:
            True if step should be run, False if it should be skipped
        """
        # Check if step is in the skip list
        if step_name in self.config.pipeline.skip:
            return False

        # Check if step is already completed
        if self.state.is_step_completed(step_name):
            logger.info(f"Step '{step_name}' already completed, skipping")
            return False

        # Check if step failed previously
        if self.state.is_step_failed(step_name):
            logger.warning(f"Step '{step_name}' failed previously, retrying")
            return True

        return True

    def run_demultiplex(self) -> Dict[str, Path]:
        """
        Run demultiplexing step.

        Returns:
            Dictionary with output paths

        Raises:
            ValueError: If demultiplex is enabled but metadata is missing
        """
        step_name = "demultiplex"

        if not self.config.demultiplex.enabled:
            self.state.skip_step(step_name, reason="Demultiplexing disabled in config")
            return {}

        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            logger.info(f"Running demultiplexing: {self.config.demultiplex.protocol}")

            if self.config.demultiplex.protocol == "ligation":
                outputs = self._run_ligation_demux()
            elif self.config.demultiplex.protocol == "standard":
                outputs = self._run_standard_demux()
            else:
                raise ValueError(
                    f"Unknown demultiplex protocol: {self.config.demultiplex.protocol}"
                )

            self.state.complete_step(step_name, outputs)
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs

        except Exception as e:
            self.state.fail_step(step_name, e)
            self._save_state()
            log_pipeline_step(step_name, "error", logger)
            raise

    def _run_ligation_demux(self) -> Dict[str, Path]:
        """Run ligation-based demultiplexing."""
        if self.config.demultiplex.metadata is None:
            raise ValueError("Metadata file required for ligation demultiplexing")

        trimmer = LigationTrimmer()
        output_dir = (
            self.config.paths.output / "01_trim" / self.config.marker.name / "demux"
        )

        # Process library
        outputs = trimmer.process_library(
            raw_reads_dir=self.config.paths.raw_data,
            library_name=self.config.marker.name,
            metadata_csv=self.config.demultiplex.metadata,
            output_base_dir=output_dir,
            forward_primer=self.config.marker.primers.forward,
            reverse_primer=self.config.marker.primers.reverse,
        )

        return {"demux_dir": output_dir, "trimmed_dir": outputs}

    def _run_standard_demux(self) -> Dict[str, Path]:
        """Run standard demultiplexing."""
        if self.config.demultiplex.metadata is None:
            raise ValueError("Metadata file required for standard demultiplexing")

        # Generate tag files
        tag_gen = TagFileGenerator()
        tag_dir = self.config.paths.output / "tags"
        tag_files = tag_gen.generate_standard_tag_files(
            metadata_csv=self.config.demultiplex.metadata, output_dir=tag_dir
        )

        # TODO: Implement standard demultiplexing workflow
        # This would use CutadaptRunner.demultiplex_by_tags()

        return {"tag_files": tag_files, "tag_dir": tag_dir}

    def run_trim(self) -> Dict[str, Path]:
        """
        Run primer trimming step.

        Returns:
            Dictionary with output paths
        """
        step_name = "trim"

        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            logger.info("Running primer trimming")

            trimmer = StandardTrimmer(
                cores=self.config.trimming.cores,
                error_rate=self.config.trimming.max_error_rate,
                min_length=self.config.trimming.min_length
            )
            
            output_dir = self.config.paths.output / "01_trim" / self.config.marker.name

            # Get list of samples from raw data directory
            samples = self._get_sample_list()
            logger.info(f"Found {len(samples)} samples to trim")

            trimmed_outputs = {}
            for sample_name in samples:
                logger.info(f"Trimming sample: {sample_name}")

                r1_input = self._find_read_file(sample_name, "R1")
                r2_input = self._find_read_file(sample_name, "R2")

                outputs = trimmer.trim_sample(
                    r1_input=r1_input,
                    r2_input=r2_input,
                    output_dir=output_dir,
                    sample_name=sample_name,
                    forward_primer=self.config.marker.primers.forward,
                    reverse_primer=self.config.marker.primers.reverse,
                    keep_untrimmed=not self.config.trimming.discard_untrimmed
                )

                trimmed_outputs[sample_name] = outputs

            outputs = {"trimmed_dir": output_dir, "samples": trimmed_outputs}

            self.state.complete_step(step_name, outputs)
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs

        except Exception as e:
            self.state.fail_step(step_name, e)
            self._save_state()
            log_pipeline_step(step_name, "error", logger)
            raise

    def run_dada2(self) -> Dict[str, Path]:
        """
        Run DADA2 processing step.

        Returns:
            Dictionary with output paths
        """
        step_name = "dada2"

        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            logger.info("Running DADA2 processing")

            # Determine trimmed reads directory
            if self.state.is_step_completed("trim"):
                trim_step = self.state.get_step("trim")
                trimmed_reads_dir = trim_step.outputs.get("trimmed_dir") if trim_step else None
                if trimmed_reads_dir is None:
                    raise ValueError("Trim step completed but no trimmed_dir in outputs")
                trimmed_reads_dir = Path(trimmed_reads_dir)
            else:
                # Use raw data if trimming was skipped
                trimmed_reads_dir = self.config.paths.raw_data

            processor = Dada2Processor(
                marker=self.config.marker.name,
                trimmed_reads_dir=trimmed_reads_dir,
                output_base_dir=self.config.paths.output,
            )

            outputs = processor.process(
                max_ee=self.config.dada2.filter.max_ee,
                trunc_q=self.config.dada2.filter.trunc_q,
                min_overlap=self.config.dada2.merge.min_overlap,
                collect_metrics=self.config.metrics.generate_plots,
            )

            self.state.complete_step(step_name, outputs)
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs

        except Exception as e:
            self.state.fail_step(step_name, e)
            self._save_state()
            log_pipeline_step(step_name, "error", logger)
            raise

    def run_taxonomy(self) -> Dict[str, Path]:
        """
        Run taxonomic assignment step.

        Returns:
            Dictionary with output paths
        """
        step_name = "taxonomy"

        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            logger.info(f"Running taxonomic assignment: {self.config.taxonomy.method}")

            # Get DADA2 outputs
            if not self.state.is_step_completed("dada2"):
                raise ValueError("DADA2 step must be completed before taxonomy")

            dada2_step = self.state.get_step("dada2")
            if dada2_step is None:
                raise ValueError("DADA2 step not found in pipeline state")
            query_fasta = dada2_step.outputs.get("query_fasta")
            asv_count_csv = dada2_step.outputs.get("seqtab_clean_t")
            if query_fasta is None or asv_count_csv is None:
                raise ValueError(
                    f"DADA2 outputs incomplete: query_fasta={query_fasta}, "
                    f"seqtab_clean_t={asv_count_csv}"
                )

            # Create taxonomic assigner
            assigner = TaxonomicAssigner(
                method=self.config.taxonomy.method,
                marker=self.config.marker.name,
                output_dir=self.config.paths.output,
            )

            # Get database config for the selected method
            db_config = self.config.taxonomy.get_database_config()

            # Prepare method-specific kwargs
            kwargs = {}
            if self.config.taxonomy.method == "dada2":
                kwargs = {
                    "rdp_db_path": db_config.all,
                    "species_db_path": db_config.species,
                }
            elif self.config.taxonomy.method == "blast":
                kwargs = {
                    "reference_fasta": db_config.fasta,
                    "threshold_species": db_config.threshold_species,
                    "threshold_genus": db_config.threshold_genus,
                    "threshold_family": db_config.threshold_family,
                }
            elif self.config.taxonomy.method == "ecotag":
                kwargs = {
                    "taxonomy_db": db_config.tree,
                    "reference_db": db_config.fasta,
                }
            elif self.config.taxonomy.method == "decipher":
                kwargs = {
                    "trained_classifier_path": db_config.trained,
                    "threshold": db_config.threshold,
                    "processors": db_config.processors,
                }

            # Run assignment
            outputs = assigner.assign_taxonomy(
                query_fasta=query_fasta, asv_count_csv=asv_count_csv, **kwargs
            )

            self.state.complete_step(step_name, outputs)
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs

        except Exception as e:
            self.state.fail_step(step_name, e)
            self._save_state()
            log_pipeline_step(step_name, "error", logger)
            raise

    def run_export(self) -> Dict[str, Path]:
        """
        Run export step (GBIF formatting).

        Returns:
            Dictionary with output paths
        """
        step_name = "export"

        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            logger.info("Running export to GBIF format")

            if not self.config.export.gbif.enabled:
                self.state.skip_step(step_name, reason="GBIF export disabled in config")
                return {}

            # Get taxonomy outputs
            if not self.state.is_step_completed("taxonomy"):
                raise ValueError("Taxonomy step must be completed before export")

            taxo_step = self.state.get_step("taxonomy")
            if taxo_step is None:
                raise ValueError("Taxonomy step not found in pipeline state")
            taxonomy_csv = taxo_step.outputs.get("final_table")

            if taxonomy_csv is None:
                raise ValueError("No taxonomy output file found")

            # Format for GBIF
            formatter = GBIFFormatter()
            output_path = (
                self.config.paths.output
                / f"{self.config.marker.name}_{self.config.taxonomy.method}_gbif.csv"
            )

            # Convert based on method
            if self.config.taxonomy.method == "dada2":
                gbif_table = formatter.from_dada2_rdp(
                    taxonomy_csv,
                    add_rank=self.config.export.gbif.add_rank,
                    add_taxon=self.config.export.gbif.add_taxon,
                )
            elif self.config.taxonomy.method == "ecotag":
                gbif_table = formatter.from_ecotag(
                    taxonomy_csv,
                    add_rank=self.config.export.gbif.add_rank,
                    add_taxon=self.config.export.gbif.add_taxon,
                )
            elif self.config.taxonomy.method == "blast":
                gbif_table = formatter.from_blast(
                    taxonomy_csv,
                    add_rank=self.config.export.gbif.add_rank,
                    add_taxon=self.config.export.gbif.add_taxon,
                )
            elif self.config.taxonomy.method == "decipher":
                gbif_table = formatter.from_decipher(
                    taxonomy_csv,
                    add_rank=self.config.export.gbif.add_rank,
                    add_taxon=self.config.export.gbif.add_taxon,
                )
            else:
                logger.warning(
                    f"GBIF export not supported for {self.config.taxonomy.method}"
                )
                gbif_table = None

            if gbif_table is not None:
                gbif_table.to_csv(output_path, index=False)
                outputs = {"gbif_csv": output_path}
            else:
                outputs = {}

            self.state.complete_step(step_name, outputs)
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs

        except Exception as e:
            self.state.fail_step(step_name, e)
            self._save_state()
            log_pipeline_step(step_name, "error", logger)
            raise

    def run(self, stop_on_error: bool = True) -> PipelineState:
        """
        Run complete pipeline.

        Args:
            stop_on_error: Whether to stop pipeline on first error (default: True)

        Returns:
            Final pipeline state

        Raises:
            Exception: If any step fails and stop_on_error=True
        """
        logger.info("=" * 80)
        logger.info(f"Starting seednap pipeline for marker: {self.config.marker.name}")
        logger.info("=" * 80)

        active_steps = [s for s in self.config.pipeline.steps if s not in self.config.pipeline.skip]
        logger.info(f"Pipeline steps: {' → '.join(active_steps)}")

        # Map step names to methods
        step_methods = {
            "demultiplex": self.run_demultiplex,
            "trim": self.run_trim,
            "dada2": self.run_dada2,
            "taxonomy": self.run_taxonomy,
            "export": self.run_export,
        }

        # Execute steps
        for step_name in active_steps:
            if step_name not in step_methods:
                logger.warning(f"Unknown step: {step_name}, skipping")
                continue

            try:
                logger.info("-" * 80)
                logger.info(f"Step: {step_name.upper()}")
                logger.info("-" * 80)

                step_method = step_methods[step_name]
                step_method()

            except Exception as e:
                logger.error(f"Step '{step_name}' failed: {e}", exc_info=True)
                if stop_on_error:
                    logger.error("Pipeline stopped due to error")
                    raise
                else:
                    logger.warning(f"Continuing pipeline despite error in '{step_name}'")

        # Mark pipeline as complete
        self.state.complete_pipeline()
        self._save_state()

        logger.info("=" * 80)
        logger.info("Pipeline completed successfully!")
        logger.info("=" * 80)

        # Print summary
        summary = self.state.get_summary()
        logger.info(f"Total duration: {summary['total_duration_seconds']:.1f}s")
        logger.info(f"Completed steps: {summary['completed']}/{summary['total_steps']}")
        if summary["failed"] > 0:
            logger.warning(f"Failed steps: {summary['failed']}")

        return self.state

    def _get_sample_list(self) -> List[str]:
        """
        Get list of sample names from raw data directory.

        Returns:
            List of sample names
        """
        import re

        raw_dir = self.config.paths.raw_data
        if not raw_dir.exists():
            raise FileNotFoundError(f"Raw data directory not found: {raw_dir}")

        # Find all R1 files
        r1_files = list(raw_dir.glob("*_R1*.fastq.gz")) + list(
            raw_dir.glob("*_R1*.fastq")
        )

        # Extract sample names (everything before _R1)
        samples = []
        for r1_file in r1_files:
            # Extract sample name using regex
            match = re.match(r"(.+?)_R[12]", r1_file.name)
            if match:
                sample_name = match.group(1)
                if sample_name not in samples:
                    samples.append(sample_name)

        return sorted(samples)

    def _find_read_file(self, sample_name: str, read: str) -> Path:
        """
        Find read file (R1 or R2) for a sample.

        Args:
            sample_name: Sample name
            read: Read number ('R1' or 'R2')

        Returns:
            Path to read file

        Raises:
            FileNotFoundError: If read file not found
        """
        raw_dir = self.config.paths.raw_data

        # Try different file name patterns
        patterns = [
            f"{sample_name}_{read}.fastq.gz",
            f"{sample_name}_{read}.fastq",
            f"{sample_name}_{read}_001.fastq.gz",
            f"{sample_name}_{read}_001.fastq",
        ]

        for pattern in patterns:
            read_file = raw_dir / pattern
            if read_file.exists():
                return read_file

        raise FileNotFoundError(
            f"Could not find {read} file for sample '{sample_name}' in {raw_dir}"
        )
