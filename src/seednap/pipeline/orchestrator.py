"""Pipeline orchestrator for end-to-end workflow execution.

This module orchestrates the complete seednap pipeline:
1. Demultiplexing (optional)
2. Primer trimming (cutadapt)
3. Denoising/clustering (DADA2 or SWARM)
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
from seednap.steps.swarm.processor import SwarmProcessor
from seednap.steps.taxonomic_assignment.assigner import TaxonomicAssigner
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

        # Remember the run-log path so the HTML report can embed the transcript.
        self._log_file = log_file

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
        (base / "02_swarm" / marker).mkdir(parents=True, exist_ok=True)
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

    def _get_trimmed_reads_dir(self) -> Path:
        """Get trimmed reads directory from trim step output or fall back to raw data."""
        if self.state.is_step_completed("trim"):
            trim_step = self.state.get_step("trim")
            trimmed_reads_dir = trim_step.outputs.get("trimmed_dir") if trim_step else None
            if trimmed_reads_dir is None:
                raise ValueError("Trim step completed but no trimmed_dir in outputs")
            return Path(trimmed_reads_dir)
        return self.config.paths.raw_data

    def run_demultiplex(self) -> Dict[str, Path]:
        """
        Run demultiplexing step.

        Returns:
            Dictionary with output paths

        Raises:
            ValueError: If demultiplex is enabled but metadata is missing
        """
        step_name = "demultiplex"

        # D3: explicit skip flag for pre-demultiplexed inputs (one FASTQ per sample
        # already in raw_data). Distinguishes "I don't want to demultiplex" from
        # "demultiplexing isn't applicable to this data".
        if self.config.demultiplex.skip:
            self.state.skip_step(
                step_name, reason="Demultiplexing skipped (raw inputs already demultiplexed)"
            )
            return {}

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

        trimmer = LigationTrimmer(
            cores=self.config.trimming.cores,
            error_rate=self.config.trimming.max_error_rate,
            min_length=self.config.trimming.min_length,
        )
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
            max_sample_failure_rate=self.config.demultiplex.max_sample_failure_rate,
        )

        return {"demux_dir": output_dir, "trimmed_dir": outputs}

    def _run_standard_demux(self) -> Dict[str, Path]:
        """Run standard demultiplexing.

        Note: this protocol is not yet implemented end-to-end. The tag-file
        generation works but the actual cutadapt demultiplex call is missing.
        Use protocol='ligation' or pre-demultiplex your data and set
        demultiplex.skip=true.
        """
        raise NotImplementedError(
            "Standard demultiplex protocol is not yet wired end-to-end. "
            "Use protocol='ligation' for ligation-based libraries, or "
            "set demultiplex.skip=true if your input is already demultiplexed."
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
                min_length=self.config.trimming.min_length,
                overlap=self.config.trimming.overlap,
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

            trimmed_reads_dir = self._get_trimmed_reads_dir()

            processor = Dada2Processor(
                marker=self.config.marker.name,
                trimmed_reads_dir=trimmed_reads_dir,
                output_base_dir=self.config.paths.output,
            )

            outputs = processor.process(
                max_ee=self.config.dada2.filter.max_ee,
                trunc_q=self.config.dada2.filter.trunc_q,
                min_overlap=self.config.dada2.merge.min_overlap,
                max_n=self.config.dada2.filter.max_n,
                rm_phix=self.config.dada2.filter.rm_phix,
                multithread=self.config.dada2.multithread,
                chimera_method=self.config.dada2.chimera.method,
                max_mismatch=self.config.dada2.merge.max_mismatch,
                pool=self.config.dada2.pool,
                min_len=self.config.dada2.filter.min_len,
                max_len=self.config.dada2.filter.max_len,
                collect_metrics=self.config.metrics.generate_plots,
            )

            self.state.complete_step(step_name, outputs)
            self._build_read_tracking_report("dada2")
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs

        except Exception as e:
            self.state.fail_step(step_name, e)
            self._save_state()
            log_pipeline_step(step_name, "error", logger)
            raise

    def run_swarm(self) -> Dict[str, Path]:
        """
        Run SWARM OTU clustering step.

        Returns:
            Dictionary with output paths
        """
        step_name = "swarm"

        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            logger.info("Running SWARM OTU clustering")

            trimmed_reads_dir = self._get_trimmed_reads_dir()

            processor = SwarmProcessor(
                marker=self.config.marker.name,
                trimmed_reads_dir=trimmed_reads_dir,
                output_base_dir=self.config.paths.output,
            )

            outputs = processor.process(
                d=self.config.swarm.clustering.d,
                fastidious=self.config.swarm.clustering.fastidious,
                boundary=self.config.swarm.clustering.boundary,
                threads=self.config.swarm.clustering.threads,
                fastq_maxdiffs=self.config.swarm.merge.fastq_maxdiffs,
                fastq_minovlen=self.config.swarm.merge.fastq_minovlen,
                allow_stagger=self.config.swarm.merge.allow_stagger,
                min_sequence_length=self.config.swarm.min_sequence_length,
                chimera_detection=self.config.swarm.chimera.method != "none",
            )

            self.state.complete_step(step_name, outputs)
            self._build_read_tracking_report("swarm")
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs

        except Exception as e:
            self.state.fail_step(step_name, e)
            self._save_state()
            log_pipeline_step(step_name, "error", logger)
            raise

    def _build_read_tracking_report(self, method: str) -> None:
        """Build the read-tracking table (+ optional HTML report) after a
        clustering step.

        Non-fatal: a reporting failure logs a ``[WARN]`` and never fails the
        run (CLAUDE.md section 4 -- the report is observational only).
        """
        if not self.config.report.read_tracking:
            return
        try:
            import pandas as pd

            from seednap.steps.report import ReadTrackingBuilder

            marker = self.config.marker.name
            out = self.config.paths.output
            report_dir = out / "04_report" / marker
            kwargs = {
                "marker": marker,
                "logs_dir": out / "logs",
                "warn_below_retention_pct": self.config.report.warn_below_retention_pct,
                "warn_step_loss_pct": self.config.report.warn_step_loss_pct,
            }
            if method == "dada2":
                kwargs["dada2_dir"] = out / "02_dada2" / marker
            elif method == "swarm":
                kwargs["swarm_otu_table"] = out / "02_swarm" / marker / "otu_table.csv"

            builder = ReadTrackingBuilder(**kwargs)
            df = builder.build()
            builder.write(report_dir, df=df)
            warns = builder.warnings(df)

            # Persist a compact summary into the step state (resume-safe).
            step = self.state.get_step(method)
            if step is not None:
                if not df.empty:
                    raw = pd.to_numeric(df["raw"], errors="coerce")
                    final = pd.to_numeric(df[builder.steps[-1]], errors="coerce")
                    pr = pd.to_numeric(df["pct_retained"], errors="coerce")
                    step.metadata["read_tracking"] = {
                        "n_samples": int(len(df)),
                        "raw_reads_total": int(raw.sum()) if raw.notna().any() else None,
                        "final_step": builder.steps[-1],
                        "final_reads_total": int(final.sum()) if final.notna().any() else None,
                        "mean_retention_pct": round(float(pr.mean()), 2) if pr.notna().any() else None,
                        "n_warnings": len(warns),
                    }
                else:
                    step.metadata["read_tracking"] = {"n_samples": 0, "n_warnings": len(warns)}
        except Exception as exc:  # noqa: BLE001 -- reporting must never fail the run
            logger.warning(
                f"[WARN] read_tracking report: expected=read-tracking table for "
                f"'{method}', got=error ({exc}), fallback=skipped (pipeline unaffected)",
            )

    def _build_html_report(self) -> None:
        """Build the full HTML run report at run end (needs taxonomy + clustering).

        Opt-in (``report.html_report``). Non-fatal: failures log a ``[WARN]`` and
        never affect the run (CLAUDE.md section 4).
        """
        if not (self.config.report.read_tracking and self.config.report.html_report):
            return
        try:
            from seednap.steps.report import HTMLReportBuilder, ReadTrackingBuilder

            marker = self.config.marker.name
            out = self.config.paths.output
            steps = set(self.config.pipeline.steps)
            kwargs = {
                "marker": marker, "logs_dir": out / "logs",
                "warn_below_retention_pct": self.config.report.warn_below_retention_pct,
                "warn_step_loss_pct": self.config.report.warn_step_loss_pct,
            }
            otu_full = None
            if "dada2" in steps:
                kwargs["dada2_dir"] = out / "02_dada2" / marker
            elif "swarm" in steps:
                kwargs["swarm_otu_table"] = out / "02_swarm" / marker / "otu_table.csv"
                otu_full = out / "02_swarm" / marker / "otu_table_full.csv"

            builder = ReadTrackingBuilder(**kwargs)
            df = builder.build()
            warns = builder.warnings(df, log=False)

            taxo_csv = None
            tstep = self.state.get_step("taxonomy")
            if tstep is not None:
                taxo_csv = tstep.outputs.get("final_table") or tstep.outputs.get("taxonomy_csv")

            # Provenance from the pipeline config (always available).
            reference_db = None
            try:
                db = self.config.taxonomy.get_database_config()
                reference_db = str(getattr(db, "fasta", None) or getattr(db, "trained", None) or "") or None
            except Exception:  # noqa: BLE001 -- provenance is best-effort, never fatal
                reference_db = None
            provenance = {
                "dataset_name": marker,
                "marker": marker,
                "primer_fwd": self.config.marker.primers.forward,
                "primer_rev": self.config.marker.primers.reverse,
                "raw_data": str(self.config.paths.raw_data),
                "reference_db": reference_db,
            }

            html_path = HTMLReportBuilder(
                marker, df, warnings=warns, steps=builder.steps,
                state=self.state.model_dump(mode="json"),
                taxonomy_csv=taxo_csv,
                otu_table_full=otu_full,
                field_metadata_csv=self.config.report.sample_metadata,
                project_metadata_csv=self.config.report.project_metadata,
                log_file=getattr(self, "_log_file", None),
                summary={
                    "warn_below_retention_pct": self.config.report.warn_below_retention_pct,
                    "subtitle": f"{len(df)} samples · marker {marker}",
                    "provenance": provenance,
                },
            ).write(out / "04_report" / marker / "report.html")
            logger.info(f"HTML run report: {html_path}")
        except Exception as exc:  # noqa: BLE001 -- reporting must never fail the run
            logger.warning(
                f"[WARN] html_report: expected=HTML report generation, "
                f"got=error ({exc}), fallback=skipped (run unaffected)",
            )

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

            # Get outputs from clustering step (DADA2 or SWARM)
            clustering_step = None
            for clustering_name in ("dada2", "swarm"):
                if self.state.is_step_completed(clustering_name):
                    clustering_step = self.state.get_step(clustering_name)
                    break

            if clustering_step is None:
                raise ValueError("DADA2 or SWARM step must be completed before taxonomy")

            query_fasta = clustering_step.outputs.get("query_fasta")
            asv_count_csv = clustering_step.outputs.get("seqtab_clean_t")
            if query_fasta is None or asv_count_csv is None:
                raise ValueError(
                    f"Clustering outputs incomplete: query_fasta={query_fasta}, "
                    f"seqtab_clean_t={asv_count_csv}"
                )
            query_fasta = Path(query_fasta)
            asv_count_csv = Path(asv_count_csv)

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
                    "multithread": self.config.dada2.multithread,
                    "bootstrap_threshold": db_config.bootstrap_threshold,
                    "contaminants": self.config.taxonomy.contaminants,
                }
            elif self.config.taxonomy.method == "blast":
                kwargs = {
                    "reference_fasta": db_config.fasta,
                    "threshold_species": db_config.threshold_species,
                    "threshold_genus": db_config.threshold_genus,
                    "threshold_family": db_config.threshold_family,
                    "threshold_order": db_config.threshold_order,
                    "threshold_class": db_config.threshold_class,
                    "top_bitscore_pct": db_config.top_bitscore_pct,
                    "contaminants": self.config.taxonomy.contaminants,
                    "perc_identity": db_config.perc_identity,
                    "qcov_hsp_perc": db_config.qcov_hsp_perc,
                    "evalue": db_config.evalue,
                    "max_target_seqs": db_config.max_target_seqs,
                    "task": db_config.task,
                }
            elif self.config.taxonomy.method == "ecotag":
                kwargs = {
                    "taxonomy_db": db_config.tree,
                    "reference_db": db_config.fasta,
                    "contaminants": self.config.taxonomy.contaminants,
                }
            elif self.config.taxonomy.method == "decipher":
                kwargs = {
                    "trained_classifier_path": db_config.trained,
                    "threshold": db_config.threshold,
                    "processors": db_config.processors,
                    "contaminants": self.config.taxonomy.contaminants,
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
            taxonomy_csv = Path(taxonomy_csv)

            # Format for GBIF
            formatter = GBIFFormatter()
            output_path = (
                self.config.paths.output
                / f"{self.config.marker.name}_{self.config.taxonomy.method}_gbif.csv"
            )

            gbif_table = formatter.from_method(
                method=self.config.taxonomy.method,
                input_path=taxonomy_csv,
                add_rank=self.config.export.gbif.add_rank,
                add_taxon=self.config.export.gbif.add_taxon,
            )

            gbif_table.to_csv(output_path, index=False)
            outputs = {"gbif_csv": output_path}

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
            "swarm": self.run_swarm,
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

        # Full HTML run report (opt-in) -- after all steps so taxonomy is available.
        self._build_html_report()

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

        # Find all R1 files (support both _R1 and .R1 naming)
        r1_patterns = [
            "*_R1*.fastq.gz", "*_R1*.fastq",
            "*.R1.fastq.gz", "*.R1.fastq",
        ]
        r1_files = []
        for pattern in r1_patterns:
            r1_files.extend(raw_dir.glob(pattern))

        # Extract sample names (everything before _R1 or .R1)
        samples = []
        for r1_file in sorted(r1_files):
            match = re.match(r"(.+?)[._]R[12]", r1_file.name)
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

        # Try different file name patterns (support both _R1 and .R1 naming)
        patterns = [
            f"{sample_name}_{read}.fastq.gz",
            f"{sample_name}_{read}.fastq",
            f"{sample_name}.{read}.fastq.gz",
            f"{sample_name}.{read}.fastq",
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
