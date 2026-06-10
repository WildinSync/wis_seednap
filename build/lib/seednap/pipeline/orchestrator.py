"""Pipeline orchestrator for end-to-end workflow execution.

This module orchestrates the complete seednap pipeline:
1. Demultiplexing (optional)
2. Primer trimming (cutadapt)
3. Denoising/clustering (DADA2 or SWARM)
4. Taxonomic assignment
5. Export to GBIF format
"""

from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union

import yaml

from seednap.__version__ import __version__ as SEEDNAP_VERSION
from seednap.config.loader import load_config
from seednap.config.models import PipelineConfig
from seednap.config.models.operational import VALID_STEPS
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
        self.config_path: Optional[Path]
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
                    f"Cannot resume: no pipeline state file at {self.state_file}. "
                    f"--resume re-reads the state JSON that a previous run of this "
                    f"marker wrote, and none exists at that path. Likely causes: this "
                    f"marker was never run, the output directory was cleared, or "
                    f"paths.output in the config points somewhere different than the "
                    f"original run. Fix: drop --resume to start a fresh run, or pass "
                    f"--state-file pointing at the existing state JSON (default "
                    f"location is <paths.output>/.<marker>_state.json)."
                )
            self.state = PipelineState.load(self.state_file)
            logger.info(f"Resuming pipeline from {self.state_file}")
            logger.info(
                f"Completed steps: {', '.join(self.state.get_completed_steps())}"
            )
            # Reproducibility: the state was written by some seednap version; if the
            # running version differs, the resumed steps may not match what already ran.
            if self.state.seednap_version != SEEDNAP_VERSION:
                logger.warning(
                    f"[WARN] resume version mismatch: "
                    f"expected={self.state.seednap_version} (version that wrote this "
                    f"run's state), got={SEEDNAP_VERSION} (running version), "
                    f"fallback=resuming anyway. Outputs from already-completed steps "
                    f"were produced by the stored version."
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

        # Build the step dispatch once and assert it covers every valid stage, so the
        # dispatch table cannot silently drift from operational.VALID_STEPS.
        self._step_methods: Dict[str, Callable[[], Dict[str, Any]]] = {
            "demultiplex": self.run_demultiplex,
            "trim": self.run_trim,
            "dada2": self.run_dada2,
            "swarm": self.run_swarm,
            "taxonomy": self.run_taxonomy,
            "clean": self.run_clean,
            "export": self.run_export,
            "report": self.run_report,
        }
        missing = set(VALID_STEPS) - set(self._step_methods)
        if missing:
            raise RuntimeError(
                f"Step dispatch is missing handlers for valid stage(s) "
                f"{sorted(missing)}; operational.VALID_STEPS and the orchestrator "
                f"dispatch have drifted. Add a run_* handler for each missing stage."
            )

        # Reproducibility: snapshot the effective merged config into the output tree and
        # record its path in the state JSON, so a run is reconstructable from its outputs.
        self._write_config_snapshot()

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

    def _write_config_snapshot(self) -> None:
        """Write the effective merged config (YAML) into the output tree and record
        its path in the state.

        Reproducibility: a run must be reconstructable from its outputs, so we persist
        the fully-merged config (defaults + the marker YAML) rather than relying on the
        original YAML still existing unchanged. The snapshot lives next to the state
        file at ``<paths.output>/.<marker>_config.snapshot.yaml`` and its path is stored
        on the state so the state JSON references it.
        """
        snapshot_path = (
            self.config.paths.output
            / f".{self.config.marker.name}_config.snapshot.yaml"
        )
        snapshot_path.parent.mkdir(parents=True, exist_ok=True)
        # mode="json" serializes Path/enum values to plain strings so safe_dump works.
        config_dump = self.config.model_dump(mode="json")
        with open(snapshot_path, "w", encoding="utf-8") as f:
            yaml.safe_dump(config_dump, f, default_flow_style=False, sort_keys=False)
        self.state.config_snapshot_path = snapshot_path
        self._save_state()
        logger.info(f"Wrote effective config snapshot to {snapshot_path}")

    def _should_run_step(self, step_name: str) -> bool:
        """
        Determine if a step should be run.

        A completed step is skipped. A previously-failed step is retried. Any
        other status (pending, or a step left RUNNING because a prior run was
        killed mid-step before fail_step could record the failure) falls through
        to True and is run/re-run. So an interrupted step re-runs on --resume.

        Args:
            step_name: Name of the step

        Returns:
            True if step should be run, False if it should be skipped
        """
        # Check if step is already completed
        if self.state.is_step_completed(step_name):
            logger.info(f"Step '{step_name}' already completed, skipping")
            return False

        # Check if step failed previously
        if self.state.is_step_failed(step_name):
            logger.warning(f"Step '{step_name}' failed previously, retrying")
            return True

        return True

    def _execute_step(
        self,
        step_name: str,
        body: Callable[[], Dict[str, Any]],
        post_complete: Optional[Callable[[], None]] = None,
    ) -> Dict[str, Any]:
        """Run one pipeline step's body with the shared state/logging scaffold.

        Behavior-preserving extraction of the boilerplate every uniform run_* method
        repeated: should-run guard (return cached outputs when skipped), start logging
        and state, try the body, on success record outputs/complete/save/log, on failure
        record fail/save/log and re-raise.

        Args:
            step_name: Pipeline stage name (also the state key).
            body: Callable returning the step's outputs dict; the only step-specific work.
            post_complete: Optional hook run AFTER complete_step but BEFORE _save_state,
                for work that must observe the completed step's recorded outputs (e.g.
                dada2/swarm manifest cross-checks). Preserves the original call ordering.

        Returns:
            The step's outputs dict (or the previously recorded outputs when skipped).
        """
        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            outputs = body()
            self.state.complete_step(step_name, outputs)
            if post_complete is not None:
                post_complete()
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs
        except Exception as e:
            self.state.fail_step(step_name, e)
            self._save_state()
            log_pipeline_step(step_name, "error", logger)
            raise

    def _get_trimmed_reads_dir(self) -> Path:
        """Get the trimmed-reads directory that dada2/swarm consume.

        Reads only the ``trim`` step's ``trimmed_dir`` output. Config validation
        forces ``trim`` to run before dada2/swarm, so that step is always present;
        the fall-back to raw_data is therefore only reached if trim is absent.

        Note: this deliberately does NOT read the demultiplex step's ``trimmed_dir``
        output. The ligation-demux path returns a ``trimmed_dir`` (see
        _run_ligation_demux), but nothing consumes it: the ``trim`` step re-runs over
        raw_data and produces the directory used here. The demux ``trimmed_dir`` is
        effectively unused in the dada2/swarm data flow.
        """
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
        # This step runs iff "demultiplex" is listed in pipeline.steps. If your raw inputs are
        # already demultiplexed (one FASTQ per sample), simply omit "demultiplex" from steps.
        def body() -> Dict[str, Any]:
            logger.info(f"Running demultiplexing: {self.config.demultiplex.protocol}")

            if self.config.demultiplex.protocol == "ligation":
                return self._run_ligation_demux()
            elif self.config.demultiplex.protocol == "standard":
                return self._run_standard_demux()
            else:
                raise ValueError(
                    f"Unknown demultiplex protocol: {self.config.demultiplex.protocol}"
                )

        return self._execute_step("demultiplex", body)

    def _run_ligation_demux(self) -> Dict[str, Path]:
        """Run ligation-based demultiplexing."""
        if self.config.demultiplex.metadata is None:
            raise ValueError(
                "Ligation demultiplexing requires a sample-tag metadata CSV, but "
                "demultiplex.metadata is not set in the config. The ligation protocol "
                "splits one multiplexed library FASTQ into per-sample files using the "
                "tag-to-sample mapping in that CSV (columns: eventID/sample, "
                "tag_demultiplex, library). Fix: add "
                "`demultiplex.metadata: /path/to/metadata.csv` to the marker YAML, or, "
                "if your reads are already demultiplexed (one FASTQ pair per sample), "
                "remove 'demultiplex' from pipeline.steps."
            )

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
        """Run standard (tag-based) demultiplexing.

        Note: this protocol is not yet implemented end-to-end. Use
        protocol='ligation', or pre-demultiplex your data and omit
        "demultiplex" from pipeline.steps.

        Raises:
            NotImplementedError: Always; the standard protocol is not wired in.
        """
        raise NotImplementedError(
            "Standard demultiplex protocol is not yet wired end-to-end. "
            "Use protocol='ligation' for ligation-based libraries, or "
            "omit 'demultiplex' from pipeline.steps if your input is already demultiplexed."
        )

        # TODO: Implement standard demultiplexing workflow
        # This would use CutadaptRunner.demultiplex_by_tags()

    def run_trim(self) -> Dict[str, Any]:
        """
        Run primer trimming step.

        Returns:
            Dictionary with the trimmed-reads directory and the per-sample
            trim outputs (under the ``"samples"`` key).
        """
        def body() -> Dict[str, Any]:
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

            trimmed_outputs: Dict[str, tuple] = {}
            for sample_name in samples:
                logger.info(f"Trimming sample: {sample_name}")

                r1_input = self._find_read_file(sample_name, "R1")
                r2_input = self._find_read_file(sample_name, "R2")

                sample_outputs = trimmer.trim_sample(
                    r1_input=r1_input,
                    r2_input=r2_input,
                    output_dir=output_dir,
                    sample_name=sample_name,
                    forward_primer=self.config.marker.primers.forward,
                    reverse_primer=self.config.marker.primers.reverse,
                    keep_untrimmed=False,
                    discard_untrimmed=self.config.trimming.discard_untrimmed,
                )

                trimmed_outputs[sample_name] = sample_outputs

            return {
                "trimmed_dir": output_dir,
                "samples": trimmed_outputs,
            }

        return self._execute_step("trim", body)

    def run_dada2(self) -> Dict[str, Path]:
        """
        Run DADA2 processing step.

        Returns:
            Dictionary with output paths
        """
        def body() -> Dict[str, Any]:
            logger.info("Running DADA2 processing")

            trimmed_reads_dir = self._get_trimmed_reads_dir()

            processor = Dada2Processor(
                marker=self.config.marker.name,
                trimmed_reads_dir=trimmed_reads_dir,
                output_base_dir=self.config.paths.output,
            )

            return processor.process(
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
                library_map=self._build_library_map(),
                collect_metrics=self.config.dada2.collect_metrics,
            )

        return self._execute_step(
            "dada2",
            body,
            post_complete=lambda: self._validate_manifest_against_abundance("dada2"),
        )

    def run_swarm(self) -> Dict[str, Path]:
        """
        Run SWARM OTU clustering step.

        Returns:
            Dictionary with output paths
        """
        def body() -> Dict[str, Any]:
            logger.info("Running SWARM OTU clustering")

            trimmed_reads_dir = self._get_trimmed_reads_dir()

            processor = SwarmProcessor(
                marker=self.config.marker.name,
                trimmed_reads_dir=trimmed_reads_dir,
                output_base_dir=self.config.paths.output,
            )

            return processor.process(
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

        return self._execute_step(
            "swarm",
            body,
            post_complete=lambda: self._validate_manifest_against_abundance("swarm"),
        )

    def _build_library_map(self) -> Optional[Path]:
        """Write a ``sample,library`` CSV for DADA2-by-library, derived from the manifest's
        seq_run_id grouping.

        Grouping source precedence: report.sample_metadata (field CSV) is preferred over
        demultiplex.metadata (lab CSV) as the primary manifest source (``src``). The lab CSV
        is passed as the manifest's extra ``lab_csv`` only when it is distinct from ``src``
        (i.e. when the field CSV was chosen as primary); if both point at the same file, no
        extra lab CSV is supplied.

        Returns the CSV path, or None when per_library is off or no grouping source exists
        (the R script then runs the standard single-batch path). A single-library grouping is
        a no-op there too, so writing the map is always safe.
        """
        if not self.config.dada2.per_library:
            return None
        field_csv = self.config.report.sample_metadata
        lab_csv = self.config.demultiplex.metadata
        src = field_csv or lab_csv
        if src is None:
            logger.warning(
                "[WARN] dada2 per_library: expected=report.sample_metadata or "
                "demultiplex.metadata for the library grouping, got=none, "
                "fallback=standard single-batch DADA2"
            )
            return None
        try:
            import pandas as pd

            from seednap.config.manifest_migrate import migrate_to_manifest

            extra_lab = Path(lab_csv) if (lab_csv and str(lab_csv) != str(src)) else None
            manifest = migrate_to_manifest(
                Path(src),
                lab_csv=extra_lab,
                project_csv=self.config.report.project_metadata,
                target_gene=self.config.marker.name,
            )
            df = pd.DataFrame(
                [{"sample": r.eventID, "library": r.seq_run_id} for r in manifest.rows]
            )
            out = self.config.paths.output / "02_dada2" / self.config.marker.name / "library_map.csv"
            out.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(out, index=False)
            logger.info(
                f"DADA2-by-library: wrote library map ({len(df)} samples, "
                f"{df['library'].nunique()} libraries) to {out}"
            )
            return out
        except Exception as exc:  # noqa: BLE001 -- never fail the run building a helper file
            logger.warning(
                f"[WARN] dada2 per_library: could not build the library map ({exc}), "
                f"fallback=standard single-batch DADA2"
            )
            return None

    def _report_dir(self) -> Path:
        """Per-marker directory for report artifacts.

        Honors ``report.output_dir`` (a per-marker subdirectory is created
        inside it); otherwise defaults to ``<paths.output>/04_report/<marker>``.
        """
        marker = self.config.marker.name
        base = self.config.report.output_dir
        if base is not None:
            return Path(base) / marker
        return self.config.paths.output / "04_report" / marker

    def _build_read_tracking_report(self, method: str) -> None:
        """Build the read-tracking table (+ optional HTML report) after a
        clustering step.

        Non-fatal: a reporting failure logs a ``[WARN]`` and never fails the
        run (the no-silent-fallbacks policy -- the report is observational only).
        """
        try:
            import pandas as pd

            from seednap.steps.report import ReadTrackingBuilder

            marker = self.config.marker.name
            out = self.config.paths.output
            report_dir = self._report_dir()
            kwargs: Dict[str, Any] = {
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
            # Run-level step summary: total reads + ASV/OTU count after each step.
            builder.write_step_summary(report_dir, summary_df=builder.step_summary(df))
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
                    # E3: per-sample counts keyed on eventID -- the substrate the FAIRe
                    # manifest's reads_* columns consume. JSON-safe (NA -> None, numpy -> int),
                    # serialized into the state JSON so it survives --resume.
                    per_sample: Dict[str, Dict[str, Optional[int]]] = {}
                    for _, r in df.iterrows():
                        rec: Dict[str, Optional[int]] = {}
                        for s in builder.steps:
                            v = pd.to_numeric(pd.Series([r[s]]), errors="coerce").iloc[0]
                            rec[s] = int(v) if pd.notna(v) else None
                        per_sample[str(r["sample"])] = rec
                    step.metadata["read_tracking_per_sample"] = per_sample
                else:
                    step.metadata["read_tracking"] = {"n_samples": 0, "n_warnings": len(warns)}
        except Exception as exc:  # noqa: BLE001 -- reporting must never fail the run
            logger.warning(
                f"[WARN] read_tracking report: expected=read-tracking table for "
                f"'{method}', got=error ({exc}), fallback=skipped (pipeline unaffected)",
            )

    def _validate_manifest_against_abundance(self, method: str) -> None:
        """Cross-check the FAIRe manifest's eventIDs against the abundance table.

        When per-sample (field) metadata is configured (``report.sample_metadata``), derive
        a manifest from it and assert its eventID set matches the abundance table's sample
        columns -- the up-front silent-ID-mismatch guard (the no-silent-fallbacks policy), catching e.g. an
        unlabelled ``Blank-PCR-3`` column. Warn-only and non-fatal: it never alters or fails
        the run.
        """
        field_csv = self.config.report.sample_metadata
        if field_csv is None:
            return
        try:
            from seednap.config.manifest import validate_against_abundance
            from seednap.config.manifest_migrate import migrate_to_manifest

            marker = self.config.marker.name
            out = self.config.paths.output
            # Sample columns live in the SWARM otu_table (sequences x samples) or the
            # transposed DADA2 seqtab (seqtab_clean_t.csv = ASVs x samples).
            if method == "dada2":
                abundance = out / "02_dada2" / marker / "seqtab_clean_t.csv"
            else:
                abundance = out / "02_swarm" / marker / "otu_table.csv"
            if not Path(field_csv).exists() or not abundance.exists():
                return
            manifest = migrate_to_manifest(
                Path(field_csv),
                project_csv=self.config.report.project_metadata,
                target_gene=marker,
            )
            validate_against_abundance(manifest, abundance)
        except Exception as exc:  # noqa: BLE001 -- validation must never fail the run
            logger.warning(
                f"[WARN] manifest validation: expected=cross-CSV eventID check for "
                f"'{method}', got=error ({exc}), fallback=skipped (pipeline unaffected)",
            )

    def _build_html_report(self) -> None:
        """Build the full HTML run report (needs taxonomy + clustering).

        Called by the ``report`` step; gated by ``report.html_report`` (default on, set
        ``false`` to write only the read-tracking table). Non-fatal: failures log a ``[WARN]``
        and never affect the run (the no-silent-fallbacks policy).
        """
        if not self.config.report.html_report:
            return
        try:
            from seednap.steps.report import HTMLReportBuilder, ReadTrackingBuilder

            marker = self.config.marker.name
            out = self.config.paths.output
            steps = set(self.config.pipeline.steps)
            kwargs: Dict[str, Any] = {
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
            step_summary_df = builder.step_summary(df)

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
                "seednap_version": self.state.seednap_version,
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
                step_summary_df=step_summary_df,
                summary={
                    "warn_below_retention_pct": self.config.report.warn_below_retention_pct,
                    "subtitle": f"{len(df)} samples · marker {marker}",
                    "provenance": provenance,
                },
            ).write(self._report_dir() / "report.html")
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
        def body() -> Dict[str, Any]:
            logger.info(f"Running taxonomic assignment: {self.config.taxonomy.method}")

            # Get outputs from clustering step (DADA2 or SWARM)
            clustering_step = None
            for clustering_name in ("dada2", "swarm"):
                if self.state.is_step_completed(clustering_name):
                    clustering_step = self.state.get_step(clustering_name)
                    break

            if clustering_step is None:
                raise ValueError(
                    "Cannot run taxonomy: it needs ASVs/OTUs from a completed feature "
                    "step, but neither 'dada2' nor 'swarm' is marked completed in this "
                    "run's state. (pipeline.steps ordering is validated at config load, "
                    "so the feature step IS configured; it simply did not finish.) This "
                    "usually means the earlier dada2/swarm step failed or was "
                    "interrupted, and you reached taxonomy via --continue-on-error or "
                    "--resume. Fix: open <paths.output>/.<marker>_state.json and check "
                    "the 'dada2' (or 'swarm') step; if status is 'failed' or missing, "
                    "re-run that step to completion first (re-run with --resume, or "
                    "without --continue-on-error so a feature-step failure stops the "
                    "pipeline) before running taxonomy."
                )

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
                    "lca_pident_delta": db_config.lca_pident_delta,
                    "lca_algorithm": db_config.lca_algorithm,
                    "lca_pid": db_config.lca_pid,
                    "lca_diff": db_config.lca_diff,
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
            return assigner.assign_taxonomy(
                query_fasta=query_fasta, asv_count_csv=asv_count_csv, **kwargs
            )

        return self._execute_step("taxonomy", body)

    def run_clean(self) -> Dict[str, Path]:
        """Control-decontamination step: clean the taxonomy table against its negative
        controls (list ``clean`` in pipeline.steps, after a feature step). Needs
        ``report.sample_metadata`` (the control-identity source). Non-fatal: a problem skips
        cleaning with a ``[WARN]`` rather than failing the run (export falls back to the
        uncleaned table)."""
        step_name = "clean"
        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        field_csv = self.config.report.sample_metadata
        if field_csv is None:
            self.state.skip_step(
                step_name,
                reason="Cleaning needs report.sample_metadata (control identity); skipped",
            )
            logger.warning(
                "[WARN] cleaning: expected=report.sample_metadata for control identity, "
                "got=none, fallback=cleaning skipped (export uses the uncleaned table)"
            )
            return {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()
        try:
            import pandas as pd

            from seednap.config.manifest import classify_control
            from seednap.config.manifest_migrate import migrate_to_manifest
            from seednap.steps.cleaning import CleaningProcessor

            taxo_step = self.state.get_step("taxonomy")
            taxonomy_csv = taxo_step.outputs.get("final_table") if taxo_step else None
            if taxonomy_csv is None:
                raise ValueError("Cleaning requires a completed taxonomy step with a final_table")
            taxonomy_csv = Path(taxonomy_csv)

            manifest = migrate_to_manifest(
                Path(field_csv),
                project_csv=self.config.report.project_metadata,
                target_gene=self.config.marker.name,
            )
            df = pd.read_csv(taxonomy_csv)
            id_col = "ASV_ID" if "ASV_ID" in df.columns else str(df.columns[0])
            manifest_ids = set(manifest.event_ids())
            sample_cols = [
                c for c in df.columns
                if c != id_col and (c in manifest_ids or classify_control(str(c)).is_control)
            ]

            cleaned, report, result = CleaningProcessor(mode=self.config.cleaning.mode).clean(
                df, manifest, id_col=id_col, sample_cols=sample_cols
            )

            cleaned_path = (
                self.config.paths.output
                / f"{self.config.marker.name}_{self.config.taxonomy.method}_cleaned.csv"
            )
            cleaned.to_csv(cleaned_path, index=False)
            report_dir = self._report_dir()
            report_dir.mkdir(parents=True, exist_ok=True)
            report_path = report_dir / "cleaning_report.csv"
            report.to_csv(report_path, index=False)

            outputs = {"cleaned_table": cleaned_path, "cleaning_report": report_path}
            step = self.state.get_step(step_name)
            if step is not None:
                step.metadata["cleaning"] = result.model_dump()
            self.state.complete_step(step_name, outputs)
            self._save_state()
            log_pipeline_step(step_name, "complete", logger)
            return outputs
        except Exception as e:
            # Cleaning is observational/optional; never fail the run over it.
            self.state.skip_step(step_name, reason=f"Cleaning error: {e}")
            self._save_state()
            logger.warning(
                f"[WARN] cleaning: expected=cleaned table, got=error ({e}), "
                f"fallback=skipped (export uses the uncleaned table)"
            )
            return {}

    def _warn_if_export_predates_clean(self, export_step: Optional[Any]) -> None:
        """Warn when an already-completed export predates a (re-)completed clean step.

        On --resume a clean step that was SKIPPED in run 1 (transient error) can
        re-run and now COMPLETE, writing a fresh cleaned table. But if export was
        already COMPLETED against the uncleaned table, _should_run_step('export')
        returns False and export is not re-run, so the GBIF CSV silently stays
        stale. Surface this per the no-silent-fallbacks policy; the user must
        re-run export to pick up the cleaned table.
        """
        if export_step is None or export_step.completed_at is None:
            return
        clean_step = self.state.get_step("clean")
        if (
            clean_step is not None
            and clean_step.completed_at is not None
            and clean_step.outputs.get("cleaned_table")
            and clean_step.completed_at > export_step.completed_at
        ):
            logger.warning(
                "[WARN] export: expected=GBIF CSV reflecting the decontaminated "
                f"table (clean completed {clean_step.completed_at.isoformat()}), "
                f"got=stale export completed earlier ({export_step.completed_at.isoformat()}) "
                "against the uncleaned table, fallback=existing export left as-is. "
                "Re-run the 'export' step (e.g. delete the 'export' entry from the "
                "state JSON, or use the standalone export command) so the GBIF CSV "
                "reflects the cleaned table."
            )

    def run_export(self) -> Dict[str, Path]:
        """
        Run export step (GBIF formatting).

        Returns:
            Dictionary with output paths
        """
        step_name = "export"

        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            # Export-staleness guard on --resume after a clean retry: if the clean
            # step completed AFTER this already-completed export, the GBIF CSV still
            # reflects the pre-clean (uncleaned) table and is silently stale.
            self._warn_if_export_predates_clean(step)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            logger.info("Running export to GBIF format")

            # Get taxonomy outputs
            if not self.state.is_step_completed("taxonomy"):
                raise ValueError(
                    "Cannot run export: GBIF formatting needs the assigned-taxonomy "
                    "table, but the 'taxonomy' step did not complete in this run (it "
                    "failed earlier and the run continued past it under "
                    "--continue-on-error). Fix the taxonomy failure first: read its "
                    "error in the run log or in the 'taxonomy' step's status/error "
                    "fields in <paths.output>/.<marker>_state.json, resolve the cause, "
                    "then re-run (use --resume to retry from the failed step). Note: a "
                    "missing or mis-ordered 'taxonomy' stage is not the cause here -- "
                    "pipeline.steps ordering is validated at config load and would have "
                    "been rejected before the run started."
                )

            taxo_step = self.state.get_step("taxonomy")
            if taxo_step is None:
                raise ValueError("Taxonomy step not found in pipeline state")
            taxonomy_csv = taxo_step.outputs.get("final_table")

            if taxonomy_csv is None:
                raise ValueError(
                    "Export cannot start: the completed taxonomy step recorded no "
                    "'final_table' output, so there is no merged taxonomy+abundance CSV "
                    "to format for GBIF. In a normal single-version run every method "
                    "(blast/dada2/ecotag/decipher) writes final_table, so the usual "
                    "cause is resuming export against a state JSON "
                    "(<paths.output>/.<marker>_state.json) written by an older seednap "
                    "that used a different output key. Fix: look for the merged table at "
                    "<paths.output>/<marker>_<method>.csv (e.g. <marker>_dada2RDP.csv, "
                    "<marker>_blast.csv) and the taxonomy log; if it is missing or the "
                    "state is stale, re-run the 'taxonomy' step (or delete the "
                    "taxonomy/export entries from the state JSON and re-run) so it "
                    "regenerates final_table before export."
                )
            taxonomy_csv = Path(taxonomy_csv)

            # Prefer the decontaminated table when the cleaning step produced one.
            clean_step = self.state.get_step("clean")
            if clean_step is not None and clean_step.outputs.get("cleaned_table"):
                cleaned = Path(clean_step.outputs["cleaned_table"])
                if cleaned.exists():
                    logger.info(f"Export using cleaned taxonomy table: {cleaned}")
                    taxonomy_csv = cleaned

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

    def run_report(self) -> Dict[str, Path]:
        """Reporting step: write the per-step read/sequence tracking table + step summary, and
        (when ``report.html_report`` is on) the self-contained HTML run report. Observational:
        a reporting failure logs a ``[WARN]`` and never fails the run (the no-silent-fallbacks policy)."""
        step_name = "report"
        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        # The read-tracking table is keyed off the feature path that produced the table.
        method = (
            "dada2" if self.state.is_step_completed("dada2")
            else "swarm" if self.state.is_step_completed("swarm")
            else None
        )
        if method is None:
            self.state.skip_step(
                step_name, reason="No completed dada2/swarm step to report on"
            )
            logger.warning(
                "[WARN] report: expected=a completed feature step (dada2 or swarm), "
                "got=none, fallback=report skipped"
            )
            self._save_state()
            return {}

        self._build_read_tracking_report(method)
        self._build_html_report()
        self.state.complete_step(step_name, {})
        self._save_state()
        log_pipeline_step(step_name, "complete", logger)
        return {}

    def run(self, stop_on_error: bool = True) -> PipelineState:
        """
        Run complete pipeline.

        With stop_on_error=True (default) the first failing step re-raises and the
        run aborts. With stop_on_error=False (--continue-on-error) a failing step is
        logged and the loop proceeds to the next step; the run then reaches the end
        and is marked complete even though one or more steps failed. In that partial-
        failure case, check the returned state (or the "Failed steps: N" warning) for
        per-step status rather than relying on the loop finishing.

        Args:
            stop_on_error: Whether to stop pipeline on first error (default: True)

        Returns:
            Final pipeline state (may contain failed steps when stop_on_error=False)

        Raises:
            Exception: If a step fails and stop_on_error=True. When stop_on_error=False,
                step failures are logged and the run continues; no exception is raised.
        """
        logger.info("=" * 80)
        logger.info(f"Starting seednap pipeline for marker: {self.config.marker.name}")
        logger.info("=" * 80)

        # pipeline.steps is the single source of truth: a stage runs iff listed, in order.
        # The order was already validated against the dependency DAG at config load.
        active_steps = list(self.config.pipeline.steps)
        logger.info(f"Pipeline steps: {' → '.join(active_steps)}")

        # Dispatch is built (and validated against VALID_STEPS) once in __init__.
        step_methods = self._step_methods

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

        # Mark pipeline as complete. NOTE: this point is reached unconditionally once
        # the step loop finishes. Under stop_on_error=False a failed step does not abort
        # the loop, so "completed" here means the loop ran to the end, not that every step
        # succeeded. The "Failed steps: N" warning below is the only signal of partial
        # failure; per-step status lives in the returned state.
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
            raise FileNotFoundError(
                f"Raw data directory not found: {raw_dir}. This is paths.raw_data in your config; "
                f"point it at the directory holding your paired-end FASTQ files (named like "
                f"<sample>_R1.fastq.gz / <sample>_R2.fastq.gz). Confirm the path exists and is "
                f"readable (ls '{raw_dir}'); a common cause is a config copied from another dataset."
            )

        # Find all R1 files (support both _R1 and .R1 naming)
        r1_patterns = [
            "*_R1*.fastq.gz", "*_R1*.fastq",
            "*.R1.fastq.gz", "*.R1.fastq",
        ]
        r1_files: List[Path] = []
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
            f"Could not find the {read} file for sample '{sample_name}' in {raw_dir}. seednap "
            f"expects paired files named like {sample_name}_{read}.fastq.gz (also accepts "
            f"{sample_name}.{read}.fastq.gz and {sample_name}_{read}_001.fastq.gz). One mate of "
            f"the pair is missing or named inconsistently; confirm both R1 and R2 exist with "
            f"matching sample names."
        )
