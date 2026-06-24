"""Pipeline orchestrator for end-to-end eDNA metabarcoding workflow execution.

This module is the top-level conductor for a single marker's run. eDNA
metabarcoding takes raw amplicon sequencing reads (PCR products of a short
"marker" gene region, e.g. teleo or MiFish for fish) and turns them into a
table of which taxa were detected in which water/soil samples. The
orchestrator runs the stages of that conversion in order, records what
happened so a run can be reconstructed, and lets a failed run resume from the
step that broke.

The stages, in pipeline order, are:
1. Demultiplexing (optional): split one multiplexed library FASTQ into
   per-sample files using sample-specific tags.
2. Primer trimming (cutadapt): remove the PCR primer sequences that flank
   every read so only the biological marker region remains.
3. Denoising/clustering: collapse millions of reads into biological features.
   DADA2 produces ASVs (Amplicon Sequence Variants, exact denoised sequences);
   SWARM produces OTUs (Operational Taxonomic Units, clusters of similar
   sequences). Both remove chimeras (artefactual sequences made of two parents).
4. Taxonomic assignment: label each ASV/OTU with a species/genus/family name
   by comparing it to a reference database (BLAST, DADA2 RDP, DECIPHER, ecotag).
5. Cleaning (optional): subtract contamination seen in negative controls.
6. Export: reshape the table into the GBIF / Darwin Core format for submission.
7. Reporting: read-tracking table and a self-contained HTML run report.

This file sits at ``src/seednap/pipeline/orchestrator.py`` and is driven by the
Click CLI; it delegates the heavy lifting to the processors under
``src/seednap/steps/`` and tracks progress via ``pipeline/state.py``.
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
        Initialize pipeline orchestrator: load config, set up logging, state, dispatch.

        Loads (or accepts) the marker configuration, creates the output directory
        tree, configures logging, and either resumes from an existing run's state
        JSON or starts a fresh run with every configured step marked pending. Also
        builds the step-name -> handler dispatch table and snapshots the effective
        merged config into the output tree for reproducibility.

        Args:
            config: Pipeline configuration, either a path (str/Path) to a marker
                YAML file to load, or an already-built PipelineConfig object.
            state_file: Path to the JSON state file tracking per-step progress. If
                None, defaults to ``<paths.output>/.<marker>_state.json``.
            resume: If True, re-read an existing state file and continue from where
                a previous run of this marker left off; if False, start fresh.

        Raises:
            FileNotFoundError: If a config path is given but the file does not exist
                (raised by load_config).
            ValueError: If resume=True but no state file exists at the resolved path.
            RuntimeError: If the step dispatch table is missing a handler for any
                stage in operational.VALID_STEPS (the two have drifted).
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
            # running version differs (or predates version stamping), the resumed steps may
            # not match what already ran.
            if self.state.seednap_version is None:
                logger.warning(
                    f"[WARN] resume version unknown: "
                    f"expected=a recorded seednap version, got=None (this state file "
                    f"predates version stamping), got_running={SEEDNAP_VERSION}, "
                    f"fallback=resuming anyway. Already-completed steps were produced by "
                    f"an unknown earlier version."
                )
            elif self.state.seednap_version != SEEDNAP_VERSION:
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
            "darwincore": self.run_darwincore,
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
        """Configure the run logger from the config's logging section.

        Wires up console and (when ``logging.file`` is set) file logging via the
        shared logging utility, naming the log file per marker and run timestamp so
        every run leaves a re-readable transcript. Also stashes the log file path on
        the instance so the HTML report can embed the transcript later.

        Args:
            None.

        Returns:
            None.
        """
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
        """Create the per-marker output directory tree (idempotent).

        Pre-creates the numbered stage directories under ``paths.output``
        (``01_trim``, ``02_dada2``, ``02_swarm``, ``03_taxo``) and the logs
        directory so downstream steps can write without checking for parents.
        Existing directories are left untouched.

        Args:
            None.

        Returns:
            None.
        """
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
        """Persist the current pipeline state to the run's state JSON.

        Writes ``self.state`` to ``self.state_file`` so progress survives an
        interruption and ``--resume`` can pick up from the last recorded step.

        Args:
            None.

        Returns:
            None.
        """
        self.state.save(self.state_file)

    def _write_config_snapshot(self) -> None:
        """Write the effective merged config (YAML) into the output tree and record
        its path in the state.

        Reproducibility: a run must be reconstructable from its outputs, so we persist
        the fully-merged config (defaults + the marker YAML) rather than relying on the
        original YAML still existing unchanged. The snapshot lives next to the state
        file at ``<paths.output>/.<marker>_config.snapshot.yaml`` and its path is stored
        on the state so the state JSON references it.

        Args:
            None.

        Returns:
            None. Side effects: writes the snapshot YAML, sets
            ``self.state.config_snapshot_path``, and saves the state.
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

        Raises:
            Exception: Re-raises whatever ``body`` raises, after recording the
                failure in the state and logging it. The state save and error log
                happen first so the failure is durable before propagation.
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

        Args:
            None.

        Returns:
            Path to the directory of primer-trimmed FASTQs that DADA2/SWARM read
            from. Equals the completed trim step's ``trimmed_dir`` output, or
            ``paths.raw_data`` only if no trim step ran.

        Raises:
            ValueError: If the trim step is marked completed but recorded no
                ``trimmed_dir`` in its outputs (a corrupt or stale state file).
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
        Run the demultiplexing step: split a multiplexed library into per-sample FASTQs.

        In a multiplexed sequencing library, many samples are pooled into one set of
        FASTQ files, each read carrying a short sample-specific tag. Demultiplexing
        reads those tags and writes one FASTQ pair per sample so the rest of the
        pipeline can treat each sample independently. The protocol
        (``demultiplex.protocol``) selects ligation- vs standard-tag handling.

        Returns:
            Dictionary of output paths. For the ligation protocol: ``demux_dir`` (the
            per-sample output directory) and ``trimmed_dir`` (the trimmer's per-sample
            result). On skip (already completed), the previously recorded outputs.

        Raises:
            ValueError: If the configured ``demultiplex.protocol`` is unknown, or (for
                the ligation protocol) if the required sample-tag metadata CSV is not
                set in the config.
            NotImplementedError: If protocol is 'standard' (not yet wired end-to-end).
        """
        # This step runs iff "demultiplex" is listed in pipeline.steps. If your raw inputs are
        # already demultiplexed (one FASTQ per sample), simply omit "demultiplex" from steps.
        def body() -> Dict[str, Any]:
            """Do this step's work and return its output paths (invoked by _execute_step).

            Returns:
                Dict mapping this step's output names to their file paths.
            """
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
        """Demultiplex a ligation-protocol library and trim primers in one pass.

        Ligation-protocol libraries carry the sample tag ligated onto the read; a
        single multiplexed FASTQ is split into per-sample files using the
        tag-to-sample mapping in the metadata CSV, and primers are trimmed at the
        same time via LigationTrimmer. A sample failing to meet the configured
        failure-rate threshold aborts the step (enforced inside the trimmer).

        Args:
            None.

        Returns:
            Dictionary with ``demux_dir`` (the per-sample output directory under
            ``01_trim/<marker>/demux``) and ``trimmed_dir`` (the trimmer's returned
            per-sample result).

        Raises:
            ValueError: If ``demultiplex.metadata`` (the sample-tag CSV) is not set.
        """
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
        Run the primer trimming step over every sample (cutadapt).

        Every read still begins/ends with the PCR primer sequences used to amplify
        the marker; these are technical, not biological, and must be removed before
        denoising or clustering or they corrupt feature inference. This step locates
        each sample's R1/R2 FASTQ pair under ``paths.raw_data`` and runs the
        StandardTrimmer (cutadapt) to strip the forward/reverse primers, writing the
        trimmed pairs under ``01_trim/<marker>``.

        Returns:
            Dictionary with ``trimmed_dir`` (the output directory consumed by
            DADA2/SWARM) and ``samples`` (a dict mapping each sample name to its
            per-sample trim output tuple). On skip, the previously recorded outputs.

        Raises:
            FileNotFoundError: If the raw data directory or a sample's R1/R2 file
                cannot be found (raised by the sample-discovery helpers).
        """
        def body() -> Dict[str, Any]:
            """Do this step's work and return its output paths (invoked by _execute_step).

            Returns:
                Dict mapping this step's output names to their file paths.
            """
            logger.info("Running primer trimming")

            trimmer = StandardTrimmer(
                cores=self.config.trimming.cores,
                error_rate=self.config.trimming.max_error_rate,
                min_length=self.config.trimming.min_length,
                overlap=self.config.trimming.overlap,
            )

            output_dir = self.config.paths.output / "01_trim" / self.config.marker.name

            # Remove trimmed reads left over from a PREVIOUS run before producing this
            # run's set. The downstream feature step (SWARM/DADA2) discovers its inputs by
            # scanning this directory, so without this a re-run that finds a different
            # sample set (e.g. raw_data was corrected, or the new raw dir holds fewer
            # samples) would leave stale per-sample FASTQs here and silently process the
            # earlier run's samples -- producing results for the wrong dataset while
            # reporting success. This runs only on a fresh trim; --resume skips the trim
            # step entirely, so cached outputs are preserved. (no-silent-fallbacks policy)
            if output_dir.exists():
                stale = (
                    list(output_dir.glob("*.fastq"))
                    + list(output_dir.glob("*.fastq.gz"))
                    + list((output_dir / "logs").glob("*_trim_pass*.txt"))
                )
                for f in stale:
                    f.unlink()
                if stale:
                    logger.info(
                        f"Cleared {len(stale)} file(s) from a previous trim run in "
                        f"{output_dir} before re-trimming, so stale samples cannot be "
                        f"reused by the downstream feature step."
                    )

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
        Run the DADA2 denoising step: turn trimmed reads into ASVs.

        DADA2 is the ASV (Amplicon Sequence Variant) path: it models per-run
        sequencing error to infer the exact biological sequences present, merges
        read pairs, builds a sample-by-sequence count table, and removes chimeras
        (artefactual sequences spliced from two real templates during PCR). All the
        filter/merge/chimera knobs come from the config's dada2 section. After the
        step completes, a post-hook cross-checks the manifest eventIDs against the
        abundance table to catch silent sample-ID mismatches.

        Returns:
            Dictionary of output paths from the DADA2 processor (includes the query
            FASTA of ASV sequences and the transposed clean count table). On skip,
            the previously recorded outputs.

        Raises:
            Exception: Re-raises any failure from the DADA2 processor or the R
                subprocess it drives (recorded in the state before propagation).
        """
        def body() -> Dict[str, Any]:
            """Do this step's work and return its output paths (invoked by _execute_step).

            Returns:
                Dict mapping this step's output names to their file paths.
            """
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
        Run the SWARM clustering step: turn trimmed reads into OTUs.

        SWARM is the OTU (Operational Taxonomic Unit) path: read pairs are merged
        and dereplicated (via vsearch), then SWARM agglomerates near-identical
        sequences into clusters using a single-linkage local threshold ``d``
        (optionally with fastidious mode), and chimeras are removed when the config
        requests it. OTUs are a coarser feature than DADA2's ASVs, so their count
        legitimately differs. The same manifest-vs-abundance post-hook runs to catch
        silent sample-ID mismatches.

        Returns:
            Dictionary of output paths from the SWARM processor (includes the query
            FASTA of OTU representative sequences and the OTU count table). On skip,
            the previously recorded outputs.

        Raises:
            Exception: Re-raises any failure from the SWARM processor or the vsearch/
                swarm subprocesses it drives (recorded in the state before propagation).
        """
        def body() -> Dict[str, Any]:
            """Do this step's work and return its output paths (invoked by _execute_step).

            Returns:
                Dict mapping this step's output names to their file paths.
            """
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

        Args:
            None.

        Returns:
            Path to the written ``library_map.csv`` (columns: ``sample``, ``library``)
            under ``02_dada2/<marker>``, or None when ``dada2.per_library`` is off, no
            grouping source is configured, or building the map fails (a ``[WARN]`` is
            logged in the latter two cases and DADA2 falls back to single-batch).
        """
        if not self.config.dada2.per_library:
            return None
        field_csv = self.config.report.sample_metadata
        lab_csv = self.config.demultiplex.metadata
        src = field_csv or lab_csv
        if src is None:
            # No metadata grouping configured. If raw_data is organized into per-library
            # subdirectories (one folder per sequencing library/run of already-demultiplexed
            # per-sample FASTQs), derive the sample->library map from the subfolder each
            # sample's R1 file lives in -- no lab metadata needed.
            subdir_map = self._library_map_from_subdirs()
            if subdir_map is not None:
                return subdir_map
            logger.warning(
                "[WARN] dada2 per_library: expected=report.sample_metadata or "
                "demultiplex.metadata (or a per-library subdirectory layout under raw_data) "
                "for the library grouping, got=none, fallback=standard single-batch DADA2"
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

    def _library_map_from_subdirs(self) -> Optional[Path]:
        """Derive a ``sample,library`` map from a per-library subdirectory layout of raw_data.

        Used for DADA2-by-library when ``dada2.per_library`` is on but no metadata grouping is
        configured. When raw_data holds no FASTQs at its top level but keeps them in
        subdirectories (one folder per sequencing library/run of already-demultiplexed
        per-sample reads), each sample's library is the immediate subfolder its R1 file lives
        in. The resulting ``sample,library`` CSV is identical in shape to the manifest-derived
        map, so the DADA2 R per-library branch consumes it unchanged -- no metadata required.

        Returns:
            Path to the written ``library_map.csv`` when at least two libraries are found,
            else None (a flat layout, a single library, or a sample name colliding across
            subdirectories all fall back to the standard single-batch DADA2 path).
        """
        import re

        import pandas as pd

        raw_dir = self.config.paths.raw_data
        if not raw_dir.exists():
            return None
        r1_patterns = ["*_R1*.fastq.gz", "*_R1*.fastq", "*.R1.fastq.gz", "*.R1.fastq"]
        # Only derive from subdirectories when the top level itself has no per-sample FASTQs
        # (mirrors sample discovery, which prefers the top level and recurses only if empty).
        if any(True for pat in r1_patterns for _ in raw_dir.glob(pat)):
            return None
        sample_lib: Dict[str, str] = {}
        for pat in r1_patterns:
            for f in raw_dir.rglob(pat):
                if f.parent == raw_dir:
                    continue
                m = re.match(r"(.+?)[._]R[12]", f.name)
                if not m:
                    continue
                sample = m.group(1)
                library = f.parent.relative_to(raw_dir).parts[0]
                if sample in sample_lib and sample_lib[sample] != library:
                    logger.warning(
                        f"[WARN] dada2 per_library: sample {sample!r} appears under two "
                        f"libraries ({sample_lib[sample]!r}, {library!r}); cannot derive a "
                        f"library grouping from subdirectories, fallback=single-batch DADA2"
                    )
                    return None
                sample_lib[sample] = library
        if len(set(sample_lib.values())) < 2:
            return None
        df = pd.DataFrame(
            [{"sample": s, "library": lib} for s, lib in sorted(sample_lib.items())]
        )
        out = self.config.paths.output / "02_dada2" / self.config.marker.name / "library_map.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
        logger.info(
            f"DADA2-by-library: derived library map from raw_data subdirectories "
            f"({len(df)} samples, {df['library'].nunique()} libraries) -> {out}"
        )
        return out

    def _report_dir(self) -> Path:
        """Per-marker directory for report artifacts.

        Honors ``report.output_dir`` (a per-marker subdirectory is created
        inside it); otherwise defaults to ``<paths.output>/04_report/<marker>``.

        Args:
            None.

        Returns:
            Path to the marker's report directory. Not created here; callers create
            it before writing.
        """
        marker = self.config.marker.name
        base = self.config.report.output_dir
        if base is not None:
            return Path(base) / marker
        return self.config.paths.output / "04_report" / marker

    def _build_read_tracking_report(self, method: str) -> None:
        """Build the per-sample read-tracking table after a feature step.

        The read-tracking table records how many reads each sample retained at each
        stage (raw -> trimmed -> denoised/clustered), so a biologist can spot a
        sample that lost most of its reads (a failed PCR, a dirty blank) before
        trusting its detections. This writes the table and a run-level step summary,
        evaluates retention warnings, and stores a compact per-sample summary into
        the step's state metadata so it survives ``--resume`` and feeds the FAIRe
        manifest's ``reads_*`` columns.

        Non-fatal: a reporting failure logs a ``[WARN]`` and never fails the
        run (the no-silent-fallbacks policy -- the report is observational only).

        Args:
            method: The feature step that produced the table, ``"dada2"`` or
                ``"swarm"``; selects which directory/table the builder reads.

        Returns:
            None. Side effects: writes the tracking table and step summary under the
            report directory and mutates the matching step's ``metadata``.
        """
        try:
            import pandas as pd

            from seednap.steps.report import ReadTrackingBuilder

            marker = self.config.marker.name
            out = self.config.paths.output
            report_dir = self._report_dir()
            kwargs: Dict[str, Any] = {
                "marker": marker,
                # Cutadapt per-sample logs are written by the trim step under
                # <output>/01_trim/<marker>/logs (see trimming_pipeline.StandardTrimmer,
                # log_dir = output_dir / "logs"). Read them from the same place so the
                # report can recover raw/trimmed counts; otherwise % retained is NA.
                "logs_dir": out / "01_trim" / marker / "logs",
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

        Args:
            method: The feature step whose abundance table to check, ``"dada2"``
                (reads ``seqtab_clean_t.csv``) or ``"swarm"`` (reads ``otu_table.csv``).

        Returns:
            None. Returns early without checking when no field metadata is
            configured or the inputs are missing; any mismatch surfaces as a
            ``[WARN]`` log inside the validator.
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

        Args:
            None.

        Returns:
            None. Side effects: writes ``report.html`` under the report directory
            when enabled; returns early (no-op) when ``report.html_report`` is off.
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
        Run the taxonomic assignment step: label each ASV/OTU with a taxon name.

        Each ASV/OTU is just an anonymous DNA sequence until it is matched against a
        reference database of known sequences to assign a taxonomy (species, genus,
        family, ...). This step picks up the query FASTA and count table from the
        completed feature step (DADA2 or SWARM), selects the method-specific
        parameters and reference database from the config (``blast``, ``dada2`` RDP,
        ``ecotag``, or ``decipher``), and delegates to the TaxonomicAssigner, which
        writes a merged taxonomy+abundance table.

        Returns:
            Dictionary of output paths from the assigner, including ``final_table``
            (the merged taxonomy+abundance CSV consumed by clean/export). On skip,
            the previously recorded outputs.

        Raises:
            ValueError: If neither dada2 nor swarm is marked completed in this run's
                state, or if the completed feature step recorded no ``query_fasta``
                or count table.
            Exception: Re-raises any failure from the assigner or its subprocesses
                (recorded in the state before propagation).
        """
        def body() -> Dict[str, Any]:
            """Do this step's work and return its output paths (invoked by _execute_step).

            Returns:
                Dict mapping this step's output names to their file paths.
            """
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
        uncleaned table).

        Negative controls (extraction blanks, PCR blanks) capture contamination that
        leaked in during lab handling; cleaning subtracts what those controls saw
        from the real samples so reported detections are not lab artefacts. The
        cleaning ``mode`` comes from the config's cleaning section.

        Returns:
            Dictionary with ``cleaned_table`` (the decontaminated CSV) and
            ``cleaning_report`` (the per-feature cleaning report CSV), or an empty
            dict when cleaning is skipped (no control metadata) or errors out. On a
            prior skip/completion, the previously recorded outputs.

        Raises:
            None. All failures are caught and recorded as a skipped step with a
            ``[WARN]``; the run is never failed over cleaning.
        """
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

        Args:
            export_step: The export step's StepState (or None). Only triggers a
                warning when it is a completed export that finished before a later
                completed clean step that produced a cleaned table.

        Returns:
            None. Side effect: logs a ``[WARN]`` when the existing export is stale
            relative to the cleaned table; otherwise a no-op.
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
        Run the export step: reshape the final table into GBIF / Darwin Core format.

        GBIF (the Global Biodiversity Information Facility) and the Darwin Core
        standard define how occurrence records must be structured for submission and
        sharing. This step reads the merged taxonomy+abundance table from the
        taxonomy step (preferring the decontaminated table if the clean step produced
        one) and writes a GBIF-formatted CSV, optionally adding rank/taxon columns.

        Returns:
            Dictionary with ``gbif_csv`` (path to the GBIF-formatted output CSV). On
            skip, the previously recorded outputs (with a staleness warning if a
            later clean step has since produced a cleaned table).

        Raises:
            ValueError: If the taxonomy step did not complete, is missing from the
                state, or recorded no ``final_table`` to format.
            Exception: Re-raises any failure from the GBIF formatter (recorded in the
                state before propagation).
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

    def run_darwincore(self) -> Dict[str, Path]:
        """DarwinCore occurrence export: build the GBIF-ready occurrence CSV in-pipeline.

        Runs the DarwinCore builder as a pipeline step: it joins the long-format export output
        (from the 'export' step) to the per-sample and per-project metadata
        (``report.sample_metadata`` / ``report.project_metadata``), fills the DarwinCore fields,
        removes control and non-target rows, and (unless ``export.darwincore.skip_enrichment``)
        enriches the higher ranks from NCBI/WoRMS. Both metadata files are required: the step
        fails fast with a clear error if either is unset, rather than emitting an occurrence
        file with blank provenance.

        Returns:
            Dictionary with ``darwincore_csv`` (path to the DarwinCore occurrence CSV), or the
            previously recorded outputs on a skip.

        Raises:
            ValueError: if the 'export' step did not complete or recorded no GBIF table, or if
                ``report.sample_metadata`` / ``report.project_metadata`` is not configured.
            Exception: re-raises any failure from the DarwinCore builder (recorded in state).
        """
        step_name = "darwincore"
        if not self._should_run_step(step_name):
            step = self.state.get_step(step_name)
            return step.outputs if step else {}

        log_pipeline_step(step_name, "start", logger)
        self.state.start_step(step_name)
        self._save_state()

        try:
            from seednap.steps.formatting.darwincore_builder import DarwinCoreBuilder

            logger.info("Building DarwinCore occurrence file")
            export_step = self.state.get_step("export")
            gbif_csv = export_step.outputs.get("gbif_csv") if export_step else None
            if gbif_csv is None:
                raise ValueError(
                    "Cannot build the DarwinCore file: the 'export' step did not complete or "
                    "recorded no 'gbif_csv' (the long-format table the DarwinCore builder joins "
                    "metadata onto). Ensure 'export' runs and completes before 'darwincore'."
                )
            sample_meta = self.config.report.sample_metadata
            project_meta = self.config.report.project_metadata
            missing = [
                name for name, val in (
                    ("report.sample_metadata", sample_meta),
                    ("report.project_metadata", project_meta),
                ) if val is None
            ]
            if missing:
                raise ValueError(
                    f"The 'darwincore' step needs per-sample and per-project metadata, but "
                    f"{', '.join(missing)} is not set. Set both to this dataset's metadata CSVs "
                    f"(they supply each occurrence's eventDate, coordinates, recorder, sequencing "
                    f"method and reference database), or remove 'darwincore' from pipeline.steps "
                    f"if you only need the long-format export."
                )

            # narrowed by the `missing` check above: both are non-None here
            assert sample_meta is not None and project_meta is not None

            # Auto-fill the reference-database and chimera-removal provenance from the run
            # config (the single source of truth) so they need not be re-entered in the
            # project metadata; a differing project value is reported by the builder.
            otu_db = None
            try:
                db = self.config.taxonomy.get_database_config()
                db_path = (
                    getattr(db, "fasta", None)
                    or getattr(db, "all", None)
                    or getattr(db, "trained", None)
                )
                if db_path:
                    otu_db = Path(db_path).name
            except Exception:  # noqa: BLE001 -- provenance is best-effort; never fail the step
                otu_db = None
            chimera_check = None
            if self.state.is_step_completed("dada2"):
                method = getattr(self.config.dada2.chimera, "method", "consensus")
                chimera_check = (
                    "not performed" if method == "none"
                    else f"removeBimeraDenovo (DADA2 {method})"
                )
            elif self.state.is_step_completed("swarm"):
                chimera_check = "uchime_denovo (VSEARCH)"

            output_path = (
                self.config.paths.output
                / f"{self.config.marker.name}_{self.config.taxonomy.method}_darwincore.csv"
            )
            builder = DarwinCoreBuilder(
                taxonomy_results_path=Path(gbif_csv),
                sample_metadata_path=Path(sample_meta),
                project_metadata_path=Path(project_meta),
                output_path=output_path,
                summarise_pcr_replicates=self.config.export.darwincore.summarise_pcr_replicates,
                skip_enrichment=self.config.export.darwincore.skip_enrichment,
                otu_db=otu_db,
                chimera_check=chimera_check,
            )
            builder.build()
            outputs: Dict[str, Any] = {"darwincore_csv": output_path}
            dropped = getattr(builder, "dropped_report_path", None)
            if dropped is not None:
                outputs["dropped_report"] = dropped

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
        a reporting failure logs a ``[WARN]`` and never fails the run (the no-silent-fallbacks policy).

        Returns:
            Empty dict. Side effects: writes the read-tracking table, step summary,
            and (when enabled) the HTML report. Marks the step skipped (with a
            ``[WARN]``) when no completed dada2/swarm step exists to report on. On a
            prior skip/completion, the previously recorded outputs.

        Raises:
            None. Reporting failures are caught downstream and logged as ``[WARN]``;
            the run is never failed over the report.
        """
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
        Discover sample names by scanning the raw data directory for R1 FASTQs.

        Globs ``paths.raw_data`` for forward-read (R1) files under several common
        naming conventions and derives each sample name as the text before the
        ``_R1``/``.R1`` marker, so the trim step knows which samples to process.

        Args:
            None.

        Returns:
            Sorted list of unique sample names found in the raw data directory.

        Raises:
            FileNotFoundError: If ``paths.raw_data`` does not exist.
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

        # Find all R1 files (support both _R1 and .R1 naming). Search the top level first;
        # if nothing is there, search subdirectories recursively, so already-demultiplexed
        # raw data organised into per-library / per-run subfolders is picked up without the
        # user having to flatten it.
        r1_patterns = [
            "*_R1*.fastq.gz", "*_R1*.fastq",
            "*.R1.fastq.gz", "*.R1.fastq",
        ]
        r1_files: List[Path] = []
        for pattern in r1_patterns:
            r1_files.extend(raw_dir.glob(pattern))
        if not r1_files:
            for pattern in r1_patterns:
                r1_files.extend(raw_dir.rglob(pattern))
            if r1_files:
                logger.info(
                    f"No FASTQs at the top level of {raw_dir}; found {len(r1_files)} R1 file(s) "
                    f"in subdirectories and will trim those (per-library/per-run layout)."
                )

        # No per-sample forward-read FASTQs anywhere (top level or subdirectories). Fail
        # loudly rather than returning an empty list and silently producing an empty run.
        # (no-silent-fallbacks policy)
        if not r1_files:
            raise FileNotFoundError(
                f"No forward-read FASTQ files found anywhere under paths.raw_data ({raw_dir}) "
                f"(searched the top level and subdirectories). Expected per-sample files named "
                f"<sample>_R1.fastq.gz / <sample>_R2.fastq.gz (or .R1/.R2). Check that "
                f"paths.raw_data points at the intended directory and that the files use that naming."
            )

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

        # Not at the top level: search subdirectories (per-library/per-run layout). Require
        # exactly one match -- if a sample name resolves to several files across subfolders
        # the choice is ambiguous, so raise rather than silently pick one.
        matches: List[Path] = []
        for pattern in patterns:
            matches.extend(raw_dir.rglob(pattern))
        unique = list(dict.fromkeys(matches))
        if len(unique) == 1:
            return unique[0]
        if len(unique) > 1:
            raise FileNotFoundError(
                f"Ambiguous {read} file for sample '{sample_name}': found {len(unique)} matching "
                f"files under {raw_dir} ({', '.join(str(p) for p in unique[:4])}). Sample names "
                f"must be unique across subdirectories; rename or separate the colliding samples."
            )

        raise FileNotFoundError(
            f"Could not find the {read} file for sample '{sample_name}' under {raw_dir} "
            f"(searched the top level and subdirectories). seednap expects paired files named "
            f"like {sample_name}_{read}.fastq.gz (also accepts {sample_name}.{read}.fastq.gz and "
            f"{sample_name}_{read}_001.fastq.gz). One mate of the pair is missing or named "
            f"inconsistently; confirm both R1 and R2 exist with matching sample names."
        )
