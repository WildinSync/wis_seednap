"""Post-processing & operational config: cleaning, logging, step order."""

from typing import List, Literal

from pydantic import Field, model_validator

from seednap.config.models.base import StrictModel

# Valid pipeline stages. `pipeline.steps` is the single source of truth for what runs and in
# what order; a stage runs iff listed. dada2 and swarm are mutually exclusive feature paths.
VALID_STEPS = ("demultiplex", "trim", "dada2", "swarm", "taxonomy", "clean", "export", "report")


# ===========================================================================
# POST-PROCESSING & OPERATIONAL: control decontamination, logging, step order
# ===========================================================================


class CleaningConfig(StrictModel):
    """Control decontamination (cleaning) of the abundance table.

    Off by default. ``mode='flag'`` annotates OTUs/ASVs found in negative controls without
    changing counts (high-consequence subtraction stays opt-in); ``mode='subtract'`` removes
    those reads from the associated samples (extraction blanks clean their extraction batch,
    PCR blanks clean the whole dataset). Control identity comes from the FAIRe manifest.
    """

    mode: Literal["flag", "subtract"] = Field(
        default="flag",
        description="'flag' annotates control OTUs without changing counts; 'subtract' removes them",
    )


class LoggingConfig(StrictModel):
    """Logging configuration."""

    level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = Field(
        default="INFO", description="Logging level"
    )
    format: Literal["simple", "detailed", "json"] = Field(
        default="detailed", description="Log format"
    )
    file: bool = Field(default=True, description="Write logs to file")
    console: bool = Field(default=True, description="Write logs to console")


class PipelineStepsConfig(StrictModel):
    """Ordered list of pipeline stages to run.

    ``steps`` is the single source of truth for what runs and in what order: a stage runs iff
    it is listed, and each stage reads its own config section for parameters (defaults if the
    section is omitted). The order is validated against stage dependencies at load time.
    """

    steps: List[str] = Field(
        default=["trim", "dada2", "taxonomy", "export", "report"],
        description=(
            "Pipeline stages to run, in order; a stage runs iff listed. "
            f"Valid stages: {', '.join(VALID_STEPS)}. dada2 and swarm are mutually exclusive."
        ),
    )

    @model_validator(mode="after")
    def _validate_step_dependencies(self) -> "PipelineStepsConfig":
        """Reject unknown stages, duplicates, both feature paths, and invalid ordering.

        Each message names the offending stage and the exact fix so an error at load time is
        self-explanatory rather than a mid-run crash.
        """
        steps = self.steps
        unknown = [s for s in steps if s not in VALID_STEPS]
        if unknown:
            raise ValueError(
                f"pipeline.steps lists unknown stage(s) {unknown}. "
                f"Valid stages are: {', '.join(VALID_STEPS)}."
            )
        dupes = sorted({s for s in steps if steps.count(s) > 1})
        if dupes:
            raise ValueError(
                f"pipeline.steps lists duplicate stage(s) {dupes}; list each stage at most once."
            )
        if "dada2" in steps and "swarm" in steps:
            raise ValueError(
                "pipeline.steps lists both 'dada2' and 'swarm', but they are mutually exclusive "
                "feature-generation paths (DADA2 ASVs vs SWARM OTUs). Keep exactly one."
            )

        pos = {s: i for i, s in enumerate(steps)}

        def requires_before(stage: str, prereq: str) -> None:
            """Require ``prereq`` to be present and earlier than ``stage`` when ``stage`` is listed."""
            if stage in pos:
                if prereq not in pos:
                    raise ValueError(
                        f"pipeline.steps: '{stage}' requires '{prereq}' earlier in the list; "
                        f"add '{prereq}' before '{stage}'."
                    )
                if pos[prereq] > pos[stage]:
                    raise ValueError(
                        f"pipeline.steps: '{stage}' must come after '{prereq}'; reorder the list."
                    )

        # demultiplex (if present) precedes trim
        if "demultiplex" in pos and "trim" in pos and pos["demultiplex"] > pos["trim"]:
            raise ValueError("pipeline.steps: 'demultiplex' must come before 'trim'.")
        # feature paths need trim earlier
        requires_before("dada2", "trim")
        requires_before("swarm", "trim")
        # stages that consume the abundance table need a feature path earlier.
        # taxonomy/clean accept whichever feature path is present (dada2 OR swarm, which are
        # mutually exclusive), so they can't use requires_before (which names one prereq).
        feature_steps = [f for f in ("dada2", "swarm") if f in pos]
        for consumer in ("taxonomy", "clean"):
            if consumer in pos:
                if not feature_steps:
                    raise ValueError(
                        f"pipeline.steps: '{consumer}' requires a feature step ('dada2' or "
                        f"'swarm') earlier in the list."
                    )
                if not any(pos[f] < pos[consumer] for f in feature_steps):
                    raise ValueError(
                        f"pipeline.steps: '{consumer}' must come after '{feature_steps[0]}'."
                    )
        # clean consumes the taxonomy-annotated table, so it must follow taxonomy.
        requires_before("clean", "taxonomy")
        # export needs taxonomy earlier
        requires_before("export", "taxonomy")
        # when cleaning is requested, export must use the cleaned table, so clean
        # must come before export; otherwise export silently emits uncleaned counts.
        if "clean" in pos and "export" in pos and pos["clean"] > pos["export"]:
            raise ValueError(
                "pipeline.steps: 'clean' must come before 'export'; otherwise export "
                "would emit the uncleaned abundance table while you asked for control "
                "decontamination. Reorder the list so 'clean' precedes 'export'."
            )
        return self
