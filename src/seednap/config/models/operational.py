"""Post-processing & operational config: cleaning, logging, step order."""

from typing import List, Literal

from pydantic import Field

from seednap.config.models.base import StrictModel


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

    enabled: bool = Field(default=False, description="Run the cleaning step (default: off)")
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
    """Pipeline steps configuration."""

    steps: List[str] = Field(
        default=["trim", "dada2", "taxonomy", "export"],
        description="Pipeline steps to execute in order",
    )
    skip: List[str] = Field(default_factory=list, description="Steps to skip")
