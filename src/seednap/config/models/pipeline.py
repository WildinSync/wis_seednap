"""Root pipeline config: composes every section above."""

from typing import Any

from pydantic import Field, model_validator

from seednap.config.models.base import StrictModel
from seednap.config.models.dada2 import Dada2Config
from seednap.config.models.input import (
    DemultiplexConfig,
    MarkerConfig,
    PathsConfig,
)
from seednap.config.models.operational import (
    CleaningConfig,
    LoggingConfig,
    PipelineStepsConfig,
)
from seednap.config.models.outputs import (
    ExportConfig,
    ReportConfig,
)
from seednap.config.models.swarm import SwarmConfig
from seednap.config.models.taxonomy import TaxonomicAssignmentConfig
from seednap.config.models.trimming import TrimmingConfig


# ===========================================================================
# ROOT: the complete pipeline config (composes every section above)
# ===========================================================================


class PipelineConfig(StrictModel):
    """Complete pipeline configuration: the root model composing every config section.

    One validated instance per marker run is the single source of truth the orchestrator
    drives from. ``marker`` and ``taxonomy`` are required; every other section has built-in
    defaults so a marker YAML need only state what differs.

    Attributes:
        marker: Marker identity and primers (required).
        paths: Raw input, output, and log directories.
        demultiplex: Demultiplexing step config.
        trimming: Primer-trimming step config.
        dada2: DADA2 (ASV) clustering path config.
        swarm: SWARM (OTU) clustering path config.
        taxonomy: Taxonomic assignment method and databases (required).
        export: GBIF / DarwinCore export config.
        report: Run-reporting config.
        cleaning: Control-decontamination config.
        logging: Run logging config.
        pipeline: Ordered list of stages to run.
    """

    marker: MarkerConfig = Field(..., description="Marker configuration")
    paths: PathsConfig = Field(default_factory=PathsConfig, description="Path configuration")
    demultiplex: DemultiplexConfig = Field(
        default_factory=DemultiplexConfig, description="Demultiplexing configuration"
    )
    trimming: TrimmingConfig = Field(
        default_factory=TrimmingConfig, description="Primer trimming configuration"
    )
    dada2: Dada2Config = Field(default_factory=Dada2Config, description="DADA2 configuration")
    swarm: SwarmConfig = Field(default_factory=SwarmConfig, description="SWARM clustering configuration")
    taxonomy: TaxonomicAssignmentConfig = Field(
        ..., description="Taxonomic assignment configuration"
    )
    export: ExportConfig = Field(
        default_factory=ExportConfig, description="Export configuration"
    )
    report: ReportConfig = Field(
        default_factory=ReportConfig, description="Run reporting configuration"
    )
    cleaning: CleaningConfig = Field(
        default_factory=CleaningConfig, description="Control decontamination configuration"
    )
    logging: LoggingConfig = Field(
        default_factory=LoggingConfig, description="Logging configuration"
    )
    pipeline: PipelineStepsConfig = Field(
        default_factory=PipelineStepsConfig, description="Pipeline steps configuration"
    )

    @model_validator(mode="after")
    def _validate_demultiplex_protocol(self) -> "PipelineConfig":
        """Reject an unrunnable demultiplexing protocol at load time, not mid-run.

        Only 'ligation' demultiplexing is implemented. Checked here (cross-field) rather than on
        DemultiplexConfig because the protocol only matters when demultiplexing actually runs
        (i.e. 'demultiplex' is listed in pipeline.steps). This catches both protocol='standard'
        (raises NotImplementedError mid-run) and the default protocol='none' (raises "Unknown
        protocol" mid-run) -- the latter is the likely slip of adding 'demultiplex' to steps
        without setting the protocol.

        Returns:
            The validated model (``self``) unchanged.

        Raises:
            ValueError: if ``demultiplex`` is in ``pipeline.steps`` while ``demultiplex.protocol``
                is anything other than the implemented ``ligation`` protocol.
        """
        if "demultiplex" in self.pipeline.steps and self.demultiplex.protocol != "ligation":
            raise ValueError(
                f"'demultiplex' is in pipeline.steps but demultiplex.protocol is "
                f"'{self.demultiplex.protocol}'. seednap currently implements only the 'ligation' "
                f"demultiplexing protocol. Either set demultiplex.protocol to 'ligation', or -- if "
                f"your reads are already demultiplexed (one sample per FASTQ pair) -- remove "
                f"'demultiplex' from pipeline.steps so the pipeline starts at trimming."
            )
        return self

    def model_post_init(self, __context: Any) -> None:
        """Create the output and log directories after the config validates.

        Side-effecting post-init: ``paths.output`` and ``paths.logs`` are created (parents
        included) as soon as a config loads, so loading is not a read-only operation.

        Args:
            __context: Pydantic post-init context object (unused).

        Raises:
            ValueError: if either directory cannot be created (e.g. a read-only mount or a
                path owned by another user); the message names the offending ``paths.*`` key.
        """
        # Create output directories if they don't exist
        for path_name in ["output", "logs"]:
            path = getattr(self.paths, path_name)
            try:
                path.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                raise ValueError(
                    f"Cannot create the paths.{path_name} directory '{path}': {e}. "
                    f"seednap creates paths.output and paths.logs when the config loads, and this "
                    f"location is not writable -- often a read-only mount, or a directory owned by "
                    f"another user (a common cause on the shared server is a config copied from a "
                    f"colleague that still points at their home). Set paths.{path_name} to a "
                    f"directory you own and can write to."
                ) from e
