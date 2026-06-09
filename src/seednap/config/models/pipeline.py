"""Root pipeline config: composes every section above."""

from typing import Any

from pydantic import Field

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
    MetricsConfig,
    ReportConfig,
)
from seednap.config.models.swarm import SwarmConfig
from seednap.config.models.taxonomy import TaxonomicAssignmentConfig
from seednap.config.models.trimming import TrimmingConfig


# ===========================================================================
# ROOT: the complete pipeline config (composes every section above)
# ===========================================================================


class PipelineConfig(StrictModel):
    """Complete pipeline configuration."""

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
    metrics: MetricsConfig = Field(
        default_factory=MetricsConfig, description="Metrics configuration"
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

    def model_post_init(self, __context: Any) -> None:
        """Post-initialization validation."""
        # Create output directories if they don't exist
        for path_name in ["output", "logs"]:
            path = getattr(self.paths, path_name)
            path.mkdir(parents=True, exist_ok=True)
