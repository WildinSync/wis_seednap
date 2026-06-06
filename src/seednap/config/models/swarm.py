"""SWARM (OTU) clustering path config."""

from typing import Literal

from pydantic import Field

from seednap.config.models.base import StrictModel


# ===========================================================================
# CLUSTERING -- SWARM (OTU) path   [used only when "swarm" is in pipeline.steps;
# mutually exclusive with the DADA2 path above]
# ===========================================================================


class SwarmMergeConfig(StrictModel):
    """vsearch read merging parameters for SWARM pipeline."""

    fastq_maxdiffs: int = Field(default=10, ge=0, description="Max differences in overlap region")
    fastq_minovlen: int = Field(default=10, ge=1, description="Min overlap length for merging")
    allow_stagger: bool = Field(default=False, description="Allow merging of staggered reads")


class SwarmClusteringConfig(StrictModel):
    """SWARM clustering algorithm parameters."""

    d: int = Field(default=1, ge=1, description="Clustering distance threshold")
    fastidious: bool = Field(default=True, description="Enable fastidious mode (refine singletons)")
    boundary: int = Field(default=3, ge=1, description="Min mass for large OTUs in fastidious mode")
    threads: int = Field(default=4, ge=1, description="Number of threads")


class SwarmChimeraConfig(StrictModel):
    """SWARM chimera detection parameters."""

    method: Literal["denovo", "none"] = Field(
        default="denovo", description="Chimera detection method"
    )


class SwarmConfig(StrictModel):
    """SWARM OTU clustering pipeline configuration."""

    merge: SwarmMergeConfig = Field(default_factory=SwarmMergeConfig)
    clustering: SwarmClusteringConfig = Field(default_factory=SwarmClusteringConfig)
    chimera: SwarmChimeraConfig = Field(default_factory=SwarmChimeraConfig)
    min_sequence_length: int = Field(
        default=20, ge=1, description="Min sequence length after merging"
    )
