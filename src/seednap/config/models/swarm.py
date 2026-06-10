"""SWARM (OTU) clustering path config."""

from typing import Literal

from pydantic import Field

from seednap.config.models.base import StrictModel


# ===========================================================================
# CLUSTERING -- SWARM (OTU) path   [used only when "swarm" is in pipeline.steps;
# mutually exclusive with the DADA2 path above]
# ===========================================================================


class SwarmMergeConfig(StrictModel):
    """vsearch paired-read merging parameters for the SWARM path.

    On the SWARM path, ``vsearch`` merges the paired reads into full amplicons before
    clustering (DADA2 does its own merging).

    Attributes:
        fastq_maxdiffs: Maximum mismatches tolerated in the read overlap.
        fastq_minovlen: Minimum overlap length (bases) required to merge a pair.
        allow_stagger: Allow merging of staggered read pairs (reads of unequal alignment).
    """

    fastq_maxdiffs: int = Field(default=10, ge=0, description="Max differences in overlap region")
    fastq_minovlen: int = Field(default=10, ge=1, description="Min overlap length for merging")
    allow_stagger: bool = Field(default=False, description="Allow merging of staggered reads")


class SwarmClusteringConfig(StrictModel):
    """SWARM clustering algorithm parameters.

    SWARM groups amplicons into OTUs by single-linkage clustering at a small local distance
    ``d`` rather than a fixed global identity cutoff, giving fine-grained, abundance-aware
    clusters.

    Attributes:
        d: Local clustering distance threshold (number of differences between linked amplicons).
        fastidious: Enable fastidious mode, which grafts low-abundance singleton clusters onto
            larger ones to reduce over-splitting.
        boundary: Minimum abundance (mass) for an OTU to count as "large" in fastidious mode.
        threads: Number of threads SWARM may use.
    """

    d: int = Field(default=1, ge=1, description="Clustering distance threshold")
    fastidious: bool = Field(default=True, description="Enable fastidious mode (refine singletons)")
    boundary: int = Field(default=3, ge=1, description="Min mass for large OTUs in fastidious mode")
    threads: int = Field(default=4, ge=1, description="Number of threads")


class SwarmChimeraConfig(StrictModel):
    """SWARM chimera detection parameters.

    Chimeras are artefactual PCR sequences formed from two unrelated templates; ``denovo``
    detection (via vsearch) screens them out before they become spurious OTUs.

    Attributes:
        method: Chimera detection method (denovo / none).
    """

    method: Literal["denovo", "none"] = Field(
        default="denovo", description="Chimera detection method"
    )


class SwarmConfig(StrictModel):
    """SWARM (OTU) clustering path configuration.

    The SWARM path merges reads with vsearch, screens chimeras, and clusters amplicons into
    operational taxonomic units (OTUs); it is the alternative to the DADA2 ASV path (the two
    are mutually exclusive). This composes the merge/clustering/chimera sub-configs.

    Attributes:
        merge: vsearch paired-read merging parameters.
        clustering: SWARM clustering algorithm parameters.
        chimera: Chimera detection parameters.
        min_sequence_length: Minimum merged-amplicon length (bases) to keep before clustering.
    """

    merge: SwarmMergeConfig = Field(default_factory=SwarmMergeConfig)
    clustering: SwarmClusteringConfig = Field(default_factory=SwarmClusteringConfig)
    chimera: SwarmChimeraConfig = Field(default_factory=SwarmChimeraConfig)
    min_sequence_length: int = Field(
        default=20, ge=1, description="Min sequence length after merging"
    )
