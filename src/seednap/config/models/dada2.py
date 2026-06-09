"""DADA2 (ASV) clustering path config."""

from typing import Literal, Optional

from pydantic import Field

from seednap.config.models.base import StrictModel


# ===========================================================================
# CLUSTERING -- DADA2 (ASV) path   [used only when "dada2" is in pipeline.steps;
# mutually exclusive with the SWARM path below]
# ===========================================================================


class Dada2FilterConfig(StrictModel):
    """DADA2 filtering parameters."""

    max_ee: float = Field(default=2.0, ge=0, description="Maximum expected errors")
    trunc_q: int = Field(default=11, ge=0, description="Truncate reads at first base with quality <= trunc_q")
    max_n: int = Field(default=0, ge=0, description="Maximum number of N bases allowed")
    rm_phix: bool = Field(default=True, description="Remove PhiX reads")
    min_len: Optional[int] = Field(default=None, ge=1, description="Minimum read length (optional)")
    max_len: Optional[int] = Field(default=None, ge=1, description="Maximum read length (optional)")


class Dada2MergeConfig(StrictModel):
    """DADA2 read merging parameters."""

    min_overlap: int = Field(default=20, ge=1, description="Minimum overlap for merging")
    max_mismatch: int = Field(default=0, ge=0, description="Maximum mismatches in overlap region")


class Dada2ChimeraConfig(StrictModel):
    """DADA2 chimera removal parameters."""

    method: Literal["consensus", "pooled", "none"] = Field(
        default="consensus", description="Chimera detection method"
    )


class Dada2Config(StrictModel):
    """DADA2 processing configuration."""

    filter: Dada2FilterConfig = Field(default_factory=Dada2FilterConfig)
    merge: Dada2MergeConfig = Field(default_factory=Dada2MergeConfig)
    chimera: Dada2ChimeraConfig = Field(default_factory=Dada2ChimeraConfig)
    pool: bool = Field(default=False, description="Pool samples for denoising")
    multithread: bool = Field(default=True, description="Use multithreading")
    collect_metrics: bool = Field(
        default=True,
        description="Collect ASV summary stats to metrics.json/csv + console (DADA2 path only)",
    )
    # DADA2-by-library: learn error models per sequencing run/library, then merge. Off by
    # default (and a no-op for single-library datasets, which is every current run). Only
    # matters when >= 2 libraries are combined into one DADA2 invocation. Requires a
    # sample->library grouping, taken from the manifest (report.sample_metadata) seq_run_id.
    per_library: bool = Field(
        default=False,
        description="Learn DADA2 error models per library/run, then merge (default: off)",
    )
