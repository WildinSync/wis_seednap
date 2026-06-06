"""Trimming config: Cutadapt two-pass primer removal."""

from typing import Literal

from pydantic import Field

from seednap.config.models.base import StrictModel


# ===========================================================================
# TRIMMING: Cutadapt two-pass primer removal
# ===========================================================================


class TrimmingConfig(StrictModel):
    """Primer trimming configuration."""

    tool: Literal["cutadapt"] = Field(default="cutadapt", description="Trimming tool to use")
    min_length: int = Field(default=20, ge=1, description="Minimum read length after trimming")
    max_error_rate: float = Field(
        default=0.1, ge=0.0, le=1.0, description="Maximum error rate for primer matching"
    )
    cores: int = Field(default=1, ge=1, description="Number of CPU cores to use")
    discard_untrimmed: bool = Field(
        default=True, description="Discard reads without detected primers"
    )
    overlap: int = Field(default=3, ge=1, description="Minimum overlap for primer detection")
