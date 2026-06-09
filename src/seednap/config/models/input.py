"""Input config: marker identity, primers, paths, demultiplexing."""

from pathlib import Path
from typing import Any, Literal, Optional

from pydantic import Field, field_validator

from seednap.config.models.base import StrictModel


# ===========================================================================
# INPUT: marker identity, primers, paths, demultiplexing
# ===========================================================================


class PrimerConfig(StrictModel):
    """Primer pair configuration."""

    forward: str = Field(..., min_length=10, description="Forward primer sequence (5' to 3')")
    reverse: str = Field(..., min_length=10, description="Reverse primer sequence (5' to 3')")

    @field_validator("forward", "reverse")
    @classmethod
    def validate_dna_sequence(cls, v: str) -> str:
        """Validate that the sequence contains only valid DNA bases including IUPAC ambiguity codes."""
        valid_bases = set("ACGTRYMKSWHBVDN")
        v_upper = v.upper()
        if not all(base in valid_bases for base in v_upper):
            invalid_bases = set(v_upper) - valid_bases
            raise ValueError(
                f"Invalid DNA sequence. Contains invalid bases: {invalid_bases}. "
                f"Valid bases are: {', '.join(sorted(valid_bases))}"
            )
        return v_upper


class MarkerConfig(StrictModel):
    """Marker-specific configuration (e.g., teleo, 16S, COI)."""

    name: str = Field(..., description="Marker name (lowercase)")
    description: Optional[str] = Field(None, description="Marker description")
    primers: PrimerConfig = Field(..., description="Primer pair configuration")


class PathsConfig(StrictModel):
    """File paths configuration."""

    raw_data: Path = Field(default=Path("data/raw"), description="Raw FASTQ data directory")
    output: Path = Field(default=Path("outputs"), description="Output directory")
    logs: Path = Field(default=Path("logs"), description="Log files directory")

    @field_validator("raw_data", "output", "logs")
    @classmethod
    def expand_path(cls, v: Path) -> Path:
        """Expand ~ and relative paths to absolute paths."""
        return v.expanduser().resolve()


class DemultiplexConfig(StrictModel):
    """Demultiplexing configuration."""

    protocol: Literal["ligation", "standard", "none"] = Field(
        default="none", description="Demultiplexing protocol type"
    )
    metadata: Optional[Path] = Field(default=None, description="Path to metadata CSV file")
    # If more than this fraction of samples fail during demultiplexing, abort.
    # Otherwise log the failures and continue. Default 0.5 = abort if >50% fail.
    max_sample_failure_rate: float = Field(
        default=0.5, ge=0.0, le=1.0,
        description="Abort demultiplex if more than this fraction of samples fail",
    )

    @field_validator("metadata")
    @classmethod
    def validate_metadata_path(cls, v: Optional[Path], info: Any) -> Optional[Path]:
        """Validate that metadata file exists if demultiplexing is enabled."""
        if v is not None:
            v = v.expanduser().resolve()
            # Note: we don't check file existence here since config might be created
            # before the file exists. Validation happens at runtime.
        return v
