"""Input config: marker identity, primers, paths, demultiplexing."""

from pathlib import Path
from typing import Any, Literal, Optional

from pydantic import Field, ValidationInfo, field_validator

from seednap.config.models.base import StrictModel
from seednap.utils.logging import get_logger

logger = get_logger(__name__)

# (field, primer) pairs already reported as inosine-converted: the config is parsed more
# than once per command (validate, then load), so warn once per primer per process.
_INOSINE_WARNED: set = set()


# ===========================================================================
# INPUT: marker identity, primers, paths, demultiplexing
# ===========================================================================


class PrimerConfig(StrictModel):
    """Forward and reverse PCR primer sequences for the marker.

    Primers are the short oligonucleotides that bracket and amplify the target barcode region;
    trimming removes them from the reads before clustering. Both are given 5' to 3' and may use
    IUPAC ambiguity codes. Inosine (``I``) is accepted and converted to ``N``.

    Attributes:
        forward: Forward primer sequence, 5' to 3' (>= 10 bases).
        reverse: Reverse primer sequence, 5' to 3' (>= 10 bases).
    """

    forward: str = Field(..., min_length=10, description="Forward primer sequence (5' to 3')")
    reverse: str = Field(..., min_length=10, description="Reverse primer sequence (5' to 3')")

    @field_validator("forward", "reverse")
    @classmethod
    def validate_dna_sequence(cls, v: str, info: ValidationInfo) -> str:
        """Validate that a primer contains only valid DNA bases (incl. IUPAC ambiguity codes).

        Inosine (``I``) is not an IUPAC code, but degenerate primers often use it
        because it pairs with any base. It is replaced by ``N``, which Cutadapt and the
        reverse-complement helper understand, and a ``[WARN]`` records the conversion.

        Args:
            v: A primer sequence string (any casing).
            info: Pydantic validation context; ``info.field_name`` names the primer.

        Returns:
            The sequence upper-cased, with any ``I`` replaced by ``N``.

        Raises:
            ValueError: if the sequence contains any character outside the IUPAC DNA alphabet
                ``ACGTRYMKSWHBVDN`` (plus ``I``); the message lists the offending bases.
        """
        valid_bases = set("ACGTRYMKSWHBVDN")
        v_upper = v.upper()
        if "I" in v_upper:
            converted = v_upper.replace("I", "N")
            if (info.field_name, v_upper) not in _INOSINE_WARNED:
                _INOSINE_WARNED.add((info.field_name, v_upper))
                logger.warning(
                    f"[WARN] primers.{info.field_name}: inosine (I) is not an IUPAC base, "
                    f"treated as N ({v_upper} -> {converted})"
                )
            v_upper = converted
        if not all(base in valid_bases for base in v_upper):
            invalid_bases = set(v_upper) - valid_bases
            raise ValueError(
                f"Invalid DNA sequence. Contains invalid bases: {invalid_bases}. "
                f"Valid bases are: {', '.join(sorted(valid_bases))} "
                f"(inosine I is also accepted and treated as N)"
            )
        return v_upper


class MarkerConfig(StrictModel):
    """Identity of the metabarcoding marker being run (e.g. teleo, 16S, COI).

    A marker is one primer-defined barcode assay; one SeeDNAP config processes one marker.

    Attributes:
        name: Marker name, lowercase (e.g. ``teleo``, ``mifish``).
        description: Optional free-text description of the marker.
        primers: The forward/reverse primer pair for this marker.
    """

    name: str = Field(..., description="Marker name (lowercase)")
    description: Optional[str] = Field(None, description="Marker description")
    primers: PrimerConfig = Field(..., description="Primer pair configuration")


class PathsConfig(StrictModel):
    """Filesystem locations for raw input, outputs, and logs.

    Attributes:
        raw_data: Directory of paired-end input FASTQ files.
        output: Directory the pipeline writes its step outputs under.
        logs: Directory for run log files.
    """

    raw_data: Path = Field(default=Path("data/raw"), description="Raw FASTQ data directory")
    output: Path = Field(default=Path("outputs"), description="Output directory")
    logs: Path = Field(default=Path("logs"), description="Log files directory")

    @field_validator("raw_data", "output", "logs")
    @classmethod
    def expand_path(cls, v: Path) -> Path:
        """Expand ``~`` and resolve a relative path to an absolute path.

        Args:
            v: A configured raw_data / output / logs path.

        Returns:
            The path with ``~`` expanded and resolved to absolute.
        """
        return v.expanduser().resolve()


class DemultiplexConfig(StrictModel):
    """Demultiplexing step config: how to split a multiplexed run into per-sample reads.

    Demultiplexing assigns reads from a pooled sequencing run to individual samples by their
    barcodes (MIDs); it runs only when ``demultiplex`` is listed in ``pipeline.steps``, before
    trimming. Modern datasets are usually already demultiplexed (one sample per FASTQ pair),
    in which case this step is omitted.

    Attributes:
        protocol: Demultiplexing protocol (ligation / standard / none); only ``ligation`` is
            implemented (enforced at the root config).
        metadata: Optional path to the metadata CSV mapping barcodes to samples.
        max_sample_failure_rate: Abort demultiplexing if more than this fraction (0.0-1.0) of
            samples fail; otherwise the failures are logged and the run continues.
    """

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
        """Normalize the metadata path (expand ~ and resolve to absolute).

        Does not check that the file exists; existence is verified at runtime
        since the config may be written before the metadata file is created.

        Args:
            v: The configured metadata path, or None.
            info: Pydantic validation context (unused; present for the validator signature).

        Returns:
            The path with ``~`` expanded and resolved to absolute, or None unchanged.
        """
        if v is not None:
            v = v.expanduser().resolve()
            # Note: we don't check file existence here since config might be created
            # before the file exists. Validation happens at runtime.
        return v
