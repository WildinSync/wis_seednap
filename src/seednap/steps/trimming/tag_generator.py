"""Generate cutadapt tag files from metadata for demultiplexing.

This module generates cutadapt adapter (tag) FASTA files from metadata CSV
files for both standard and ligation-based demultiplexing.
"""

import logging
from pathlib import Path
from typing import Union

import pandas as pd

from seednap.utils.sequences import reverse_complement

logger = logging.getLogger(__name__)


class TagFileGenerator:
    """Generate cutadapt tag files for demultiplexing.

    Creates FASTA files with tag sequences for cutadapt demultiplexing, for
    both standard (per-run) and ligation-based (per-library) layouts.
    """

    def __init__(self, min_overlap: int = 8):
        """
        Initialize tag file generator.

        Args:
            min_overlap: Minimum overlap required for tag matching (default: 8)
        """
        self.min_overlap = min_overlap

    def _format_tag_sequence(self, tag: str) -> str:
        """
        Format tag sequence for cutadapt file adapter specification.

        The format tells cutadapt to match either the tag or its reverse complement
        with a specified minimum overlap.

        Args:
            tag: Tag sequence (DNA)

        Returns:
            Formatted tag string: "TAG;min_overlap=N...RC_TAG;min_overlap=N"
        """
        tag = tag.upper()
        tag_rc = reverse_complement(tag).lower()

        return f"{tag};min_overlap={self.min_overlap}...{tag_rc};min_overlap={self.min_overlap}"

    # IUPAC nucleotide alphabet accepted in tag sequences.
    _IUPAC_DNA = set("ACGTRYSWKMBDHVN")

    def _format_validated_tag(
        self,
        value: object,
        metadata_csv: Union[str, Path],
        sample_name: object,
        group_label: object,
        group_field: str,
        tag_col: str,
    ) -> str:
        """
        Validate one tag cell and format it for cutadapt.

        Raises a row-context error instead of letting an empty/NaN/non-DNA tag
        surface as a bare ``AttributeError`` (e.g. ``float`` has no ``upper``)
        or get silently turned into a corrupt cutadapt adapter.

        Args:
            value: Raw tag value from the metadata cell.
            metadata_csv: Source metadata CSV path (for the error message).
            sample_name: Sample identifier for the offending row.
            group_label: Run/library identifier for the offending row.
            group_field: Name of the grouping field ('run' or 'library').
            tag_col: Original tag column header (for the error message).

        Returns:
            Formatted cutadapt tag string.

        Raises:
            ValueError: If the tag is empty/NaN or contains non-DNA characters.
        """
        # Treat NaN/None and non-string cells as invalid (a NaN float would
        # otherwise raise a context-free AttributeError on .upper()).
        if not isinstance(value, str) or pd.isna(value):
            tag = ""
        else:
            tag = value.strip()
        if not tag or not set(tag.upper()) <= self._IUPAC_DNA:
            raise ValueError(
                f"Invalid tag in {metadata_csv} for sample '{sample_name}' "
                f"({group_field} '{group_label}'): the '{tag_col}' column value "
                f"is {value!r}, which is not a valid DNA/IUPAC sequence. Every "
                f"tag must be a non-empty string of ACGT/IUPAC characters; empty "
                f"cells, whitespace, and non-ACGT characters are not allowed "
                f"(they would otherwise produce a corrupt or degenerate cutadapt "
                f"adapter). Fix this row in the metadata CSV and re-run."
            )
        return self._format_tag_sequence(tag)

    def _write_fasta(self, df: pd.DataFrame, output_path: Union[str, Path]) -> None:
        """
        Write DataFrame to FASTA file for cutadapt.

        Args:
            df: DataFrame with columns 'sample_name' and 'tag_formatted'
            output_path: Path to output FASTA file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            for _, row in df.iterrows():
                f.write(f">{row['sample_name']}\n")
                f.write(f"{row['tag_formatted']}\n")

        logger.info(f"Wrote tag file with {len(df)} samples to {output_path}")

    def generate_standard_tag_files(
        self,
        metadata_csv: Union[str, Path],
        output_dir: Union[str, Path],
        sample_col: str = "sample_name",
        tag_col: str = "tag",
        run_col: str = "run",
    ) -> dict:
        """
        Generate tag files for standard (non-ligation) demultiplexing.

        Creates one tag file per sequencing run/library.

        Expected metadata columns:
        - sample_name (or specified by sample_col): Sample identifier
        - tag (or specified by tag_col): Tag sequence
        - run (or specified by run_col): Sequencing run/library identifier

        Args:
            metadata_csv: Path to metadata CSV file
            output_dir: Directory for output tag files
            sample_col: Name of sample name column (default: 'sample_name')
            tag_col: Name of tag column (default: 'tag')
            run_col: Name of run column (default: 'run')

        Returns:
            Dictionary mapping run names to output file paths

        Raises:
            FileNotFoundError: If metadata CSV doesn't exist
            ValueError: If required columns are missing
        """
        # Read metadata
        metadata_csv = Path(metadata_csv)
        if not metadata_csv.exists():
            raise FileNotFoundError(
                f"Metadata CSV not found: {metadata_csv}\n\n"
                "Ligation demultiplexing needs the sample/tag/library metadata "
                "file. Check the path you set: either demultiplex.metadata in "
                "your marker config (for `seednap run-pipeline`) or the "
                "METADATA_CSV argument to `seednap demultiplex`. It must point "
                "at an existing CSV containing the eventID, tag_demultiplex, "
                "and library columns."
            )

        df = pd.read_csv(metadata_csv)

        # Validate columns
        required_cols = [sample_col, tag_col, run_col]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(
                f"Metadata CSV {metadata_csv} is missing required column(s): "
                f"{missing_cols}. Standard demultiplexing needs columns: "
                f"{required_cols} (sample id, tag sequence, run/library). "
                f"Found columns: {list(df.columns)}. "
                f"Rename your CSV headers to match."
            )

        logger.info(f"Loaded metadata with {len(df)} samples from {metadata_csv}")

        # Standardize column names
        df = df.rename(
            columns={sample_col: "sample_name", tag_col: "tag", run_col: "run"}
        )

        # Format tags
        df["tag_formatted"] = df.apply(
            lambda row: self._format_validated_tag(
                row["tag"],
                metadata_csv,
                row["sample_name"],
                row["run"],
                "run",
                tag_col,
            ),
            axis=1,
        )

        # Split by run/library
        output_dir = Path(output_dir)
        output_files = {}

        for run_name, run_df in df.groupby("run"):
            output_path = output_dir / f"{run_name}.fasta"
            self._write_fasta(run_df[["sample_name", "tag_formatted"]], output_path)
            output_files[run_name] = output_path

        logger.info(f"Generated {len(output_files)} tag files in {output_dir}")
        return output_files

    def generate_ligation_tag_files(
        self,
        metadata_csv: Union[str, Path],
        output_dir: Union[str, Path],
        sample_col: str = "eventID",
        tag_col: str = "tag_demultiplex",
        library_col: str = "library",
    ) -> dict:
        """
        Generate tag files for ligation-based demultiplexing.

        Creates one tag file per library.

        Expected metadata columns:
        - eventID (or specified by sample_col): Sample/event identifier
        - tag_demultiplex (or specified by tag_col): Tag sequence
        - library (or specified by library_col): Library identifier

        Args:
            metadata_csv: Path to metadata CSV file
            output_dir: Directory for output tag files (default: outputs/00_demultiplex_ligation/cutadapt_tags/)
            sample_col: Name of sample column (default: 'eventID')
            tag_col: Name of tag column (default: 'tag_demultiplex')
            library_col: Name of library column (default: 'library')

        Returns:
            Dictionary mapping library names to output file paths

        Raises:
            FileNotFoundError: If metadata CSV doesn't exist
            ValueError: If required columns are missing
        """
        # Read metadata
        metadata_csv = Path(metadata_csv)
        if not metadata_csv.exists():
            raise FileNotFoundError(
                f"Metadata CSV not found: {metadata_csv}\n\n"
                "Ligation demultiplexing needs the sample/tag/library metadata "
                "file. Check the path you set: either demultiplex.metadata in "
                "your marker config (for `seednap run-pipeline`) or the "
                "METADATA_CSV argument to `seednap demultiplex`. It must point "
                "at an existing CSV containing the eventID, tag_demultiplex, "
                "and library columns."
            )

        df = pd.read_csv(metadata_csv)

        # Validate columns
        required_cols = [sample_col, tag_col, library_col]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(
                f"Metadata CSV {metadata_csv} is missing required column(s): "
                f"{missing_cols}. Ligation demultiplexing needs columns: "
                f"{required_cols} (sample id, tag sequence, library). "
                f"Found columns: {list(df.columns)}. "
                f"Rename your CSV headers to match."
            )

        logger.info(f"Loaded ligation metadata with {len(df)} samples from {metadata_csv}")

        # Standardize column names
        df = df.rename(
            columns={
                sample_col: "sample_name",
                tag_col: "tag",
                library_col: "library",
            }
        )

        # Format tags
        df["tag_formatted"] = df.apply(
            lambda row: self._format_validated_tag(
                row["tag"],
                metadata_csv,
                row["sample_name"],
                row["library"],
                "library",
                tag_col,
            ),
            axis=1,
        )

        # Split by library
        output_dir = Path(output_dir)
        output_files = {}

        for library_name, library_df in df.groupby("library"):
            output_path = output_dir / f"{library_name}.fasta"
            self._write_fasta(
                library_df[["sample_name", "tag_formatted"]], output_path
            )
            output_files[library_name] = output_path

        logger.info(f"Generated {len(output_files)} ligation tag files in {output_dir}")
        return output_files
