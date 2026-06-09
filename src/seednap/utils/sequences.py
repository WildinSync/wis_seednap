"""Sequence manipulation utilities."""

from pathlib import Path
from typing import Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def reverse_complement(sequence: str) -> str:
    """
    Calculate reverse complement of a DNA sequence.

    Handles IUPAC ambiguity codes correctly.

    Args:
        sequence: DNA sequence string (can include IUPAC ambiguity codes)

    Returns:
        Reverse complement of the input sequence

    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("ATCGRYMKSWHBVDN")
        'NHBVDWSMKRYCGAT'
    """
    seq = Seq(sequence.upper())
    return str(seq.reverse_complement())


def df_to_fasta(
    df: pd.DataFrame,
    output_path: Union[str, Path],
    id_col: str = "id",
    seq_col: str = "sequence",
    description_col: Union[str, None] = None,
) -> None:
    """
    Convert DataFrame to FASTA file.

    Args:
        df: DataFrame containing sequences
        output_path: Path to output FASTA file
        id_col: Name of column containing sequence IDs (default: 'id')
        seq_col: Name of column containing sequences (default: 'sequence')
        description_col: Optional column name for sequence descriptions

    Raises:
        ValueError: If required columns are missing from DataFrame
        IOError: If output file cannot be written

    Examples:
        >>> df = pd.DataFrame({'id': ['seq1', 'seq2'], 'sequence': ['ATCG', 'GCTA']})
        >>> df_to_fasta(df, 'output.fasta')
    """
    # Validate columns
    if id_col not in df.columns:
        raise ValueError(f"ID column '{id_col}' not found in DataFrame")
    if seq_col not in df.columns:
        raise ValueError(f"Sequence column '{seq_col}' not found in DataFrame")

    # Convert to Path
    output_path = Path(output_path)

    # Create output directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create SeqRecord objects
    records = []
    for _, row in df.iterrows():
        seq_id = str(row[id_col])
        sequence = str(row[seq_col])

        # Get description if column specified
        description = ""
        if description_col and description_col in df.columns:
            description = str(row[description_col])

        record = SeqRecord(Seq(sequence), id=seq_id, description=description)
        records.append(record)

    # Write to FASTA
    try:
        with open(output_path, "w") as f:
            SeqIO.write(records, f, "fasta")
    except IOError as e:
        raise IOError(f"Failed to write FASTA file to {output_path}: {e}") from e


def fasta_to_df(fasta_path: Union[str, Path], include_description: bool = False) -> pd.DataFrame:
    """
    Read FASTA file into DataFrame.

    Args:
        fasta_path: Path to FASTA file
        include_description: Whether to include description column (default: False)

    Returns:
        DataFrame with columns: 'id', 'sequence', and optionally 'description'

    Raises:
        FileNotFoundError: If FASTA file does not exist
        ValueError: If FASTA file is empty or invalid

    Examples:
        >>> df = fasta_to_df('sequences.fasta')
        >>> df.columns
        Index(['id', 'sequence'], dtype='object')
    """
    # Convert to Path
    fasta_path = Path(fasta_path)

    # Check file exists
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    # Parse FASTA
    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
    except Exception as e:
        raise ValueError(
            f"Failed to parse FASTA file {fasta_path}: {e}. "
            f"This file is not valid FASTA: a FASTA file has records that begin "
            f"with a '>' header line followed by sequence lines. The path most "
            f"likely points at the wrong file (e.g. a CSV, a plain-text or HTML "
            f"file, or a still-compressed .fasta.gz), or at a non-FASTA artifact "
            f"produced by an earlier pipeline step. Confirm this path is the "
            f"ASV/OTU sequence FASTA from the DADA2/SWARM step (decompress it "
            f"first if it is gzipped), then re-run."
        ) from e

    # Check file not empty
    if len(records) == 0:
        raise ValueError(f"FASTA file is empty: {fasta_path}")

    # Convert to DataFrame
    data = {"id": [record.id for record in records], "sequence": [str(record.seq) for record in records]}

    if include_description:
        data["description"] = [record.description for record in records]

    return pd.DataFrame(data)
