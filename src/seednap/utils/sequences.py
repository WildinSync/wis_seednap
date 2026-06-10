"""Sequence manipulation utilities.

Small Biopython-backed helpers shared across the pipeline for DNA sequence
handling and FASTA <-> DataFrame conversion. ASV/OTU sequences move between
tabular form (counts and taxonomy tables) and FASTA (the input expected by
vsearch, swarm, and blastn), so these helpers bridge the two and provide
reverse-complementing for primer/orientation handling. FASTA I/O is delegated
to Biopython's SeqIO rather than a hand-rolled parser. Lives in seednap/utils/.
"""

from pathlib import Path
from typing import Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def reverse_complement(sequence: str) -> str:
    """
    Calculate the reverse complement of a DNA sequence.

    The reverse complement is the sequence read on the opposite strand in the
    5'->3' direction (each base swapped for its pairing partner and the order
    reversed); needed when a marker's primer or read can arrive in either
    orientation. Uppercases the input and handles IUPAC ambiguity codes (e.g.
    R, Y, N) correctly via Biopython.

    Args:
        sequence: DNA sequence string (may include IUPAC ambiguity codes);
            case-insensitive (uppercased internally).

    Returns:
        The reverse complement of the input sequence, in uppercase.

    Raises:
        TypeError: If `sequence` is not a string (it has no .upper()).

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
    Write a DataFrame of sequences out to a FASTA file.

    Each row becomes one FASTA record (header = id_col, body = seq_col), the
    form the downstream tools (vsearch, swarm, blastn) read. The output
    directory is created if missing.

    Args:
        df: DataFrame containing sequences, one per row.
        output_path: Path to the FASTA file to write.
        id_col: Name of the column holding sequence IDs (default: 'id').
        seq_col: Name of the column holding sequences (default: 'sequence').
        description_col: Optional column whose value becomes each record's
            FASTA description; ignored if None or absent from the DataFrame.

    Returns:
        None. The FASTA file is written to `output_path` as a side effect.

    Raises:
        ValueError: If `id_col` or `seq_col` is missing from the DataFrame.
        IOError: If the output file cannot be written.

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
    Read a FASTA file into a DataFrame.

    The inverse of df_to_fasta: parses ASV/OTU (or reference) sequence records
    into tabular form for merging with count/taxonomy tables. Empty or
    non-FASTA input is rejected with a descriptive error rather than returning
    an empty frame.

    Args:
        fasta_path: Path to the FASTA file to read.
        include_description: Whether to add a 'description' column carrying
            each record's full header line (default: False).

    Returns:
        DataFrame with one row per record and columns 'id', 'sequence', and
        optionally 'description'.

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
