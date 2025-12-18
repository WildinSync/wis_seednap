"""Sequence manipulation utilities."""

from Bio.Seq import Seq


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
        'NHVBDWSMKRYCGAT'
    """
    seq = Seq(sequence.upper())
    return str(seq.reverse_complement())
