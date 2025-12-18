"""Tests for sequence utilities."""

import pytest

from seednap.utils.sequences import reverse_complement


class TestReverseComplement:
    """Tests for reverse_complement function."""

    def test_simple_sequence(self) -> None:
        """Test reverse complement of simple DNA sequence."""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("GGGG") == "CCCC"

    def test_longer_sequence(self) -> None:
        """Test reverse complement of longer sequence."""
        seq = "ACACCGCCCGTCACTCT"
        expected = "AGAGTGACGGGCGGTGT"
        assert reverse_complement(seq) == expected

    def test_ambiguous_bases(self) -> None:
        """Test reverse complement with IUPAC ambiguity codes."""
        # R (A or G) -> Y (T or C)
        # Y (C or T) -> R (G or A)
        # M (A or C) -> K (G or T)
        # K (G or T) -> M (A or C)
        # S (G or C) -> S (C or G)
        # W (A or T) -> W (T or A)
        # H (A, C, or T) -> D (A, G, or T)
        # B (C, G, or T) -> V (A, C, or G)
        # V (A, C, or G) -> B (C, G, or T)
        # D (A, G, or T) -> H (A, C, or T)
        # N (any) -> N (any)

        assert reverse_complement("RYMKSWHBVDN") == "NHBVDWSMKRY"

    def test_case_insensitive(self) -> None:
        """Test that lowercase input is handled correctly."""
        assert reverse_complement("atcg") == "CGAT"
        assert reverse_complement("AtCg") == "CGAT"

    def test_palindrome(self) -> None:
        """Test palindromic sequences."""
        # GATC is its own reverse complement
        assert reverse_complement("GATC") == "GATC"

    def test_empty_sequence(self) -> None:
        """Test empty sequence returns empty string."""
        assert reverse_complement("") == ""

    def test_single_base(self) -> None:
        """Test single base reverse complement."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("G") == "C"
        assert reverse_complement("C") == "G"

    @pytest.mark.parametrize(
        "input_seq,expected",
        [
            ("ACGT", "ACGT"),  # Reverse of reverse
            ("TGCA", "TGCA"),
            ("AAATTT", "AAATTT"),
        ],
    )
    def test_double_reverse_complement(self, input_seq: str, expected: str) -> None:
        """Test that reverse complement of reverse complement returns original."""
        rc = reverse_complement(input_seq)
        rc_rc = reverse_complement(rc)
        assert rc_rc == expected
