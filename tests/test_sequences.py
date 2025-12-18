"""Tests for sequence utilities."""

from pathlib import Path

import pandas as pd
import pytest

from seednap.utils.sequences import df_to_fasta, fasta_to_df, reverse_complement


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


class TestDfToFasta:
    """Tests for df_to_fasta function."""

    def test_basic_conversion(self, temp_dir: Path) -> None:
        """Test basic DataFrame to FASTA conversion."""
        df = pd.DataFrame({"id": ["seq1", "seq2", "seq3"], "sequence": ["ATCG", "GCTA", "TTAA"]})

        output_path = temp_dir / "test.fasta"
        df_to_fasta(df, output_path)

        assert output_path.exists()

        # Read back and verify
        with open(output_path) as f:
            content = f.read()

        assert ">seq1" in content
        assert "ATCG" in content
        assert ">seq2" in content
        assert "GCTA" in content
        assert ">seq3" in content
        assert "TTAA" in content

    def test_with_description(self, temp_dir: Path) -> None:
        """Test FASTA conversion with description column."""
        df = pd.DataFrame({
            "id": ["seq1"],
            "sequence": ["ATCG"],
            "description": ["Test sequence"]
        })

        output_path = temp_dir / "test_desc.fasta"
        df_to_fasta(df, output_path, description_col="description")

        assert output_path.exists()

        with open(output_path) as f:
            content = f.read()

        assert ">seq1 Test sequence" in content

    def test_custom_column_names(self, temp_dir: Path) -> None:
        """Test with custom column names."""
        df = pd.DataFrame({"asv_id": ["asv1"], "seq": ["ATCG"]})

        output_path = temp_dir / "custom.fasta"
        df_to_fasta(df, output_path, id_col="asv_id", seq_col="seq")

        assert output_path.exists()

    def test_missing_id_column(self, temp_dir: Path) -> None:
        """Test error when ID column is missing."""
        df = pd.DataFrame({"sequence": ["ATCG"]})

        output_path = temp_dir / "test.fasta"

        with pytest.raises(ValueError) as exc_info:
            df_to_fasta(df, output_path)

        assert "ID column 'id' not found" in str(exc_info.value)

    def test_missing_sequence_column(self, temp_dir: Path) -> None:
        """Test error when sequence column is missing."""
        df = pd.DataFrame({"id": ["seq1"]})

        output_path = temp_dir / "test.fasta"

        with pytest.raises(ValueError) as exc_info:
            df_to_fasta(df, output_path)

        assert "Sequence column 'sequence' not found" in str(exc_info.value)

    def test_creates_output_directory(self, temp_dir: Path) -> None:
        """Test that output directory is created if it doesn't exist."""
        output_path = temp_dir / "subdir" / "test.fasta"
        df = pd.DataFrame({"id": ["seq1"], "sequence": ["ATCG"]})

        df_to_fasta(df, output_path)

        assert output_path.exists()
        assert output_path.parent.exists()


class TestFastaToDf:
    """Tests for fasta_to_df function."""

    def test_basic_conversion(self, sample_fasta: Path) -> None:
        """Test basic FASTA to DataFrame conversion."""
        df = fasta_to_df(sample_fasta)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert list(df.columns) == ["id", "sequence"]
        assert df["id"].tolist() == ["seq1", "seq2"]
        assert "ACGTACGTACGTACGT" in df["sequence"].values
        assert "TGCATGCATGCATGCA" in df["sequence"].values

    def test_with_description(self, temp_dir: Path) -> None:
        """Test reading FASTA with descriptions."""
        fasta_path = temp_dir / "with_desc.fasta"
        with open(fasta_path, "w") as f:
            f.write(">seq1 This is a description\n")
            f.write("ATCG\n")
            f.write(">seq2 Another description\n")
            f.write("GCTA\n")

        df = fasta_to_df(fasta_path, include_description=True)

        assert "description" in df.columns
        assert len(df) == 2
        assert "seq1 This is a description" in df["description"].values

    def test_file_not_found(self, temp_dir: Path) -> None:
        """Test error when FASTA file doesn't exist."""
        non_existent = temp_dir / "doesnt_exist.fasta"

        with pytest.raises(FileNotFoundError) as exc_info:
            fasta_to_df(non_existent)

        assert "not found" in str(exc_info.value).lower()

    def test_empty_fasta(self, temp_dir: Path) -> None:
        """Test error when FASTA file is empty."""
        empty_fasta = temp_dir / "empty.fasta"
        empty_fasta.touch()

        with pytest.raises(ValueError) as exc_info:
            fasta_to_df(empty_fasta)

        assert "empty" in str(exc_info.value).lower()

    def test_roundtrip(self, temp_dir: Path) -> None:
        """Test that df_to_fasta and fasta_to_df are inverse operations."""
        original_df = pd.DataFrame({
            "id": ["seq1", "seq2", "seq3"],
            "sequence": ["ATCG", "GCTA", "TTAA"]
        })

        fasta_path = temp_dir / "roundtrip.fasta"

        # Write to FASTA
        df_to_fasta(original_df, fasta_path)

        # Read back
        result_df = fasta_to_df(fasta_path)

        # Compare
        pd.testing.assert_frame_equal(original_df, result_df)
