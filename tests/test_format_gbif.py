"""Tests for GBIF formatting functions."""

from pathlib import Path

import pandas as pd
import pytest

from seednap.steps.format_gbif import (
    add_rank_dada,
    add_taxon_dada,
    format_dada2_to_gbif,
    format_ecotag_to_gbif,
)


class TestAddRankDada:
    """Tests for add_rank_dada function."""

    def test_species_level(self) -> None:
        """Test rank determination for species-level assignment."""
        df = pd.DataFrame({
            "kingdom": ["Animalia"],
            "phylum": ["Chordata"],
            "class": ["Actinopteri"],
            "order": ["Perciformes"],
            "family": ["Sparidae"],
            "genus": ["Diplodus"],
            "species": ["Diplodus_sargus"],
        })

        result = add_rank_dada(df)

        assert result["rank"].iloc[0] == "species"
        assert result["species"].iloc[0] == "Diplodus_sargus"

    def test_genus_level_with_ambiguous_species(self) -> None:
        """Test rank determination when species has '/' indicating ambiguity."""
        df = pd.DataFrame({
            "genus": ["Diplodus"],
            "species": ["Diplodus_sargus/Diplodus_vulgaris"],
            "family": ["Sparidae"],
            "order": ["Perciformes"],
            "class": ["Actinopteri"],
            "phylum": ["Chordata"],
            "kingdom": ["Animalia"],
        })

        result = add_rank_dada(df)

        assert result["rank"].iloc[0] == "genus"
        # Species should be set to NA when rank is genus and contains "/"
        assert pd.isna(result["species"].iloc[0])

    def test_family_level(self) -> None:
        """Test rank determination for family-level assignment."""
        df = pd.DataFrame({
            "family": ["Sparidae"],
            "genus": [pd.NA],
            "species": [pd.NA],
            "order": ["Perciformes"],
            "class": ["Actinopteri"],
            "phylum": ["Chordata"],
            "kingdom": ["Animalia"],
        })

        result = add_rank_dada(df)

        assert result["rank"].iloc[0] == "family"

    def test_higher_level(self) -> None:
        """Test rank determination for higher taxonomic levels."""
        df = pd.DataFrame({
            "order": ["Perciformes"],
            "family": [pd.NA],
            "genus": [pd.NA],
            "species": [pd.NA],
            "class": ["Actinopteri"],
            "phylum": ["Chordata"],
            "kingdom": ["Animalia"],
        })

        result = add_rank_dada(df)

        assert result["rank"].iloc[0] == "higher"


class TestAddTaxonDada:
    """Tests for add_taxon_dada function."""

    def test_species_level_taxon(self) -> None:
        """Test taxon extraction at species level."""
        df = pd.DataFrame({
            "rank": ["species"],
            "species": ["Diplodus_sargus"],
            "genus": ["Diplodus"],
            "family": ["Sparidae"],
        })

        result = add_taxon_dada(df)

        assert result["taxon"].iloc[0] == "Diplodus_sargus"

    def test_genus_level_taxon(self) -> None:
        """Test taxon extraction at genus level."""
        df = pd.DataFrame({"rank": ["genus"], "genus": ["Diplodus"], "family": ["Sparidae"]})

        result = add_taxon_dada(df)

        assert result["taxon"].iloc[0] == "Diplodus"

    def test_family_level_taxon(self) -> None:
        """Test taxon extraction at family level."""
        df = pd.DataFrame({
            "rank": ["family"],
            "family": ["Sparidae"],
            "order": ["Perciformes"],
        })

        result = add_taxon_dada(df)

        assert result["taxon"].iloc[0] == "Sparidae"

    def test_higher_level_taxon(self) -> None:
        """Test taxon extraction at higher taxonomic levels."""
        df = pd.DataFrame({
            "rank": ["higher"],
            "order": ["Perciformes"],
            "class": ["Actinopteri"],
            "phylum": ["Chordata"],
        })

        result = add_taxon_dada(df)

        assert result["taxon"].iloc[0] == "Perciformes"


class TestFormatDada2ToGbif:
    """Tests for format_dada2_to_gbif function."""

    @pytest.fixture
    def sample_dada2_output(self, temp_dir: Path) -> Path:
        """Create a sample DADA2 output CSV file."""
        df = pd.DataFrame({
            "kingdom": ["Animalia", "Animalia"],
            "phylum": ["Chordata", "Chordata"],
            "class": ["Actinopteri", "Actinopteri"],
            "order": ["Perciformes", "Perciformes"],
            "family": ["Sparidae", "Sparidae"],
            "genus": ["Diplodus", "Diplodus"],
            "species": ["Diplodus_sargus", "Diplodus_sargus/Diplodus_vulgaris"],
            "sequence": ["ATCGATCGATCG", "GCTAGCTAGCTA"],
            "sample1": [100, 0],
            "sample2": [0, 50],
            "sample3": [25, 30],
        })

        csv_path = temp_dir / "dada2_output.csv"
        df.to_csv(csv_path, index=False)
        return csv_path

    def test_basic_conversion(self, sample_dada2_output: Path, temp_dir: Path) -> None:
        """Test basic DADA2 to GBIF conversion."""
        output_path = temp_dir / "gbif_output.csv"

        result = format_dada2_to_gbif(sample_dada2_output, output_path)

        # Check output file was created
        assert output_path.exists()

        # Check DataFrame structure
        assert "eventID" in result.columns
        assert "rank" in result.columns
        assert "taxon" in result.columns
        assert "nb_reads" in result.columns

        # Check all zero counts are filtered out
        assert (result["nb_reads"] > 0).all()

        # Check total records (sample1: 2 non-zero, sample2: 1, sample3: 2) = 5 records
        assert len(result) == 4  # 100, 50, 25, 30 are non-zero

    def test_rank_assignment(self, sample_dada2_output: Path) -> None:
        """Test that ranks are correctly assigned."""
        result = format_dada2_to_gbif(sample_dada2_output)

        # First species has no "/" so should be species level
        species_level = result[result["sequence"] == "ATCGATCGATCG"]
        assert (species_level["rank"] == "species").all()

        # Second species has "/" so should be genus level
        genus_level = result[result["sequence"] == "GCTAGCTAGCTA"]
        assert (genus_level["rank"] == "genus").all()

    def test_taxon_extraction(self, sample_dada2_output: Path) -> None:
        """Test that taxon column contains correct values."""
        result = format_dada2_to_gbif(sample_dada2_output)

        # Species level should have species name as taxon
        species_level = result[(result["sequence"] == "ATCGATCGATCG") & (result["eventID"] == "sample1")]
        assert species_level["taxon"].iloc[0] == "Diplodus_sargus"

        # Genus level should have genus name as taxon
        genus_level = result[(result["sequence"] == "GCTAGCTAGCTA") & (result["eventID"] == "sample2")]
        assert genus_level["taxon"].iloc[0] == "Diplodus"

    def test_eventid_renaming(self, sample_dada2_output: Path) -> None:
        """Test that filter_code is renamed to eventID."""
        result = format_dada2_to_gbif(sample_dada2_output)

        assert "eventID" in result.columns
        assert "filter_code" not in result.columns
        assert set(result["eventID"].unique()) == {"sample1", "sample2", "sample3"}

    def test_removes_x_column(self, temp_dir: Path) -> None:
        """Test that R's X index column is removed."""
        df = pd.DataFrame({
            "X": [0, 1],
            "kingdom": ["Animalia", "Animalia"],
            "phylum": ["Chordata", "Chordata"],
            "class": ["A", "A"],
            "order": ["O", "O"],
            "family": ["F", "F"],
            "genus": ["G", "G"],
            "species": ["S1", "S2"],
            "sequence": ["ATCG", "GCTA"],
            "sample1": [10, 20],
        })

        csv_path = temp_dir / "with_x.csv"
        df.to_csv(csv_path, index=False)

        result = format_dada2_to_gbif(csv_path)

        assert "X" not in result.columns

    def test_file_not_found(self, temp_dir: Path) -> None:
        """Test error handling for missing file."""
        non_existent = temp_dir / "doesnt_exist.csv"

        with pytest.raises(FileNotFoundError):
            format_dada2_to_gbif(non_existent)

    def test_missing_required_columns(self, temp_dir: Path) -> None:
        """Test error handling for missing required columns."""
        df = pd.DataFrame({"kingdom": ["Animalia"], "sample1": [10]})

        csv_path = temp_dir / "incomplete.csv"
        df.to_csv(csv_path, index=False)

        with pytest.raises(ValueError) as exc_info:
            format_dada2_to_gbif(csv_path)

        assert "Missing required columns" in str(exc_info.value)


class TestFormatEcotagToGbif:
    """Tests for format_ecotag_to_gbif function."""

    @pytest.fixture
    def sample_ecotag_output(self, temp_dir: Path) -> Path:
        """Create a sample Ecotag output CSV file."""
        df = pd.DataFrame({
            "id": ["seq1", "seq2"],
            "definition": ["def1", "def2"],
            "count": [100, 50],
            "scientific_name": ["Name1", "Name2"],
            "order_name": ["Perciformes", "Perciformes"],
            "family_name": ["Sparidae", "Sparidae"],
            "genus_name": ["Diplodus", "Diplodus"],
            "species_name": ["Diplodus_sargus", "Diplodus_vulgaris"],
            "sequence": ["ATCGATCGATCG", "GCTAGCTAGCTA"],
            "rank": ["species", "species"],
            "best_identity": [0.99, 0.98],
            "sample1": [100, 0],
            "sample2": [0, 50],
        })

        csv_path = temp_dir / "ecotag_output.csv"
        df.to_csv(csv_path, index=False)
        return csv_path

    def test_column_renaming(self, sample_ecotag_output: Path) -> None:
        """Test that Ecotag columns are renamed to standard names."""
        result = format_ecotag_to_gbif(sample_ecotag_output)

        assert "family" in result.columns
        assert "genus" in result.columns
        assert "species" in result.columns
        assert "order" in result.columns
        assert "family_name" not in result.columns
        assert "genus_name" not in result.columns

    def test_adds_placeholder_columns(self, sample_ecotag_output: Path) -> None:
        """Test that kingdom, phylum, class columns are added."""
        result = format_ecotag_to_gbif(sample_ecotag_output)

        assert "kingdom" in result.columns
        assert "phylum" in result.columns
        assert "class" in result.columns

    def test_removes_ecotag_metadata(self, sample_ecotag_output: Path) -> None:
        """Test that Ecotag-specific metadata columns are removed."""
        result = format_ecotag_to_gbif(sample_ecotag_output)

        assert "id" not in result.columns
        assert "definition" not in result.columns
        assert "count" not in result.columns
        assert "scientific_name" not in result.columns
        assert "best_identity" not in result.columns

    def test_basic_conversion(self, sample_ecotag_output: Path, temp_dir: Path) -> None:
        """Test basic Ecotag to GBIF conversion."""
        output_path = temp_dir / "ecotag_gbif.csv"

        result = format_ecotag_to_gbif(sample_ecotag_output, output_path)

        # Check output file was created
        assert output_path.exists()

        # Check structure
        assert "eventID" in result.columns
        assert "nb_reads" in result.columns

        # Check zero counts filtered
        assert (result["nb_reads"] > 0).all()

        # Check number of records (2 non-zero values)
        assert len(result) == 2
