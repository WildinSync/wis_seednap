"""Tests for BLAST taxonomic assignment."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from seednap.steps.taxonomic_assignment.blast import (
    BlastLCAResolver,
    BlastOutputFormatter,
    BlastPhyloFilter,
    BlastTaxonomicAssigner,
)
from seednap.steps.taxonomic_assignment.blast_runner import BlastDatabaseError, BlastRunner


class TestBlastOutputFormatter:
    """Tests for BlastOutputFormatter class."""

    @pytest.fixture
    def sample_reference_fasta(self, temp_dir: Path) -> Path:
        """Create a sample reference FASTA with phylogeny in headers."""
        fasta_path = temp_dir / "ref_db.fasta"
        with open(fasta_path, "w") as f:
            f.write(">seq1\tAnimalia;Chordata;Actinopteri;Perciformes;Sparidae;Diplodus;Diplodus_sargus\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write(">seq2\tAnimalia;Chordata;Actinopteri;Perciformes;Sparidae;Diplodus;Diplodus_vulgaris\n")
            f.write("GCTAGCTAGCTAGCTA\n")
        return fasta_path

    @pytest.fixture
    def sample_blast_output(self, temp_dir: Path) -> Path:
        """Create a sample BLAST output TSV."""
        tsv_path = temp_dir / "blast_output.tsv"
        # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq
        with open(tsv_path, "w") as f:
            f.write("ASV1\tseq1\t98.5\t100\t1\t0\t1\t100\t1\t100\t1e-50\t200\tATCG\tATCG\n")
            f.write("ASV2\tseq2\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-40\t180\tGCTA\tGCTA\n")
        return tsv_path

    def test_initialization(self, sample_reference_fasta: Path) -> None:
        """Test formatter initialization."""
        formatter = BlastOutputFormatter(sample_reference_fasta)
        assert formatter.reference_fasta == sample_reference_fasta
        assert len(formatter._phylo_dict) == 2
        assert "seq1" in formatter._phylo_dict

    def test_initialization_file_not_found(self, temp_dir: Path) -> None:
        """Test error when reference FASTA doesn't exist."""
        non_existent = temp_dir / "doesnt_exist.fasta"
        with pytest.raises(FileNotFoundError):
            BlastOutputFormatter(non_existent)

    def test_format_blast_output(
        self, sample_reference_fasta: Path, sample_blast_output: Path
    ) -> None:
        """Test formatting BLAST output with phylogeny."""
        formatter = BlastOutputFormatter(sample_reference_fasta)
        result = formatter.format_blast_output(sample_blast_output)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert "qseqid" in result.columns
        assert "kingdom" in result.columns
        assert "species" in result.columns
        assert "blast_rank" in result.columns

        # Check phylogeny was extracted
        assert result.iloc[0]["kingdom"] == "Animalia"
        assert result.iloc[0]["genus"] == "Diplodus"
        assert result.iloc[0]["species"] == "Diplodus_sargus"

    def test_blast_rank_assignment(
        self, sample_reference_fasta: Path, temp_dir: Path
    ) -> None:
        """Test that blast_rank is assigned correctly."""
        # Create TSV with multiple hits for same query
        tsv_path = temp_dir / "multi_hit.tsv"
        with open(tsv_path, "w") as f:
            f.write("ASV1\tseq1\t98.5\t100\t1\t0\t1\t100\t1\t100\t1e-50\t200\tATCG\tATCG\n")
            f.write("ASV1\tseq2\t97.0\t100\t3\t0\t1\t100\t1\t100\t1e-45\t190\tATCG\tATCG\n")

        formatter = BlastOutputFormatter(sample_reference_fasta)
        result = formatter.format_blast_output(tsv_path)

        # Check blast_rank
        asv1_hits = result[result["qseqid"] == "ASV1"]
        assert len(asv1_hits) == 2
        assert list(asv1_hits["blast_rank"].values) == [1, 2]


class TestBlastPhyloFilter:
    """Tests for BlastPhyloFilter class."""

    def test_filter_by_thresholds(self) -> None:
        """Test filtering by percent identity thresholds."""
        df = pd.DataFrame({
            "pident": [99.0, 97.0, 95.0, 85.0],
            "species": ["S1", "S2", "S3", "S4"],
            "genus": ["G1", "G2", "G3", "G4"],
            "family": ["F1", "F2", "F3", "F4"],
        })

        filter = BlastPhyloFilter(
            threshold_species=98.0, threshold_genus=96.0, threshold_family=86.5
        )
        result = filter.filter(df)

        # Species: need >= 98.0
        assert pd.notna(result.iloc[0]["species"])  # 99.0 passes
        assert pd.isna(result.iloc[1]["species"])  # 97.0 fails
        assert pd.isna(result.iloc[2]["species"])  # 95.0 fails

        # Genus: need >= 96.0
        assert pd.notna(result.iloc[1]["genus"])  # 97.0 passes
        assert pd.isna(result.iloc[2]["genus"])  # 95.0 fails

        # Family: need >= 86.5
        assert pd.isna(result.iloc[3]["family"])  # 85.0 fails

    def test_custom_thresholds(self) -> None:
        """Test with custom thresholds."""
        df = pd.DataFrame({
            "pident": [95.0],
            "species": ["S1"],
            "genus": ["G1"],
            "family": ["F1"],
        })

        filter = BlastPhyloFilter(threshold_species=90.0, threshold_genus=90.0, threshold_family=90.0)
        result = filter.filter(df)

        # All should pass with 95.0 pident
        assert pd.notna(result.iloc[0]["species"])
        assert pd.notna(result.iloc[0]["genus"])
        assert pd.notna(result.iloc[0]["family"])


class TestBlastLCAResolver:
    """Tests for BlastLCAResolver class."""

    def test_single_hit_no_ambiguity(self) -> None:
        """Test with single hit (no ambiguity)."""
        df = pd.DataFrame({
            "qseqid": ["ASV1"],
            "bitscore": [200],
            "kingdom": ["Animalia"],
            "phylum": ["Chordata"],
            "class": ["Actinopteri"],
            "order": ["Perciformes"],
            "family": ["Sparidae"],
            "genus": ["Diplodus"],
            "species": ["Diplodus_sargus"],
        })

        resolver = BlastLCAResolver()
        result = resolver.resolve_ambiguous_hits(df)

        assert len(result) == 1
        assert result.iloc[0]["keep_for_analysis"] == True  # noqa: E712

    def test_multiple_hits_same_phylogeny(self) -> None:
        """Test multiple hits with same phylogeny (no conflict)."""
        df = pd.DataFrame({
            "qseqid": ["ASV1", "ASV1"],
            "bitscore": [200, 200],
            "kingdom": ["Animalia", "Animalia"],
            "phylum": ["Chordata", "Chordata"],
            "class": ["Actinopteri", "Actinopteri"],
            "order": ["Perciformes", "Perciformes"],
            "family": ["Sparidae", "Sparidae"],
            "genus": ["Diplodus", "Diplodus"],
            "species": ["Diplodus_sargus", "Diplodus_sargus"],
        })

        resolver = BlastLCAResolver()
        result = resolver.resolve_ambiguous_hits(df)

        # Both should be kept since they agree
        kept = result[result["keep_for_analysis"] == True]  # noqa: E712
        assert len(kept) == 2

    def test_ambiguous_hits_lca_resolution(self) -> None:
        """Test LCA resolution for conflicting hits."""
        df = pd.DataFrame({
            "qseqid": ["ASV1", "ASV1"],
            "bitscore": [200, 200],  # Same score = ambiguous
            "kingdom": ["Animalia", "Animalia"],
            "phylum": ["Chordata", "Chordata"],
            "class": ["Actinopteri", "Actinopteri"],
            "order": ["Perciformes", "Perciformes"],
            "family": ["Sparidae", "Sparidae"],
            "genus": ["Diplodus", "Diplodus"],
            "species": ["Diplodus_sargus", "Diplodus_vulgaris"],  # Different species
        })

        resolver = BlastLCAResolver()
        result = resolver.resolve_ambiguous_hits(df)

        # Should create LCA row
        kept = result[result["keep_for_analysis"] == True]  # noqa: E712
        assert len(kept) == 1

        # LCA should have genus but not species
        lca_row = kept.iloc[0]
        assert lca_row["genus"] == "Diplodus"
        assert pd.isna(lca_row["species"])  # Species differs, so set to None

    def test_different_bitscores_no_lca(self) -> None:
        """Test that different bitscores don't trigger LCA."""
        df = pd.DataFrame({
            "qseqid": ["ASV1", "ASV1"],
            "bitscore": [200, 190],  # Different scores
            "kingdom": ["Animalia", "Animalia"],
            "phylum": ["Chordata", "Chordata"],
            "class": ["Actinopteri", "Actinopteri"],
            "order": ["Perciformes", "Perciformes"],
            "family": ["Sparidae", "Sparidae"],
            "genus": ["Diplodus", "Diplodus"],
            "species": ["Diplodus_sargus", "Diplodus_vulgaris"],
        })

        resolver = BlastLCAResolver()
        result = resolver.resolve_ambiguous_hits(df)

        # Only best hit should be kept
        kept = result[result["keep_for_analysis"] == True]  # noqa: E712
        assert len(kept) == 1
        assert kept.iloc[0]["species"] == "Diplodus_sargus"


class TestBlastRunner:
    """Tests for BlastRunner class."""

    def test_initialization(self) -> None:
        """Test runner initialization with parameters."""
        runner = BlastRunner(perc_identity=95.0, evalue=1e-30)
        assert runner.perc_identity == 95.0
        assert runner.evalue == 1e-30
        assert runner.qcov_hsp_perc == 80.0  # default

    def test_check_blast_db_exists_false(self, temp_dir: Path) -> None:
        """Test database check when files don't exist."""
        fasta_path = temp_dir / "test.fasta"
        fasta_path.touch()

        runner = BlastRunner()
        assert runner.check_blast_db_exists(fasta_path) is False

    def test_check_blast_db_exists_true(self, temp_dir: Path) -> None:
        """Test database check when files exist."""
        fasta_path = temp_dir / "test.fasta"
        fasta_path.touch()

        # Create fake database files
        for ext in [".nhr", ".nin", ".nsq"]:
            (temp_dir / f"test.fasta{ext}").touch()

        runner = BlastRunner()
        assert runner.check_blast_db_exists(fasta_path) is True

    @patch("subprocess.run")
    def test_create_blast_db_success(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test successful database creation."""
        fasta_path = temp_dir / "test.fasta"
        fasta_path.touch()

        mock_run.return_value = MagicMock(stdout="Database created", returncode=0)

        runner = BlastRunner()
        runner.create_blast_db(fasta_path)

        # Check makeblastdb was called
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert "makeblastdb" in args
        assert "-dbtype" in args
        assert "nucl" in args

    @patch("subprocess.run")
    def test_create_blast_db_file_not_found(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test error when FASTA file doesn't exist."""
        fasta_path = temp_dir / "doesnt_exist.fasta"

        runner = BlastRunner()
        with pytest.raises(FileNotFoundError):
            runner.create_blast_db(fasta_path)

    @patch("subprocess.run")
    def test_run_blastn_success(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test successful BLAST search."""
        query_fasta = temp_dir / "query.fasta"
        db_fasta = temp_dir / "db.fasta"
        output_tsv = temp_dir / "output.tsv"

        # Create files
        query_fasta.touch()
        db_fasta.touch()

        # Create fake database files
        for ext in [".nhr", ".nin", ".nsq"]:
            (temp_dir / f"db.fasta{ext}").touch()

        mock_run.return_value = MagicMock(stdout="BLAST completed", returncode=0)

        runner = BlastRunner(perc_identity=85.0)
        runner.run_blastn(query_fasta, db_fasta, output_tsv)

        # Check blastn was called
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert "blastn" in args
        assert "-perc_identity" in args
        assert "85.0" in args


class TestBlastTaxonomicAssigner:
    """Integration tests for complete BLAST taxonomic assignment."""

    @pytest.fixture
    def complete_test_setup(self, temp_dir: Path) -> dict:
        """Create complete test setup with all required files."""
        # Reference FASTA
        ref_fasta = temp_dir / "ref.fasta"
        with open(ref_fasta, "w") as f:
            f.write(">seq1\tAnimalia;Chordata;Actinopteri;Perciformes;Sparidae;Diplodus;Diplodus_sargus\n")
            f.write("ATCGATCGATCG\n")

        # BLAST output
        blast_tsv = temp_dir / "blast.tsv"
        with open(blast_tsv, "w") as f:
            f.write("ASV_1\tseq1\t98.5\t12\t0\t0\t1\t12\t1\t12\t1e-20\t50\tATCG\tATCG\n")

        # ASV count table
        asv_count = temp_dir / "counts.csv"
        pd.DataFrame({
            "ATCGATCGATCG": [100]  # Sequence as column
        }, index=["sample1"]).to_csv(asv_count)

        # ASV FASTA
        asv_fasta = temp_dir / "asvs.fasta"
        with open(asv_fasta, "w") as f:
            f.write(">ASV_1\n")
            f.write("ATCGATCGATCG\n")

        return {
            "ref_fasta": ref_fasta,
            "blast_tsv": blast_tsv,
            "asv_count": asv_count,
            "asv_fasta": asv_fasta,
        }

    def test_complete_assignment(self, complete_test_setup: dict, temp_dir: Path) -> None:
        """Test complete taxonomic assignment workflow."""
        assigner = BlastTaxonomicAssigner(
            reference_fasta=complete_test_setup["ref_fasta"], threshold_species=95.0
        )

        output_path = temp_dir / "final.csv"
        result = assigner.assign_taxonomy(
            blast_tsv=complete_test_setup["blast_tsv"],
            asv_count_csv=complete_test_setup["asv_count"],
            asv_fasta=complete_test_setup["asv_fasta"],
            output_path=output_path,
        )

        # Check result
        assert len(result) > 0
        assert "ASV_ID" in result.columns
        assert "species" in result.columns
        assert "pident" in result.columns

        # Check output file was created
        assert output_path.exists()
