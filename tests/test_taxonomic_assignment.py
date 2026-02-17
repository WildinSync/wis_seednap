"""Tests for taxonomic assignment orchestration (Phase 5)."""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from seednap.steps.taxonomic_assignment import (
    DecipherError,
    DecipherRunner,
    EcotagError,
    EcotagRunner,
    TaxonomicAssigner,
    TaxonomyMethod,
)


class TestEcotagRunner:
    """Tests for EcotagRunner class."""

    @patch("subprocess.run")
    def test_initialization(self, mock_run: MagicMock) -> None:
        """Test runner initialization and OBITools check."""
        # Mock successful command checks
        mock_run.return_value = MagicMock(stdout="", stderr="", returncode=0)

        runner = EcotagRunner(timeout=1800)
        assert runner.timeout == 1800

        # Check that all three commands were verified
        assert mock_run.call_count == 3  # ecotag, obiannotate, obitab

    @patch("subprocess.run")
    def test_initialization_missing_obitools(self, mock_run: MagicMock) -> None:
        """Test error when OBITools is not installed."""
        mock_run.side_effect = FileNotFoundError("ecotag not found")

        with pytest.raises(EcotagError, match="OBITools command 'ecotag' not found"):
            EcotagRunner()

    @patch("subprocess.run")
    def test_run_ecotag(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test ecotag execution."""
        # Mock OBITools check
        mock_run.return_value = MagicMock(stdout="", stderr="", returncode=0)
        runner = EcotagRunner()

        # Create input files
        query_fasta = temp_dir / "query.fasta"
        query_fasta.write_text(">ASV1\nATCG\n")

        taxonomy_db = temp_dir / "taxonomy.db"
        taxonomy_db.write_text("dummy taxonomy")

        reference_db = temp_dir / "reference.db"
        reference_db.write_text("dummy reference")

        output_fasta = temp_dir / "output.fasta"

        # Mock ecotag output
        mock_run.return_value = MagicMock(
            stdout=">ASV1|taxonomy=Animalia;Chordata\nATCG\n", stderr="", returncode=0
        )

        result = runner.run_ecotag(
            query_fasta=query_fasta,
            taxonomy_db=taxonomy_db,
            reference_db=reference_db,
            output_fasta=output_fasta,
        )

        assert result == output_fasta
        assert output_fasta.exists()

        # Check command was called correctly
        call_args = [call[0][0] for call in mock_run.call_args_list]
        # Find the ecotag call (not the version checks)
        ecotag_calls = [args for args in call_args if "ecotag" in args and "-t" in args]
        assert len(ecotag_calls) > 0

    @patch("subprocess.run")
    def test_clean_annotations(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test annotation cleaning with obiannotate."""
        # Mock OBITools check
        mock_run.return_value = MagicMock(stdout="", stderr="", returncode=0)
        runner = EcotagRunner()

        # Create input file
        input_fasta = temp_dir / "input.fasta"
        input_fasta.write_text(">ASV1|tag1=value1|tag2=value2\nATCG\n")

        output_fasta = temp_dir / "cleaned.fasta"

        # Mock obiannotate output
        mock_run.return_value = MagicMock(stdout=">ASV1|taxonomy=X\nATCG\n", returncode=0)

        result = runner.clean_annotations(
            input_fasta=input_fasta,
            output_fasta=output_fasta,
            tags_to_delete=["tag1", "tag2"],
        )

        assert result == output_fasta
        assert output_fasta.exists()

    @patch("subprocess.run")
    def test_convert_to_table(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test FASTA to table conversion with obitab."""
        # Mock OBITools check
        mock_run.return_value = MagicMock(stdout="", stderr="", returncode=0)
        runner = EcotagRunner()

        # Create input file
        input_fasta = temp_dir / "input.fasta"
        input_fasta.write_text(">ASV1|taxonomy=X\nATCG\n")

        output_tsv = temp_dir / "output.tsv"

        # Mock obitab output
        mock_run.return_value = MagicMock(
            stdout="sequence\ttaxonomy\nATCG\tX\n", returncode=0
        )

        result = runner.convert_to_table(input_fasta=input_fasta, output_tsv=output_tsv)

        assert result == output_tsv
        assert output_tsv.exists()

    @patch("subprocess.run")
    def test_run_complete_workflow(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test complete ecotag workflow."""
        # Mock all OBITools commands
        mock_run.return_value = MagicMock(stdout="", stderr="", returncode=0)
        runner = EcotagRunner()

        # Create inputs
        query_fasta = temp_dir / "query.fasta"
        query_fasta.write_text(">ASV1\nATCG\n")

        taxonomy_db = temp_dir / "taxonomy.db"
        taxonomy_db.write_text("dummy")

        reference_db = temp_dir / "reference.db"
        reference_db.write_text("dummy")

        output_dir = temp_dir / "ecotag_out"

        # Mock outputs for each step
        def mock_side_effect(*args, **kwargs):
            cmd = args[0]
            if "ecotag" in cmd and "-t" in cmd:
                return MagicMock(stdout=">ASV1|tax=X\nATCG\n", returncode=0)
            elif "obiannotate" in cmd:
                return MagicMock(stdout=">ASV1|tax=X\nATCG\n", returncode=0)
            elif "obitab" in cmd:
                return MagicMock(stdout="sequence\ttax\nATCG\tX\n", returncode=0)
            else:
                return MagicMock(stdout="", returncode=0)

        mock_run.side_effect = mock_side_effect

        outputs = runner.run_complete_workflow(
            query_fasta=query_fasta,
            taxonomy_db=taxonomy_db,
            reference_db=reference_db,
            output_dir=output_dir,
            marker="test",
        )

        assert "ecotag_fasta" in outputs
        assert "cleaned_fasta" in outputs
        assert "taxonomy_tsv" in outputs
        assert outputs["taxonomy_tsv"].exists()

    def test_link_with_abundance_table(self, temp_dir: Path) -> None:
        """Test linking ecotag taxonomy with abundance table."""
        # Create taxonomy TSV
        taxonomy_tsv = temp_dir / "taxonomy.tsv"
        pd.DataFrame({
            "sequence": ["ATCG", "GCTA"],
            "kingdom": ["Animalia", "Plantae"],
            "species": ["Fish", "Plant"],
        }).to_csv(taxonomy_tsv, sep="\t", index=False)

        # Create abundance CSV
        abundance_csv = temp_dir / "abundance.csv"
        pd.DataFrame({
            "sample1": [100, 50],
            "sample2": [200, 75],
        }, index=["ATCG", "GCTA"]).to_csv(abundance_csv)

        output_csv = temp_dir / "complete.csv"

        # Need to create a runner just for the method (no subprocess calls)
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout="", returncode=0)
            runner = EcotagRunner()

        result = runner.link_with_abundance_table(
            taxonomy_tsv=taxonomy_tsv,
            abundance_csv=abundance_csv,
            output_csv=output_csv,
        )

        assert result.exists()
        df = pd.read_csv(result)
        assert "sequence" in df.columns
        assert "kingdom" in df.columns
        assert "sample1" in df.columns


class TestDecipherRunner:
    """Tests for DecipherRunner class."""

    @patch("subprocess.run")
    def test_initialization(self, mock_run: MagicMock) -> None:
        """Test runner initialization and R/DECIPHER check."""
        # Mock Rscript version check
        mock_run.return_value = MagicMock(stdout="", stderr="R version 4.2.0", returncode=0)

        # Mock DECIPHER package check
        def mock_side_effect(*args, **kwargs):
            if "--version" in args[0]:
                return MagicMock(stdout="", stderr="R version", returncode=0)
            elif "-e" in args[0]:
                return MagicMock(stdout="2.26.0", stderr="", returncode=0)
            return MagicMock(stdout="", returncode=0)

        mock_run.side_effect = mock_side_effect

        runner = DecipherRunner(timeout=3600)
        assert runner.timeout == 3600

    @patch("subprocess.run")
    def test_initialization_no_decipher(self, mock_run: MagicMock) -> None:
        """Test error when DECIPHER is not installed."""
        # Mock Rscript OK
        def mock_side_effect(*args, **kwargs):
            if "--version" in args[0]:
                return MagicMock(stdout="", stderr="R version", returncode=0)
            elif "-e" in args[0]:
                raise subprocess.CalledProcessError(1, "Rscript", stderr="Package not installed")
            return MagicMock(stdout="", returncode=0)

        mock_run.side_effect = mock_side_effect

        with pytest.raises(DecipherError, match="DECIPHER R package not installed"):
            DecipherRunner()

    @patch("subprocess.run")
    def test_run_decipher_assignment(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test DECIPHER taxonomic assignment."""
        # Mock R checks
        def mock_side_effect(*args, **kwargs):
            if "--version" in args[0]:
                return MagicMock(stdout="", stderr="R version", returncode=0)
            elif "-e" in args[0] and "packageVersion" in str(args):
                return MagicMock(stdout="2.26.0", returncode=0)
            else:
                # R script execution
                return MagicMock(stdout="DECIPHER completed", returncode=0)

        mock_run.side_effect = mock_side_effect

        runner = DecipherRunner()

        # Create required files
        output_dir = temp_dir
        marker_dir = output_dir / "02_dada2" / "test"
        marker_dir.mkdir(parents=True)

        seqtab_rds = marker_dir / "seqtab_clean.rds"
        seqtab_rds.write_text("dummy rds")

        trained_classifier = temp_dir / "trained.rds"
        trained_classifier.write_text("dummy classifier")

        # Create expected output files
        taxonomy_csv = marker_dir / "taxo_assigned_decipher.csv"
        taxonomy_csv.write_text("sequence,kingdom\nATCG,Animalia\n")

        complete_csv = output_dir / "test_decipher.csv"
        complete_csv.write_text("sequence,kingdom,sample1\nATCG,Animalia,100\n")

        outputs = runner.run_decipher_assignment(
            marker="test",
            output_dir=output_dir,
            trained_classifier_path=trained_classifier,
        )

        assert "taxonomy" in outputs
        assert "final_table" in outputs


class TestTaxonomicAssigner:
    """Tests for TaxonomicAssigner unified interface."""

    def test_initialization(self, temp_dir: Path) -> None:
        """Test assigner initialization."""
        assigner = TaxonomicAssigner(
            method="blast", marker="teleo", output_dir=temp_dir
        )

        assert assigner.method == TaxonomyMethod.BLAST
        assert assigner.marker == "teleo"
        assert assigner.taxo_dir.exists()

    def test_initialization_with_enum(self, temp_dir: Path) -> None:
        """Test initialization with TaxonomyMethod enum."""
        assigner = TaxonomicAssigner(
            method=TaxonomyMethod.DADA2, marker="amph", output_dir=temp_dir
        )

        assert assigner.method == TaxonomyMethod.DADA2
        assert assigner.marker == "amph"

    def test_initialization_invalid_method(self, temp_dir: Path) -> None:
        """Test error with invalid method."""
        with pytest.raises(ValueError):
            TaxonomicAssigner(method="invalid_method", marker="test", output_dir=temp_dir)

    def test_get_method_requirements(self) -> None:
        """Test getting method requirements."""
        reqs = TaxonomicAssigner.get_method_requirements(TaxonomyMethod.BLAST)
        assert "reference_fasta" in reqs

        reqs = TaxonomicAssigner.get_method_requirements("dada2")
        assert "rdp_db_path" in reqs
        assert "species_db_path" in reqs

    @patch("seednap.steps.taxonomic_assignment.assigner.BlastRunner")
    @patch("seednap.steps.taxonomic_assignment.assigner.BlastTaxonomicAssigner")
    def test_assign_blast(
        self, mock_assigner: MagicMock, mock_runner: MagicMock, temp_dir: Path
    ) -> None:
        """Test BLAST assignment through unified interface."""
        # Create input files
        query_fasta = temp_dir / "query.fasta"
        query_fasta.write_text(">ASV1\nATCG\n")

        asv_count = temp_dir / "counts.csv"
        asv_count.write_text("sequence,sample1\nATCG,100\n")

        reference_fasta = temp_dir / "ref.fasta"
        reference_fasta.write_text(">ref1|Animalia;Chordata\nATCG\n")

        # Mock BLAST runner
        mock_runner_instance = MagicMock()
        mock_runner_instance.check_blast_db_exists.return_value = True
        mock_runner.return_value = mock_runner_instance

        # Mock BLAST assigner
        mock_assigner_instance = MagicMock()
        mock_assigner_instance.assign_taxonomy.return_value = pd.DataFrame({
            "ASV_ID": ["ASV1"],
            "kingdom": ["Animalia"],
        })
        mock_assigner.return_value = mock_assigner_instance

        # Run assignment
        assigner = TaxonomicAssigner(method="blast", marker="test", output_dir=temp_dir)
        outputs = assigner.assign_taxonomy(
            query_fasta=query_fasta,
            asv_count_csv=asv_count,
            reference_fasta=reference_fasta,
        )

        assert "final_table" in outputs
        mock_runner_instance.run_blastn.assert_called_once()
        mock_assigner_instance.assign_taxonomy.assert_called_once()

    def test_assign_blast_missing_reference(self, temp_dir: Path) -> None:
        """Test error when BLAST reference is missing."""
        query_fasta = temp_dir / "query.fasta"
        query_fasta.write_text(">ASV1\nATCG\n")

        asv_count = temp_dir / "counts.csv"
        asv_count.write_text("sequence,sample1\nATCG,100\n")

        assigner = TaxonomicAssigner(method="blast", marker="test", output_dir=temp_dir)

        with pytest.raises(ValueError, match="reference_fasta is required"):
            assigner.assign_taxonomy(
                query_fasta=query_fasta,
                asv_count_csv=asv_count,
                # Missing reference_fasta!
            )

    @patch("seednap.steps.taxonomic_assignment.dada2_taxonomy_runner.Dada2TaxonomyRunner")
    def test_assign_dada2(self, mock_runner: MagicMock, temp_dir: Path) -> None:
        """Test DADA2 assignment through unified interface."""
        query_fasta = temp_dir / "query.fasta"
        query_fasta.write_text(">ASV1\nATCG\n")

        asv_count = temp_dir / "counts.csv"
        asv_count.write_text("sequence,sample1\nATCG,100\n")

        rdp_db = temp_dir / "rdp.fasta"
        rdp_db.write_text(">tax1\nATCG\n")

        species_db = temp_dir / "species.fasta"
        species_db.write_text(">species1\nATCG\n")

        # Mock DADA2 runner
        mock_runner_instance = MagicMock()
        mock_runner_instance.run_dada2_taxonomy.return_value = {
            "taxonomy": temp_dir / "taxonomy.csv",
            "final_table": temp_dir / "complete.csv",
        }
        mock_runner.return_value = mock_runner_instance

        assigner = TaxonomicAssigner(method="dada2", marker="test", output_dir=temp_dir)
        outputs = assigner.assign_taxonomy(
            query_fasta=query_fasta,
            asv_count_csv=asv_count,
            rdp_db_path=rdp_db,
            species_db_path=species_db,
        )

        assert "taxonomy" in outputs
        assert "final_table" in outputs
        mock_runner_instance.run_dada2_taxonomy.assert_called_once()

    @patch("seednap.steps.taxonomic_assignment.assigner.EcotagRunner")
    def test_assign_ecotag(self, mock_runner: MagicMock, temp_dir: Path) -> None:
        """Test ecotag assignment through unified interface."""
        query_fasta = temp_dir / "query.fasta"
        query_fasta.write_text(">ASV1\nATCG\n")

        asv_count = temp_dir / "counts.csv"
        asv_count.write_text("sequence,sample1\nATCG,100\n")

        taxonomy_db = temp_dir / "taxonomy.db"
        taxonomy_db.write_text("dummy")

        reference_db = temp_dir / "reference.db"
        reference_db.write_text("dummy")

        # Mock ecotag runner
        mock_runner_instance = MagicMock()
        mock_runner_instance.run_complete_workflow.return_value = {
            "ecotag_fasta": temp_dir / "ecotag.fasta",
            "cleaned_fasta": temp_dir / "cleaned.fasta",
            "taxonomy_tsv": temp_dir / "taxonomy.tsv",
        }
        mock_runner_instance.link_with_abundance_table.return_value = temp_dir / "complete.csv"
        mock_runner.return_value = mock_runner_instance

        assigner = TaxonomicAssigner(method="ecotag", marker="test", output_dir=temp_dir)
        outputs = assigner.assign_taxonomy(
            query_fasta=query_fasta,
            asv_count_csv=asv_count,
            taxonomy_db=taxonomy_db,
            reference_db=reference_db,
        )

        assert "final_table" in outputs
        mock_runner_instance.run_complete_workflow.assert_called_once()
        mock_runner_instance.link_with_abundance_table.assert_called_once()

    @patch("seednap.steps.taxonomic_assignment.assigner.DecipherRunner")
    def test_assign_decipher(self, mock_runner: MagicMock, temp_dir: Path) -> None:
        """Test DECIPHER assignment through unified interface."""
        query_fasta = temp_dir / "query.fasta"
        query_fasta.write_text(">ASV1\nATCG\n")

        asv_count = temp_dir / "counts.csv"
        asv_count.write_text("sequence,sample1\nATCG,100\n")

        trained_classifier = temp_dir / "trained.rds"
        trained_classifier.write_text("dummy")

        # Mock DECIPHER runner
        mock_runner_instance = MagicMock()
        mock_runner_instance.run_decipher_assignment.return_value = {
            "taxonomy": temp_dir / "taxonomy.csv",
            "final_table": temp_dir / "complete.csv",
        }
        mock_runner.return_value = mock_runner_instance

        assigner = TaxonomicAssigner(method="decipher", marker="test", output_dir=temp_dir)
        outputs = assigner.assign_taxonomy(
            query_fasta=query_fasta,
            asv_count_csv=asv_count,
            trained_classifier_path=trained_classifier,
        )

        assert "taxonomy" in outputs
        assert "final_table" in outputs
        mock_runner_instance.run_decipher_assignment.assert_called_once()


class TestTaxonomyMethod:
    """Tests for TaxonomyMethod enum."""

    def test_enum_values(self) -> None:
        """Test enum values."""
        assert TaxonomyMethod.BLAST.value == "blast"
        assert TaxonomyMethod.DADA2.value == "dada2"
        assert TaxonomyMethod.ECOTAG.value == "ecotag"
        assert TaxonomyMethod.DECIPHER.value == "decipher"

    def test_enum_from_string(self) -> None:
        """Test creating enum from string."""
        method = TaxonomyMethod("blast")
        assert method == TaxonomyMethod.BLAST

        method = TaxonomyMethod("dada2")
        assert method == TaxonomyMethod.DADA2
