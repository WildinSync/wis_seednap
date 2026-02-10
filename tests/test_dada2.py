"""Tests for DADA2 processing and metrics collection."""

import gzip
import json
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from seednap.steps.dada2 import (
    ASVMetrics,
    Dada2Error,
    Dada2Processor,
    Dada2Runner,
    MetricsCollector,
    ReadMetrics,
    SampleMetrics,
)


class TestDada2Runner:
    """Tests for Dada2Runner class."""

    @patch("subprocess.run")
    def test_initialization(self, mock_run: MagicMock) -> None:
        """Test runner initialization and R check."""
        # Mock successful Rscript check
        mock_run.return_value = MagicMock(
            stdout="", stderr="R scripting front-end version 4.2.0", returncode=0
        )

        runner = Dada2Runner(timeout=3600)
        assert runner.timeout == 3600

        # Check that Rscript --version was called
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert args[0] == "Rscript"
        assert args[1] == "--version"

    @patch("subprocess.run")
    def test_initialization_no_r(self, mock_run: MagicMock) -> None:
        """Test error when R is not installed."""
        mock_run.side_effect = FileNotFoundError("Rscript not found")

        with pytest.raises(Dada2Error, match="Rscript not found"):
            Dada2Runner()

    @patch("subprocess.run")
    def test_run_r_script_success(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test successful R script execution."""
        # Mock Rscript check
        mock_run.return_value = MagicMock(stdout="", stderr="R version 4.2.0", returncode=0)
        runner = Dada2Runner()

        # Create dummy R script
        script = temp_dir / "test.R"
        script.write_text("cat('Hello from R')")

        # Mock script execution
        mock_run.return_value = MagicMock(
            stdout="Hello from R", stderr="", returncode=0
        )

        result = runner._run_r_script(script_path=script, args=["arg1", "arg2"])

        assert result == "Hello from R"
        # Check command structure
        call_args = mock_run.call_args[0][0]
        assert call_args[0] == "Rscript"
        assert str(script) in call_args
        assert "arg1" in call_args
        assert "arg2" in call_args

    @patch("subprocess.run")
    def test_run_r_script_failure(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test R script execution failure."""
        # Mock Rscript check
        mock_run.return_value = MagicMock(stdout="", stderr="R version 4.2.0", returncode=0)
        runner = Dada2Runner()

        script = temp_dir / "failing_script.R"
        script.write_text("stop('Error!')")

        # Mock script failure
        mock_run.side_effect = subprocess.CalledProcessError(
            1, "Rscript", stderr="Error in script"
        )

        with pytest.raises(Dada2Error, match="R script failed"):
            runner._run_r_script(script_path=script, args=[])

    @patch("subprocess.run")
    def test_run_r_script_with_log(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test R script execution with log file."""
        # Mock Rscript check
        mock_run.return_value = MagicMock(stdout="", stderr="R version 4.2.0", returncode=0)
        runner = Dada2Runner()

        script = temp_dir / "test.R"
        script.write_text("cat('Test')")
        log_file = temp_dir / "test.log"

        # Mock script execution
        mock_run.return_value = MagicMock(
            stdout="Test output", stderr="Test stderr", returncode=0
        )

        runner._run_r_script(script_path=script, args=[], log_file=log_file)

        # Check log file was created
        assert log_file.exists()
        log_content = log_file.read_text()
        assert "Test output" in log_content
        assert "Test stderr" in log_content

    @patch("subprocess.run")
    def test_check_dada2_packages(self, mock_run: MagicMock) -> None:
        """Test checking DADA2 package versions."""
        # Mock Rscript check
        mock_run.return_value = MagicMock(stdout="", stderr="R version 4.2.0", returncode=0)
        runner = Dada2Runner()

        # Mock package check
        package_output = (
            "dada2:1.26.0\n"
            "Biostrings:2.66.0\n"
            "DECIPHER:2.26.0\n"
            "dplyr:1.1.0\n"
            "ggplot2:3.4.0\n"
            "patchwork:1.1.2\n"
        )
        mock_run.return_value = MagicMock(stdout=package_output, stderr="", returncode=0)

        versions = runner.check_dada2_packages()

        assert versions["dada2"] == "1.26.0"
        assert versions["Biostrings"] == "2.66.0"
        assert versions["DECIPHER"] == "2.26.0"

    @patch("subprocess.run")
    def test_check_dada2_packages_missing(self, mock_run: MagicMock) -> None:
        """Test error when packages are missing."""
        # Mock Rscript check
        mock_run.return_value = MagicMock(stdout="", stderr="R version 4.2.0", returncode=0)
        runner = Dada2Runner()

        # Mock missing package
        package_output = "dada2:NOT_INSTALLED\n"
        mock_run.return_value = MagicMock(stdout=package_output, stderr="", returncode=0)

        with pytest.raises(Dada2Error, match="Required R packages not installed"):
            runner.check_dada2_packages()


class TestReadMetrics:
    """Tests for ReadMetrics dataclass."""

    def test_initialization(self) -> None:
        """Test metrics initialization."""
        metrics = ReadMetrics()
        assert metrics.raw_reads == 0
        assert metrics.trimmed_reads == 0
        assert metrics.filtered_reads == 0

    def test_to_dict(self) -> None:
        """Test conversion to dictionary."""
        metrics = ReadMetrics(raw_reads=1000, trimmed_reads=900, filtered_reads=800)
        d = metrics.to_dict()

        assert d["raw_reads"] == 1000
        assert d["trimmed_reads"] == 900
        assert d["filtered_reads"] == 800

    def test_retention_rates(self) -> None:
        """Test retention rate calculation."""
        metrics = ReadMetrics(
            raw_reads=1000,
            trimmed_reads=900,
            filtered_reads=800,
            denoised_reads=750,
            merged_reads=700,
            non_chimeric_reads=650,
        )

        rates = metrics.get_retention_rates()

        assert rates["trimming"] == 90.0
        assert rates["filtering"] == 80.0
        assert rates["denoising"] == 75.0
        assert rates["merging"] == 70.0
        assert rates["chimera_removal"] == 65.0

    def test_retention_rates_zero_raw(self) -> None:
        """Test retention rates with zero raw reads."""
        metrics = ReadMetrics()
        rates = metrics.get_retention_rates()
        assert rates == {}


class TestASVMetrics:
    """Tests for ASVMetrics dataclass."""

    def test_initialization(self) -> None:
        """Test ASV metrics initialization."""
        metrics = ASVMetrics()
        assert metrics.num_asvs == 0
        assert metrics.num_samples == 0

    def test_to_dict(self) -> None:
        """Test conversion to dictionary."""
        metrics = ASVMetrics(
            num_asvs=100,
            num_samples=50,
            total_abundance=10000,
            min_length=200,
            max_length=300,
            mean_length=250.5,
        )

        d = metrics.to_dict()
        assert d["num_asvs"] == 100
        assert d["mean_length"] == 250.5


class TestSampleMetrics:
    """Tests for SampleMetrics dataclass."""

    def test_initialization(self) -> None:
        """Test sample metrics initialization."""
        metrics = SampleMetrics(sample_name="sample1")
        assert metrics.sample_name == "sample1"
        assert metrics.num_asvs == 0

    def test_to_dict(self) -> None:
        """Test conversion to dictionary."""
        read_metrics = ReadMetrics(raw_reads=1000, trimmed_reads=900)
        metrics = SampleMetrics(
            sample_name="sample1", reads=read_metrics, num_asvs=50, total_abundance=5000
        )

        d = metrics.to_dict()
        assert d["sample_name"] == "sample1"
        assert d["num_asvs"] == 50
        assert d["reads"]["raw_reads"] == 1000


class TestMetricsCollector:
    """Tests for MetricsCollector class."""

    def test_initialization(self, temp_dir: Path) -> None:
        """Test collector initialization."""
        collector = MetricsCollector(marker="teleo", output_dir=temp_dir)

        assert collector.marker == "teleo"
        assert collector.metrics_dir.exists()
        assert collector.metrics_dir.name == "metrics"

    def test_count_fastq_reads(self, temp_dir: Path) -> None:
        """Test FASTQ read counting."""
        # Create test FASTQ file
        fastq = temp_dir / "test.fastq"
        fastq.write_text(
            "@read1\n"
            "ATCG\n"
            "+\n"
            "IIII\n"
            "@read2\n"
            "GCTA\n"
            "+\n"
            "IIII\n"
        )

        collector = MetricsCollector(marker="test", output_dir=temp_dir)
        count = collector.count_fastq_reads(fastq)

        assert count == 2

    def test_count_fastq_reads_gzipped(self, temp_dir: Path) -> None:
        """Test counting reads in gzipped FASTQ."""
        fastq_gz = temp_dir / "test.fastq.gz"

        with gzip.open(fastq_gz, "wt") as f:
            f.write("@read1\nATCG\n+\nIIII\n")
            f.write("@read2\nGCTA\n+\nIIII\n")

        collector = MetricsCollector(marker="test", output_dir=temp_dir)
        count = collector.count_fastq_reads(fastq_gz)

        assert count == 2

    def test_count_fastq_reads_missing_file(self, temp_dir: Path) -> None:
        """Test counting reads from nonexistent file."""
        collector = MetricsCollector(marker="test", output_dir=temp_dir)
        count = collector.count_fastq_reads(temp_dir / "nonexistent.fastq")

        assert count == 0

    def test_collect_asv_metrics(self, temp_dir: Path) -> None:
        """Test ASV metrics collection."""
        # Create sequence table (transposed)
        seqtab = temp_dir / "seqtab_clean_t.csv"
        df = pd.DataFrame({
            "sample1": [100, 50, 25],
            "sample2": [200, 75, 30],
        }, index=["ATCGATCG", "GCTAGCTA", "TTAATTAA"]) 
        df.to_csv(seqtab)

        # Create correspondence file
        corresp = temp_dir / "corresp_seq.csv"
        corresp_df = pd.DataFrame({
            "ASV_n": ["ASV1", "ASV2", "ASV3"],
            "sequence": ["ATCGATCG", "GCTAGCTA", "TTAATTAA"],
        })
        corresp_df.to_csv(corresp, index=False)

        collector = MetricsCollector(marker="test", output_dir=temp_dir)
        collector.collect_asv_metrics(seqtab_path=seqtab, corresp_seq_path=corresp)

        assert collector.asv_metrics.num_asvs == 3
        assert collector.asv_metrics.num_samples == 2
        assert collector.asv_metrics.total_abundance == 480
        assert collector.asv_metrics.min_length == 8
        assert collector.asv_metrics.max_length == 8

    def test_generate_summary_report(self, temp_dir: Path) -> None:
        """Test summary report generation."""
        collector = MetricsCollector(marker="teleo", output_dir=temp_dir)
        collector.read_metrics.raw_reads = 10000
        collector.read_metrics.trimmed_reads = 9000
        collector.read_metrics.filtered_reads = 8000
        collector.read_metrics.non_chimeric_reads = 7000
        collector.asv_metrics.num_asvs = 150
        collector.asv_metrics.num_samples = 20

        report = collector.generate_summary_report()

        assert "TELEO" in report
        assert "10,000" in report
        assert "9,000" in report
        assert "150" in report
        assert "90.0%" in report  # Trimming retention rate

    def test_export_to_json(self, temp_dir: Path) -> None:
        """Test metrics export to JSON."""
        collector = MetricsCollector(marker="test", output_dir=temp_dir)
        collector.read_metrics.raw_reads = 1000
        collector.asv_metrics.num_asvs = 50

        json_path = collector.export_to_json()

        assert json_path.exists()
        with open(json_path) as f:
            data = json.load(f)

        assert data["marker"] == "test"
        assert data["read_metrics"]["raw_reads"] == 1000
        assert data["asv_metrics"]["num_asvs"] == 50

    def test_export_to_csv(self, temp_dir: Path) -> None:
        """Test metrics export to CSV."""
        collector = MetricsCollector(marker="test", output_dir=temp_dir)
        collector.read_metrics.raw_reads = 1000
        collector.asv_metrics.num_asvs = 50

        csv_path = collector.export_to_csv()

        assert csv_path.exists()
        df = pd.read_csv(csv_path)

        assert "metric" in df.columns
        assert "value" in df.columns
        assert "category" in df.columns
        assert len(df) > 0


class TestDada2Processor:
    """Tests for Dada2Processor class."""

    def test_initialization(self, temp_dir: Path) -> None:
        """Test processor initialization."""
        # Create trimmed reads directory
        trimmed_dir = temp_dir / "trimmed"
        trimmed_dir.mkdir()

        processor = Dada2Processor(
            marker="teleo",
            trimmed_reads_dir=trimmed_dir,
            output_base_dir=temp_dir,
        )

        assert processor.marker == "teleo"
        assert processor.trimmed_reads_dir == trimmed_dir
        assert processor.output_dir.exists()

    def test_initialization_missing_dir(self, temp_dir: Path) -> None:
        """Test error when trimmed reads directory doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Trimmed reads directory not found"):
            Dada2Processor(
                marker="test",
                trimmed_reads_dir=temp_dir / "nonexistent",
                output_base_dir=temp_dir,
            )

    @patch("seednap.steps.dada2.dada2_runner.Dada2Runner.check_dada2_packages")
    @patch("seednap.steps.dada2.dada2_runner.Dada2Runner.run_dada2_process")
    @patch("subprocess.run")
    def test_process(
        self,
        mock_subprocess: MagicMock,
        mock_run_dada2: MagicMock,
        mock_check_packages: MagicMock,
        temp_dir: Path,
    ) -> None:
        """Test complete DADA2 processing."""
        # Setup
        trimmed_dir = temp_dir / "trimmed"
        trimmed_dir.mkdir()

        # Mock Rscript check
        mock_subprocess.return_value = MagicMock(
            stdout="", stderr="R version 4.2.0", returncode=0
        )

        # Mock package check
        mock_check_packages.return_value = {
            "dada2": "1.26.0",
            "Biostrings": "2.66.0",
            "DECIPHER": "2.26.0",
            "dplyr": "1.1.0",
            "ggplot2": "3.4.0",
            "patchwork": "1.1.2",
        }

        # Create expected output files
        output_dir = temp_dir / "02_dada2" / "teleo"
        output_dir.mkdir(parents=True)

        seqtab_clean_t = output_dir / "seqtab_clean_t.csv"
        pd.DataFrame({
            "sample1": [100],
        }, index=["ATCG"]).to_csv(seqtab_clean_t)

        corresp_seq = output_dir / "corresp_seq.csv"
        pd.DataFrame({
            "ASV_n": ["ASV1"],
        }, index=["ATCG"]).to_csv(corresp_seq, index=False)

        # Mock DADA2 process return
        mock_run_dada2.return_value = {
            "seqtab": output_dir / "seqtab.rds",
            "seqtab_clean": output_dir / "seqtab_clean.csv",
            "seqtab_clean_rds": output_dir / "seqtab_clean.rds",
            "seqtab_clean_t": seqtab_clean_t,
            "query_fasta": output_dir / "query.fasta",
            "corresp_seq": corresp_seq,
            "metrics_dir": output_dir / "QC",
        }

        # Run processor
        processor = Dada2Processor(
            marker="teleo",
            trimmed_reads_dir=trimmed_dir,
            output_base_dir=temp_dir,
        )

        outputs = processor.process(collect_metrics=True)

        # Verify
        assert "seqtab_clean_t" in outputs
        mock_check_packages.assert_called_once()
        mock_run_dada2.assert_called_once()

    @patch("seednap.steps.taxonomic_assignment.dada2_taxonomy_runner.Dada2TaxonomyRunner.run_dada2_taxonomy")
    @patch("subprocess.run")
    def test_assign_taxonomy(
        self, mock_subprocess: MagicMock, mock_run_taxonomy: MagicMock, temp_dir: Path
    ) -> None:
        """Test taxonomic assignment."""
        # Setup
        trimmed_dir = temp_dir / "trimmed"
        trimmed_dir.mkdir()

        output_dir = temp_dir / "02_dada2" / "test"
        output_dir.mkdir(parents=True)

        # Create required seqtab file
        seqtab_rds = output_dir / "seqtab_clean.rds"
        seqtab_rds.write_text("dummy rds file")

        # Mock Rscript check
        mock_subprocess.return_value = MagicMock(
            stdout="", stderr="R version 4.2.0", returncode=0
        )

        # Mock taxonomy return
        mock_run_taxonomy.return_value = {
            "taxonomy": output_dir / "taxonomy_dada2RDP.csv",
            "complete": temp_dir / "test_dada2RDP.csv",
        }

        # Create DB files
        rdp_db = temp_dir / "rdp.fasta"
        rdp_db.write_text(">seq1\nATCG\n")
        species_db = temp_dir / "species.fasta"
        species_db.write_text(">seq1\nATCG\n")

        # Run
        processor = Dada2Processor(
            marker="test",
            trimmed_reads_dir=trimmed_dir,
            output_base_dir=temp_dir,
        )

        outputs = processor.assign_taxonomy(rdp_db_path=rdp_db, species_db_path=species_db)

        # Verify
        assert "taxonomy" in outputs
        assert "complete" in outputs
        mock_run_taxonomy.assert_called_once()

    def test_assign_taxonomy_missing_seqtab(self, temp_dir: Path) -> None:
        """Test error when sequence table is missing."""
        trimmed_dir = temp_dir / "trimmed"
        trimmed_dir.mkdir()

        rdp_db = temp_dir / "rdp.fasta"
        rdp_db.write_text(">seq1\nATCG\n")
        species_db = temp_dir / "species.fasta"
        species_db.write_text(">seq1\nATCG\n")

        processor = Dada2Processor(
            marker="test",
            trimmed_reads_dir=trimmed_dir,
            output_base_dir=temp_dir,
        )

        with pytest.raises(FileNotFoundError, match="Sequence table not found"):
            processor.assign_taxonomy(rdp_db_path=rdp_db, species_db_path=species_db)
