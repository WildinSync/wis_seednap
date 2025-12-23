"""Tests for primer trimming and demultiplexing."""

import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from seednap.steps.trimming import (
    CutadaptError,
    CutadaptRunner,
    LigationTrimmer,
    StandardTrimmer,
    TagFileGenerator,
)


class TestCutadaptRunner:
    """Tests for CutadaptRunner class."""

    def test_initialization(self) -> None:
        """Test runner initialization with parameters."""
        runner = CutadaptRunner(cores=4, error_rate=0.15, min_length=30)
        assert runner.cores == 4
        assert runner.error_rate == 0.15
        assert runner.min_length == 30
        assert runner.no_indels is False

    def test_initialization_defaults(self) -> None:
        """Test runner initialization with defaults."""
        runner = CutadaptRunner()
        assert runner.cores == 1
        assert runner.error_rate == 0.1
        assert runner.min_length == 20
        assert runner.min_overlap == 3

    @patch("subprocess.run")
    def test_trim_primers_single_end(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test single-end primer trimming."""
        # Create input file
        r1_input = temp_dir / "test_R1.fastq"
        r1_input.write_text("@read1\nATCG\n+\nIIII\n")

        # Expected output
        r1_output = temp_dir / "trimmed_R1.fastq"

        # Mock successful cutadapt
        mock_run.return_value = MagicMock(stdout="Cutadapt output", returncode=0)

        runner = CutadaptRunner()
        result = runner.trim_primers(
            r1_input=r1_input,
            r1_output=r1_output,
            forward_primer="ATCG",
        )

        # Check command was called
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert "cutadapt" in args
        assert "-g" in args
        assert "ATCG" in args

    @patch("subprocess.run")
    def test_trim_primers_paired_end(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test paired-end primer trimming."""
        # Create input files
        r1_input = temp_dir / "test_R1.fastq.gz"
        r2_input = temp_dir / "test_R2.fastq.gz"

        with gzip.open(r1_input, "wt") as f:
            f.write("@read1\nATCG\n+\nIIII\n")
        with gzip.open(r2_input, "wt") as f:
            f.write("@read1\nGCTA\n+\nIIII\n")

        # Expected outputs
        r1_output = temp_dir / "trimmed_R1.fastq"
        r2_output = temp_dir / "trimmed_R2.fastq"

        # Mock successful cutadapt
        mock_run.return_value = MagicMock(stdout="Cutadapt output", returncode=0)

        runner = CutadaptRunner()
        result = runner.trim_primers(
            r1_input=r1_input,
            r1_output=r1_output,
            r2_input=r2_input,
            r2_output=r2_output,
            forward_primer="ATCG",
            reverse_primer="GCTA",
        )

        # Check command structure
        args = mock_run.call_args[0][0]
        assert "-p" in args  # Paired-end flag
        assert args.count("-g") == 1  # Forward primer for R1
        assert args.count("-G") == 1  # Reverse primer for R2

    @patch("subprocess.run")
    def test_trim_with_untrimmed_output(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test saving untrimmed reads."""
        r1_input = temp_dir / "test_R1.fastq"
        r1_input.write_text("@read1\nATCG\n+\nIIII\n")

        r1_output = temp_dir / "trimmed_R1.fastq"
        untrimmed_r1 = temp_dir / "untrimmed_R1.fastq"

        mock_run.return_value = MagicMock(stdout="Cutadapt output", returncode=0)

        runner = CutadaptRunner()
        runner.trim_primers(
            r1_input=r1_input,
            r1_output=r1_output,
            forward_primer="ATCG",
            untrimmed_r1=untrimmed_r1,
        )

        args = mock_run.call_args[0][0]
        assert "--untrimmed-output" in args
        assert str(untrimmed_r1) in args

    @patch("subprocess.run")
    def test_demultiplex_by_tags(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test tag-based demultiplexing."""
        # Create input files
        r1_input = temp_dir / "lib_R1.fastq.gz"
        r2_input = temp_dir / "lib_R2.fastq.gz"
        tag_file = temp_dir / "tags.fasta"

        with gzip.open(r1_input, "wt") as f:
            f.write("@read1\nATCG\n+\nIIII\n")
        with gzip.open(r2_input, "wt") as f:
            f.write("@read1\nGCTA\n+\nIIII\n")

        tag_file.write_text(">sample1\nACGT;min_overlap=8...tgca;min_overlap=8\n")

        output_dir = temp_dir / "demux"

        mock_run.return_value = MagicMock(stdout="Cutadapt output", returncode=0)

        runner = CutadaptRunner()
        runner.demultiplex_by_tags(
            r1_input=r1_input,
            r2_input=r2_input,
            tag_file=tag_file,
            output_dir=output_dir,
        )

        # Check command structure
        args = mock_run.call_args[0][0]
        assert f"file:{tag_file}" in args
        assert "-e" in args
        assert "0.0" in args  # Exact matching
        assert "--no-indels" in args
        assert "--discard-untrimmed" in args

    @patch("subprocess.run")
    def test_detect_primers_no_trim(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test primer detection without trimming."""
        r1_input = temp_dir / "test_R1.fastq.gz"
        r2_input = temp_dir / "test_R2.fastq.gz"

        with gzip.open(r1_input, "wt") as f:
            f.write("@read1\nATCG\n+\nIIII\n")
        with gzip.open(r2_input, "wt") as f:
            f.write("@read1\nGCTA\n+\nIIII\n")

        r1_output = temp_dir / "detected_R1.fastq.gz"
        r2_output = temp_dir / "detected_R2.fastq.gz"

        mock_run.return_value = MagicMock(stdout="Cutadapt output", returncode=0)

        runner = CutadaptRunner()
        runner.detect_primers_no_trim(
            r1_input=r1_input,
            r1_output=r1_output,
            r2_input=r2_input,
            r2_output=r2_output,
            adapter_5p_r1="^ATCG...TAGC",
            adapter_5p_r2="^GCTA...CGAT",
        )

        args = mock_run.call_args[0][0]
        assert "--action=none" in args
        assert "--no-indels" in args

    @patch("subprocess.run")
    def test_command_failure(self, mock_run: MagicMock, temp_dir: Path) -> None:
        """Test error handling when cutadapt fails."""
        r1_input = temp_dir / "test_R1.fastq"
        r1_input.write_text("@read1\nATCG\n+\nIIII\n")
        r1_output = temp_dir / "trimmed_R1.fastq"

        # Mock failed cutadapt
        import subprocess

        mock_run.side_effect = subprocess.CalledProcessError(
            1, "cutadapt", stderr="Error: invalid adapter"
        )

        runner = CutadaptRunner()
        with pytest.raises(CutadaptError, match="cutadapt failed"):
            runner.trim_primers(r1_input=r1_input, r1_output=r1_output, forward_primer="ATCG")

    def test_missing_input_file(self, temp_dir: Path) -> None:
        """Test error when input file doesn't exist."""
        runner = CutadaptRunner()

        with pytest.raises(FileNotFoundError):
            runner.trim_primers(
                r1_input=temp_dir / "nonexistent.fastq",
                r1_output=temp_dir / "output.fastq",
                forward_primer="ATCG",
            )


class TestTagFileGenerator:
    """Tests for TagFileGenerator class."""

    def test_initialization(self) -> None:
        """Test generator initialization."""
        generator = TagFileGenerator(min_overlap=10)
        assert generator.min_overlap == 10

    def test_format_tag_sequence(self) -> None:
        """Test tag sequence formatting."""
        generator = TagFileGenerator(min_overlap=8)
        result = generator._format_tag_sequence("ATCG")

        assert result.startswith("ATCG;min_overlap=8...")
        assert "cgat;min_overlap=8" in result  # Reverse complement lowercased

    def test_generate_standard_tag_files(self, temp_dir: Path) -> None:
        """Test standard tag file generation."""
        # Create metadata CSV
        metadata = temp_dir / "metadata.csv"
        df = pd.DataFrame({
            "sample_name": ["sample1", "sample2", "sample3"],
            "tag": ["ACGT", "TGCA", "GGCC"],
            "run": ["run1", "run1", "run2"],
        })
        df.to_csv(metadata, index=False)

        output_dir = temp_dir / "tags"

        generator = TagFileGenerator()
        result = generator.generate_standard_tag_files(
            metadata_csv=metadata, output_dir=output_dir
        )

        # Check outputs
        assert "run1" in result
        assert "run2" in result
        assert result["run1"].exists()
        assert result["run2"].exists()

        # Check run1 file content
        run1_content = result["run1"].read_text()
        assert ">sample1" in run1_content
        assert ">sample2" in run1_content
        assert "ACGT;min_overlap=8" in run1_content

        # Check run2 file content
        run2_content = result["run2"].read_text()
        assert ">sample3" in run2_content

    def test_generate_ligation_tag_files(self, temp_dir: Path) -> None:
        """Test ligation tag file generation."""
        # Create metadata CSV
        metadata = temp_dir / "ligation_metadata.csv"
        df = pd.DataFrame({
            "eventID": ["event1", "event2"],
            "tag_demultiplex": ["ACGT", "TGCA"],
            "library": ["lib1", "lib1"],
        })
        df.to_csv(metadata, index=False)

        output_dir = temp_dir / "ligation_tags"

        generator = TagFileGenerator()
        result = generator.generate_ligation_tag_files(
            metadata_csv=metadata, output_dir=output_dir
        )

        # Check outputs
        assert "lib1" in result
        assert result["lib1"].exists()

        # Check file content
        content = result["lib1"].read_text()
        assert ">event1" in content
        assert ">event2" in content

    def test_missing_columns_error(self, temp_dir: Path) -> None:
        """Test error when required columns are missing."""
        metadata = temp_dir / "bad_metadata.csv"
        df = pd.DataFrame({"wrong_column": ["value"]})
        df.to_csv(metadata, index=False)

        generator = TagFileGenerator()
        with pytest.raises(ValueError, match="Missing required columns"):
            generator.generate_standard_tag_files(metadata_csv=metadata, output_dir=temp_dir)


class TestStandardTrimmer:
    """Tests for StandardTrimmer class."""

    def test_initialization(self) -> None:
        """Test trimmer initialization."""
        trimmer = StandardTrimmer(cores=4, error_rate=0.15, min_length=30)
        assert trimmer.cutadapt.cores == 4
        assert trimmer.cutadapt.error_rate == 0.15
        assert trimmer.cutadapt.min_length == 30

    @patch("seednap.steps.trimming.cutadapt_runner.CutadaptRunner.trim_primers")
    def test_trim_sample(self, mock_trim: MagicMock, temp_dir: Path) -> None:
        """Test single sample trimming."""
        # Create input files
        r1_input = temp_dir / "sample_R1.fastq.gz"
        r2_input = temp_dir / "sample_R2.fastq.gz"

        with gzip.open(r1_input, "wt") as f:
            f.write("@read1\nATCG\n+\nIIII\n")
        with gzip.open(r2_input, "wt") as f:
            f.write("@read1\nGCTA\n+\nIIII\n")

        output_dir = temp_dir / "trimmed"

        # Mock cutadapt calls - create the output files it expects
        def mock_trim_side_effect(*args, **kwargs):
            # Create the output files that trim_primers would create
            r1_out = kwargs.get("r1_output")
            r2_out = kwargs.get("r2_output")
            if r1_out:
                r1_out.parent.mkdir(parents=True, exist_ok=True)
                r1_out.write_text("@read1\nATCG\n+\nIIII\n")
            if r2_out:
                r2_out.parent.mkdir(parents=True, exist_ok=True)
                r2_out.write_text("@read1\nGCTA\n+\nIIII\n")
            return "Cutadapt output"

        mock_trim.side_effect = mock_trim_side_effect

        trimmer = StandardTrimmer()
        r1_out, r2_out = trimmer.trim_sample(
            r1_input=r1_input,
            r2_input=r2_input,
            output_dir=output_dir,
            sample_name="sample",
            forward_primer="ATCG",
            reverse_primer="GCTA",
        )

        # Check outputs
        assert r1_out == output_dir / "sample.R1.fastq"
        assert r2_out == output_dir / "sample.R2.fastq"

        # Check cutadapt was called twice (two-pass)
        assert mock_trim.call_count == 2

    @patch("seednap.steps.trimming.cutadapt_runner.CutadaptRunner.trim_primers")
    def test_trim_directory(self, mock_trim: MagicMock, temp_dir: Path) -> None:
        """Test trimming all samples in directory."""
        # Create raw reads directory
        raw_dir = temp_dir / "raw"
        raw_dir.mkdir()

        # Create sample files
        for sample in ["sample1", "sample2"]:
            r1 = raw_dir / f"{sample}_R1.fastq.gz"
            r2 = raw_dir / f"{sample}_R2.fastq.gz"

            with gzip.open(r1, "wt") as f:
                f.write("@read1\nATCG\n+\nIIII\n")
            with gzip.open(r2, "wt") as f:
                f.write("@read1\nGCTA\n+\nIIII\n")

        output_dir = temp_dir / "trimmed"

        # Mock cutadapt - create output files
        def mock_trim_side_effect(*args, **kwargs):
            r1_out = kwargs.get("r1_output")
            r2_out = kwargs.get("r2_output")
            if r1_out:
                r1_out.parent.mkdir(parents=True, exist_ok=True)
                r1_out.write_text("@read1\nATCG\n+\nIIII\n")
            if r2_out:
                r2_out.parent.mkdir(parents=True, exist_ok=True)
                r2_out.write_text("@read1\nGCTA\n+\nIIII\n")
            return "Cutadapt output"

        mock_trim.side_effect = mock_trim_side_effect

        trimmer = StandardTrimmer()
        results = trimmer.trim_directory(
            raw_reads_dir=raw_dir,
            output_dir=output_dir,
            forward_primer="ATCG",
            reverse_primer="GCTA",
        )

        # Check results
        assert len(results) == 2
        # 2 samples × 2 passes each = 4 calls
        assert mock_trim.call_count == 4


class TestLigationTrimmer:
    """Tests for LigationTrimmer class."""

    def test_initialization(self) -> None:
        """Test trimmer initialization."""
        trimmer = LigationTrimmer(cores=4, error_rate=0.2, min_tag_overlap=10)
        assert trimmer.cutadapt.cores == 4
        assert trimmer.cutadapt.error_rate == 0.2
        assert trimmer.tag_generator.min_overlap == 10

    def test_merge_gzip_files(self, temp_dir: Path) -> None:
        """Test merging gzipped files."""
        # Create input files
        file1 = temp_dir / "file1.gz"
        file2 = temp_dir / "file2.gz"

        with gzip.open(file1, "wt") as f:
            f.write("Content from file 1\n")
        with gzip.open(file2, "wt") as f:
            f.write("Content from file 2\n")

        output = temp_dir / "merged.gz"

        # Merge
        LigationTrimmer._merge_gzip_files([file1, file2], output)

        # Check output
        assert output.exists()
        with gzip.open(output, "rt") as f:
            content = f.read()
            assert "Content from file 1" in content
            assert "Content from file 2" in content

    @patch("seednap.steps.trimming.cutadapt_runner.CutadaptRunner.demultiplex_by_tags")
    @patch("seednap.steps.trimming.cutadapt_runner.CutadaptRunner.detect_primers_no_trim")
    @patch("seednap.steps.trimming.tag_generator.TagFileGenerator.generate_ligation_tag_files")
    def test_process_library(
        self,
        mock_tag_gen: MagicMock,
        mock_detect: MagicMock,
        mock_demux: MagicMock,
        temp_dir: Path,
    ) -> None:
        """Test complete ligation library processing."""
        # Create raw library files
        raw_dir = temp_dir / "raw"
        raw_dir.mkdir()

        r1 = raw_dir / "lib1_R1.fastq.gz"
        r2 = raw_dir / "lib1_R2.fastq.gz"

        with gzip.open(r1, "wt") as f:
            f.write("@read1\nATCG\n+\nIIII\n")
        with gzip.open(r2, "wt") as f:
            f.write("@read1\nGCTA\n+\nIIII\n")

        # Create metadata
        metadata = temp_dir / "metadata.csv"
        pd.DataFrame({
            "eventID": ["sample1"],
            "tag_demultiplex": ["ACGT"],
            "library": ["lib1"],
        }).to_csv(metadata, index=False)

        # Mock tag generation
        tag_file = temp_dir / "tags" / "lib1.fasta"
        tag_file.parent.mkdir()
        tag_file.write_text(">sample1\nACGT\n")
        mock_tag_gen.return_value = {"lib1": tag_file}

        # Mock demultiplexing - create expected demux files
        demux_dir = temp_dir / "output" / "00_demultiplex_ligation" / "demultiplex"
        demux_dir.mkdir(parents=True)

        for sample in ["sample1"]:
            for read in ["R1", "R2"]:
                demux_file = demux_dir / f"{sample}.{read}.fastq.gz"
                with gzip.open(demux_file, "wt") as f:
                    f.write("@read1\nATCG\n+\nIIII\n")

        # Mock primer detection - create expected detected files
        def create_detected_files(*args, **kwargs):
            r1_out = kwargs.get("r1_output")
            r2_out = kwargs.get("r2_output")
            r1_out.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(r1_out, "wt") as f:
                f.write("@read1\nATCG\n+\nIIII\n")
            with gzip.open(r2_out, "wt") as f:
                f.write("@read1\nGCTA\n+\nIIII\n")

        mock_detect.side_effect = create_detected_files

        output_dir = temp_dir / "output"

        trimmer = LigationTrimmer()
        realigned_dir = trimmer.process_library(
            raw_reads_dir=raw_dir,
            library_name="lib1",
            metadata_csv=metadata,
            output_base_dir=output_dir,
            forward_primer="ATCG",
            reverse_primer="GCTA",
            gunzip_output=False,  # Skip gunzip for test
        )

        # Check results
        assert realigned_dir.exists()
        assert (realigned_dir / "sample1.R1.fastq.gz").exists()
        assert (realigned_dir / "sample1.R2.fastq.gz").exists()

        # Check method calls
        mock_tag_gen.assert_called_once()
        mock_demux.assert_called_once()
        # 2 rounds of detection for 1 sample = 2 calls
        assert mock_detect.call_count == 2
