"""Cutadapt subprocess wrapper for primer trimming and demultiplexing.

This module provides a Python wrapper around the cutadapt command-line tool,
handling primer trimming, demultiplexing, and adapter removal for eDNA metabarcoding.
"""

import logging
import subprocess
from pathlib import Path
from typing import List, Optional, Union

logger = logging.getLogger(__name__)


class CutadaptError(Exception):
    """Exception raised for cutadapt command errors."""

    pass


class CutadaptRunner:
    """Run cutadapt commands for primer trimming and demultiplexing.

    This class wraps the cutadapt command-line tool and provides methods for:
    - Standard primer trimming (two-pass approach)
    - Tag-based demultiplexing
    - Primer detection with various filtering options
    """

    def __init__(
        self,
        cores: int = 1,
        error_rate: float = 0.1,
        min_length: int = 20,
        min_overlap: int = 3,
        no_indels: bool = False,
        timeout: int = 7200,
    ):
        """
        Initialize cutadapt runner.

        Args:
            cores: Number of CPU cores to use (default: 1)
            error_rate: Maximum allowed error rate (default: 0.1 = 10%)
            min_length: Minimum read length after trimming (default: 20)
            min_overlap: Minimum overlap between read and adapter (default: 3)
            no_indels: Forbid insertions/deletions in adapters (default: False)
            timeout: Command timeout in seconds (default: 7200 = 2 hours)
        """
        self.cores = cores
        self.error_rate = error_rate
        self.min_length = min_length
        self.min_overlap = min_overlap
        self.no_indels = no_indels
        self.timeout = timeout

    def _build_base_command(self) -> List[str]:
        """Build base cutadapt command with common parameters."""
        cmd = [
            "cutadapt",
            "-j",
            str(self.cores),
            "-e",
            str(self.error_rate),
            "-m",
            str(self.min_length),
            "-O",
            str(self.min_overlap),
        ]

        if self.no_indels:
            cmd.append("--no-indels")

        return cmd

    def _run_command(self, cmd: List[str], log_file: Optional[Union[str, Path]] = None) -> str:
        """
        Execute cutadapt command.

        Args:
            cmd: Command list to execute
            log_file: Optional path to log file for stdout/stderr

        Returns:
            stdout from command

        Raises:
            CutadaptError: If command fails
        """
        logger.info(f"Running cutadapt: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True, timeout=self.timeout
            )

            # Write to log file if specified
            if log_file:
                log_path = Path(log_file)
                log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(log_path, "a") as f:
                    f.write(result.stdout)
                    f.write(result.stderr)

            logger.debug(f"cutadapt completed successfully")
            return result.stdout

        except subprocess.CalledProcessError as e:
            error_msg = f"cutadapt failed: {e.stderr}"
            logger.error(error_msg)
            raise CutadaptError(error_msg) from e

        except subprocess.TimeoutExpired as e:
            error_msg = f"cutadapt timed out after {self.timeout} seconds"
            logger.error(error_msg)
            raise CutadaptError(error_msg) from e

        except FileNotFoundError as e:
            error_msg = "cutadapt command not found. Is cutadapt installed?"
            logger.error(error_msg)
            raise CutadaptError(error_msg) from e

    def trim_primers(
        self,
        r1_input: Union[str, Path],
        r1_output: Union[str, Path],
        r2_input: Optional[Union[str, Path]] = None,
        r2_output: Optional[Union[str, Path]] = None,
        forward_primer: Optional[str] = None,
        reverse_primer: Optional[str] = None,
        adapter_5p_r1: Optional[str] = None,
        adapter_3p_r1: Optional[str] = None,
        adapter_5p_r2: Optional[str] = None,
        adapter_3p_r2: Optional[str] = None,
        untrimmed_r1: Optional[Union[str, Path]] = None,
        untrimmed_r2: Optional[Union[str, Path]] = None,
        discard_untrimmed: bool = False,
        log_file: Optional[Union[str, Path]] = None,
    ) -> str:
        """
        Trim adapters/primers from reads.

        Supports both single-end and paired-end reads. Adapters can be specified
        either as simple primer sequences or with specific orientations (5' or 3').

        Args:
            r1_input: Input R1 FASTQ file (can be gzipped)
            r1_output: Output R1 FASTQ file
            r2_input: Input R2 FASTQ file (optional for paired-end)
            r2_output: Output R2 FASTQ file (optional for paired-end)
            forward_primer: Forward primer sequence (added as -g and -G)
            reverse_primer: Reverse primer sequence (added as -g and -G)
            adapter_5p_r1: 5' adapter for R1 (-g)
            adapter_3p_r1: 3' adapter for R1 (-a)
            adapter_5p_r2: 5' adapter for R2 (-G)
            adapter_3p_r2: 3' adapter for R2 (-A)
            untrimmed_r1: Save untrimmed R1 reads to file
            untrimmed_r2: Save untrimmed R2 reads to file
            discard_untrimmed: Discard reads without adapters
            log_file: Path to log file

        Returns:
            stdout from cutadapt

        Raises:
            ValueError: If inputs are invalid
            CutadaptError: If cutadapt command fails
        """
        # Validate inputs
        r1_input = Path(r1_input)
        if not r1_input.exists():
            raise FileNotFoundError(f"R1 input file not found: {r1_input}")

        paired_end = r2_input is not None
        if paired_end:
            r2_input = Path(r2_input)
            if not r2_input.exists():
                raise FileNotFoundError(f"R2 input file not found: {r2_input}")
            if not r2_output:
                raise ValueError("r2_output required for paired-end reads")

        # Build command
        cmd = self._build_base_command()

        # Add adapters
        if forward_primer:
            cmd.extend(["-g", forward_primer])
            if paired_end:
                cmd.extend(["-G", reverse_primer or forward_primer])

        if reverse_primer and not forward_primer:
            cmd.extend(["-g", reverse_primer])
            if paired_end:
                cmd.extend(["-G", reverse_primer])

        if adapter_5p_r1:
            cmd.extend(["-g", adapter_5p_r1])
        if adapter_3p_r1:
            cmd.extend(["-a", adapter_3p_r1])
        if paired_end:
            if adapter_5p_r2:
                cmd.extend(["-G", adapter_5p_r2])
            if adapter_3p_r2:
                cmd.extend(["-A", adapter_3p_r2])

        # Add untrimmed output
        if untrimmed_r1:
            Path(untrimmed_r1).parent.mkdir(parents=True, exist_ok=True)
            cmd.extend(["--untrimmed-output", str(untrimmed_r1)])
        if untrimmed_r2 and paired_end:
            Path(untrimmed_r2).parent.mkdir(parents=True, exist_ok=True)
            cmd.extend(["--untrimmed-paired-output", str(untrimmed_r2)])

        # Discard untrimmed
        if discard_untrimmed:
            cmd.append("--discard-untrimmed")

        # Add output files
        r1_output = Path(r1_output)
        r1_output.parent.mkdir(parents=True, exist_ok=True)
        cmd.extend(["-o", str(r1_output)])

        if paired_end:
            r2_output = Path(r2_output)
            r2_output.parent.mkdir(parents=True, exist_ok=True)
            cmd.extend(["-p", str(r2_output)])

        # Add input files
        cmd.append(str(r1_input))
        if paired_end:
            cmd.append(str(r2_input))

        return self._run_command(cmd, log_file)

    def demultiplex_by_tags(
        self,
        r1_input: Union[str, Path],
        r2_input: Union[str, Path],
        tag_file: Union[str, Path],
        output_dir: Union[str, Path],
        discard_untrimmed: bool = True,
        log_file: Optional[Union[str, Path]] = None,
    ) -> str:
        """
        Demultiplex reads by tags using file-based adapter specification.

        Tag file should be in FASTA format with names matching sample names.
        Output files will be named {sample_name}.R1.fastq.gz and {sample_name}.R2.fastq.gz.

        Args:
            r1_input: Input R1 FASTQ file
            r2_input: Input R2 FASTQ file
            tag_file: FASTA file with tag sequences (format: >sample_name\\nTAGSEQ)
            output_dir: Directory for demultiplexed output files
            discard_untrimmed: Discard reads without matching tags (default: True)
            log_file: Path to log file

        Returns:
            stdout from cutadapt

        Raises:
            FileNotFoundError: If input files don't exist
            CutadaptError: If cutadapt command fails
        """
        # Validate inputs
        r1_input = Path(r1_input)
        r2_input = Path(r2_input)
        tag_file = Path(tag_file)

        if not r1_input.exists():
            raise FileNotFoundError(f"R1 input not found: {r1_input}")
        if not r2_input.exists():
            raise FileNotFoundError(f"R2 input not found: {r2_input}")
        if not tag_file.exists():
            raise FileNotFoundError(f"Tag file not found: {tag_file}")

        # Create output directory
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Build command - use exact error rate and no-indels for demultiplexing
        cmd = [
            "cutadapt",
            "-j",
            str(self.cores),
            "-e",
            "0.0",  # Exact matching for tags
            "--no-indels",
            "-a",
            f"file:{tag_file}",
            "-A",
            f"file:{tag_file}",
            "-o",
            str(output_dir / "{name}.R1.fastq.gz"),
            "-p",
            str(output_dir / "{name}.R2.fastq.gz"),
            str(r1_input),
            str(r2_input),
        ]

        if discard_untrimmed:
            cmd.append("--discard-untrimmed")

        return self._run_command(cmd, log_file)

    def detect_primers_no_trim(
        self,
        r1_input: Union[str, Path],
        r1_output: Union[str, Path],
        r2_input: Union[str, Path],
        r2_output: Union[str, Path],
        adapter_5p_r1: str,
        adapter_5p_r2: str,
        discard_untrimmed: bool = True,
        log_file: Optional[Union[str, Path]] = None,
    ) -> str:
        """
        Detect primers without trimming them (action=none).

        This is useful for filtering reads that have primers in the expected positions
        while keeping the full-length sequences (e.g., for ligation-based libraries).

        Args:
            r1_input: Input R1 FASTQ file
            r1_output: Output R1 FASTQ file
            r2_input: Input R2 FASTQ file
            r2_output: Output R2 FASTQ file
            adapter_5p_r1: 5' adapter pattern for R1 (e.g., "^PRIMER...RC_PRIMER")
            adapter_5p_r2: 5' adapter pattern for R2 (e.g., "^PRIMER...RC_PRIMER")
            discard_untrimmed: Discard reads without primers (default: True)
            log_file: Path to log file

        Returns:
            stdout from cutadapt

        Raises:
            FileNotFoundError: If input files don't exist
            CutadaptError: If cutadapt command fails
        """
        # Validate inputs
        r1_input = Path(r1_input)
        r2_input = Path(r2_input)

        if not r1_input.exists():
            raise FileNotFoundError(f"R1 input not found: {r1_input}")
        if not r2_input.exists():
            raise FileNotFoundError(f"R2 input not found: {r2_input}")

        # Create output directories
        r1_output = Path(r1_output)
        r2_output = Path(r2_output)
        r1_output.parent.mkdir(parents=True, exist_ok=True)
        r2_output.parent.mkdir(parents=True, exist_ok=True)

        # Build command
        cmd = [
            "cutadapt",
            "-j",
            str(self.cores),
            "--action=none",  # Don't trim, just filter
            "-e",
            str(self.error_rate),
            "--no-indels",
            "-m",
            str(self.min_length),
            "-g",
            adapter_5p_r1,
            "-G" if adapter_5p_r2 else "-A",
            adapter_5p_r2,
            "-o",
            str(r1_output),
            "-p",
            str(r2_output),
            str(r1_input),
            str(r2_input),
        ]

        if discard_untrimmed:
            cmd.append("--discard-untrimmed")

        return self._run_command(cmd, log_file)
