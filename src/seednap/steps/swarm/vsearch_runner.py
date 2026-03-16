"""vsearch wrapper for read merging, dereplication, sorting, and chimera detection.

Provides Python methods around vsearch CLI commands used in the SWARM pipeline.
"""

import logging
import re
import subprocess
from pathlib import Path
from typing import Optional, Tuple, Union

from seednap.utils.subprocess import run_subprocess

logger = logging.getLogger(__name__)


class VsearchError(Exception):
    """Exception raised for vsearch errors."""

    pass


class VsearchRunner:
    """
    Wrapper around vsearch CLI for SWARM pipeline operations.

    Handles read merging, dereplication, abundance sorting,
    and de novo chimera detection.
    """

    def __init__(self, timeout: int = 3600):
        """
        Initialize vsearch runner.

        Args:
            timeout: Command timeout in seconds (default: 1 hour)

        Raises:
            VsearchError: If vsearch is not installed
        """
        self.timeout = timeout
        self.version = self._check_availability()

    @staticmethod
    def _parse_version(text: str) -> Tuple[int, ...]:
        """Extract (major, minor, patch) from vsearch --version output."""
        m = re.search(r"vsearch\s+v?(\d+)\.(\d+)\.(\d+)", text)
        if m:
            return tuple(int(x) for x in m.groups())
        return (0, 0, 0)

    def _check_availability(self) -> Tuple[int, ...]:
        """Check that vsearch is installed and return its version tuple."""
        try:
            result = subprocess.run(
                ["vsearch", "--version"],
                capture_output=True, text=True, timeout=10,
            )
            # vsearch prints version to stderr
            version_text = result.stderr or result.stdout or ""
            version = self._parse_version(version_text)
            logger.info(f"Detected vsearch version: {'.'.join(str(v) for v in version)}")
            return version
        except FileNotFoundError:
            raise VsearchError("vsearch is not installed or not on PATH")

    def merge_pairs(
        self,
        r1: Union[str, Path],
        r2: Union[str, Path],
        output: Union[str, Path],
        *,
        fastq_maxdiffs: int = 10,
        fastq_minovlen: int = 10,
        allow_stagger: bool = False,
        fastq_minmergelen: int = 0,
        fastq_maxns: Optional[int] = None,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Path:
        """
        Merge paired-end reads with vsearch.

        Args:
            r1: Path to R1 FASTQ file
            r2: Path to R2 FASTQ file
            output: Path to merged FASTQ output
            fastq_maxdiffs: Max differences in overlap region
            fastq_minovlen: Min overlap length for merging
            allow_stagger: Allow merging of staggered reads
            fastq_minmergelen: Min merged read length (0 = no filter)
            fastq_maxns: Max number of N bases allowed (None = no filter)
            log_file: Optional log file path

        Returns:
            Path to merged FASTQ file

        Raises:
            VsearchError: If merging fails
        """
        r1, r2, output = Path(r1), Path(r2), Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "vsearch",
            "--fastq_mergepairs", str(r1),
            "--reverse", str(r2),
            "--fastqout", str(output),
            "--fastq_maxdiffs", str(fastq_maxdiffs),
            "--fastq_minovlen", str(fastq_minovlen),
        ]

        if allow_stagger:
            cmd.append("--fastq_allowmergestagger")

        if fastq_minmergelen > 0:
            cmd.extend(["--fastq_minmergelen", str(fastq_minmergelen)])

        if fastq_maxns is not None:
            cmd.extend(["--fastq_maxns", str(fastq_maxns)])

        run_subprocess(cmd, timeout=self.timeout, log_file=log_file, error_class=VsearchError)

        logger.info(f"Merged pairs → {output}")
        return output

    def dereplicate(
        self,
        input_fasta: Union[str, Path],
        output_fasta: Union[str, Path],
        *,
        min_unique_size: int = 1,
        sizein: bool = False,
        relabel_sha1: bool = False,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Path:
        """
        Dereplicate sequences and annotate with abundance.

        Args:
            input_fasta: Input FASTA/FASTQ file
            output_fasta: Output dereplicated FASTA with ;size=N; annotations
            min_unique_size: Minimum abundance to keep a sequence
            sizein: Read ;size=N; annotations from input (for summing abundances)
            relabel_sha1: Relabel sequences with SHA1 hash of sequence content
            log_file: Optional log file path

        Returns:
            Path to dereplicated FASTA file

        Raises:
            VsearchError: If dereplication fails
        """
        input_fasta, output_fasta = Path(input_fasta), Path(output_fasta)
        output_fasta.parent.mkdir(parents=True, exist_ok=True)

        # vsearch >= 2.28 rejects FASTQ input with --derep_fulllength;
        # use --fastx_uniques instead (available since vsearch 2.17).
        # For older versions, --derep_fulllength accepts both FASTA and FASTQ.
        input_str = str(input_fasta)
        is_fastq = input_str.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))
        use_fastx_uniques = is_fastq and self.version >= (2, 28, 0)

        derep_cmd = "--fastx_uniques" if use_fastx_uniques else "--derep_fulllength"
        # --fastx_uniques uses --fastaout/--fastqout; --derep_fulllength uses --output
        output_flag = "--fastaout" if use_fastx_uniques else "--output"

        cmd = [
            "vsearch",
            derep_cmd, input_str,
            output_flag, str(output_fasta),
            "--sizeout",
            "--fasta_width", "0",
            "--minuniquesize", str(min_unique_size),
        ]

        if sizein:
            cmd.append("--sizein")

        if relabel_sha1:
            cmd.append("--relabel_sha1")

        run_subprocess(cmd, timeout=self.timeout, log_file=log_file, error_class=VsearchError)

        logger.info(f"Dereplicated → {output_fasta}")
        return output_fasta

    def sort_by_size(
        self,
        input_fasta: Union[str, Path],
        output_fasta: Union[str, Path],
        *,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Path:
        """
        Sort FASTA sequences by decreasing abundance.

        Args:
            input_fasta: Input FASTA with ;size=N; annotations
            output_fasta: Output sorted FASTA
            log_file: Optional log file path

        Returns:
            Path to sorted FASTA file

        Raises:
            VsearchError: If sorting fails
        """
        input_fasta, output_fasta = Path(input_fasta), Path(output_fasta)
        output_fasta.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "vsearch",
            "--sortbysize", str(input_fasta),
            "--output", str(output_fasta),
            "--fasta_width", "0",
        ]

        run_subprocess(cmd, timeout=self.timeout, log_file=log_file, error_class=VsearchError)

        logger.info(f"Sorted by size → {output_fasta}")
        return output_fasta

    def chimera_denovo(
        self,
        input_fasta: Union[str, Path],
        output_uchime: Union[str, Path],
        *,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Path:
        """
        Run de novo chimera detection with UCHIME.

        Args:
            input_fasta: Input FASTA (sorted by abundance)
            output_uchime: Output UCHIME results file
            log_file: Optional log file path

        Returns:
            Path to UCHIME output file

        Raises:
            VsearchError: If chimera detection fails
        """
        input_fasta, output_uchime = Path(input_fasta), Path(output_uchime)
        output_uchime.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "vsearch",
            "--uchime_denovo", str(input_fasta),
            "--uchimeout", str(output_uchime),
        ]

        run_subprocess(cmd, timeout=self.timeout, log_file=log_file, error_class=VsearchError)

        logger.info(f"Chimera detection → {output_uchime}")
        return output_uchime
