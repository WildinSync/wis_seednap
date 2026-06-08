"""High-level primer trimming workflows for eDNA metabarcoding.

This module provides orchestration classes for complete trimming workflows:
- StandardTrimmer: Two-pass primer trimming for standard libraries
- LigationTrimmer: Tag demultiplexing + primer detection for ligation-based libraries
"""

import gzip
import logging
import shutil
from pathlib import Path
from typing import List, Union

from seednap.steps.trimming.cutadapt_runner import CutadaptRunner
from seednap.steps.trimming.tag_generator import TagFileGenerator
from seednap.utils.sequences import reverse_complement

logger = logging.getLogger(__name__)


class StandardTrimmer:
    """Two-pass primer trimming workflow for standard (non-ligation) libraries.

    The two-pass approach:
    1. Pass 1: Trim forward/reverse primers from 5' ends (-g/-G)
    2. Pass 2: Trim reverse complement primers from 3' ends (-a/-A)
    """

    def __init__(
        self,
        cores: int = 1,
        error_rate: float = 0.1,
        min_length: int = 20,
        overlap: int = 3,
    ):
        """
        Initialize standard trimmer.

        Args:
            cores: Number of CPU cores for cutadapt
            error_rate: Maximum allowed error rate (default: 0.1)
            min_length: Minimum read length after trimming (default: 20)
            overlap: Minimum overlap for primer detection (default: 3)
        """
        self.cutadapt = CutadaptRunner(
            cores=cores, error_rate=error_rate, min_length=min_length,
            min_overlap=overlap,
        )

    def trim_sample(
        self,
        r1_input: Union[str, Path],
        r2_input: Union[str, Path],
        output_dir: Union[str, Path],
        sample_name: str,
        forward_primer: str,
        reverse_primer: str,
        keep_untrimmed: bool = False,
        discard_untrimmed: bool = True,
    ) -> tuple:
        """
        Perform two-pass primer trimming on a single sample.

        Args:
            r1_input: Input R1 FASTQ file
            r2_input: Input R2 FASTQ file
            output_dir: Output directory
            sample_name: Sample name (for output file naming)
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence
            keep_untrimmed: Route reads lacking the 5' primer to a side file for
                inspection (default: False). When True they are written aside and
                removed from the main output (overrides ``discard_untrimmed``).
            discard_untrimmed: When True (default), pass-1 cutadapt drops reads in
                which the 5' primer was not found (``--discard-untrimmed``). When
                False, such reads are kept in the output. This is the config knob
                ``trimming.discard_untrimmed``; passing it here is what makes it
                take effect (it was previously inert -- untrimmed reads flowed
                through regardless).

        Returns:
            Tuple of (r1_output_path, r2_output_path)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Calculate reverse complements for 3' trimming
        fwd_rc = reverse_complement(forward_primer)
        rev_rc = reverse_complement(reverse_primer)

        logger.info(f"Starting two-pass trimming for sample: {sample_name}")

        # Temporary files for pass 1
        r1_temp = output_dir / f"{sample_name}.R1_TEMPORARY.fastq"
        r2_temp = output_dir / f"{sample_name}.R2_TEMPORARY.fastq"

        # Final output files
        r1_final = output_dir / f"{sample_name}.R1.fastq"
        r2_final = output_dir / f"{sample_name}.R2.fastq"

        # Untrimmed files (if keeping)
        untrimmed_r1 = output_dir / f"untrimmed_{sample_name}.R1.fastq" if keep_untrimmed else None
        untrimmed_r2 = output_dir / f"untrimmed_{sample_name}.R2.fastq" if keep_untrimmed else None

        # Pass 1: Trim 5' primers (-g/-G). Reads without the 5' primer are routed to a
        # side file when keep_untrimmed, else discarded when discard_untrimmed, else kept.
        do_discard = discard_untrimmed and not keep_untrimmed
        if not keep_untrimmed and not discard_untrimmed:
            logger.warning(
                f"[WARN] trimming {sample_name}: expected=primer on every read, "
                f"got=discard_untrimmed=False, fallback=reads lacking the 5' primer are "
                f"KEPT in the trimmed output"
            )
        logger.info(f"Pass 1: Trimming 5' primers for {sample_name}")
        self.cutadapt.trim_primers(
            r1_input=r1_input,
            r1_output=r1_temp,
            r2_input=r2_input,
            r2_output=r2_temp,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            untrimmed_r1=untrimmed_r1,
            untrimmed_r2=untrimmed_r2,
            discard_untrimmed=do_discard,
            log_file=output_dir.parent.parent / "logs" / f"{sample_name}_trim_pass1.txt",
        )

        # Pass 2: Trim 3' primers (-a/-A) on reverse complements
        logger.info(f"Pass 2: Trimming 3' primers for {sample_name}")
        self.cutadapt.trim_primers(
            r1_input=r1_temp,
            r1_output=r1_final,
            r2_input=r2_temp,
            r2_output=r2_final,
            adapter_3p_r1=rev_rc,
            adapter_3p_r2=fwd_rc,
            log_file=output_dir.parent.parent / "logs" / f"{sample_name}_trim_pass2.txt",
        )

        # Clean up temporary files
        r1_temp.unlink()
        r2_temp.unlink()

        if not keep_untrimmed and untrimmed_r1 and untrimmed_r2:
            if untrimmed_r1.exists():
                untrimmed_r1.unlink()
            if untrimmed_r2.exists():
                untrimmed_r2.unlink()

        logger.info(f"Completed two-pass trimming for {sample_name}")
        return (r1_final, r2_final)

    def trim_directory(
        self,
        raw_reads_dir: Union[str, Path],
        output_dir: Union[str, Path],
        forward_primer: str,
        reverse_primer: str,
        keep_untrimmed: bool = False,
        discard_untrimmed: bool = True,
    ) -> List[tuple]:
        """
        Trim all samples in a directory.

        Automatically detects paired R1/R2 files and processes them.

        Args:
            raw_reads_dir: Directory containing raw FASTQ files
            output_dir: Output directory for trimmed reads
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence
            keep_untrimmed: Save untrimmed reads (default: False)

        Returns:
            List of (r1_output, r2_output) tuples for each sample
        """
        raw_reads_dir = Path(raw_reads_dir)
        if not raw_reads_dir.exists():
            raise FileNotFoundError(f"Raw reads directory not found: {raw_reads_dir}")

        # Find all R1 files
        r1_files = sorted(raw_reads_dir.glob("*_R1.fastq.gz"))
        if not r1_files:
            r1_files = sorted(raw_reads_dir.glob("*_R1.fastq"))

        if not r1_files:
            logger.warning(f"No R1 FASTQ files found in {raw_reads_dir}")
            return []

        logger.info(f"Found {len(r1_files)} samples to trim")

        results = []
        for r1_file in r1_files:
            # Get sample name and corresponding R2 file
            sample_name = r1_file.name.replace("_R1.fastq.gz", "").replace("_R1.fastq", "")

            # Find R2 file
            r2_pattern = f"{sample_name}_R2.fastq*"
            r2_files = list(raw_reads_dir.glob(r2_pattern))

            if not r2_files:
                logger.warning(f"No R2 file found for {sample_name}, skipping")
                continue

            r2_file = r2_files[0]

            # Trim sample
            result = self.trim_sample(  # discard_untrimmed threaded below
                r1_input=r1_file,
                r2_input=r2_file,
                output_dir=output_dir,
                sample_name=sample_name,
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                keep_untrimmed=keep_untrimmed,
                discard_untrimmed=discard_untrimmed,
            )
            results.append(result)

        logger.info(f"Completed trimming {len(results)} samples")
        return results


class LigationTrimmer:
    """Complete workflow for ligation-based library demultiplexing and trimming.

    The ligation demultiplexing workflow:
    1. Generate tag files from metadata
    2. Demultiplex by tags
    3. Detect primers (expected orientation)
    4. Detect primers (reverse orientation)
    5. Merge and realign reads
    """

    def __init__(
        self,
        cores: int = 1,
        error_rate: float = 0.1,
        min_length: int = 20,
        min_tag_overlap: int = 8,
    ):
        """
        Initialize ligation trimmer.

        Args:
            cores: Number of CPU cores for cutadapt
            error_rate: Maximum allowed error rate for primer detection
            min_length: Minimum read length
            min_tag_overlap: Minimum overlap for tag matching (default: 8)
        """
        self.cutadapt = CutadaptRunner(
            cores=cores,
            error_rate=error_rate,
            min_length=min_length,
            no_indels=True,  # For tag matching
        )
        self.tag_generator = TagFileGenerator(min_overlap=min_tag_overlap)
        self.cores = cores

    def process_library(
        self,
        raw_reads_dir: Union[str, Path],
        library_name: str,
        metadata_csv: Union[str, Path],
        output_base_dir: Union[str, Path],
        forward_primer: str,
        reverse_primer: str,
        gunzip_output: bool = True,
        max_sample_failure_rate: float = 0.5,
    ) -> Path:
        """
        Process a single ligation-based library through complete workflow.

        Per-sample errors (cutadapt failure, missing tag match, etc.) are
        collected and logged; the library aborts only if more than
        `max_sample_failure_rate` of samples fail. This prevents one bad
        sample from killing an entire 200-sample library.

        Args:
            raw_reads_dir: Directory with raw library FASTQ files
            library_name: Library identifier (matches filename prefix)
            metadata_csv: Metadata CSV with sample/tag/library columns
            output_base_dir: Base output directory
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence
            gunzip_output: Gunzip final output files (default: True)
            max_sample_failure_rate: Abort if more than this fraction of samples
                fail. Default 0.5 (50%). Set to 1.0 to never abort.

        Returns:
            Path to realigned output directory

        Raises:
            FileNotFoundError: If input files not found
            ValueError: If metadata is invalid or too many samples fail
        """
        logger.info(f"Processing ligation library: {library_name}")

        raw_reads_dir = Path(raw_reads_dir)
        output_base_dir = Path(output_base_dir)

        # Find library FASTQ files
        r1_matches = list(raw_reads_dir.glob(f"{library_name}*_R1.fastq.gz"))
        r2_matches = list(raw_reads_dir.glob(f"{library_name}*_R2.fastq.gz"))

        if not r1_matches or not r2_matches:
            raise FileNotFoundError(
                f"Could not find R1/R2 files for library {library_name} in {raw_reads_dir}"
            )

        r1_file = r1_matches[0]
        r2_file = r2_matches[0]

        # Step 1: Generate tag files
        logger.info("Step 1: Generating tag files from metadata")
        tag_dir = output_base_dir / "00_demultiplex_ligation" / "cutadapt_tags"
        tag_files = self.tag_generator.generate_ligation_tag_files(
            metadata_csv=metadata_csv, output_dir=tag_dir
        )

        if library_name not in tag_files:
            raise ValueError(f"Library {library_name} not found in metadata")

        tag_file = tag_files[library_name]

        # Step 2: Demultiplex by tags
        logger.info("Step 2: Demultiplexing by tags")
        demux_dir = output_base_dir / "00_demultiplex_ligation" / "demultiplex"
        self.cutadapt.demultiplex_by_tags(
            r1_input=r1_file,
            r2_input=r2_file,
            tag_file=tag_file,
            output_dir=demux_dir,
            discard_untrimmed=True,
        )

        # Get list of demultiplexed samples
        samples = sorted(
            set(
                [
                    f.name.replace(".R1.fastq.gz", "").replace(".R2.fastq.gz", "")
                    for f in demux_dir.glob("*.R*.fastq.gz")
                    if "unknown" not in f.name.lower()
                ]
            )
        )

        logger.info(f"Demultiplexed {len(samples)} samples")

        # Calculate primer patterns for detection
        fwd_rc = reverse_complement(forward_primer)
        rev_rc = reverse_complement(reverse_primer)

        # Pattern for expected orientation
        pattern_r1_expected = f"^{forward_primer}...{rev_rc}"
        pattern_r2_expected = f"^{reverse_primer}...{fwd_rc}"

        # Pattern for reverse orientation
        pattern_r1_reverse = f"^{reverse_primer}...{fwd_rc}"
        pattern_r2_reverse = f"^{forward_primer}...{rev_rc}"

        # Step 3: Detect primers in expected orientation (round 1)
        logger.info("Step 3: Detecting primers (expected orientation)")
        primer_detect_dir = output_base_dir / "00_demultiplex_ligation" / "primer_detection"
        primer_detect_dir.mkdir(parents=True, exist_ok=True)

        # Track per-sample failures so one bad sample doesn't kill the library.
        failed_samples: List[str] = []
        for sample in samples:
            try:
                self.cutadapt.detect_primers_no_trim(
                    r1_input=demux_dir / f"{sample}.R1.fastq.gz",
                    r1_output=primer_detect_dir / f"trim_round1_{sample}.R1.fastq.gz",
                    r2_input=demux_dir / f"{sample}.R2.fastq.gz",
                    r2_output=primer_detect_dir / f"trim_round1_{sample}.R2.fastq.gz",
                    adapter_5p_r1=pattern_r1_expected,
                    adapter_5p_r2=pattern_r2_expected,
                    discard_untrimmed=True,
                )
            except Exception as e:
                logger.warning(
                    f"Step 3 (primer detect, expected orientation) failed for "
                    f"sample '{sample}': {e}"
                )
                failed_samples.append(sample)

        # Step 4: Detect primers in reverse orientation (round 2)
        logger.info("Step 4: Detecting primers (reverse orientation)")
        for sample in samples:
            if sample in failed_samples:
                continue  # already failed in step 3
            try:
                self.cutadapt.detect_primers_no_trim(
                    r1_input=demux_dir / f"{sample}.R1.fastq.gz",
                    r1_output=primer_detect_dir / f"trim_round2_{sample}.R1.fastq.gz",
                    r2_input=demux_dir / f"{sample}.R2.fastq.gz",
                    r2_output=primer_detect_dir / f"trim_round2_{sample}.R2.fastq.gz",
                    adapter_5p_r1=pattern_r1_reverse,
                    adapter_5p_r2=pattern_r2_reverse,
                    discard_untrimmed=True,
                )
            except Exception as e:
                logger.warning(
                    f"Step 4 (primer detect, reverse orientation) failed for "
                    f"sample '{sample}': {e}"
                )
                failed_samples.append(sample)

        # Bail out early if too many samples failed.
        if samples and (len(failed_samples) / len(samples)) > max_sample_failure_rate:
            raise ValueError(
                f"Demultiplexing failed for {len(failed_samples)} of "
                f"{len(samples)} samples ({100 * len(failed_samples) / len(samples):.0f}%) "
                f"in library '{library_name}', exceeding the "
                f"max_sample_failure_rate of {max_sample_failure_rate:.0%}. "
                f"Failed samples (first 10): {failed_samples[:10]}"
            )

        # Step 5: Merge and realign reads (skip failed samples)
        logger.info("Step 5: Merging and realigning reads")
        realigned_dir = output_base_dir / "00_demultiplex_ligation" / "realigned"
        realigned_dir.mkdir(parents=True, exist_ok=True)

        for sample in samples:
            if sample in failed_samples:
                continue
            try:
                # Merge R1: round1_R1 + round2_R2 (swapped!)
                self._merge_gzip_files(
                    [
                        primer_detect_dir / f"trim_round1_{sample}.R1.fastq.gz",
                        primer_detect_dir / f"trim_round2_{sample}.R2.fastq.gz",
                    ],
                    realigned_dir / f"{sample}.R1.fastq.gz",
                )

                # Merge R2: round1_R2 + round2_R1 (swapped!)
                self._merge_gzip_files(
                    [
                        primer_detect_dir / f"trim_round1_{sample}.R2.fastq.gz",
                        primer_detect_dir / f"trim_round2_{sample}.R1.fastq.gz",
                    ],
                    realigned_dir / f"{sample}.R2.fastq.gz",
                )
            except Exception as e:
                logger.warning(
                    f"Step 5 (merge realigned) failed for sample '{sample}': {e}"
                )
                failed_samples.append(sample)

        if failed_samples:
            logger.warning(
                f"Demultiplex completed with {len(failed_samples)} failed sample(s) "
                f"out of {len(samples)} for library '{library_name}': "
                f"{failed_samples[:10]}"
                f"{'...' if len(failed_samples) > 10 else ''}"
            )

        # Step 6: Gunzip if requested
        if gunzip_output:
            logger.info("Step 6: Gunzipping output files")
            for gz_file in realigned_dir.glob("*.fastq.gz"):
                output_file = gz_file.with_suffix("")
                with gzip.open(gz_file, "rb") as f_in:
                    with open(output_file, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                gz_file.unlink()

        logger.info(f"Completed ligation library processing: {library_name}")
        return realigned_dir

    @staticmethod
    def _merge_gzip_files(input_files: List[Path], output_file: Path) -> None:
        """
        Merge multiple gzipped files into one.

        Args:
            input_files: List of input gzipped files
            output_file: Output gzipped file
        """
        with gzip.open(output_file, "wb") as f_out:
            for input_file in input_files:
                if input_file.exists():
                    with gzip.open(input_file, "rb") as f_in:
                        shutil.copyfileobj(f_in, f_out)
