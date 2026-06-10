"""SWARM OTU clustering workflow orchestration.

Coordinates the full SWARM pipeline: vsearch merge → dereplicate →
SWARM cluster → chimera detection → OTU table building.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple, Union

from seednap.steps.swarm.otu_table_builder import OtuTableBuilder
from seednap.steps.swarm.swarm_runner import SwarmClusterer
from seednap.steps.swarm.vsearch_runner import VsearchRunner

logger = logging.getLogger(__name__)


class SwarmProcessor:
    """
    Orchestrate complete SWARM OTU clustering workflow.

    Coordinates vsearch (merging, dereplication, sorting, chimera detection)
    and SWARM (clustering) to produce an OTU table from trimmed paired-end reads.
    """

    def __init__(
        self,
        marker: str,
        trimmed_reads_dir: Union[str, Path],
        output_base_dir: Union[str, Path],
        timeout: int = 3600,
    ):
        """
        Initialize SWARM processor.

        Args:
            marker: Marker name (e.g., 'teleo', 'amph')
            trimmed_reads_dir: Directory with primer-trimmed FASTQ files
            output_base_dir: Base output directory
            timeout: Timeout for external commands in seconds (default: 1 hour)
        """
        self.marker = marker.lower()
        self.trimmed_reads_dir = Path(trimmed_reads_dir)
        self.output_base_dir = Path(output_base_dir)

        if not self.trimmed_reads_dir.exists():
            raise FileNotFoundError(
                f"Trimmed reads directory not found: {self.trimmed_reads_dir}. "
                "SWARM clustering reads its input FASTQs from this directory. When the "
                "pipeline runs SWARM without a preceding trim step, this path comes from "
                "paths.raw_data in the marker YAML; otherwise it is the trim step's output. "
                "The likely cause is a wrong/typo'd paths.raw_data, or the directory was "
                "never created because the trim step did not run for this marker. Fix: "
                "confirm paths.raw_data points at an existing directory of R1/R2 FASTQ "
                "pairs, or add the trim step before swarm in pipeline.steps (it writes "
                "trimmed pairs to outputs/01_trim/<marker>/). If invoking the swarm CLI "
                "directly, pass an existing directory as the TRIMMED_READS_DIR positional "
                "argument."
            )

        self.vsearch = VsearchRunner(timeout=timeout)
        self.swarm = SwarmClusterer(timeout=timeout)
        self.table_builder = OtuTableBuilder()

        self.output_dir = self.output_base_dir / "02_swarm" / self.marker
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Initialized SWARM processor for marker: {self.marker}")

    def process(
        self,
        *,
        d: int = 1,
        fastidious: bool = True,
        boundary: int = 3,
        threads: int = 4,
        fastq_maxdiffs: int = 10,
        fastq_minovlen: int = 10,
        allow_stagger: bool = False,
        min_sequence_length: int = 20,
        chimera_detection: bool = True,
    ) -> Dict[str, Path]:
        """
        Run complete SWARM clustering workflow.

        Steps:
        1. Find trimmed R1/R2 pairs
        2. Merge pairs per sample (vsearch)
        3. Dereplicate per sample (vsearch)
        4. Combine and globally dereplicate
        5. SWARM clustering
        6. Sort representatives by abundance
        7. Chimera detection (optional)
        8. Build OTU contingency table
        9. Write normalized outputs for taxonomy

        Args:
            d: SWARM distance threshold (default: 1)
            fastidious: Enable fastidious mode (default: True)
            boundary: Fastidious boundary value (default: 3)
            threads: Number of threads (default: 4)
            fastq_maxdiffs: Max differences in overlap for merging (default: 10)
            fastq_minovlen: Min overlap length for merging (default: 10)
            allow_stagger: Allow staggered read merging (default: False)
            min_sequence_length: Min sequence length after merging (default: 20)
            chimera_detection: Run chimera detection (default: True)

        Returns:
            Dictionary with output paths (keys match DADA2 output convention):
            - query_fasta: Representative sequences FASTA (for taxonomy)
            - seqtab_clean_t: OTU abundance table CSV (for taxonomy)
            - otu_table_full: Full OTU table with metadata
            - merged_dir: Directory with merged reads
        """
        logger.info(f"Starting SWARM workflow for {self.marker}")
        logger.info(f"Parameters: d={d}, fastidious={fastidious}, threads={threads}")

        # Step 1: Find sample pairs
        sample_pairs = self._find_sample_pairs(self.trimmed_reads_dir)
        logger.info(f"Found {len(sample_pairs)} samples to process")

        if not sample_pairs:
            raise FileNotFoundError(
                f"No R1/R2 FASTQ pairs found in {self.trimmed_reads_dir}. SWARM needs "
                "paired-end FASTQ files named so each R1 has a matching R2. No filenames "
                "in this directory match a supported pattern, or R1 files are present "
                "without a matching R2. Supported patterns: {sample}.R1.fastq, "
                "{sample}_R1.fastq, {sample}_R1_001.fastq and their .gz variants, each "
                "with the same name using R2. Note: the .fq / .fq.gz extension and "
                "SRA-style _1/_2 naming are NOT matched here. Fix: point this at the "
                "directory that holds the paired reads (the trim step's 01_trim/<marker>/ "
                "output, or your pre-trimmed reads), and rename files to one of the "
                "patterns above if needed. If you expected pairs to be here, check the run "
                "log for per-sample 'No R2 found for ... skipping' warnings, which flag "
                "orphaned R1 files."
            )

        # Step 2 & 3: Merge and dereplicate per sample
        merged_dir = self.output_dir / "merged"
        derep_dir = self.output_dir / "dereplicated"
        log_dir = self.output_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)

        sample_fastas = []
        skipped_samples = []
        for sample_name, r1_path, r2_path in sample_pairs:
            logger.info(f"Processing sample: {sample_name}")

            # Merge (filter N bases with --fastq_maxns 0)
            merged_path = merged_dir / f"{sample_name}.merged.fastq"
            self.vsearch.merge_pairs(
                r1=r1_path,
                r2=r2_path,
                output=merged_path,
                fastq_maxdiffs=fastq_maxdiffs,
                fastq_minovlen=fastq_minovlen,
                allow_stagger=allow_stagger,
                fastq_minmergelen=min_sequence_length,
                fastq_maxns=0,
                log_file=log_dir / f"{sample_name}_merge.log",
            )

            # Skip empty merged files (blanks/negative controls)
            if not merged_path.exists() or merged_path.stat().st_size == 0:
                logger.warning(
                    f"Sample {sample_name}: merged file is empty, skipping"
                )
                skipped_samples.append(sample_name)
                continue

            # Dereplicate with SHA1 relabeling so same sequence
            # gets the same ID across all samples
            derep_path = derep_dir / f"{sample_name}.fasta"
            self.vsearch.dereplicate(
                input_fasta=merged_path,
                output_fasta=derep_path,
                min_unique_size=1,
                relabel_sha1=True,
                log_file=log_dir / f"{sample_name}_derep.log",
            )

            # Skip if dereplication produced empty output
            if not derep_path.exists() or derep_path.stat().st_size == 0:
                logger.warning(
                    f"Sample {sample_name}: no sequences after dereplication, skipping"
                )
                skipped_samples.append(sample_name)
                continue

            sample_fastas.append(derep_path)

        if skipped_samples:
            logger.info(
                f"Skipped {len(skipped_samples)} empty samples: "
                f"{', '.join(skipped_samples[:10])}"
                f"{'...' if len(skipped_samples) > 10 else ''}"
            )

        if not sample_fastas:
            raise ValueError(
                "All samples produced empty output after read merging, so the SWARM path has no "
                "sequences to cluster. Usually the paired reads do not overlap (the amplicon is "
                "longer than read1 + read2) or the primers were not trimmed. Check the trim-step "
                "output is non-empty, and the swarm.merge settings (fastq_minovlen, "
                "fastq_maxdiffs) for your amplicon length."
            )

        # Step 4: Combine and globally dereplicate (with --sizein to sum abundances)
        combined_path = self.output_dir / "combined.fasta"
        self._combine_fastas(sample_fastas, combined_path)

        all_uniq_path = self.output_dir / "all.uniq.fasta"
        self.vsearch.dereplicate(
            input_fasta=combined_path,
            output_fasta=all_uniq_path,
            min_unique_size=1,
            sizein=True,
            log_file=log_dir / "global_derep.log",
        )

        # Step 5: SWARM clustering
        swarm_outputs = self.swarm.cluster(
            input_fasta=all_uniq_path,
            output_dir=self.output_dir,
            d=d,
            fastidious=fastidious,
            boundary=boundary,
            threads=threads,
            log_file=log_dir / "swarm.log",
        )

        # Step 6: Sort representatives by abundance
        sorted_reps = self.output_dir / "cluster_representatives.sorted.fasta"
        self.vsearch.sort_by_size(
            input_fasta=swarm_outputs["representatives"],
            output_fasta=sorted_reps,
            log_file=log_dir / "sort.log",
        )

        # Step 7: Chimera detection
        uchime_path = None
        if chimera_detection:
            uchime_path = self.output_dir / "cluster_representatives.uchime"
            self.vsearch.chimera_denovo(
                input_fasta=sorted_reps,
                output_uchime=uchime_path,
                log_file=log_dir / "chimera.log",
            )

        # Step 8: Build OTU table
        logger.info("Building OTU contingency table...")
        otu_table = self.table_builder.build(
            representatives_fasta=sorted_reps,
            stats_file=swarm_outputs["stats_file"],
            swarm_file=swarm_outputs["swarm_file"],
            uchime_file=uchime_path,
            sample_fastas=sample_fastas,
        )

        # Step 9: Write outputs
        otu_table_full_path = self.output_dir / "otu_table_full.csv"
        otu_table.to_csv(otu_table_full_path, index=False)
        logger.info(f"Full OTU table → {otu_table_full_path}")

        query_fasta_path = self.output_dir / "query.fasta"
        abundance_csv_path = self.output_dir / "otu_table.csv"
        self.table_builder.to_taxonomy_input(
            otu_table=otu_table,
            query_fasta_path=query_fasta_path,
            abundance_csv_path=abundance_csv_path,
        )

        logger.info(f"SWARM workflow completed for {self.marker}")

        return {
            "query_fasta": query_fasta_path,
            "seqtab_clean_t": abundance_csv_path,
            "otu_table_full": otu_table_full_path,
            "merged_dir": merged_dir,
        }

    @staticmethod
    def _find_sample_pairs(
        trimmed_dir: Path,
    ) -> List[Tuple[str, Path, Path]]:
        """
        Find R1/R2 FASTQ pairs in trimmed reads directory.

        Supports patterns: {sample}.R1.fastq, {sample}_R1.fastq,
        {sample}_R1_001.fastq, and .gz variants.

        Args:
            trimmed_dir: Directory containing trimmed FASTQ files

        Returns:
            List of (sample_name, r1_path, r2_path) tuples
        """
        r1_patterns = ["*.R1.fastq", "*_R1.fastq", "*_R1_001.fastq",
                        "*.R1.fastq.gz", "*_R1.fastq.gz", "*_R1_001.fastq.gz"]

        r1_files: List[Path] = []
        for pattern in r1_patterns:
            r1_files.extend(trimmed_dir.glob(pattern))

        pairs = []
        seen = set()
        for r1 in sorted(r1_files):
            # Extract sample name
            match = re.match(r"(.+?)[._]R1", r1.name)
            if not match:
                continue
            sample_name = match.group(1)
            if sample_name in seen:
                continue
            seen.add(sample_name)

            # Find corresponding R2. Rewrite ONLY the matched R1 read-token, not
            # an arbitrary "R1" substring: the sample-name prefix may itself
            # contain "R1" (e.g. MR12_R1.fastq, R1B-site_R1.fastq), and a blind
            # .replace("R1", "R2", 1) would corrupt the prefix and silently drop
            # the sample (or pair it with the wrong R2). match.end() sits just
            # after the token's "R1", so match.end()-2 is the start of that "R1".
            token_start = match.end() - 2
            r2_name = r1.name[:token_start] + "R2" + r1.name[match.end():]
            r2 = r1.parent / r2_name
            if r2.exists():
                pairs.append((sample_name, r1, r2))
            else:
                logger.warning(f"No R2 found for {r1.name}, skipping sample {sample_name}")

        return pairs

    @staticmethod
    def _combine_fastas(fasta_list: List[Path], output_path: Path) -> Path:
        """
        Concatenate multiple FASTA files into one.

        Args:
            fasta_list: List of FASTA file paths
            output_path: Output combined FASTA path

        Returns:
            Path to combined FASTA
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as out:
            for fasta_path in fasta_list:
                with open(fasta_path) as f:
                    out.write(f.read())

        logger.debug(f"Combined {len(fasta_list)} FASTA files → {output_path}")
        return output_path
