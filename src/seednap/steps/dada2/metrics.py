"""Metrics collection and reporting for eDNA pipeline.

This module provides classes for tracking and reporting metrics throughout
the seednap pipeline, including read counts, quality statistics, and ASV metrics.
"""

import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ReadMetrics:
    """Metrics for read counts at different pipeline stages."""

    raw_reads: int = 0
    trimmed_reads: int = 0
    filtered_reads: int = 0
    denoised_reads: int = 0
    merged_reads: int = 0
    non_chimeric_reads: int = 0

    def to_dict(self) -> Dict[str, int]:
        """Convert to dictionary."""
        return asdict(self)

    def get_retention_rates(self) -> Dict[str, float]:
        """
        Calculate retention rates at each step (as percentage of raw reads).

        Returns:
            Dictionary mapping step names to retention percentages
        """
        if self.raw_reads == 0:
            return {}

        return {
            "trimming": (self.trimmed_reads / self.raw_reads) * 100,
            "filtering": (self.filtered_reads / self.raw_reads) * 100,
            "denoising": (self.denoised_reads / self.raw_reads) * 100,
            "merging": (self.merged_reads / self.raw_reads) * 100,
            "chimera_removal": (self.non_chimeric_reads / self.raw_reads) * 100,
        }


@dataclass
class ASVMetrics:
    """Metrics for Amplicon Sequence Variants (ASVs)."""

    num_asvs: int = 0
    num_samples: int = 0
    total_abundance: int = 0
    min_length: int = 0
    max_length: int = 0
    mean_length: float = 0.0
    median_abundance: float = 0.0

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)


@dataclass
class SampleMetrics:
    """Metrics for a single sample."""

    sample_name: str
    reads: ReadMetrics = field(default_factory=ReadMetrics)
    num_asvs: int = 0
    total_abundance: int = 0

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            "sample_name": self.sample_name,
            "reads": self.reads.to_dict(),
            "num_asvs": self.num_asvs,
            "total_abundance": self.total_abundance,
        }


class MetricsCollector:
    """
    Collect and track metrics throughout the seednap pipeline.

    This class provides methods to:
    - Track read counts at each pipeline stage
    - Calculate retention rates
    - Compute ASV statistics
    - Generate summary reports
    - Export metrics to JSON/CSV
    """

    def __init__(self, marker: str, output_dir: Union[str, Path]):
        """
        Initialize metrics collector.

        Args:
            marker: Marker name (e.g., 'teleo', 'amph')
            output_dir: Base output directory for metrics
        """
        self.marker = marker
        self.output_dir = Path(output_dir)
        self.metrics_dir = self.output_dir / "02_dada2" / marker / "metrics"
        self.metrics_dir.mkdir(parents=True, exist_ok=True)

        self.read_metrics = ReadMetrics()
        self.asv_metrics = ASVMetrics()
        self.sample_metrics: Dict[str, SampleMetrics] = {}

    def count_fastq_reads(self, fastq_path: Union[str, Path]) -> int:
        """
        Count reads in a FASTQ file.

        Args:
            fastq_path: Path to FASTQ file (can be gzipped)

        Returns:
            Number of reads
        """
        import gzip

        fastq_path = Path(fastq_path)
        if not fastq_path.exists():
            logger.warning(f"FASTQ file not found: {fastq_path}")
            return 0

        try:
            if fastq_path.suffix == ".gz":
                opener = gzip.open
                mode = "rt"
            else:
                opener = open
                mode = "r"

            count = 0
            with opener(fastq_path, mode) as f:
                for i, line in enumerate(f):
                    if i % 4 == 0:  # FASTQ format: every 4th line is a read header
                        count += 1
            return count
        except Exception as e:
            logger.error(f"Error counting reads in {fastq_path}: {e}")
            return 0

    def collect_trimming_metrics(self, raw_dir: Union[str, Path], trimmed_dir: Union[str, Path]) -> None:
        """
        Collect metrics from primer trimming step.

        Args:
            raw_dir: Directory with raw FASTQ files
            trimmed_dir: Directory with trimmed FASTQ files
        """
        logger.info("Collecting trimming metrics")

        raw_dir = Path(raw_dir)
        trimmed_dir = Path(trimmed_dir)

        # Count raw reads
        raw_r1_files = sorted(raw_dir.glob("*_R1.fastq*"))
        for r1_file in raw_r1_files:
            self.read_metrics.raw_reads += self.count_fastq_reads(r1_file)

        # Count trimmed reads
        trimmed_r1_files = sorted(trimmed_dir.glob("*_R1.fastq*"))
        for r1_file in trimmed_r1_files:
            self.read_metrics.trimmed_reads += self.count_fastq_reads(r1_file)

        logger.info(
            f"Trimming: {self.read_metrics.raw_reads} raw reads → "
            f"{self.read_metrics.trimmed_reads} trimmed reads "
            f"({self.read_metrics.trimmed_reads / max(self.read_metrics.raw_reads, 1) * 100:.1f}%)"
        )

    def collect_filtering_metrics(self, filtered_dir: Union[str, Path]) -> None:
        """
        Collect metrics from DADA2 filtering step.

        Args:
            filtered_dir: Directory with filtered FASTQ files
        """
        logger.info("Collecting filtering metrics")

        filtered_dir = Path(filtered_dir)
        filtered_r1_files = sorted(filtered_dir.glob("*_R1.fastq*"))

        for r1_file in filtered_r1_files:
            self.read_metrics.filtered_reads += self.count_fastq_reads(r1_file)

        logger.info(
            f"Filtering: {self.read_metrics.trimmed_reads} trimmed reads → "
            f"{self.read_metrics.filtered_reads} filtered reads "
            f"({self.read_metrics.filtered_reads / max(self.read_metrics.trimmed_reads, 1) * 100:.1f}%)"
        )

    def collect_asv_metrics(
        self, seqtab_path: Union[str, Path], corresp_seq_path: Optional[Union[str, Path]] = None
    ) -> None:
        """
        Collect metrics from ASV table.

        Args:
            seqtab_path: Path to sequence table CSV (transposed)
            corresp_seq_path: Optional path to ASV correspondence CSV
        """
        logger.info("Collecting ASV metrics")

        seqtab_path = Path(seqtab_path)
        if not seqtab_path.exists():
            logger.warning(f"Sequence table not found: {seqtab_path}")
            return

        # Read sequence table
        seqtab = pd.read_csv(seqtab_path, index_col=0)

        # Calculate metrics
        self.asv_metrics.num_asvs = len(seqtab)
        self.asv_metrics.num_samples = len(seqtab.columns)
        self.asv_metrics.total_abundance = seqtab.values.sum()

        # Sequence lengths (from index if available)
        if corresp_seq_path:
            corresp = pd.read_csv(corresp_seq_path)
            if "sequence" in corresp.columns:
                lengths = corresp["sequence"].str.len()
                self.asv_metrics.min_length = int(lengths.min())
                self.asv_metrics.max_length = int(lengths.max())
                self.asv_metrics.mean_length = float(lengths.mean())

        # Median abundance per ASV
        abundances = seqtab.sum(axis=1)
        self.asv_metrics.median_abundance = float(abundances.median())

        # Total reads after chimera removal
        self.read_metrics.non_chimeric_reads = int(seqtab.values.sum())

        logger.info(
            f"ASVs: {self.asv_metrics.num_asvs} ASVs across {self.asv_metrics.num_samples} samples, "
            f"total abundance: {self.asv_metrics.total_abundance:,}"
        )

    def generate_summary_report(self) -> str:
        """
        Generate a human-readable summary report.

        Returns:
            Formatted summary string
        """
        retention = self.read_metrics.get_retention_rates()

        report = [
            f"\n{'='*80}",
            f"Pipeline Metrics Summary - {self.marker.upper()}",
            f"{'='*80}",
            "",
            "Read Counts:",
            f"  Raw reads:             {self.read_metrics.raw_reads:>12,}",
            f"  After trimming:        {self.read_metrics.trimmed_reads:>12,} ({retention.get('trimming', 0):.1f}%)",
            f"  After filtering:       {self.read_metrics.filtered_reads:>12,} ({retention.get('filtering', 0):.1f}%)",
            f"  After merging:         {self.read_metrics.merged_reads:>12,} ({retention.get('merging', 0):.1f}%)",
            f"  After chimera removal: {self.read_metrics.non_chimeric_reads:>12,} ({retention.get('chimera_removal', 0):.1f}%)",
            "",
            "ASV Statistics:",
            f"  Number of ASVs:        {self.asv_metrics.num_asvs:>12,}",
            f"  Number of samples:     {self.asv_metrics.num_samples:>12,}",
            f"  Total abundance:       {self.asv_metrics.total_abundance:>12,}",
            f"  Median ASV abundance:  {self.asv_metrics.median_abundance:>12,.0f}",
        ]

        if self.asv_metrics.mean_length > 0:
            report.extend(
                [
                    f"  Min sequence length:   {self.asv_metrics.min_length:>12,} bp",
                    f"  Max sequence length:   {self.asv_metrics.max_length:>12,} bp",
                    f"  Mean sequence length:  {self.asv_metrics.mean_length:>12,.1f} bp",
                ]
            )

        report.append(f"{'='*80}\n")
        return "\n".join(report)

    def export_to_json(self, output_path: Optional[Union[str, Path]] = None) -> Path:
        """
        Export metrics to JSON file.

        Args:
            output_path: Optional output path (default: metrics_dir/metrics.json)

        Returns:
            Path to exported JSON file
        """
        if output_path is None:
            output_path = self.metrics_dir / "metrics.json"

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        metrics_dict = {
            "marker": self.marker,
            "read_metrics": self.read_metrics.to_dict(),
            "retention_rates": self.read_metrics.get_retention_rates(),
            "asv_metrics": self.asv_metrics.to_dict(),
        }

        with open(output_path, "w") as f:
            json.dump(metrics_dict, f, indent=2)

        logger.info(f"Exported metrics to {output_path}")
        return output_path

    def export_to_csv(self, output_path: Optional[Union[str, Path]] = None) -> Path:
        """
        Export metrics to CSV file.

        Args:
            output_path: Optional output path (default: metrics_dir/metrics.csv)

        Returns:
            Path to exported CSV file
        """
        if output_path is None:
            output_path = self.metrics_dir / "metrics.csv"

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Combine all metrics into DataFrame
        data = []
        data.append({"metric": "marker", "value": self.marker, "category": "general"})

        # Read metrics
        for key, value in self.read_metrics.to_dict().items():
            data.append({"metric": key, "value": value, "category": "reads"})

        # Retention rates
        for key, value in self.read_metrics.get_retention_rates().items():
            data.append({"metric": f"{key}_retention_pct", "value": f"{value:.2f}", "category": "retention"})

        # ASV metrics
        for key, value in self.asv_metrics.to_dict().items():
            data.append({"metric": key, "value": value, "category": "asvs"})

        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)

        logger.info(f"Exported metrics to {output_path}")
        return output_path
