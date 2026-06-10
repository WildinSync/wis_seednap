"""DADA2 ASV/sequence metrics for the eDNA pipeline.

Computes ASV statistics from the sequence table and exports a summary (JSON/CSV). Per-step
read counts live in the run report (read_tracking.csv / step_summary.csv), not here.
"""

import json
import logging
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
from pandas.errors import EmptyDataError

from seednap.errors import SeednapError

logger = logging.getLogger(__name__)


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


class MetricsCollector:
    """
    Collect and report DADA2 ASV/sequence statistics.

    This class computes ASV statistics from the sequence table, generates a summary report,
    and exports metrics to JSON/CSV. Per-step read counts are tracked separately in the run
    report (read_tracking.csv / step_summary.csv).
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

        self.asv_metrics = ASVMetrics()

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

        dada2_dir = self.output_dir / "02_dada2" / self.marker
        processing_log = dada2_dir / "dada2_processing.log"
        track_reads = dada2_dir / "track_reads.csv"

        # Read sequence table. A header-only or empty file (0 ASVs surviving
        # filtering/merging/chimera removal) is the realistic failure here; turn the
        # raw pandas error into a DADA2-context message rather than a bare traceback.
        try:
            seqtab = pd.read_csv(seqtab_path, index_col=0)
        except EmptyDataError as exc:
            raise SeednapError(
                f"Could not compute ASV metrics for marker '{self.marker}': the DADA2 "
                f"sequence table {seqtab_path} is empty (0 bytes / no columns)",
                why=(
                    "The DADA2 R step completed and this metrics summary ran afterwards, but the "
                    "sequence table has no content. If the file is truncated rather than "
                    "header-only, the R run was interrupted mid-write; otherwise every sequence "
                    "was dropped during filtering, failed to merge, or was removed as a chimera, "
                    "leaving zero ASVs."
                ),
                fix=(
                    f"Inspect {processing_log} and the per-sample read-tracking table "
                    f"({track_reads}) to see where reads were lost, then either re-run the dada2 "
                    f"step (if the table was truncated) or loosen dada2.filter (max_ee, trunc_q, "
                    f"min_len/max_len), dada2.merge (min_overlap), or dada2.chimera in the marker "
                    f"YAML."
                ),
            ) from exc

        if len(seqtab) == 0:
            raise SeednapError(
                f"Could not compute ASV metrics for marker '{self.marker}': the DADA2 run "
                f"produced zero ASVs (the sequence table {seqtab_path} has no rows)",
                why=(
                    "The DADA2 R step itself completed successfully and wrote its tables; only "
                    "this metrics summary failed. Zero ASVs usually means every sequence was "
                    "dropped during filtering, failed to merge, or was removed as a chimera."
                ),
                fix=(
                    f"Inspect {processing_log} and the per-sample read-tracking table "
                    f"({track_reads}) to see where reads were lost, then consider loosening "
                    f"dada2.filter (max_ee, trunc_q, min_len/max_len), dada2.merge (min_overlap), "
                    f"or dada2.chimera in the marker YAML. If {seqtab_path} is truncated rather "
                    f"than header-only, the R run was interrupted mid-write; re-run the dada2 step."
                ),
            )

        # Calculate metrics
        self.asv_metrics.num_asvs = len(seqtab)
        self.asv_metrics.num_samples = len(seqtab.columns)
        self.asv_metrics.total_abundance = int(seqtab.values.sum())

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
        report = [
            f"\n{'='*80}",
            f"Pipeline Metrics Summary - {self.marker.upper()}",
            f"{'='*80}",
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

        report.append("")
        report.append("Per-step read counts: see read_tracking.csv / step_summary.csv (04_report/).")
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
        data: List[Dict[str, object]] = []
        data.append({"metric": "marker", "value": self.marker, "category": "general"})

        # ASV metrics
        for key, value in self.asv_metrics.to_dict().items():
            data.append({"metric": key, "value": value, "category": "asvs"})

        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)

        logger.info(f"Exported metrics to {output_path}")
        return output_path
