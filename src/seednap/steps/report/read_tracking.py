"""Per-sample read/sequence tracking across pipeline steps.

Assembles the classic eDNA "read tracking" table -- how many reads/sequences
survive each step -- from artifacts the pipeline already writes:

- **raw / trimmed** come from the per-sample Cutadapt logs
  (``<sample>_trim_pass1.txt`` / ``_trim_pass2.txt``);
- **DADA2 path** (``filtered -> denoised -> merged -> nonchim``) comes from the
  ``track_reads.csv`` emitted by ``seednap/scripts/dada2_process.R``;
- **SWARM path** (``clustered``) comes from per-sample column sums of
  ``otu_table.csv``.

The reported chain adapts to the method (``DADA2_STEPS`` vs ``SWARM_STEPS``).

Design rule (the no-silent-fallbacks policy): a count that cannot be measured is recorded
as *absent* (``pandas.NA``), never as a silent ``0``. "Absent" and "genuinely
zero" are distinguished, and an absent count raises a ``[WARN]`` so a broken
measurement is never mistaken for real data loss.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Union, cast

import pandas as pd

from seednap.utils.logging import get_logger

logger = get_logger(__name__)

DADA2_STEPS = ["raw", "trimmed", "filtered", "denoised", "merged", "nonchim"]
SWARM_STEPS = ["raw", "trimmed", "clustered"]

# Cutadapt summary lines (numbers carry thousands separators, e.g. 705,447).
_RE_PROCESSED = re.compile(r"Total read pairs processed:\s*([\d,]+)")
_RE_WRITTEN = re.compile(r"Pairs written \(passing filters\):\s*([\d,]+)")


def _parse_int(text: str) -> int:
    """Parse an integer that may carry thousands separators (e.g. ``705,447``)."""
    return int(text.replace(",", ""))


def _first_match(path: Path, pattern: re.Pattern) -> Optional[int]:
    """Return the first integer matched by ``pattern`` in ``path``, or None."""
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                m = pattern.search(line)
                if m:
                    return _parse_int(m.group(1))
    except OSError as exc:
        logger.warning(
            f"[WARN] read_tracking: expected=readable Cutadapt log, "
            f"got=unreadable ({path.name}: {exc}), fallback=absent",
        )
        return None
    return None


class ReadTrackingBuilder:
    """Build the per-sample read-tracking table from on-disk artifacts.

    Pass ``dada2_dir`` for the DADA2 chain or ``swarm_otu_table`` for the SWARM
    chain; with neither, only raw/trimmed are reported.
    """

    def __init__(
        self,
        marker: str,
        logs_dir: Union[str, Path],
        dada2_dir: Optional[Union[str, Path]] = None,
        swarm_otu_table: Optional[Union[str, Path]] = None,
        warn_below_retention_pct: float = 30.0,
        warn_step_loss_pct: float = 70.0,
    ) -> None:
        self.marker = marker
        self.logs_dir = Path(logs_dir)
        self.dada2_dir = Path(dada2_dir) if dada2_dir else None
        self.swarm_otu_table = Path(swarm_otu_table) if swarm_otu_table else None
        self.warn_below_retention_pct = warn_below_retention_pct
        self.warn_step_loss_pct = warn_step_loss_pct
        if self.dada2_dir is not None:
            self.steps = DADA2_STEPS
        elif self.swarm_otu_table is not None:
            self.steps = SWARM_STEPS
        else:
            self.steps = ["raw", "trimmed"]

    # ------------------------------------------------------------------
    # Count sources
    # ------------------------------------------------------------------

    def _trim_counts(self) -> Dict[str, Dict[str, Optional[int]]]:
        """Per-sample raw (pass1 processed) and trimmed (pass2 written) counts."""
        counts: Dict[str, Dict[str, Optional[int]]] = {}
        if not self.logs_dir.is_dir():
            logger.warning(
                f"[WARN] read_tracking: expected=trim log dir, "
                f"got=missing ({self.logs_dir}), fallback=no raw/trimmed counts",
            )
            return counts
        for pass1 in sorted(self.logs_dir.glob("*_trim_pass1.txt")):
            sample = pass1.name[: -len("_trim_pass1.txt")]
            pass2 = self.logs_dir / f"{sample}_trim_pass2.txt"
            raw = _first_match(pass1, _RE_PROCESSED)
            trimmed = _first_match(pass2, _RE_WRITTEN) if pass2.exists() else None
            if trimmed is None and pass2.exists():
                logger.warning(
                    f"[WARN] read_tracking: expected='Pairs written' in {pass2.name}, "
                    f"got=not found, fallback=absent",
                )
            counts[sample] = {"raw": raw, "trimmed": trimmed}
        return counts

    def _dada2_counts(self) -> pd.DataFrame:
        """Read the DADA2 ``track_reads.csv`` (filtered/denoised/merged/nonchim)."""
        # Only called when dada2_dir is set (see build()); narrow for the type checker.
        track = cast(Path, self.dada2_dir) / "track_reads.csv"
        if not track.exists():
            logger.warning(
                f"[WARN] read_tracking: expected=DADA2 track_reads.csv, "
                f"got=missing ({track}), fallback=no per-step DADA2 counts",
            )
            return pd.DataFrame()
        try:
            df = pd.read_csv(track)
        except pd.errors.EmptyDataError:
            logger.warning(
                f"[WARN] read_tracking: expected=DADA2 track_reads.csv with data, "
                f"got=empty ({track}), fallback=no per-step DADA2 counts",
            )
            return pd.DataFrame()
        if "sample" not in df.columns:
            logger.warning(
                f"[WARN] read_tracking: expected='sample' column in {track.name}, "
                f"got={list(df.columns)}, fallback=skip DADA2 counts",
            )
            return pd.DataFrame()
        return df.set_index("sample")

    def _swarm_counts(self) -> Dict[str, int]:
        """Per-sample 'clustered' reads = column sums of ``otu_table.csv``."""
        # Only called when swarm_otu_table is set (see build()); narrow for the checker.
        table = cast(Path, self.swarm_otu_table)
        if not table.exists():
            logger.warning(
                f"[WARN] read_tracking: expected=SWARM otu_table.csv, "
                f"got=missing ({table}), fallback=no clustered counts",
            )
            return {}
        try:
            df = pd.read_csv(table)
        except (pd.errors.EmptyDataError, OSError) as exc:
            logger.warning(
                f"[WARN] read_tracking: expected=readable otu_table.csv, "
                f"got=unreadable ({table}: {exc}), fallback=no clustered counts",
            )
            return {}
        # First column is the sequence/OTU id; the rest are per-sample counts.
        sample_cols = [c for c in df.columns[1:]]
        sums = df[sample_cols].apply(pd.to_numeric, errors="coerce").sum()
        return {str(k): int(v) for k, v in sums.items()}

    # ------------------------------------------------------------------
    # Assembly
    # ------------------------------------------------------------------

    def build(self) -> pd.DataFrame:
        """Assemble the per-sample tracking table; missing counts stay ``NA``."""
        trim = self._trim_counts()
        dada = self._dada2_counts() if self.dada2_dir is not None else pd.DataFrame()
        swarm = self._swarm_counts() if self.swarm_otu_table is not None else {}

        samples = sorted(set(trim) | set(dada.index.astype(str)) | set(swarm))
        if not samples:
            logger.warning(
                "[WARN] read_tracking: expected=samples from logs/track, "
                "got=none, fallback=empty table",
            )
            return pd.DataFrame(columns=["sample", *self.steps, "pct_retained"])

        rows: List[Dict[str, object]] = []
        for sample in samples:
            row: Dict[str, object] = {"sample": sample}
            t = trim.get(sample, {})
            row["raw"] = t.get("raw")
            row["trimmed"] = t.get("trimmed")
            if self.dada2_dir is not None:
                if sample in dada.index:
                    d = dada.loc[sample]
                    # DADA2 'input' == reads handed to filterAndTrim (the trimmed
                    # reads); use it only as a fallback if the trim log was absent.
                    if row["trimmed"] is None and "input" in d and pd.notna(d["input"]):
                        row["trimmed"] = int(d["input"])
                    for step in ("filtered", "denoised", "merged", "nonchim"):
                        row[step] = int(d[step]) if step in d and pd.notna(d[step]) else pd.NA
                else:
                    for step in ("filtered", "denoised", "merged", "nonchim"):
                        row[step] = pd.NA
            elif self.swarm_otu_table is not None:
                row["clustered"] = swarm.get(sample, pd.NA)
            rows.append(row)

        df = pd.DataFrame(rows, columns=["sample", *self.steps])
        # % of raw reads surviving to the final step of this method.
        raw = pd.to_numeric(df["raw"], errors="coerce")
        final = pd.to_numeric(df[self.steps[-1]], errors="coerce")
        df["pct_retained"] = (final / raw * 100).round(2)
        return df

    # ------------------------------------------------------------------
    # Warnings + output
    # ------------------------------------------------------------------

    def warnings(self, df: pd.DataFrame, log: bool = True) -> List[str]:
        """Data-loss + measurement warnings.

        When ``log`` is true (pipeline runs), each is emitted as a ``[WARN]`` to
        the configured logger so it lands in the run log (the no-silent-fallbacks policy).
        The standalone CLI passes ``log=False`` to avoid flooding the console --
        the same messages appear in the HTML report's "Notable events".
        """
        msgs: List[str] = []
        for _, r in df.iterrows():
            sample = r["sample"]
            absent = [s for s in self.steps if pd.isna(r[s])]
            if absent:
                msgs.append(
                    f"[WARN] read_tracking {sample}: expected=counts for "
                    f"{absent}, got=absent (not measured), fallback=NA"
                )
            pr = r["pct_retained"]
            if pd.notna(pr) and pr < self.warn_below_retention_pct:
                msgs.append(
                    f"[WARN] read_tracking {sample}: low overall retention "
                    f"{pr:.1f}% < {self.warn_below_retention_pct:.0f}% "
                    f"(raw={r['raw']} -> {self.steps[-1]}={r[self.steps[-1]]})"
                )
            for a, b in zip(self.steps, self.steps[1:]):
                va, vb = r[a], r[b]
                if pd.notna(va) and pd.notna(vb) and va > 0:
                    loss = (1 - vb / va) * 100
                    if loss > self.warn_step_loss_pct:
                        msgs.append(
                            f"[WARN] read_tracking {sample}: {a}->{b} dropped "
                            f"{loss:.1f}% ({va} -> {vb}), "
                            f"threshold {self.warn_step_loss_pct:.0f}%"
                        )
        if log:
            for m in msgs:
                logger.warning(m)
        return msgs

    def write(
        self, output_dir: Union[str, Path], df: Optional[pd.DataFrame] = None
    ) -> Dict[str, Path]:
        """Write ``read_tracking.csv`` and a human-readable ``.txt``.

        Pass a pre-built ``df`` to avoid recomputation. Warnings are NOT logged
        here -- call :meth:`warnings` once in the caller.
        """
        if df is None:
            df = self.build()
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        csv_path = out_dir / "read_tracking.csv"
        txt_path = out_dir / "read_tracking.txt"
        df.to_csv(csv_path, index=False)
        txt_path.write_text(self._render_text(df), encoding="utf-8")
        logger.info(f"Wrote read-tracking table: {csv_path}")
        return {"read_tracking_csv": csv_path, "read_tracking_txt": txt_path}

    def _render_text(self, df: pd.DataFrame) -> str:
        """Aligned plain-text table for the ``.txt`` artifact."""
        if df.empty:
            return f"Read tracking ({self.marker}): no samples found.\n"
        display = df.copy()
        for col in self.steps:
            display[col] = display[col].apply(
                lambda v: "NA" if pd.isna(v) else f"{int(v):,}"
            )
        display["pct_retained"] = display["pct_retained"].apply(
            lambda v: "NA" if pd.isna(v) else f"{v:.1f}%"
        )
        return f"Read tracking -- {self.marker}\n" + str(display.to_string(index=False)) + "\n"

    # ------------------------------------------------------------------
    # Step summary: run-level reads + feature (ASV/OTU) counts per step
    # ------------------------------------------------------------------
    def _feature_counts(self) -> Dict[str, Optional[int]]:
        """Per-step feature counts (number of ASVs/OTUs), only where a feature table exists.

        DADA2: read the run-level ``feature_counts.csv`` written by ``dada2_process.R``
        (the merged and non-chimeric ASV counts). SWARM: the number of OTUs (rows of the OTU
        table) at the ``clustered`` step. Read-level steps (raw/trimmed/filtered/denoised) have
        no feature table and are simply absent from the returned dict (reported as NA upstream).
        A missing or unreadable source raises a ``[WARN]`` rather than guessing (section 4).
        """
        counts: Dict[str, Optional[int]] = {}
        if self.dada2_dir is not None:
            fc = self.dada2_dir / "feature_counts.csv"
            if not fc.exists():
                logger.warning(
                    f"[WARN] step_summary: expected=DADA2 feature_counts.csv, "
                    f"got=missing ({fc}), fallback=no per-step ASV counts",
                )
                return counts
            try:
                df = pd.read_csv(fc)
                for _, r in df.iterrows():
                    counts[str(r["step"])] = int(r["n_features"])
            except (pd.errors.EmptyDataError, OSError, KeyError, ValueError) as exc:
                logger.warning(
                    f"[WARN] step_summary: expected=readable feature_counts.csv, "
                    f"got=unreadable ({fc}: {exc}), fallback=no per-step ASV counts",
                )
        elif self.swarm_otu_table is not None:
            if not self.swarm_otu_table.exists():
                logger.warning(
                    f"[WARN] step_summary: expected=OTU table for the OTU count, "
                    f"got=missing ({self.swarm_otu_table}), fallback=no OTU count",
                )
            else:
                try:
                    counts["clustered"] = int(len(pd.read_csv(self.swarm_otu_table)))
                except (pd.errors.EmptyDataError, OSError) as exc:
                    logger.warning(
                        f"[WARN] step_summary: expected=readable OTU table, "
                        f"got=unreadable ({self.swarm_otu_table}: {exc}), fallback=no OTU count",
                    )
        return counts

    def step_summary(self, tracking_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """Run-level summary: total reads and feature count after each pipeline step.

        Columns: ``step``; ``total_reads`` (sum of the per-sample read-tracking counts, NA if
        no sample had a measurable count at that step); ``n_features`` (number of ASVs for the
        DADA2 path / OTUs for SWARM, populated only at the stages where a feature table exists
        -- merged and nonchim for DADA2, clustered for SWARM -- and NA at the read-level steps).

        If a step is measured for some samples but NA (unmeasured) for others, the total is the
        sum over the measured samples and a ``[WARN]`` names the step and the unmeasured samples,
        so an incomplete run total is never reported silently (the no-silent-fallbacks policy).
        """
        if tracking_df is None:
            tracking_df = self.build()
        features = self._feature_counts()
        rows: List[Dict[str, object]] = []
        for step in self.steps:
            if step in tracking_df.columns:
                vals = pd.to_numeric(tracking_df[step], errors="coerce")
                if not vals.notna().any():
                    total: object = pd.NA
                else:
                    total = int(vals.sum())  # skipna: sum over the measured samples
                    if vals.isna().any():
                        missing = tracking_df.loc[vals.isna(), "sample"].astype(str).tolist()
                        logger.warning(
                            f"[WARN] step_summary: step '{step}' total_reads={total:,} is summed "
                            f"over {int(vals.notna().sum())}/{len(vals)} samples; "
                            f"{len(missing)} unmeasured (NA): {missing[:10]}"
                            f"{' ...' if len(missing) > 10 else ''} -- run total may be incomplete",
                        )
            else:
                total = pd.NA
            rows.append(
                {"step": step, "total_reads": total, "n_features": features.get(step, pd.NA)}
            )
        return pd.DataFrame(rows, columns=["step", "total_reads", "n_features"])

    def write_step_summary(
        self, output_dir: Union[str, Path], summary_df: Optional[pd.DataFrame] = None
    ) -> Path:
        """Write ``step_summary.csv`` (run-level reads + feature counts after each step)."""
        if summary_df is None:
            summary_df = self.step_summary()
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        csv_path = out_dir / "step_summary.csv"
        summary_df.to_csv(csv_path, index=False)
        logger.info(f"Wrote step summary: {csv_path}")
        return csv_path
