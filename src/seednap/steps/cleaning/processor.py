"""Presence-based control decontamination of an abundance table.

The lab's validated control-cleaning standard:

* an **extraction blank** (``Blank-ext`` / neg_cont_type "extraction negative") cleans the
  biological samples that share its ``extraction_ID``;
* a **PCR blank** (``Blank-PCR`` / "PCR negative", which carries no ``extraction_ID``) cleans
  the whole dataset.

Cleaning is **presence-based**: any OTU/ASV that has reads in an applicable control is removed
from (``mode="subtract"``) or flagged in (``mode="flag"``, the default) the associated samples.
Counts are never altered in flag mode; both modes emit a per-sample report
(``reads_before/after``, ``n_otus_removed``, ``n_reads_removed``, ``driving_controls``).

Control identity comes from the FAIRe manifest (``samp_category`` / ``neg_cont_type`` /
``extraction_ID``); a control column present in the abundance table but absent from the
manifest (e.g. an unlabelled ``Blank-PCR-3``) is classified by name via
:func:`~seednap.config.manifest.classify_control` and warned about. No silent fallbacks
(the no-silent-fallbacks policy): zero controls, an extraction blank matching no sample, and orphan
control columns all emit ``[WARN]``; nothing is removed without being counted and reported.

Statistical modes (decontam / microDecon prevalence) are a documented future option; this
module implements the presence-based standard only.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from seednap.config.manifest import SampleManifest, classify_control

logger = logging.getLogger(__name__)

_CLEANING_MODES = ("flag", "subtract")
# OTU-level flag column added to the output marking OTUs seen in any negative control.
CONTROL_FLAG_COL = "in_negative_control"


class CleaningResult(BaseModel):
    """Summary of a cleaning run (the per-sample table is returned alongside, as a DataFrame)."""

    model_config = ConfigDict(extra="forbid")

    mode: str
    n_controls: int = Field(description="Negative-control columns found in the abundance table")
    n_samples: int = Field(description="Biological sample columns cleaned/flagged")
    n_otus_flagged: int = Field(description="OTUs present in at least one negative control")
    total_reads_removed: int = Field(description="Reads removed across all samples (0 in flag mode)")


def _is_extraction_neg(neg_cont_type: Optional[str]) -> bool:
    """Return True if ``neg_cont_type`` names an extraction negative (extraction blank)."""
    return neg_cont_type is not None and "extraction" in neg_cont_type.lower()


class CleaningProcessor:
    """Decontaminate an abundance table against its negative controls."""

    def __init__(self, mode: str = "flag"):
        if mode not in _CLEANING_MODES:
            raise ValueError(f"mode must be one of {_CLEANING_MODES}; got {mode!r}")
        self.mode = mode

    def clean(
        self,
        abundance: pd.DataFrame,
        manifest: SampleManifest,
        *,
        id_col: str,
        sample_cols: Optional[List[str]] = None,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, CleaningResult]:
        """Clean ``abundance`` (one row per OTU, numeric sample columns) using ``manifest``.

        Args:
            abundance: OTU x sample table; ``id_col`` is the OTU identifier column.
            manifest: provides per-eventID samp_category / neg_cont_type / extraction_ID.
            id_col: name of the OTU identifier column (e.g. "sequence" or "ASV_ID").
            sample_cols: explicit sample columns. Required for a taxonomy table, which has
                numeric *non-sample* columns (e.g. ``pident``); when omitted, every numeric
                column other than ``id_col`` is treated as a sample (correct for a pure
                OTU/ASV count matrix).

        Returns:
            (cleaned_df, report_df, result). ``cleaned_df`` is a copy with a
            ``in_negative_control`` OTU flag column and, in subtract mode, control reads
            zeroed in the associated samples. ``report_df`` has one row per biological sample.
        """
        if id_col not in abundance.columns:
            raise ValueError(f"id_col {id_col!r} not in abundance columns")

        df = abundance.copy()
        by_event = {r.eventID: r for r in manifest.rows}

        if sample_cols is None:
            # No explicit list: numeric columns are counts; taxonomy/meta columns are strings.
            sample_cols = [
                c for c in df.columns
                if c != id_col and pd.api.types.is_numeric_dtype(df[c])
            ]
        else:
            sample_cols = [c for c in sample_cols if c in df.columns and c != id_col]

        control_cols: Dict[str, Dict[str, Optional[str]]] = {}
        bio_cols: List[str] = []
        for c in sample_cols:
            row = by_event.get(c)
            if row is not None and row.samp_category == "negative control":
                control_cols[c] = {"neg_cont_type": row.neg_cont_type, "extraction_id": row.extraction_ID}
            elif row is not None and row.is_control:
                # A positive control / PCR standard deliberately contains target species; using
                # it for decontamination would erase legitimate reads. Exclude it, but record it.
                bio_cols.append(c)
                logger.warning(
                    f"[WARN] cleaning: expected=negative control for decontamination, "
                    f"got={c!r} is samp_category={row.samp_category!r} (not a negative control), "
                    f"fallback=not used as a decontamination control"
                )
            elif row is not None:
                bio_cols.append(c)
            else:
                cls = classify_control(c)
                if cls.neg_cont_type is not None and cls.samp_category == "negative control":
                    control_cols[c] = {"neg_cont_type": cls.neg_cont_type, "extraction_id": None}
                    logger.warning(
                        f"[WARN] cleaning: expected control {c!r} in the manifest, got=absent, "
                        f"fallback=classified by name as {cls.neg_cont_type!r}"
                    )
                elif cls.is_control:
                    # Positive control / PCR standard classified by name: not a decontamination
                    # control. Keep its reads (treated as a biological column) but record it.
                    bio_cols.append(c)
                    logger.warning(
                        f"[WARN] cleaning: column {c!r} is absent from the manifest and "
                        f"classifies as samp_category={cls.samp_category!r} (not a negative "
                        f"control), fallback=not used as a decontamination control"
                    )
                elif cls.warn_reason:
                    # A control-looking name that matched no known control pattern; surface
                    # the specific reason (the no-silent-fallbacks policy) rather than the
                    # weaker generic message.
                    bio_cols.append(c)
                    logger.warning(
                        f"[WARN] cleaning: column {c!r} is absent from the manifest; "
                        f"treated as a biological sample (whole-dataset controls apply). "
                        f"{cls.warn_reason}"
                    )
                else:
                    bio_cols.append(c)
                    logger.warning(
                        f"[WARN] cleaning: sample column {c!r} is absent from the manifest; "
                        f"treated as a biological sample (whole-dataset controls apply)"
                    )

        if not control_cols:
            logger.warning(
                "[WARN] cleaning: expected=at least one negative control column, got=none, "
                "fallback=no cleaning performed (contamination cannot be assessed)"
            )

        # Partition controls: extraction-scoped (have an extraction_ID) vs whole-dataset.
        ext_controls: Dict[str, str] = {}
        whole_ds_controls: List[str] = []
        for c, info in control_cols.items():
            eid = info.get("extraction_id")
            if _is_extraction_neg(info.get("neg_cont_type")) and eid:
                ext_controls[c] = str(eid)
            else:
                whole_ds_controls.append(c)

        # Warn for extraction blanks whose batch matches no biological sample.
        bio_ext = {c: (by_event[c].extraction_ID if c in by_event else None) for c in bio_cols}
        sample_eids = {e for e in bio_ext.values() if e}
        for c, eid in ext_controls.items():
            if eid not in sample_eids:
                logger.warning(
                    f"[WARN] cleaning: extraction blank {c!r} extraction_ID={eid!r} matches no "
                    f"biological sample; it cleans nothing"
                )

        # OTU-level flag: present in ANY control (annotation; never alters counts).
        all_controls = list(control_cols)
        if all_controls:
            df[CONTROL_FLAG_COL] = (df[all_controls] > 0).any(axis=1)
        else:
            df[CONTROL_FLAG_COL] = False
        n_otus_flagged = int(df[CONTROL_FLAG_COL].sum())

        report_rows: List[Dict[str, object]] = []
        total_removed = 0
        for s in bio_cols:
            applicable = list(whole_ds_controls) + [
                c for c, eid in ext_controls.items() if bio_ext.get(s) == eid
            ]
            before = int(df[s].sum())
            if applicable:
                control_present = (df[applicable] > 0).any(axis=1)
                removed_mask = control_present & (df[s] > 0)
                n_otus = int(removed_mask.sum())
                n_reads = int(df.loc[removed_mask, s].sum())
                if self.mode == "subtract" and n_otus:
                    df.loc[removed_mask, s] = 0
                    total_removed += n_reads
            else:
                n_otus = n_reads = 0
            after = int(df[s].sum())
            report_rows.append({
                "eventID": s,
                "reads_before": before,
                "reads_after": after,
                "n_otus_removed": n_otus,
                "n_reads_removed": n_reads,
                "driving_controls": "|".join(applicable),
            })

        report = pd.DataFrame(
            report_rows,
            columns=["eventID", "reads_before", "reads_after", "n_otus_removed",
                     "n_reads_removed", "driving_controls"],
        )
        result = CleaningResult(
            mode=self.mode,
            n_controls=len(control_cols),
            n_samples=len(bio_cols),
            n_otus_flagged=n_otus_flagged,
            total_reads_removed=total_removed,
        )
        logger.info(
            f"cleaning ({self.mode}): {result.n_controls} control(s), {result.n_samples} "
            f"sample(s), {n_otus_flagged} OTU(s) flagged, {total_removed} reads removed"
        )
        return df, report, result
