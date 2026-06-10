"""FAIRe-anchored sample manifest: one canonical, validated record per sample-library.

SeeDNAP today carries sample information across three loosely-specified CSVs (a demux
``metadata_lab_*.csv``, a per-sample ``metadata_field_*.csv``, and a one-row
``metadata_proj_*.csv``) with inconsistent dialects, no cross-validation, and a
library/run grouping column that is dropped after trimming. This module defines the
single canonical manifest those CSVs migrate into: one row per sample-library, anchored
on the **FAIRe (FAIR eDNA) v1.0.2** ``sampleMetadata`` + ``experimentRunMetadata`` term
names (which are themselves assembled from MIxS and Darwin Core), with DwC ``eventID``
kept as the canonical sample key (FAIRe ``samp_name`` is an alias equal to it).

What lives here:

* :class:`SampleManifestRow` -- the strict (``extra="forbid"``) canonical row model.
* :class:`SampleManifest` -- a validated collection of rows plus convenience accessors.
* :func:`load_manifest` -- read + validate a canonical manifest CSV.
* :func:`validate_against_abundance` -- the up-front cross-CSV ``eventID`` check (the
  silent-ID-mismatch guard, the no-silent-fallbacks policy).
* :func:`classify_control` -- the single source of truth for control classification by
  name pattern, a strict superset of the legacy ``blank|CNEG|CMET|CEXT`` regex.

Deriving a manifest from today's messy CSVs (BOM, dotted dates, header-casing, control
prefixes) lives in :mod:`seednap.config.manifest_migrate`.

Design notes (corrections from the manifest cross-check, ``03-manifest-verification``):

* ``samp_collect_method`` uses the MIxS canonical spelling; the GBIF DNA-derived-data
  extension spells it ``samp_collec_method`` (no "t"). Any GBIF export must translate.
* ``assay_name`` is FAIRe-Mandatory for metabarcoding; we model it Optional and let the
  migrator populate it from the marker, a deliberate, documented downgrade.
* Control identity is the two-field FAIRe model ``samp_category`` +
  ``neg_cont_type``/``pos_cont_type`` -- there is no single ``control_type`` field.
* Missing values (empty, ``NA``, and the INSDC tokens FAIRe mandates on control rows,
  e.g. ``not applicable: control sample``) all normalize to ``None`` here; a typed,
  ISO-8601 canonical form is enforced for ``eventDate``.

No silent fallbacks (the no-silent-fallbacks policy): a missing required field raises with the offending
``eventID``/file named; ambiguity is surfaced by the migrator as ``[WARN]``, never guessed.
"""

from __future__ import annotations

import logging
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

logger = logging.getLogger(__name__)

# Tokens that mean "no value" in the lab CSVs, including the INSDC missing-value
# vocabulary FAIRe mandates for Mandatory fields on control rows. Compared
# case-insensitively against stripped cell values.
MISSING_VALUE_TOKENS = frozenset(
    {
        "",
        "na",
        "n/a",
        "nan",
        "none",
        "null",
        "#n/a",
        "not applicable",
        "not applicable: control sample",
        "not applicable: sample group",
        "not collected",
        "not collected: not available",
        "not provided",
        "missing",
    }
)

# FAIRe samp_category controlled vocabulary (slots/samp_category.yaml). "other:" is a
# FAIRe convention: a literal "other:" optionally followed by free text.
SAMP_CATEGORIES = ("sample", "negative control", "positive control", "PCR standard")

# MIxS neg_cont_type controlled vocabulary (mixs/0001321).
NEG_CONT_TYPES = (
    "site negative",
    "field negative",
    "process negative",
    "extraction negative",
    "PCR negative",
)


def _is_other(value: str) -> bool:
    """Test whether a value uses the FAIRe ``other:`` free-text escape hatch.

    FAIRe controlled-vocabulary slots permit a value outside the fixed list by
    prefixing it with the literal ``other:`` followed by free text; this is how a
    real-world control type that is not in the standard vocabulary is recorded.

    Args:
        value: A candidate vocabulary value (leading/trailing whitespace tolerated).

    Returns:
        True if the value begins (case-insensitively) with ``other:``.
    """
    return value.strip().lower().startswith("other:")


# --------------------------------------------------------------------------- #
# Control classification (single source of truth)
# --------------------------------------------------------------------------- #
class ControlClass(BaseModel):
    """Result of classifying a sample name into the FAIRe two-field control model."""

    model_config = ConfigDict(extra="forbid")

    samp_category: str
    neg_cont_type: Optional[str] = None
    pos_cont_type: Optional[str] = None
    rule: str = Field(description="Which classification rule matched (for logging/provenance)")
    warn_reason: Optional[str] = Field(
        default=None,
        description="Set when the classification is an inference or an ambiguous "
        "control-looking name; the caller must emit a [WARN] (the no-silent-fallbacks policy)",
    )

    @property
    def is_control(self) -> bool:
        """Report whether this classification is a control rather than a field sample.

        Any non-``sample`` FAIRe category (negative control, positive control, or PCR
        standard) is a control: a deliberately non-biological well used to detect
        contamination or to calibrate, not a collected environmental specimen.

        Returns:
            True for any category other than ``sample``.
        """
        return self.samp_category != "sample"

    @property
    def is_pcr_blank(self) -> bool:
        """Report whether this control is a PCR-stage (rather than extraction-stage) blank.

        A PCR blank (no-template control added at amplification) legitimately has no
        extraction batch, so a null ``extraction_ID`` is expected for it and it is scoped
        to the whole dataset rather than to one extraction batch.

        Returns:
            True if the matched classification rule was a PCR-blank rule
            (``blank-pcr``, ``pcr-nc``, or ``cpcr``).
        """
        return self.rule in ("blank-pcr", "pcr-nc", "cpcr")


# Ordered (regex, ControlClass-kwargs) rules. First match wins. Patterns are matched
# case-insensitively against the *stripped* sample name. They tolerate the real-world
# separators and replicate/run suffixes seen across eras (Blank-ext-2, Blank_PCR-1,
# Blank-ext-2run2, CNEG01_03-MB1123A4, EXT_NC, PCR_NC, water).
_CONTROL_RULES: Tuple[Tuple[re.Pattern, Dict[str, Any]], ...] = (
    (re.compile(r"^blank[\s_-]*ext", re.I),
     dict(samp_category="negative control", neg_cont_type="extraction negative", rule="blank-ext")),
    (re.compile(r"^blank[\s_-]*pcr", re.I),
     dict(samp_category="negative control", neg_cont_type="PCR negative", rule="blank-pcr")),
    (re.compile(r"^pcr[\s_-]*nc\b", re.I),
     dict(samp_category="negative control", neg_cont_type="PCR negative", rule="pcr-nc")),
    (re.compile(r"^ext[\s_-]*nc\b", re.I),
     dict(samp_category="negative control", neg_cont_type="extraction negative", rule="ext-nc")),
    (re.compile(r"^cpcr", re.I),
     dict(samp_category="negative control", neg_cont_type="PCR negative", rule="cpcr")),
    (re.compile(r"^cneg", re.I),
     dict(samp_category="negative control", neg_cont_type="PCR negative", rule="cneg")),
    (re.compile(r"^cext", re.I),
     dict(samp_category="negative control", neg_cont_type="extraction negative", rule="cext")),
    # CMET: lab-internal, most likely a process/filtration ("Controle Methode") blank, but
    # this is an inference, not a documented standard -> classify but WARN.
    (re.compile(r"^cmet", re.I),
     dict(samp_category="negative control", neg_cont_type="process negative", rule="cmet",
          warn_reason="CMET mapped to 'process negative' by inference (unconfirmed lab convention)")),
    (re.compile(r"^(cpos|pos[\s_-]*c|mock)", re.I),
     dict(samp_category="positive control", pos_cont_type="other: positive/mock control",
          rule="positive",
          warn_reason="positive/mock control mapped to a generic pos_cont_type; confirm the exact type")),
    (re.compile(r"^water\b", re.I),
     dict(samp_category="negative control", neg_cont_type="other: water control", rule="water",
          warn_reason="'water' control mapped to neg_cont_type 'other:'; confirm the exact control type")),
    # Bare "blank" with no ext/pcr qualifier is genuinely ambiguous.
    (re.compile(r"^blank", re.I),
     dict(samp_category="negative control", neg_cont_type="other: unspecified blank", rule="blank-bare",
          warn_reason="bare 'blank' prefix is ambiguous (not ext/pcr); classified as 'other:'")),
)

# A name that *looks* like a control but matches no rule above must not be silently
# treated as a biological sample (it would inject contamination reads into the dataset).
_CONTROL_LIKE = re.compile(r"(blank|neg|ctrl|control|\bnc\b|_nc|nc_|water)", re.I)


def classify_control(sample_name: str) -> ControlClass:
    """Classify a sample name into the FAIRe two-field control model.

    The single source of truth for control identity. A strict superset of the legacy
    ``blank|CNEG|CMET|CEXT`` regex: it additionally catches ``CPCR``, ``EXT_NC``,
    ``PCR_NC``, ``water`` and the underscore/space-separated and run-suffixed forms.
    Controls (blanks, negative/positive controls) must be distinguished from biological
    samples so their reads are treated as contamination signal, not as detected taxa.

    Args:
        sample_name: The raw sample/library name as it appears in the lab CSVs (matched
            case-insensitively against its stripped form; None is tolerated).

    Returns:
        A :class:`ControlClass` carrying ``samp_category``, the matching
        ``neg_cont_type``/``pos_cont_type`` (if any), the ``rule`` that matched, and an
        optional ``warn_reason``. When ``warn_reason`` is set, or when a control-looking
        name fails to classify (``rule == 'unclassified-control-like'``), the caller must
        emit a ``[WARN]`` so the assumption is recorded rather than silently applied.
    """
    name = (sample_name or "").strip()
    for pattern, kwargs in _CONTROL_RULES:
        if pattern.match(name):
            return ControlClass(**kwargs)
    if _CONTROL_LIKE.search(name):
        return ControlClass(
            samp_category="sample",
            rule="unclassified-control-like",
            warn_reason=f"sample name {name!r} looks like a control but matched no known "
            f"control pattern; left as a biological 'sample' -- verify it is not a control",
        )
    return ControlClass(samp_category="sample", rule="sample")


# --------------------------------------------------------------------------- #
# Canonical manifest row
# --------------------------------------------------------------------------- #
_ISO_DATE = re.compile(r"^\d{4}(?:-(\d{2})(?:-(\d{2}))?)?$")


class SampleManifestRow(BaseModel):
    """One canonical manifest row: a single sample in a single sequencing run/library.

    FAIRe-anchored, strict (``extra="forbid"``): an unknown column name is a typo and
    errors at load (the strict-validation policy). The same biological sample sequenced in two runs is
    two rows (FAIRe ``experimentRunMetadata`` granularity), which is what enables
    DADA2-by-library and per-run control association.
    """

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    # --- identity & event (Darwin Core) ---
    eventID: str = Field(..., min_length=1, description="Canonical sample key (DwC eventID; FAIRe samp_name alias)")
    parentEventID: Optional[str] = Field(default=None, description="DwC parentEventID (site/event grouping)")
    eventDate: Optional[str] = Field(default=None, description="DwC eventDate, ISO-8601 (YYYY[-MM[-DD]])")
    decimalLatitude: Optional[float] = Field(default=None, ge=-90, le=90, allow_inf_nan=False, description="DwC decimalLatitude (WGS84)")
    decimalLongitude: Optional[float] = Field(default=None, ge=-180, le=180, allow_inf_nan=False, description="DwC decimalLongitude (WGS84)")
    geodeticDatum: Optional[str] = Field(default=None, description="DwC geodeticDatum (e.g. EPSG:4326)")
    # No ge=0: marine datasets encode depth as a negative elevation below the surface
    # (e.g. -840 m), so depth may legitimately be negative or positive. inf/nan rejected.
    minimumDepthInMeters: Optional[float] = Field(default=None, allow_inf_nan=False, description="DwC minimumDepthInMeters")
    maximumDepthInMeters: Optional[float] = Field(default=None, allow_inf_nan=False, description="DwC maximumDepthInMeters")
    materialSampleID: Optional[str] = Field(default=None, description="DwC materialSampleID")
    samplingProtocol: Optional[str] = Field(default=None, description="DwC samplingProtocol")

    # --- environment (MIxS) ---
    geo_loc_name: Optional[str] = Field(default=None, description="MIxS geo_loc_name (INSDC country[:region])")
    env_broad_scale: Optional[str] = Field(default=None, description="MIxS env_broad_scale (ENVO biome)")
    env_local_scale: Optional[str] = Field(default=None, description="MIxS env_local_scale (ENVO feature)")
    env_medium: Optional[str] = Field(default=None, description="MIxS env_medium (ENVO material)")
    samp_collect_method: Optional[str] = Field(
        default=None,
        description="MIxS samp_collect_method (canonical spelling; GBIF extension uses "
        "'samp_collec_method' -- translate at export)",
    )
    samp_size: Optional[float] = Field(default=None, ge=0, allow_inf_nan=False, description="MIxS samp_size (amount collected; from 'volume')")
    samp_size_unit: Optional[str] = Field(default=None, description="FAIRe samp_size_unit (mL|L|...)")
    size_frac: Optional[float] = Field(default=None, ge=0, allow_inf_nan=False, description="MIxS size_frac (filter pore size, um)")
    samp_vol_we_dna_ext: Optional[float] = Field(default=None, ge=0, allow_inf_nan=False, description="MIxS samp_vol_we_dna_ext")
    sop: Optional[str] = Field(default=None, description="MIxS sop (protocol URL/DOI)")

    # --- assay / library / run (FAIRe) ---
    target_gene: Optional[str] = Field(default=None, description="MIxS target_gene (marker, e.g. 12S/16S/COI)")
    assay_name: Optional[str] = Field(
        default=None,
        description="FAIRe assay_name (FAIRe-Mandatory; modelled Optional here and populated "
        "from the marker by the migrator -- a documented downgrade)",
    )
    assay_type: str = Field(default="metabarcoding", description="FAIRe assay_type (metabarcoding|targeted)")
    seq_run_id: str = Field(
        ..., min_length=1,
        description="FAIRe seq_run_id: the DADA2-by-library error-model batch key. Required; "
        "the migrator synthesises a single value for one-library datasets with a [WARN] "
        "(never a silent single-batch assumption).",
    )
    lib_id: Optional[str] = Field(default=None, description="FAIRe lib_id (finer per-library grouping)")
    pcr_plate_id: Optional[str] = Field(default=None, description="FAIRe pcr_plate_id")
    mid_forward: Optional[str] = Field(default=None, description="FAIRe mid_forward (forward demux barcode)")
    mid_reverse: Optional[str] = Field(default=None, description="FAIRe mid_reverse (reverse demux barcode)")

    # --- controls (FAIRe two-field model) ---
    samp_category: str = Field(..., description="FAIRe samp_category (sample|negative control|positive control|PCR standard|other:)")
    neg_cont_type: Optional[str] = Field(default=None, description="MIxS neg_cont_type (required if negative control)")
    pos_cont_type: Optional[str] = Field(default=None, description="MIxS pos_cont_type (required if positive control)")
    rel_cont_id: Optional[str] = Field(default=None, description="FAIRe rel_cont_id (|-separated eventIDs of related controls)")
    extraction_ID: Optional[str] = Field(
        default=None,
        description="SeeDNAP-internal extraction batch key (EnviDat distribute-blank key); "
        "the control-association key for Blank-ext. PCR blanks legitimately have none.",
    )

    # --- QC read accounting (pipeline-written, keyed on eventID) ---
    reads_raw: Optional[int] = Field(default=None, ge=0, description="Raw read pairs (pipeline-written)")
    reads_trimmed: Optional[int] = Field(default=None, ge=0, description="After primer trimming (pipeline-written)")
    reads_filtered: Optional[int] = Field(default=None, ge=0, description="After DADA2 filterAndTrim (pipeline-written)")
    reads_denoised: Optional[int] = Field(default=None, ge=0, description="After DADA2 denoising (pipeline-written)")
    reads_merged: Optional[int] = Field(default=None, ge=0, description="After pair merging (pipeline-written)")
    reads_nonchimeric: Optional[int] = Field(default=None, ge=0, description="After chimera removal (pipeline-written)")

    @model_validator(mode="before")
    @classmethod
    def _normalise_missing(cls, data: Any) -> Any:
        """Map empty/NA/INSDC missing tokens to None; strip strings.

        Runs before typing so that ``"NA"`` in a numeric column becomes ``None`` (not a
        string that would later coerce to NaN and skip range checks) and so that the INSDC
        missing-value tokens FAIRe mandates on control rows (e.g.
        ``not applicable: control sample``) are tolerated.

        Args:
            data: The raw input passed to the model; only a dict (one CSV record) is
                rewritten, anything else is returned untouched.

        Returns:
            For a dict input, a new dict where string cells are stripped and any
            missing-value token (compared case-insensitively against
            :data:`MISSING_VALUE_TOKENS`) is replaced with None; otherwise the input
            unchanged.
        """
        if not isinstance(data, dict):
            return data
        cleaned: Dict[str, Any] = {}
        for key, value in data.items():
            if isinstance(value, str):
                stripped = value.strip()
                if stripped.lower() in MISSING_VALUE_TOKENS:
                    cleaned[key] = None
                else:
                    cleaned[key] = stripped
            else:
                cleaned[key] = value
        return cleaned

    @field_validator("samp_category")
    @classmethod
    def _validate_samp_category(cls, v: str) -> str:
        """Validate the FAIRe ``samp_category`` against its controlled vocabulary.

        ``samp_category`` is the top-level FAIRe distinction between a biological sample
        and the various control types; an unrecognised value here would silently mislabel
        a control as a sample (or vice versa).

        Args:
            v: The proposed samp_category value.

        Returns:
            The value unchanged when it is one of :data:`SAMP_CATEGORIES` or an
            ``other:``-prefixed value.

        Raises:
            ValueError: The value is neither in the controlled vocabulary nor an
                ``other:`` value.
        """
        if v in SAMP_CATEGORIES or _is_other(v):
            return v
        raise ValueError(
            f"invalid samp_category {v!r}; expected one of {SAMP_CATEGORIES} or an 'other:' value"
        )

    @field_validator("neg_cont_type")
    @classmethod
    def _validate_neg_cont_type(cls, v: Optional[str]) -> Optional[str]:
        """Validate the MIxS ``neg_cont_type`` against its controlled vocabulary.

        ``neg_cont_type`` records which negative control a blank is (field, extraction,
        PCR, etc.), which determines what contamination it can detect and which samples
        it should be associated with.

        Args:
            v: The proposed neg_cont_type value, or None if absent.

        Returns:
            The value unchanged when it is None, one of :data:`NEG_CONT_TYPES`, or an
            ``other:``-prefixed value.

        Raises:
            ValueError: A non-None value that is neither in the controlled vocabulary nor
                an ``other:`` value.
        """
        if v is None or v in NEG_CONT_TYPES or _is_other(v):
            return v
        raise ValueError(
            f"invalid neg_cont_type {v!r}; expected one of {NEG_CONT_TYPES} or an 'other:' value"
        )

    @field_validator("eventDate")
    @classmethod
    def _validate_event_date(cls, v: Optional[str]) -> Optional[str]:
        """Validate that ``eventDate`` is a real ISO-8601 (partial) calendar date.

        Canonical manifests carry ISO-8601 dates only (the migrator normalises legacy
        dotted forms). Reject anything else loudly rather than risk a silent mis-parse,
        because the collection date flows through to GBIF. The format regex only checks
        digit shape, so the month/day range and full-date calendar validity are checked
        explicitly here: an out-of-range value like ``2024-13`` or ``2024-00-45``, or an
        impossible date like Feb 30, is a silent data-corruption path into GBIF and must
        be rejected.

        Args:
            v: The proposed eventDate, expected as ``YYYY``, ``YYYY-MM`` or
                ``YYYY-MM-DD``, or None if absent.

        Returns:
            The value unchanged when it is None or a valid (partial) ISO-8601 date.

        Raises:
            ValueError: The value is non-ISO-8601 in shape, has an out-of-range month
                (not 01-12) or day (not 01-31), or is a calendar-impossible full date.
        """
        if v is None:
            return v
        m = _ISO_DATE.match(v)
        if not m:
            raise ValueError(
                f"eventDate {v!r} is not ISO-8601 (YYYY[-MM[-DD]]). Legacy dotted dates must be "
                f"normalised by the migrator; a hand-authored manifest must use ISO dates."
            )
        month, day = m.group(1), m.group(2)
        if month is not None and not (1 <= int(month) <= 12):
            raise ValueError(
                f"eventDate {v!r} has an out-of-range month {month!r} (expected 01-12)."
            )
        if day is not None and not (1 <= int(day) <= 31):
            raise ValueError(
                f"eventDate {v!r} has an out-of-range day {day!r} (expected 01-31)."
            )
        # Full YYYY-MM-DD: reject calendar-impossible dates (e.g. Feb 30, Apr 31).
        if month is not None and day is not None:
            from datetime import date as _date

            try:
                _date(int(v[:4]), int(month), int(day))
            except ValueError as exc:
                raise ValueError(
                    f"eventDate {v!r} is not a real calendar date ({exc})."
                ) from exc
        return v

    @model_validator(mode="after")
    def _validate_conditional_requirements(self) -> "SampleManifestRow":
        """Enforce FAIRe Mandatory-if rules tying a control category to its control type.

        FAIRe requires that a negative control declare which negative it is
        (``neg_cont_type``) and a positive control declare its ``pos_cont_type``; without
        these a control cannot be correctly associated with the samples it guards.
        Note: eventDate is *not* hard-required here. A manifest can legitimately be built
        at the demux stage (identity + library + tag) before field metadata supplies dates;
        a missing eventDate on a biological sample is surfaced as a loud
        :meth:`SampleManifest.check_completeness` ``[WARN]``, not a constructor failure.

        Returns:
            The validated model instance (``self``), unchanged.

        Raises:
            ValueError: A negative control is missing ``neg_cont_type``, or a positive
                control is missing ``pos_cont_type`` (message names the offending
                eventID).
        """
        if self.samp_category == "negative control" and not self.neg_cont_type:
            raise ValueError(
                f"negative control {self.eventID!r} is missing neg_cont_type "
                f"(FAIRe Mandatory-if)"
            )
        if self.samp_category == "positive control" and not self.pos_cont_type:
            raise ValueError(
                f"positive control {self.eventID!r} is missing pos_cont_type "
                f"(FAIRe Mandatory-if)"
            )
        return self

    @property
    def is_control(self) -> bool:
        """Report whether this row is a control rather than a biological field sample.

        Any non-``sample`` FAIRe ``samp_category`` (negative control, positive control,
        or PCR standard) is a control well used to detect contamination or calibrate,
        not a collected environmental specimen.

        Returns:
            True for any ``samp_category`` other than ``sample``.
        """
        return self.samp_category != "sample"


# --------------------------------------------------------------------------- #
# Manifest collection
# --------------------------------------------------------------------------- #
# Canonical column order for serialisation = model field declaration order.
MANIFEST_COLUMNS: Tuple[str, ...] = tuple(SampleManifestRow.model_fields.keys())


class SampleManifest:
    """A validated collection of :class:`SampleManifestRow`, plus provenance + accessors."""

    def __init__(self, rows: List[SampleManifestRow], source: Optional[Path] = None) -> None:
        """Hold validated rows and the optional source path they came from.

        Args:
            rows: The already-validated manifest rows, one per sample-library.
            source: The CSV path the rows were loaded from, kept for provenance; None
                when the manifest was assembled in memory (e.g. by the migrator).
        """
        self.rows = rows
        self.source = source

    def __len__(self) -> int:
        """Return the number of rows (sample-libraries) in the manifest.

        Returns:
            The row count.
        """
        return len(self.rows)

    # -- accessors --------------------------------------------------------- #
    def event_ids(self) -> List[str]:
        """Return the eventID (canonical sample key) of every row, in row order.

        Returns:
            A list of eventIDs; may contain repeats if a sample appears in more than
            one sequencing run.
        """
        return [r.eventID for r in self.rows]

    def controls(self) -> List[SampleManifestRow]:
        """Return the control rows (any non-``sample`` samp_category).

        Controls are negative/positive controls and PCR standards: non-biological wells
        used to assess contamination and calibration.

        Returns:
            The subset of rows that are controls, in row order.
        """
        return [r for r in self.rows if r.is_control]

    def biological_samples(self) -> List[SampleManifestRow]:
        """Return the biological (non-control) sample rows.

        Returns:
            The subset of rows whose ``samp_category`` is ``sample``, in row order.
        """
        return [r for r in self.rows if not r.is_control]

    def seq_run_ids(self) -> List[str]:
        """Return the distinct seq_run_id values in first-seen order.

        ``seq_run_id`` is the DADA2-by-library error-model batch key: samples sharing a
        run are denoised together, so this is the set of error-model batches.

        Returns:
            The unique seq_run_id values, ordered by first appearance.
        """
        seen: Dict[str, None] = {}
        for r in self.rows:
            seen.setdefault(r.seq_run_id, None)
        return list(seen)

    # -- serialisation ----------------------------------------------------- #
    def to_dataframe(self) -> pd.DataFrame:
        """Serialise the manifest to a DataFrame in canonical column order.

        Returns:
            A DataFrame with one row per sample-library and columns equal to
            :data:`MANIFEST_COLUMNS` (the model field declaration order). Missing
            optional values appear as None.
        """
        df = pd.DataFrame([r.model_dump() for r in self.rows], columns=list(MANIFEST_COLUMNS))
        return df

    def to_csv(self, path: Path) -> Path:
        """Write the manifest to ``path`` as CSV, creating parent directories.

        Args:
            path: Destination CSV path; parent directories are created if absent.

        Returns:
            The path written to (the same ``path``, as a :class:`Path`).
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.to_dataframe().to_csv(path, index=False)
        logger.info(f"Wrote sample manifest ({len(self.rows)} rows) to {path}")
        return path

    def check_controls(self) -> None:
        """Emit a ``[WARN]`` when the manifest carries no control samples.

        A run with no negative controls cannot have its contamination assessed, so the
        absence is logged loudly rather than passing silently (the no-silent-fallbacks
        policy).

        Returns:
            None. Side effect only: logs a ``[WARN]`` if no controls are present.
        """
        if not self.controls():
            logger.warning(
                "[WARN] manifest: expected=at least one negative/positive control sample, "
                "got=none detected, fallback=proceeding (contamination cannot be assessed "
                "for this run; the no-silent-fallbacks policy)"
            )

    def check_completeness(self) -> None:
        """Emit a ``[WARN]`` for biological samples missing ``eventDate``.

        eventDate (collection date) is required before GBIF export but may legitimately
        be absent in a manifest built at the demux stage, before field metadata is
        joined; this surfaces the gap loudly without blocking construction. Up to five
        offending eventIDs are named, with a count of any remainder.

        Returns:
            None. Side effect only: logs a ``[WARN]`` if any biological sample lacks an
            eventDate.
        """
        missing = [r.eventID for r in self.biological_samples() if not r.eventDate]
        if missing:
            shown = missing[:5]
            more = "" if len(missing) <= 5 else f", +{len(missing) - 5} more"
            logger.warning(
                f"[WARN] manifest: expected=eventDate on every biological sample, "
                f"got={len(missing)} without one ({shown}{more}), fallback=kept "
                f"(eventDate is required before GBIF export)"
            )


def load_manifest(path: Path) -> SampleManifest:
    """Read and strictly validate a canonical manifest CSV.

    Reads with ``utf-8-sig`` (tolerates a BOM) and ``keep_default_na=False`` so the
    row model -- not pandas -- owns null handling. Every row is validated through
    :class:`SampleManifestRow`; all row errors are aggregated into one message naming the
    offending rows. Duplicate ``(eventID, seq_run_id)`` keys raise (one row per
    sample-library, since a sample appears once per sequencing run). Emits a ``[WARN]`` if
    the manifest has no controls and a ``[WARN]`` if any biological sample lacks an
    eventDate.

    Args:
        path: Path to the canonical manifest CSV.

    Returns:
        A :class:`SampleManifest` of validated rows, with ``source`` set to ``path``.

    Raises:
        FileNotFoundError: the manifest does not exist.
        ValueError: the manifest is empty, has duplicate ``(eventID, seq_run_id)`` keys,
            or any row fails validation (the message lists every bad row).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")

    df = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    if df.empty:
        raise ValueError(f"Manifest is empty (no data rows): {path}")

    rows: List[SampleManifestRow] = []
    errors: List[str] = []
    for i, record in enumerate(df.to_dict(orient="records")):
        try:
            rows.append(SampleManifestRow(**record))
        except Exception as exc:  # pydantic ValidationError or ValueError
            ident = record.get("eventID") or f"<row {i + 2}>"  # +2: header + 1-based
            errors.append(f"  - row {i + 2} ({ident}): {exc}")

    if errors:
        raise ValueError(
            f"Manifest validation failed for {path} ({len(errors)} bad row(s)):\n"
            + "\n".join(errors)
        )

    # One row per sample-library: (eventID, seq_run_id) must be unique.
    key_counts = Counter((r.eventID, r.seq_run_id) for r in rows)
    dupes = sorted(k for k, n in key_counts.items() if n > 1)
    if dupes:
        raise ValueError(
            f"Manifest {path} has duplicate (eventID, seq_run_id) keys: {dupes}. "
            f"A sample appears once per sequencing run."
        )

    manifest = SampleManifest(rows=rows, source=path)
    manifest.check_controls()
    manifest.check_completeness()
    return manifest


# --------------------------------------------------------------------------- #
# Cross-CSV eventID validation (E2)
# --------------------------------------------------------------------------- #
class AbundanceValidationResult(BaseModel):
    """Outcome of validating a manifest's eventIDs against an abundance table."""

    model_config = ConfigDict(extra="forbid")

    abundance_samples: List[str] = Field(description="Sample columns found in the abundance table")
    orphan_abundance_columns: List[str] = Field(
        description="Sample columns in the abundance table with NO manifest row "
        "(the dangerous silent-ID-mismatch case, e.g. an unlisted Blank-PCR-3)"
    )
    manifest_extra_rows: List[str] = Field(
        description="Manifest eventIDs absent from the abundance table (harmless: dropped by "
        "the pipeline, e.g. a control filtered to nothing)"
    )

    @property
    def ok(self) -> bool:
        """Report whether every abundance sample column has a matching manifest row.

        An orphan column (a sample in the OTU/abundance table with no manifest row) is
        the dangerous silent-ID-mismatch case; its absence means the cross-check passed.

        Returns:
            True when there are no orphan abundance columns.
        """
        return not self.orphan_abundance_columns


# Non-sample metadata columns commonly present in SeeDNAP abundance/OTU tables; these are
# not biological sample columns and must be excluded before the eventID comparison.
_ABUNDANCE_META_COLUMNS = frozenset(
    {
        "", "otu", "otu_id", "asv", "asv_id", "id", "seed", "amplicon", "sequence", "seq",
        "total", "total_reads", "nb_reads", "nb_reads_total", "size", "length", "chimera",
        "cluster", "n_samples", "spread", "abundance", "identity", "taxonomy",
        "kingdom", "phylum", "class", "order", "family", "genus", "species", "taxon", "rank",
    }
)


def _abundance_sample_columns(
    abundance_csv: Path,
    id_column: Optional[str] = None,
    known_event_ids: Optional[set] = None,
) -> List[str]:
    """Return the per-sample columns of an abundance/OTU table.

    SeeDNAP OTU tables are sequences x samples (one row per OTU/ASV, one column per
    sample plus OTU-level metadata such as taxonomy or total read count). The sample
    columns are every column that is not OTU metadata. A column whose name matches a
    manifest ``eventID`` is ALWAYS kept (a real sample literally named ``total``/``order``
    must not be silently dropped by the metadata denylist, and if such a name is also a
    meta token a ``[WARN]`` is logged). Reads only the header (``nrows=0``).

    Args:
        abundance_csv: Path to the abundance/OTU CSV; read with ``utf-8-sig``.
        id_column: Name of the OTU/ASV identifier column to exclude, if known; when None,
            identifier columns are excluded via the metadata denylist instead.
        known_event_ids: The manifest eventIDs; any column matching one of these is kept
            as a sample column even if it collides with a metadata token.

    Returns:
        The list of column names judged to be per-sample columns, in header order.
    """
    known = known_event_ids or set()
    header = pd.read_csv(abundance_csv, nrows=0, encoding="utf-8-sig")
    cols = [str(c).strip() for c in header.columns]
    samples = []
    for c in cols:
        if id_column is not None and c == id_column:
            continue
        if c.lower() in _ABUNDANCE_META_COLUMNS:
            if c in known:
                logger.warning(
                    f"[WARN] manifest_vs_abundance: abundance column {c!r} matches a manifest "
                    f"eventID but also an OTU-table metadata token; kept as a sample column"
                )
            else:
                continue
        samples.append(c)
    return samples


def validate_against_abundance(
    manifest: SampleManifest,
    abundance_csv: Path,
    *,
    id_column: Optional[str] = None,
    raise_on_orphan: bool = False,
) -> AbundanceValidationResult:
    """Cross-check the manifest's ``eventID`` set against an abundance table's sample columns.

    This is the up-front silent-ID-mismatch guard (the no-silent-fallbacks policy): an
    abundance column with no manifest row (an *orphan*, e.g. a ``Blank-PCR-3`` present in
    the OTU table but absent from the field metadata) would otherwise be silently treated
    as an unknown field sample, injecting contamination reads into the dataset. Orphans
    are reported and ``[WARN]``-logged (or raised if ``raise_on_orphan``); manifest rows
    absent from the abundance table are reported separately as harmless (the pipeline may
    legitimately drop a control to zero reads).

    Args:
        manifest: The validated manifest whose eventID set is the source of truth.
        abundance_csv: Path to the abundance/OTU CSV to check against (read with
            ``utf-8-sig``).
        id_column: Name of the OTU/ASV identifier column to exclude from the sample
            columns, if known.
        raise_on_orphan: When True, raise instead of only warning if any orphan column
            is found.

    Returns:
        An :class:`AbundanceValidationResult` with the detected sample columns, the
        orphan columns (abundance samples missing from the manifest), and the extra
        manifest rows (eventIDs absent from the abundance table).

    Raises:
        FileNotFoundError: the abundance table does not exist.
        ValueError: ``raise_on_orphan`` is True and one or more orphan columns are found.
    """
    abundance_csv = Path(abundance_csv)
    if not abundance_csv.exists():
        raise FileNotFoundError(f"Abundance table not found: {abundance_csv}")

    manifest_ids = set(manifest.event_ids())
    samples = _abundance_sample_columns(abundance_csv, id_column=id_column, known_event_ids=manifest_ids)
    sample_set = set(samples)

    orphans = sorted(sample_set - manifest_ids)
    extra = sorted(manifest_ids - sample_set)

    if orphans:
        msg = (
            f"[WARN] manifest_vs_abundance: expected=every abundance sample column has a "
            f"manifest row, got={len(orphans)} orphan column(s) with no manifest row "
            f"({orphans}), fallback=" + ("raising" if raise_on_orphan else "reported, not dropped")
        )
        if raise_on_orphan:
            logger.error(msg)
            raise ValueError(
                f"Abundance table {abundance_csv} has sample columns absent from the manifest: "
                f"{orphans}. Add them to the manifest (with the correct samp_category) or "
                f"correct the abundance table."
            )
        logger.warning(msg)
    if extra:
        logger.info(
            f"manifest_vs_abundance: {len(extra)} manifest eventID(s) absent from the "
            f"abundance table (likely dropped by the pipeline): {extra}"
        )

    return AbundanceValidationResult(
        abundance_samples=samples,
        orphan_abundance_columns=orphans,
        manifest_extra_rows=extra,
    )
