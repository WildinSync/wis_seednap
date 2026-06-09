"""Derive a canonical :class:`~seednap.config.manifest.SampleManifest` from today's CSVs.

This is the bridge from SeeDNAP's three existing, messy metadata CSVs to the single
canonical manifest. It is deliberately separate from :mod:`seednap.config.manifest` (which
owns the stable canonical schema) because all the era-specific heuristics live here:

* read with ``utf-8-sig`` (a BOM glued to ``eventID`` would otherwise break every join),
* canonicalise header casing and strip unit suffixes (``Site_names`` vs ``site_names``,
  ``conductivity [µS]``/``[μS]`` -- two distinct Unicode mu codepoints),
* normalise ``eventDate`` per-file, distinguishing ``YYYY.MM.DD`` from ``DD.MM.YYYY``
  (the same dataset can use different orders across markers) and **raising** on a genuinely
  ambiguous order rather than guessing,
* classify controls by name pattern across eras via the single
  :func:`~seednap.config.manifest.classify_control`,
* synthesise ``seq_run_id`` for modern pre-demultiplexed datasets (one library per
  dataset) with a ``[WARN]`` -- never a silent single-batch assumption,
* recover ``seq_run_id`` (``library``) and ``mid_forward`` (the ``tag_demultiplex`` family,
  including the legacy ``tag_demltiplex`` typo) from a legacy demux ``metadata_lab_*.csv``.

Every dropped column, ambiguous control, missing/orphan ``extraction_ID``, defaulted
``seq_run_id`` and zero-control dataset is surfaced as a ``[WARN]`` (the no-silent-fallbacks policy); a
missing sample key or an unresolvable date order raises naming the file.

Environmental measurement columns with no FAIRe slot (pH, conductivity, weather, lab codes
EVE/EVS, transect-end coordinates, ...) are intentionally not carried into the manifest:
they are not consumed by GBIF export or by the grouping/control/QC logic the manifest
exists to serve, and the source CSV is preserved. They are listed in a ``[WARN]`` so the
drop is on the record, never silent.
"""

from __future__ import annotations

import logging
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from seednap.config.manifest import (
    MISSING_VALUE_TOKENS,
    SampleManifest,
    SampleManifestRow,
    classify_control,
)

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------- #
# Header canonicalisation
# --------------------------------------------------------------------------- #
# A trailing bracketed unit, e.g. "conductivity [µS]" / "River width [m]".
_UNIT_SUFFIX = re.compile(r"\s*\[[^\]]*\]\s*$")


def _canon_header(raw: str) -> str:
    """Fold a raw field-metadata header to a lowercase match key.

    Strips a BOM, a trailing ``[unit]`` suffix, collapses internal whitespace, and
    lower-cases. Used only to *match* a header to a canonical manifest field; the
    canonical (camelCase) name is taken from :data:`_FIELD_ALIASES`.
    """
    h = raw.lstrip("﻿").strip()
    h = _UNIT_SUFFIX.sub("", h)
    h = re.sub(r"\s+", " ", h)
    return h.lower()


def _cell(value: object) -> Optional[str]:
    """Strip a raw cell and map missing-value tokens (``NA``, empty, INSDC tokens) to None,
    so the migrator's bookkeeping never treats ``NA`` as a real value."""
    s = str(value).strip()
    return None if s.lower() in MISSING_VALUE_TOKENS else s


# Lowercased match-key -> canonical manifest field. Field metadata columns only.
_FIELD_ALIASES: Dict[str, str] = {
    "eventid": "eventID",
    "samp_name": "eventID",
    "parenteventid": "parentEventID",
    "eventdate": "eventDate",
    "decimallatitude": "decimalLatitude",
    "decimallongitude": "decimalLongitude",
    "geodeticdatum": "geodeticDatum",
    "minimumdepthinmeters": "minimumDepthInMeters",
    "maximumdepthinmeters": "maximumDepthInMeters",
    "depth": "maximumDepthInMeters",
    "materialsampleid": "materialSampleID",
    "samplingprotocol": "samplingProtocol",
    "geo_loc_name": "geo_loc_name",
    "env_broad_scale": "env_broad_scale",
    "env_local_scale": "env_local_scale",
    "env_medium": "env_medium",
    "samp_collect_method": "samp_collect_method",
    "samp_size": "samp_size",
    "volume": "samp_size",  # legacy alias already used by the GBIF builder
    "samp_size_unit": "samp_size_unit",
    "size_frac": "size_frac",
    "samp_vol_we_dna_ext": "samp_vol_we_dna_ext",
    "sop": "sop",
    "target_gene": "target_gene",
    "extraction_id": "extraction_ID",
    # library / run grouping + demux barcodes (present in demux lab CSVs, absent in modern
    # pre-demultiplexed field metadata). Recognised so a lab CSV migrates directly.
    "seq_run_id": "seq_run_id",
    "library": "seq_run_id",
    "lib_id": "lib_id",
    "pcr_plate_id": "pcr_plate_id",
    "mid_forward": "mid_forward",
    "tag": "mid_forward",
    "tag_demultiplex": "mid_forward",
    "tag_demultiplex_f": "mid_forward",
    "tag_demultiplex_i7": "mid_forward",
    "tag_demltiplex": "mid_forward",  # documented legacy typo
    "mid_reverse": "mid_reverse",
    "tag_demultiplex_r": "mid_reverse",
    "tag_demultiplex_i5": "mid_reverse",
}

# Recognised lab columns with no FAIRe slot: dropped, but expected (a single summary
# [WARN] lists them). Anything dropped that is NOT here gets an individual louder [WARN]
# because it may signal a shifted/garbled file (e.g. ogooue's stray "Test" column).
_KNOWN_UNMAPPED = frozenset(
    {
        "decimallatitude_end", "decimallongitude_end", "start_time", "duration",
        "institution", "device", "filter", "laboratory", "ecosystem", "body",
        "polygon_id", "area_basin", "eve", "evs", "site_names", "conductivity",
        "temperature", "rain occurrence before sampling", "habitat", "weather", "ph",
        "river width", "river depth",
    }
)


def _build_header_map(
    columns: List[str],
) -> Tuple[Dict[str, str], List[str], List[str], List[Tuple[str, str, str]]]:
    """Map raw field-metadata headers to canonical fields.

    Returns ``(canonical_field -> raw_column, dropped_known, dropped_unexpected, collisions)``.
    A ``pcr_code_*`` column is recognised (carried in ``dropped_known``) regardless of its
    marker token. If two raw columns canonicalise to the same field the first wins and the
    second is recorded in ``collisions`` (as ``(dropped_raw, kept_raw, canonical)``) so the
    caller can ``[WARN]`` -- never a silent discard (the no-silent-fallbacks policy).
    """
    field_to_raw: Dict[str, str] = {}
    dropped_known: List[str] = []
    dropped_unexpected: List[str] = []
    collisions: List[Tuple[str, str, str]] = []
    for raw in columns:
        key = _canon_header(raw)
        if key in _FIELD_ALIASES:
            canonical = _FIELD_ALIASES[key]
            if canonical in field_to_raw:
                collisions.append((raw, field_to_raw[canonical], canonical))
            else:
                field_to_raw[canonical] = raw
        elif key in _KNOWN_UNMAPPED or key.startswith("pcr_code_") or key.startswith("pcr_primer_"):
            dropped_known.append(raw)
        else:
            dropped_unexpected.append(raw)
    return field_to_raw, dropped_known, dropped_unexpected, collisions


# --------------------------------------------------------------------------- #
# Date normalisation (the silent-corruption guard)
# --------------------------------------------------------------------------- #
_ISO_FULL = re.compile(r"^(\d{4})-(\d{2})-(\d{2})$")
_ISO_YM = re.compile(r"^(\d{4})-(\d{2})$")
_YEAR = re.compile(r"^(\d{4})$")
# A separator-delimited 3-part date; separator may be "." or "/". A trailing clock time
# (e.g. " 00:00:00") is stripped first. Two-digit years are intentionally NOT matched here
# (they are unsafe to disambiguate -- see below).
_TIME_SUFFIX = re.compile(r"[ T]\d{1,2}:\d{2}(:\d{2})?\s*$")
_THREE_PART = re.compile(r"^(\d{1,4})[./](\d{1,2})[./](\d{1,4})$")
_TWO_PART = re.compile(r"^(\d{1,4})[./](\d{1,4})$")


_DATE_ORDERS = {
    # order -> (year_index, month_index, day_index) into the 3 positional fields
    "YMD": (0, 1, 2),
    "DMY": (2, 1, 0),
    "MDY": (2, 0, 1),
}


def normalise_event_dates(
    values: List[str], *, context: str = "eventDate", order: Optional[str] = None
) -> Dict[str, str]:
    """Build a verbatim-token -> ISO-8601 map for one file's eventDate column.

    Handles ISO (``YYYY-MM-DD``), dotted/slashed 3-part dates in any of YYYY.MM.DD,
    YYYY.DD.MM, DD.MM.YYYY and MM.DD.YYYY orders, an optional trailing clock time, and
    partial year / year-month forms. The 4-digit year fixes which end is the year; the
    day/month order of the remaining two fields is then resolved **per file** from any
    unambiguous token (a field > 12 is the day). Missing-value tokens (``NA``, empty, INSDC
    tokens) are skipped, not parsed.

    Args:
        order: optional explicit field order ``YMD``/``DMY``/``MDY`` to use when the dates
            are genuinely ambiguous (every token has day and month both <= 12). When given,
            it is applied to every 3-part token and the choice is logged; when ``None`` the
            order is auto-detected and an ambiguous file *raises* rather than being guessed.

    Raises -- never guesses -- when a file mixes year-first and year-last tokens, when the
    day/month order is contradictory, when it is ambiguous and no ``order`` was supplied,
    when a year is only two digits, or when a field is out of range.

    Returns:
        Map from each distinct verbatim token to its ISO-8601 form.
    """
    if order is not None and order.upper() not in _DATE_ORDERS:
        raise ValueError(f"{context}: unknown date order {order!r}; expected one of {list(_DATE_ORDERS)}")
    tokens = sorted({v.strip() for v in values if v and v.strip()})
    # Drop missing-value tokens up front so a mixed valid/NA column is not rejected.
    tokens = [t for t in tokens if t.lower() not in MISSING_VALUE_TOKENS]
    if not tokens:
        return {}

    simple: Dict[str, str] = {}                       # already-ISO / partial -> ISO
    three: List[Tuple[str, int, int, int]] = []       # (tok, f0, f1, f2) raw 3-part fields
    year_positions: set = set()                       # 0 (year-first) or 2 (year-last)
    unparseable: List[str] = []

    for tok in tokens:
        core = _TIME_SUFFIX.sub("", tok).strip()
        if _ISO_FULL.match(core) or _YEAR.match(core):
            simple[tok] = core
        elif (m := _ISO_YM.match(core)):
            simple[tok] = f"{int(m.group(1)):04d}-{int(m.group(2)):02d}"
        elif (m := _THREE_PART.match(core)):
            f0, f1, f2 = int(m.group(1)), int(m.group(2)), int(m.group(3))
            yp = 0 if len(m.group(1)) == 4 else (2 if len(m.group(3)) == 4 else None)
            if yp is None:
                unparseable.append(tok)  # 2-digit year -> unsafe, refuse
            else:
                year_positions.add(yp)
                three.append((tok, f0, f1, f2))
        elif (m := _TWO_PART.match(core)):
            a, b = m.group(1), m.group(2)
            if len(a) == 4:
                simple[tok] = f"{int(a):04d}-{int(b):02d}"      # YYYY.MM
            elif len(b) == 4:
                simple[tok] = f"{int(b):04d}-{int(a):02d}"      # MM.YYYY
            else:
                unparseable.append(tok)
        else:
            unparseable.append(tok)

    if unparseable:
        raise ValueError(
            f"{context}: unrecognised/ambiguous date format(s) {unparseable[:5]}; expected "
            f"ISO-8601 (YYYY-MM-DD) or a dotted/slashed date with a 4-digit year."
        )

    mapping: Dict[str, str] = dict(simple)
    if not three:
        return mapping

    # Explicit user-supplied order: apply directly (the choice is logged, not guessed).
    if order is not None:
        yi, mi, di = _DATE_ORDERS[order.upper()]
        for tok, *fields in three:
            year, month, day = fields[yi], fields[mi], fields[di]
            if len(str(year)) != 4 or not (1 <= month <= 12 and 1 <= day <= 31):
                raise ValueError(
                    f"{context}: token {tok!r} is inconsistent with the supplied order {order!r}"
                )
            mapping[tok] = f"{year:04d}-{month:02d}-{day:02d}"
        logger.warning(
            f"[WARN] {context}: expected=ISO-8601 dates, got=non-ISO dates, "
            f"fallback=parsed with the user-supplied order {order.upper()} and normalised to ISO"
        )
        return mapping

    if len(year_positions) > 1:
        raise ValueError(
            f"{context}: file mixes year-first and year-last dates "
            f"(e.g. {three[0][0]!r}) -- refusing to guess; normalise the source file."
        )
    year_pos = next(iter(year_positions))

    # The two non-year fields, in positional order (a, b).
    def _nonyear(f0: int, f1: int, f2: int) -> Tuple[int, int]:
        """Return the two non-year positional fields (a, b), dropping whichever holds the year."""
        return (f1, f2) if year_pos == 0 else (f0, f1)

    day_first_evidence = any(_nonyear(f0, f1, f2)[0] > 12 for _, f0, f1, f2 in three)   # a is day
    month_first_evidence = any(_nonyear(f0, f1, f2)[1] > 12 for _, f0, f1, f2 in three)  # b is day
    if day_first_evidence and month_first_evidence:
        raise ValueError(
            f"{context}: dates give contradictory day/month order (e.g. {three[0][0]!r}) -- "
            f"refusing to guess."
        )
    if not day_first_evidence and not month_first_evidence:
        raise ValueError(
            f"{context}: dates like {three[0][0]!r} are ambiguous (day and month both <= 12) "
            f"-- cannot tell the day/month order. Normalise the source file."
        )
    day_first = day_first_evidence  # exactly one is true

    for tok, f0, f1, f2 in three:
        year = f0 if year_pos == 0 else f2
        a, b = _nonyear(f0, f1, f2)
        day, month = (a, b) if day_first else (b, a)
        if not (1 <= month <= 12 and 1 <= day <= 31):
            raise ValueError(f"{context}: invalid month/day in {tok!r}")
        mapping[tok] = f"{year:04d}-{month:02d}-{day:02d}"

    order = ("YYYY." if year_pos == 0 else "") + ("DD.MM" if day_first else "MM.DD") + \
            ("" if year_pos == 0 else ".YYYY")
    logger.warning(
        f"[WARN] {context}: expected=ISO-8601 dates, got=non-ISO dates, "
        f"fallback=parsed as {order} (resolved from unambiguous tokens) and normalised to ISO"
    )
    return mapping


# --------------------------------------------------------------------------- #
# Numeric value cleaning (robust, never silent)
# --------------------------------------------------------------------------- #
# manifest numeric field -> (low, high) bound for a sanity range check (None = unbounded).
_NUMERIC_RANGES: Dict[str, Tuple[Optional[float], Optional[float]]] = {
    "decimalLatitude": (-90.0, 90.0),
    "decimalLongitude": (-180.0, 180.0),
    "minimumDepthInMeters": (None, None),   # marine depths may be negative
    "maximumDepthInMeters": (None, None),
    "samp_size": (0.0, None),
    "size_frac": (0.0, None),
    "samp_vol_we_dna_ext": (0.0, None),
}
_LEADING_NUMBER = re.compile(r"^(-?\d+(?:\.\d+)?)\s*\S.*$")  # "17 L" -> "17"


def _clean_numeric(field: str, value: str, event_id: str) -> Optional[str]:
    """Coerce one numeric cell to a clean float-string, or None, never crashing the row.

    Strips a trailing unit ("17 L" -> "17"), range-checks, and on any failure returns None
    with a ``[WARN]`` rather than letting one bad cell abort the whole dataset's migration.
    """
    v = (value or "").strip()
    if not v or v.lower() in MISSING_VALUE_TOKENS:
        return None
    lo, hi = _NUMERIC_RANGES.get(field, (None, None))
    cleaned = v
    try:
        f = float(v)
    except ValueError:
        m = _LEADING_NUMBER.match(v)
        if not m:
            logger.warning(
                f"[WARN] manifest_migrate: expected=numeric {field} for {event_id!r}, "
                f"got={v!r}, fallback=null"
            )
            return None
        cleaned = m.group(1)
        f = float(cleaned)
        logger.warning(
            f"[WARN] manifest_migrate: {field} for {event_id!r}: stripped non-numeric text "
            f"from {v!r}, using {cleaned}"
        )
    if not math.isfinite(f):
        logger.warning(
            f"[WARN] manifest_migrate: expected=a finite {field} for {event_id!r}, "
            f"got={v!r} (inf/nan), fallback=null"
        )
        return None
    if (lo is not None and f < lo) or (hi is not None and f > hi):
        logger.warning(
            f"[WARN] manifest_migrate: expected={field} in [{lo}, {hi}] for {event_id!r}, "
            f"got={v!r} (out of range), fallback=null"
        )
        return None
    return cleaned


# --------------------------------------------------------------------------- #
# Legacy demux lab metadata (tag / library recovery)
# --------------------------------------------------------------------------- #
_TAG_FORWARD = ("tag_demultiplex", "tag_demultiplex_f", "tag_demultiplex_forward",
                "tag_demultiplex_i7", "mid_forward")
_TAG_REVERSE = ("tag_demultiplex_r", "tag_demultiplex_reverse", "tag_demultiplex_i5", "mid_reverse")
_TAG_TYPO = "tag_demltiplex"  # documented legacy misspelling
# A column that *looks* like a demux barcode column (so an unrecognised variant warns
# rather than being silently ignored).
_TAG_LIKE = re.compile(r"(tag|mid|barcode|demult|demlt|index|i7|i5)", re.I)


def _resolve_lab_columns(columns: List[str]) -> Dict[str, Optional[str]]:
    """Locate the eventID, library, and forward/reverse tag columns in a demux lab CSV.

    Unrecognised tag-looking columns emit a ``[WARN]`` rather than being silently dropped.
    """
    by_key = {_canon_header(c): c for c in columns}
    fwd = next((by_key[k] for k in _TAG_FORWARD if k in by_key), None)
    if fwd is None and _TAG_TYPO in by_key:
        fwd = by_key[_TAG_TYPO]
        logger.warning(
            f"[WARN] manifest_migrate: expected=a 'tag_demultiplex' barcode column, "
            f"got=the legacy typo {_TAG_TYPO!r}, fallback=using it as mid_forward"
        )
    rev = next((by_key[k] for k in _TAG_REVERSE if k in by_key), None)
    resolved = {by_key.get("eventid"), by_key.get("library"), fwd, rev}
    for key, col in by_key.items():
        if col in resolved:
            continue
        if _TAG_LIKE.search(key) and not key.startswith("pcr_primer"):
            logger.warning(
                f"[WARN] manifest_migrate: expected=a recognised demux barcode column, "
                f"got=unresolved tag-like column {col!r}, fallback=ignored (no mid_forward/"
                f"mid_reverse from it). Add it to the tag-column aliases if it is a barcode."
            )
    return {
        "eventID": by_key.get("eventid"),
        "library": by_key.get("library"),
        "mid_forward": fwd,
        "mid_reverse": rev,
    }


def _load_lab_metadata(lab_csv: Path) -> Dict[str, Dict[str, Optional[str]]]:
    """Read a demux lab CSV into ``{eventID: {seq_run_id, mid_forward, mid_reverse}}``."""
    df = pd.read_csv(lab_csv, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    cols = _resolve_lab_columns(list(df.columns))
    if cols["eventID"] is None:
        raise ValueError(f"Demux lab metadata {lab_csv} has no eventID column")
    if cols["library"] is None:
        raise ValueError(
            f"Demux lab metadata {lab_csv} has no 'library' column -- cannot recover the "
            f"seq_run_id grouping. (A field-metadata CSV with no library is the "
            f"pre-demultiplexed case; pass it as the field metadata instead.)"
        )
    out: Dict[str, Dict[str, Optional[str]]] = {}
    for rec in df.to_dict(orient="records"):
        ev = _cell(rec[cols["eventID"]])
        if not ev:
            continue
        out[ev] = {
            "seq_run_id": _cell(rec[cols["library"]]),
            "mid_forward": _cell(rec[cols["mid_forward"]]) if cols["mid_forward"] else None,
            "mid_reverse": _cell(rec[cols["mid_reverse"]]) if cols["mid_reverse"] else None,
        }
    return out


# --------------------------------------------------------------------------- #
# Project metadata
# --------------------------------------------------------------------------- #
def _read_marker(project_csv: Path) -> Optional[str]:
    """Pull the marker/target_gene from a one-row project metadata CSV."""
    df = pd.read_csv(project_csv, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    by_key = { _canon_header(c): c for c in df.columns }
    col = by_key.get("marker") or by_key.get("target_gene")
    if col is None or df.empty:
        logger.warning(
            f"[WARN] manifest_migrate: expected=a 'marker' column in {project_csv}, "
            f"got=none, fallback=target_gene/assay_name left unset"
        )
        return None
    val = str(df.iloc[0][col]).strip()
    return val or None


# --------------------------------------------------------------------------- #
# Migration
# --------------------------------------------------------------------------- #
def _dataset_from_filename(field_csv: Path) -> str:
    """``metadata_field_<dataset>.csv`` -> ``<dataset>`` (else the file stem)."""
    stem = field_csv.stem
    for prefix in ("metadata_field_", "metadata_"):
        if stem.startswith(prefix):
            return stem[len(prefix):]
    return stem


def migrate_to_manifest(
    field_csv: Path,
    *,
    project_csv: Optional[Path] = None,
    lab_csv: Optional[Path] = None,
    seq_run_id: Optional[str] = None,
    target_gene: Optional[str] = None,
    dataset: Optional[str] = None,
    date_order: Optional[str] = None,
) -> SampleManifest:
    """Derive a canonical :class:`SampleManifest` from today's field/project/lab CSVs.

    Args:
        field_csv: per-sample field metadata (``metadata_field_*.csv``). Required.
        project_csv: one-row project metadata (``metadata_proj_*.csv``); supplies the
            marker -> ``target_gene``/``assay_name``.
        lab_csv: legacy demux metadata (``metadata_lab_*.csv``); supplies ``seq_run_id``
            (``library``) and ``mid_forward``/``mid_reverse`` (the tag columns).
        seq_run_id: explicit run id for the whole dataset (overrides the lab/derived value).
        target_gene: explicit marker (overrides the project CSV).
        dataset: dataset label used to derive a default ``seq_run_id`` when none is
            available (modern pre-demultiplexed data); defaults to the field filename.

    Returns:
        A validated :class:`SampleManifest`.

    Raises:
        FileNotFoundError / ValueError: missing file, missing sample key, unresolvable
            date order, or any row failing the canonical model (errors aggregated).
    """
    field_csv = Path(field_csv)
    if not field_csv.exists():
        raise FileNotFoundError(f"Field metadata not found: {field_csv}")
    dataset = dataset or _dataset_from_filename(field_csv)

    df = pd.read_csv(field_csv, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    if df.empty:
        raise ValueError(f"Field metadata is empty: {field_csv}")

    field_to_raw, dropped_known, dropped_unexpected, collisions = _build_header_map(list(df.columns))
    if "eventID" not in field_to_raw:
        raise ValueError(
            f"Field metadata {field_csv} has no eventID/samp_name column "
            f"(after BOM/casing normalisation of {list(df.columns)[:5]}...)."
        )
    if dropped_known:
        logger.warning(
            f"[WARN] manifest_migrate: {len(dropped_known)} field column(s) carry no FAIRe slot "
            f"and were not mapped into the manifest (kept in the source CSV): {dropped_known}"
        )
    for dropped_raw, kept_raw, canonical in collisions:
        logger.warning(
            f"[WARN] manifest_migrate: expected=one source column per manifest field, got=both "
            f"{kept_raw!r} and {dropped_raw!r} map to {canonical!r}, fallback=using {kept_raw!r} "
            f"and dropping {dropped_raw!r} (resolve the duplicate in {field_csv.name})"
        )
    for col in dropped_unexpected:
        logger.warning(
            f"[WARN] manifest_migrate: expected=a recognised field-metadata column, "
            f"got=unexpected column {col!r} in {field_csv.name}, fallback=dropped "
            f"(verify this is not a shifted/garbled column)"
        )

    # marker -> target_gene / assay_name
    marker = target_gene
    if marker is None and project_csv is not None:
        marker = _read_marker(Path(project_csv))

    # seq_run_id / tag source
    lab_index: Dict[str, Dict[str, Optional[str]]] = {}
    if lab_csv is not None:
        lab_index = _load_lab_metadata(Path(lab_csv))

    default_run_id = seq_run_id or (f"{dataset}_{marker}" if marker else dataset)
    used_default_run = False
    # Whether a grouping source exists at all: a library/seq_run_id column in the field CSV,
    # or a demux lab CSV. If neither, synthesising one run id for the whole dataset is the
    # legitimate pre-demultiplexed case; if one exists, a row falling back lost ITS grouping.
    has_run_col = "seq_run_id" in field_to_raw
    have_lab = bool(lab_index)
    rows_defaulted_grouping = 0

    # eventDate normalisation (per file)
    ev_raw_col = field_to_raw.get("eventDate")
    date_map: Dict[str, str] = {}
    if ev_raw_col is not None:
        date_map = normalise_event_dates(
            [str(v) for v in df[ev_raw_col].tolist()],
            context=f"{field_csv.name}:eventDate",
            order=date_order,
        )

    rows: List[SampleManifestRow] = []
    errors: List[str] = []
    extraction_of_samples: set = set()
    extraction_of_blanks: set = set()

    for i, rec in enumerate(df.to_dict(orient="records")):
        # Null missing-value tokens here too (not just inside the strict model), so the
        # migrator's own bookkeeping never treats "NA" as a real value -- e.g. an extraction
        # blank with extraction_ID="NA" must read as None, not register a fabricated batch.
        raw = {field: _cell(rec[col]) for field, col in field_to_raw.items()}
        event_id = raw.get("eventID")
        if not event_id:
            logger.warning(
                f"[WARN] manifest_migrate: expected=an eventID, got=missing on row "
                f"{i + 2} of {field_csv.name}, fallback=row skipped"
            )
            continue

        cls = classify_control(event_id)
        if cls.warn_reason:
            logger.warning(f"[WARN] manifest_migrate: control classification: {cls.warn_reason}")

        kwargs: Dict[str, Any] = dict(raw)
        # eventDate -> ISO
        ev_date = raw.get("eventDate")
        if ev_date:
            kwargs["eventDate"] = date_map.get(ev_date, ev_date)
        # numeric fields: strip units, range-check, null-with-WARN on failure (never crash)
        for nfield in _NUMERIC_RANGES:
            nval = raw.get(nfield)
            if nval:
                kwargs[nfield] = _clean_numeric(nfield, nval, event_id)

        # control identity (two-field FAIRe model)
        kwargs["samp_category"] = cls.samp_category
        kwargs["neg_cont_type"] = cls.neg_cont_type
        kwargs["pos_cont_type"] = cls.pos_cont_type

        # marker / assay
        if marker:
            kwargs.setdefault("target_gene", marker)
            kwargs["assay_name"] = marker
        kwargs["assay_type"] = "metabarcoding"

        # seq_run_id + tags precedence: a library/tag column in the field CSV itself wins,
        # then a separate demux lab CSV joined on eventID, then the synthesised dataset run id.
        lab_rec = lab_index.get(event_id)
        if raw.get("seq_run_id"):
            pass  # already mapped from a library/seq_run_id column in the field CSV
        elif lab_rec and lab_rec.get("seq_run_id"):
            kwargs["seq_run_id"] = lab_rec["seq_run_id"]
            if lab_rec.get("mid_forward") and not raw.get("mid_forward"):
                kwargs["mid_forward"] = lab_rec["mid_forward"]
            if lab_rec.get("mid_reverse") and not raw.get("mid_reverse"):
                kwargs["mid_reverse"] = lab_rec["mid_reverse"]
        else:
            kwargs["seq_run_id"] = default_run_id
            used_default_run = True
            if has_run_col or have_lab:
                # A grouping source exists but THIS row has no value in it -> it is silently
                # batched with the default unless we say so. Name the offending eventID.
                rows_defaulted_grouping += 1
                logger.warning(
                    f"[WARN] manifest_migrate: expected=a library/seq_run_id grouping for "
                    f"{event_id!r}, got=none (empty grouping cell or eventID absent from the "
                    f"lab CSV), fallback=assigned the dataset default seq_run_id="
                    f"{default_run_id!r} -- this row joins a different DADA2 batch; verify"
                )

        # extraction_ID expectations (the control-association key)
        extraction = raw.get("extraction_ID")
        if cls.rule == "blank-ext":
            if extraction:
                extraction_of_blanks.add(extraction)
            else:
                logger.warning(
                    f"[WARN] manifest_migrate: expected=an extraction_ID for extraction blank "
                    f"{event_id!r}, got=none, fallback=blank cannot be tied to a batch"
                )
        elif cls.is_control and cls.is_pcr_blank:
            pass  # PCR blanks legitimately have no extraction_ID (whole-dataset scope)
        elif not cls.is_control:
            if extraction:
                extraction_of_samples.add(extraction)
            # a biological sample with no extraction_ID is unusual but not fatal -> note it
            elif "extraction_ID" in field_to_raw:
                logger.info(
                    f"manifest_migrate: sample {event_id!r} has no extraction_ID "
                    f"(control association by extraction batch will be unavailable for it)"
                )

        try:
            rows.append(SampleManifestRow(**kwargs))
        except Exception as exc:
            errors.append(f"  - row {i + 2} ({event_id}): {str(exc).splitlines()[-1] if str(exc) else exc}")

    if errors:
        raise ValueError(
            f"Migration of {field_csv} produced {len(errors)} invalid row(s):\n"
            + "\n".join(errors)
        )

    if used_default_run and not has_run_col and not have_lab:
        # The legitimate pre-demultiplexed case: no grouping source anywhere.
        logger.warning(
            f"[WARN] manifest_migrate: expected=an explicit sequencing-run/library grouping, "
            f"got=none (no 'library'/'seq_run_id' column and no demux lab CSV), "
            f"fallback=synthesised seq_run_id={default_run_id!r} for the whole dataset "
            f"(one library per dataset). Set --seq-run-id or provide a demux lab CSV to override."
        )
    elif rows_defaulted_grouping:
        # A grouping source existed but some rows had no value -> already warned per-row above.
        logger.warning(
            f"[WARN] manifest_migrate: {rows_defaulted_grouping} row(s) had no library/run "
            f"grouping value and were assigned the dataset default seq_run_id={default_run_id!r} "
            f"(see the per-row warnings above)."
        )

    # extraction-batch orphan check: a Blank-ext batch with no biological sample, or a
    # sample whose extraction batch has no extraction blank.
    orphan_blank_batches = extraction_of_blanks - extraction_of_samples
    if orphan_blank_batches:
        logger.warning(
            f"[WARN] manifest_migrate: expected=each extraction blank's extraction_ID to match "
            f"a biological sample, got=orphan blank batch(es) {sorted(orphan_blank_batches)}, "
            f"fallback=kept (these blanks will clean no sample)"
        )
    orphan_sample_batches = extraction_of_samples - extraction_of_blanks
    if orphan_sample_batches:
        logger.info(
            f"manifest_migrate: extraction batch(es) {sorted(orphan_sample_batches)} have "
            f"samples but no extraction blank (cleaning for them falls back to PCR/whole-dataset)"
        )

    # One row per sample-library. Demux tag-map CSVs repeat an eventID across tag rows, so
    # surface duplicate keys loudly (load_manifest rejects them) rather than emit a manifest
    # that silently violates the invariant.
    dup_keys = sorted(k for k, n in Counter((r.eventID, r.seq_run_id) for r in rows).items() if n > 1)
    if dup_keys:
        logger.warning(
            f"[WARN] manifest_migrate: expected=one row per (eventID, seq_run_id), "
            f"got={len(dup_keys)} duplicated key(s) (e.g. {dup_keys[:3]}), fallback=kept all rows "
            f"(a demux tag-map source must be collapsed to one row per sample-library)"
        )

    manifest = SampleManifest(rows=rows, source=field_csv)
    manifest.check_controls()
    manifest.check_completeness()
    return manifest
