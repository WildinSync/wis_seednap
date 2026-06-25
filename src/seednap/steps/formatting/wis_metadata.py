"""Bridge the WIS database into the GBIF/DarwinCore export's metadata CSVs.

The DarwinCore export (``create-gbif`` / the ``darwincore`` pipeline step) joins the
long-format taxonomy table to two metadata tables: per-sample field metadata
(``report.sample_metadata``) and a single project row (``report.project_metadata``).
Historically those came from hand-authored CSVs. This module lets the per-sample field
metadata come instead from the **WIS database** (the normalized PostgreSQL/PostGIS schema
built by ``wis_database_creator``), so an export can be driven from the database of record
rather than a spreadsheet.

It deliberately produces the **same CSVs the builder already consumes** -- it does not change
``DarwinCoreBuilder`` at all. The CLI command ``seednap wis-metadata`` runs it; you then point
``report.sample_metadata`` / ``report.project_metadata`` at the generated files.

Design:

- The transform layer (``build_sample_metadata_df`` / ``build_project_metadata_df`` and the
  small mappers) is pure Python over plain dict rows, so it is fully unit-testable without a
  database and is where the env_medium mapping, date formatting and column assembly live.
- The query layer (``WisMetadataExporter``) is a thin PostgreSQL/PostGIS reader. SQLAlchemy and a
  PostgreSQL driver are an **optional** dependency (``pip install 'seednap[wis]'``); the import is
  lazy and, when missing, raises a clear actionable error rather than failing obscurely.

env_medium: the WIS schema stores the environmental medium as a controlled ``sample_type`` code
(``models/sample_type_catalog.py``), not as an ENVO ontology term. This bridge maps those codes
to the small env_medium vocabulary ``DarwinCoreBuilder`` recognises (which it in turn maps to
ENVO). A WIS sample_type with no aquatic / soil / sediment analogue (air, blood, honey, deadwood)
is intentionally left unmapped: the bridge passes the raw label through with a ``[WARN]`` so the
builder rejects it loudly rather than silently publishing a mislabelled record
(the no-silent-fallbacks policy).
"""

import logging
from datetime import date, datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# WIS sample_type 2-letter code -> the env_medium vocabulary DarwinCoreBuilder recognises
# (water / soil / river / marine / sediment, mapped there to ENVO). Codes are from
# wis_database_creator models/sample_type_catalog.py (SAMPLE_TYPES). Codes with no aquatic /
# soil / sediment analogue (AI air, BL blood, HN honey, DW deadwood) are intentionally absent:
# _env_medium_from_sample_type passes their raw label through with a [WARN] so the builder fails
# loudly rather than mislabelling a GBIF record.
WIS_SAMPLE_TYPE_TO_ENV_MEDIUM: Dict[str, str] = {
    "FW": "water",  # Freshwater
    "MA": "marine",  # Marine water
    "SE": "sediment",  # Sediment
    "SO": "soil",  # Soil
    "SU": "water",  # Surface (water): assumed aquatic, the common eDNA case
}

# Raw English labels for the WIS sample_type codes (same catalog), used only to make the
# passthrough [WARN] human-readable for an unmapped medium.
_WIS_SAMPLE_TYPE_LABELS: Dict[str, str] = {
    "SU": "Surface",
    "SO": "Soil",
    "FW": "Freshwater",
    "MA": "Marine water",
    "AI": "Air",
    "BL": "Blood",
    "HN": "Honey",
    "DW": "Deadwood",
    "SE": "Sediment",
}

# eventID source column: the pipeline's per-sample names (the taxonomy table's eventIDs, derived
# from the FASTQ filenames) must correspond to one of these WIS sample identifiers. ``sample_id``
# is the short operator-entered code (e.g. DAR2025005); ``material_sample_id`` is the computed
# canonical key. Selectable on the CLI so it can match a site's FASTQ naming without a code change.
EVENT_ID_FIELDS: Tuple[str, ...] = ("sample_id", "material_sample_id")


def _env_medium_from_sample_type(code: object) -> str:
    """Map a WIS ``sample_type`` code to a DarwinCoreBuilder env_medium label.

    Args:
        code: The ``sample_metadata.sample_type_code`` value (2-letter code; NaN/None tolerated).

    Returns:
        One of the builder's recognised env_medium labels (``water``/``marine``/``sediment``/
        ``soil``) for a mapped code; ``""`` for a missing code; or the raw English label for an
        unmapped code (after emitting a ``[WARN]``), so the builder's env_medium validation
        surfaces it rather than the bridge silently mislabelling the sample.
    """
    if code is None or (isinstance(code, float) and pd.isna(code)):
        return ""
    key = str(code).strip().upper()
    mapped = WIS_SAMPLE_TYPE_TO_ENV_MEDIUM.get(key)
    if mapped is not None:
        return mapped
    label = _WIS_SAMPLE_TYPE_LABELS.get(key, key)
    logger.warning(
        f"[WARN] wis-metadata: expected=a WIS sample_type code with a known env_medium "
        f"mapping, got={key!r} ({label!r}), fallback=writing the raw label so the DarwinCore "
        f"builder fails loudly rather than mislabelling the record. Add an entry to "
        f"WIS_SAMPLE_TYPE_TO_ENV_MEDIUM (and a term to the builder's _ENVO_TERMS) if this "
        f"medium should be published to GBIF.",
    )
    return label


def _format_event_date(value: object) -> str:
    """Format a WIS ``event_date`` as the builder's ``yyyy.mm.dd`` eventDate.

    Args:
        value: A ``date``/``datetime``, an ISO date string, or NaN/None.

    Returns:
        ``"YYYY.MM.DD"`` for a parseable date, or ``""`` when the value is missing or
        unparseable (an unparseable value also logs a ``[WARN]`` so it is not silently blanked).
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, (date, datetime)):
        return f"{value.year:04d}.{value.month:02d}.{value.day:02d}"
    ts = pd.to_datetime(value, errors="coerce")
    if pd.isna(ts):
        logger.warning(
            f"[WARN] wis-metadata: expected=a parseable event_date, got={value!r}, "
            f"fallback=blank eventDate for this sample.",
        )
        return ""
    return str(ts.strftime("%Y.%m.%d"))


def build_sample_metadata_df(
    rows: Sequence[Dict[str, Any]], event_id_field: str = "sample_id"
) -> pd.DataFrame:
    """Assemble the per-sample metadata CSV (builder contract) from WIS sample rows.

    Args:
        rows: Sample rows as plain dicts, each with the keys the query layer selects:
            ``sample_id``, ``material_sample_id``, ``event_date``, ``sample_depth``,
            ``sample_size_frac``, ``sample_size``, ``sample_type_code``, ``lat``, ``lon``.
        event_id_field: Which WIS identifier becomes the ``eventID`` (must be one of
            :data:`EVENT_ID_FIELDS`); defaults to ``sample_id``.

    Returns:
        A DataFrame with exactly the columns ``DarwinCoreBuilder`` reads from the sample
        metadata: ``eventID``, ``eventDate``, ``env_medium`` (builder vocabulary),
        ``decimalLatitude``, ``decimalLongitude``, ``depth``, ``size_frac``, ``samp_size``.

    Raises:
        ValueError: If ``event_id_field`` is not one of :data:`EVENT_ID_FIELDS`.
    """
    if event_id_field not in EVENT_ID_FIELDS:
        raise ValueError(
            f"event_id_field must be one of {EVENT_ID_FIELDS}, got {event_id_field!r}."
        )
    columns = [
        "eventID",
        "eventDate",
        "env_medium",
        "decimalLatitude",
        "decimalLongitude",
        "depth",
        "size_frac",
        "samp_size",
    ]
    records: List[Dict[str, Any]] = []
    for r in rows:
        records.append(
            {
                "eventID": "" if r.get(event_id_field) is None else str(r.get(event_id_field)),
                "eventDate": _format_event_date(r.get("event_date")),
                "env_medium": _env_medium_from_sample_type(r.get("sample_type_code")),
                "decimalLatitude": r.get("lat"),
                "decimalLongitude": r.get("lon"),
                "depth": r.get("sample_depth"),
                "size_frac": r.get("sample_size_frac"),
                "samp_size": r.get("sample_size"),
            }
        )
    return pd.DataFrame.from_records(records, columns=columns)


def build_project_metadata_df(
    *,
    marker: str,
    recorded_by: str,
    identification_remarks: str,
    identification_references: str,
    seq_meth: str = "",
    otu_seq_comp_appr: str = "",
) -> pd.DataFrame:
    """Assemble the single-row project metadata CSV (builder contract).

    The reference-database (``otu_db``) and chimera-removal (``chimera_check``) provenance are
    intentionally omitted: the ``darwincore`` pipeline step fills those from the run config (the
    single source of truth). The remaining required project fields are operator-supplied because
    the WIS schema does not store an ``identificationRemarks`` / ``identificationReferences``
    equivalent.

    Args:
        marker: Marker name, matching an entry in the bundled ``primers_list.csv`` (e.g. ``teleo``).
        recorded_by: Data contributor / recorder (DwC ``recordedBy``).
        identification_remarks: Free-text note on the identification method.
        identification_references: Citation(s) for the reference DB / identification method.
        seq_meth: Optional sequencing-method description (DwC ``seq_meth``).
        otu_seq_comp_appr: Optional OTU/ASV sequence-comparison approach.

    Returns:
        A one-row DataFrame with the builder's project-metadata columns.
    """
    return pd.DataFrame(
        [
            {
                "marker": marker,
                "recordedby": recorded_by,
                "identificationRemarks": identification_remarks,
                "identificationReferences": identification_references,
                "seqmet": seq_meth,
                "otu_seq_comp_appr": otu_seq_comp_appr,
            }
        ]
    )


# SQL for the per-sample field metadata. PostGIS: coordinates live only in gis_point.geom
# (SRID 4326), so decimal lat/lon come from ST_Y/ST_X of the COORDINATE-typed point joined on
# material_sample_id. Column/table names are from wis_database_creator's SQLAlchemy models.
_SAMPLE_QUERY = """
    SELECT s.sample_id          AS sample_id,
           s.material_sample_id AS material_sample_id,
           s.event_date         AS event_date,
           s.sample_depth       AS sample_depth,
           s.sample_size_frac   AS sample_size_frac,
           s.sample_size        AS sample_size,
           s.sample_type_code   AS sample_type_code,
           ST_Y(p.geom)         AS lat,
           ST_X(p.geom)         AS lon
    FROM sample_metadata s
    LEFT JOIN gis_point p
           ON p.material_sample_id = s.material_sample_id
          AND p.gis_point_location_type = 'COORDINATE'
    {where}
    ORDER BY s.sample_id
"""


class WisMetadataExporter:
    """Read per-sample field metadata from the WIS PostgreSQL/PostGIS database.

    Thin query layer over the WIS schema. SQLAlchemy and a PostgreSQL driver are an optional
    dependency (``pip install 'seednap[wis]'``); the import is lazy so the rest of the pipeline
    never requires them, and a missing install raises a clear, actionable error.
    """

    def __init__(self, database_url: str) -> None:
        """Store the connection string.

        Args:
            database_url: SQLAlchemy URL for the WIS database, e.g.
                ``postgresql://user:pass@host:5432/wis`` (matches wis_database_creator's
                ``config/settings.py`` connection string).
        """
        self.database_url = database_url

    @staticmethod
    def _require_sqlalchemy() -> Any:
        """Import SQLAlchemy, or raise a clear install hint if the optional extra is absent.

        Returns:
            The imported ``sqlalchemy`` module.

        Raises:
            RuntimeError: If SQLAlchemy is not installed, with the ``pip install 'seednap[wis]'``
                remedy (the no-silent-fallbacks policy: fail loudly with the fix, not obscurely).
        """
        try:
            import sqlalchemy
        except ImportError as exc:
            raise RuntimeError(
                "The WIS metadata bridge needs SQLAlchemy and a PostgreSQL driver, which are an "
                "optional dependency. Install them with:  pip install 'seednap[wis]'  (this adds "
                "sqlalchemy and psycopg2). They are intentionally not part of the core pipeline "
                f"so a run that does not use the database stays dependency-light. (import error: {exc})"
            ) from exc
        return sqlalchemy

    def fetch_sample_rows(
        self, monitoring: Optional[str] = None, mission: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """Query the per-sample field metadata, optionally filtered to a site / campaign.

        Args:
            monitoring: Restrict to one ``monitoring_id`` (the WIS site / long-term project),
                or None for no site filter.
            mission: Restrict to one ``mission_id`` (the WIS sampling campaign), or None.

        Returns:
            One plain dict per sample with the keys the transform layer expects (see
            :func:`build_sample_metadata_df`). ``lat``/``lon`` are None for a sample with no
            COORDINATE point.
        """
        sqlalchemy = self._require_sqlalchemy()
        clauses, params = [], {}
        if monitoring:
            clauses.append("s.monitoring_id = :monitoring")
            params["monitoring"] = monitoring
        if mission:
            clauses.append("s.mission_id = :mission")
            params["mission"] = mission
        where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
        sql = sqlalchemy.text(_SAMPLE_QUERY.format(where=where))

        engine = sqlalchemy.create_engine(self.database_url, pool_pre_ping=True)
        try:
            with engine.connect() as conn:
                result = conn.execute(sql, params)
                return [dict(row._mapping) for row in result]
        finally:
            engine.dispose()

    def export(
        self,
        *,
        output_dir: Path,
        marker: str,
        recorded_by: str,
        identification_remarks: str,
        identification_references: str,
        monitoring: Optional[str] = None,
        mission: Optional[str] = None,
        event_id_field: str = "sample_id",
        seq_meth: str = "",
        otu_seq_comp_appr: str = "",
    ) -> Tuple[Path, Path]:
        """Write ``<marker>_sample_metadata.csv`` and ``<marker>_project_metadata.csv``.

        Args:
            output_dir: Directory to write the two CSVs into (created if absent).
            marker: Marker name for the project row and the output filenames.
            recorded_by: DwC ``recordedBy`` for the project row.
            identification_remarks: Identification-method note for the project row.
            identification_references: Reference-DB / method citation for the project row.
            monitoring: Optional ``monitoring_id`` site filter.
            mission: Optional ``mission_id`` campaign filter.
            event_id_field: Which WIS identifier becomes ``eventID`` (see :data:`EVENT_ID_FIELDS`).
            seq_meth: Optional sequencing-method description.
            otu_seq_comp_appr: Optional OTU/ASV comparison-approach description.

        Returns:
            ``(sample_csv_path, project_csv_path)``.

        Raises:
            ValueError: If the query returns no samples for the given filter (an empty export
                would silently produce a GBIF file with no field metadata).
            RuntimeError: If the optional SQLAlchemy dependency is not installed.
        """
        rows = self.fetch_sample_rows(monitoring=monitoring, mission=mission)
        if not rows:
            scope = monitoring or mission or "the whole database"
            raise ValueError(
                f"No samples found in the WIS database for {scope!r}. Check the --monitoring / "
                f"--mission selector and the database URL; refusing to write empty metadata that "
                f"would leave every GBIF occurrence without a date/location."
            )

        sample_df = build_sample_metadata_df(rows, event_id_field=event_id_field)
        project_df = build_project_metadata_df(
            marker=marker,
            recorded_by=recorded_by,
            identification_remarks=identification_remarks,
            identification_references=identification_references,
            seq_meth=seq_meth,
            otu_seq_comp_appr=otu_seq_comp_appr,
        )

        output_dir.mkdir(parents=True, exist_ok=True)
        sample_csv = output_dir / f"{marker}_sample_metadata.csv"
        project_csv = output_dir / f"{marker}_project_metadata.csv"
        sample_df.to_csv(sample_csv, index=False)
        project_df.to_csv(project_csv, index=False)
        logger.info(
            f"wis-metadata: wrote {len(sample_df)} sample row(s) to {sample_csv} and the "
            f"project row to {project_csv}"
        )
        return sample_csv, project_csv
