"""Build DarwinCore-compliant GBIF occurrence CSVs from eDNA pipeline outputs."""

import hashlib
import logging
import re
from importlib import resources
from pathlib import Path
from typing import Union

import pandas as pd

from seednap.steps.formatting.non_target_filter import NonTargetFilter
from seednap.steps.formatting.taxonomy_enricher import TaxonomyEnricher

logger = logging.getLogger(__name__)

# ENVO environment medium mapping. Add new terms here -- unrecognised values
# now raise (G3) instead of silently defaulting to water.
_ENVO_TERMS = {
    "water": "liquid water [ENVO_00002006]",
    "soil": "soil [ENVO:00001998]",
    "river": "liquid water [ENVO_00002006]",
    "marine": "sea water [ENVO_00002149]",
    "sediment": "sediment [ENVO_00002007]",
}

# Required GBIF DwC eDNA-Occurrence fields. Empty values fail G1 validation.
_DWC_REQUIRED_FIELDS = (
    "occurrenceID",
    "eventID",
    "basisOfRecord",
    "target_gene",
    "pcr_primer_forward",
    "pcr_primer_reverse",
    "otu_db",
)


def _load_template(filename: str) -> pd.DataFrame:
    """Load a bundled CSV template from the seednap.data.templates package."""
    ref = resources.files("seednap.data.templates").joinpath(filename)
    with resources.as_file(ref) as path:
        return pd.read_csv(path)


class DarwinCoreBuilder:
    """
    Transform taxonomy results + metadata into a DarwinCore-compliant CSV.

    This replaces the R ``create_df_gbif`` project. It takes three input CSVs
    (taxonomy results, sample metadata, project metadata) and produces a single
    output CSV with all 38+ DarwinCore columns populated.
    """

    def __init__(
        self,
        taxonomy_results_path: Union[str, Path],
        sample_metadata_path: Union[str, Path],
        project_metadata_path: Union[str, Path],
        output_path: Union[str, Path],
        summarise_pcr_replicates: bool = False,
        skip_enrichment: bool = False,
    ) -> None:
        """Store input/output paths and build options.

        Args:
            taxonomy_results_path: CSV of taxonomy results (long format).
            sample_metadata_path: CSV of per-sample metadata.
            project_metadata_path: CSV of project-level metadata.
            output_path: Destination path for the DarwinCore CSV.
            summarise_pcr_replicates: If True, collapse PCR replicate suffixes
                and sum their reads before building the output.
            skip_enrichment: If True, skip NCBI/WORMS kingdom/phylum enrichment.
        """
        self.taxonomy_results_path = Path(taxonomy_results_path)
        self.sample_metadata_path = Path(sample_metadata_path)
        self.project_metadata_path = Path(project_metadata_path)
        self.output_path = Path(output_path)
        self.summarise_pcr_replicates = summarise_pcr_replicates
        self.skip_enrichment = skip_enrichment

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def build(self) -> Path:
        """
        Run the full DarwinCore build pipeline.

        Returns:
            Path to the written output CSV.

        Raises:
            FileNotFoundError: If any input file does not exist.
            ValueError: If required columns are missing or dates are invalid.
        """
        # Validate inputs exist
        for p in (
            self.taxonomy_results_path,
            self.sample_metadata_path,
            self.project_metadata_path,
        ):
            if not p.exists():
                raise FileNotFoundError(f"Input file not found: {p}")

        # Load CSVs
        results = pd.read_csv(self.taxonomy_results_path)
        sample_meta = pd.read_csv(self.sample_metadata_path)
        project_meta = pd.read_csv(self.project_metadata_path)

        # Legacy column renames
        if "filter_code" in results.columns:
            results = results.rename(columns={"filter_code": "eventID"})
        if "volume" in sample_meta.columns:
            sample_meta = sample_meta.rename(columns={"volume": "samp_size"})

        # G3: validate input metadata BEFORE doing any work
        self._validate_sample_metadata(sample_meta, self.sample_metadata_path)
        self._validate_project_metadata(project_meta, self.project_metadata_path)

        # The taxonomy results table (format-gbif output) must carry the long-format columns we
        # read directly below; validate up-front so a missing one is a clear error rather than a
        # raw pandas KeyError mid-build.
        results_required = ("eventID", "taxon", "nb_reads")
        missing_results = [c for c in results_required if c not in results.columns]
        if missing_results:
            raise ValueError(
                f"Taxonomy results table (TAXONOMY_RESULTS) is missing required column(s) "
                f"{missing_results}. create-gbif expects the long-format output of "
                f"`seednap format-gbif` (eventID, taxon, nb_reads, plus class/order/family/"
                f"genus/species). Run format-gbif first and pass its output here."
            )

        # PCR replicate summarisation
        if self.summarise_pcr_replicates:
            results = self._summarise_pcr_replicates(results)
            logger.info("PCR replicates summarised per sample")

        # Remove controls
        results = self._remove_controls(results)

        # Validate dates
        self._validate_dates(sample_meta["eventDate"])

        # Load bundled templates
        primers = _load_template("primers_list.csv")

        # Marker info -- G1: hard-fail if missing rather than silently writing
        # blank target_gene / primer columns into the GBIF submission. Match
        # case-insensitively so "Teleo" matches "teleo" in the primers list.
        marker = str(project_meta["marker"].iloc[0])
        info_marker = primers[primers["name"].str.lower() == marker.lower()]
        if info_marker.empty:
            available = sorted(primers["name"].unique().tolist())
            raise ValueError(
                f"Marker '{marker}' not found in primers_list.csv. "
                f"GBIF submission requires target_gene and primer info; "
                f"please add an entry for this marker. "
                f"Available markers: {available}"
            )

        # G2: surface contaminant flags before any filtering so they survive to output.
        if "is_contaminant_candidate" in results.columns:
            n_contam = int(results["is_contaminant_candidate"].sum())
            if n_contam > 0:
                logger.info(
                    f"{n_contam} OTU rows are flagged as contaminant candidates "
                    f"(propagated to GBIF output as `contamination_flag`)"
                )

        # Filter non-target taxa
        results = NonTargetFilter().filter(results, marker)

        # Sum reads per eventID
        results = self._sum_reads(results)

        # Map env_medium to ENVO terms
        sample_meta["env_medium_envo"] = sample_meta["env_medium"].map(
            self._map_env_medium
        )

        # Join sample metadata
        merged = results.merge(sample_meta, on="eventID", how="left")

        # Add project-level fields
        merged["recordedby"] = project_meta["recordedby"].iloc[0]

        # G5: stable occurrenceID = marker:eventID:sha256(sequence)[:8]
        # Independent of run order and unique by sequence content. The previous
        # scheme used a per-run taxon_seqindex which broke GBIF dataset versioning.
        merged["occurrenceID"] = self._create_occurrence_id(merged, marker)

        # Build the DarwinCore output DataFrame
        out = pd.DataFrame()

        # Taxonomy columns — scientificName needed for enrichment
        out["scientificName"] = merged["taxon"]
        out["kingdom"] = ""
        out["phylum"] = ""
        out["class"] = merged.get("class", "")
        out["order"] = merged.get("order", "")
        out["family"] = merged.get("family", "")
        out["genus"] = merged.get("genus", "")

        # Enrich kingdom/phylum via NCBI/WORMS
        if not self.skip_enrichment:
            enricher = TaxonomyEnricher()
            out = enricher.enrich(out)

        # Occurrence / event columns
        out["occurrenceID"] = merged["occurrenceID"]
        out["basisOfRecord"] = "MaterialSample"
        out["eventID"] = merged["eventID"]
        out["eventDate"] = merged["eventDate"]
        out["recordedBy"] = merged["recordedby"]
        out["organismQuantity"] = merged["nb_reads"]
        out["organismQuantityType"] = "DNA sequence reads"
        out["sampleSizeValue"] = merged["nb_reads_total"]
        out["sampleSizeUnit"] = "DNA sequence reads"
        out["identificationRemarks"] = project_meta["identificationRemarks"].iloc[0]
        out["identificationReferences"] = project_meta[
            "identificationReferences"
        ].iloc[0]

        # Location columns
        out["decimalLatitude"] = merged.get("decimalLatitude", "")
        out["decimalLongitude"] = merged.get("decimalLongitude", "")
        out["geodeticDatum"] = "EPSG:4326"
        out["maximumDepthInMeters"] = merged.get("depth", "")

        # Sample columns
        out["size_frac"] = merged.get("size_frac", "")
        out["samp_size"] = merged.get("samp_size", "")

        # Taxonomic rank
        out["TaxonRank"] = merged.get("rank", "")

        # Sequence
        if "sequence" not in merged.columns:
            raise ValueError(
                f"Taxonomy results CSV '{self.taxonomy_results_path}' has no 'sequence' "
                f"column, which is required to fill the DarwinCore 'DNA_sequence' field. "
                f"Provide the long-format taxonomy table from the 'format-gbif' step, which "
                f"always emits a lowercase 'sequence' column (it renames a capital-S "
                f"'Sequence' if needed); a table that lacks it cannot be turned into a GBIF "
                f"occurrence record."
            )
        out["DNA_sequence"] = merged["sequence"].str.upper()

        # Primer / marker columns
        if not info_marker.empty:
            row = info_marker.iloc[0]
            out["target_gene"] = row.get("target_gene", "")
            out["pcr_primer_forward"] = row.get("pcr_primer_forward", "")
            out["pcr_primer_reverse"] = row.get("pcr_primer_reverse", "")
            out["pcr_primer_name_forward"] = row.get("pcr_primer_name_forward", "")
            out["pcr_primer_name_reverse"] = row.get("pcr_primer_name_reverse", "")
            out["pcr_primer_reference"] = row.get("pcr_primer_reference", "")
        else:
            for col in (
                "target_gene",
                "pcr_primer_forward",
                "pcr_primer_reverse",
                "pcr_primer_name_forward",
                "pcr_primer_name_reverse",
                "pcr_primer_reference",
            ):
                out[col] = ""

        # Environment columns
        out["env_medium"] = merged.get("env_medium_envo", "")

        # Sequencing / bioinformatics columns
        out["lib_layout"] = "paired"
        out["seq_meth"] = project_meta.get("seqmet", pd.Series([""])).iloc[0]
        out["otu_seq_comp_appr"] = project_meta.get(
            "otu_seq_comp_appr", pd.Series([""])
        ).iloc[0]
        out["otu_db"] = project_meta.get("otu_db", pd.Series([""])).iloc[0]
        out["chimera_check"] = project_meta.get(
            "chimera_check", pd.Series([""])
        ).iloc[0]

        # G2: propagate contamination flag from upstream taxonomy to GBIF output.
        if "is_contaminant_candidate" in merged.columns:
            out["contamination_flag"] = merged["is_contaminant_candidate"].fillna(False).astype(bool)
        else:
            out["contamination_flag"] = False

        # G1: validate required DwC fields are populated BEFORE writing.
        self._check_required_fields(out)

        # Write output
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(self.output_path, index=False)
        logger.info(f"DarwinCore CSV written to {self.output_path}")

        return self.output_path

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _validate_dates(series: pd.Series) -> None:
        """Raise ValueError if any date doesn't match yyyy[.mm[.dd]]."""
        pattern = re.compile(r"^\d{4}(\.\d{2}){0,2}$")
        invalid = series.dropna().apply(lambda d: not pattern.match(str(d)))
        if invalid.any():
            bad = series[invalid].unique().tolist()
            raise ValueError(
                f"Invalid date format(s): {bad}. "
                "Expected yyyy, yyyy.mm, or yyyy.mm.dd."
            )

    @staticmethod
    def _summarise_pcr_replicates(df: pd.DataFrame) -> pd.DataFrame:
        """Strip PCR replicate suffix (_XX) from eventID and sum reads."""
        df = df.copy()
        df["sampleID"] = df["eventID"].str.replace(
            r"_\d{2}$", "", regex=True
        )

        # is_contaminant_candidate is per-sequence, so keeping it as a group key carries the
        # contaminant flag through the replicate-summing groupby (it would otherwise be dropped,
        # silently zeroing contamination_flag downstream on the --summarise-pcr path).
        group_cols = [
            c
            for c in [
                "kingdom", "phylum", "class", "order", "family",
                "genus", "species", "taxon", "rank", "sequence", "sampleID",
                "is_contaminant_candidate",
            ]
            if c in df.columns
        ]

        df = df.groupby(group_cols, as_index=False, dropna=False).agg(
            nb_reads=("nb_reads", "sum")
        )
        df["eventID"] = df["sampleID"]
        df = df.drop(columns=["sampleID"])
        return df

    @staticmethod
    def _remove_controls(df: pd.DataFrame) -> pd.DataFrame:
        """Drop control rows by eventID, using the canonical control classifier.

        Uses ``config.manifest.classify_control`` -- the FAIRe-anchored single
        source of truth for control identity -- rather than a separate ad-hoc
        regex. The previous local regex (``blank|CNEG|CMET|CEXT``) was a strict
        SUBSET of classify_control's patterns, so controls named ``CPCR*``,
        ``EXT_NC``, ``PCR_NC`` or ``water`` passed through and were written into
        the GBIF occurrence CSV as genuine biological records -- a silent
        injection of negative/positive-control reads into a GBIF submission.

        Any name that classify_control flags as control-LOOKING but cannot
        classify to a known rule (``rule == 'unclassified-control-like'``, e.g.
        an underscore-suffixed ``water_001`` that the canonical patterns do not
        resolve) is left in the output (it stays a biological sample), but a
        ``[WARN]`` is emitted naming the eventID so a possible leaked control is
        visible before GBIF submission rather than silent (no-silent-fallbacks
        policy).
        """
        from seednap.config.manifest import classify_control

        event_ids = df["eventID"].astype(str)
        classes = {eid: classify_control(eid) for eid in event_ids.unique()}

        control_ids = {eid for eid, cls in classes.items() if cls.is_control}
        unclassified = {
            eid: cls.warn_reason
            for eid, cls in classes.items()
            if cls.rule == "unclassified-control-like"
        }

        for eid, reason in unclassified.items():
            print(
                f"[WARN] _remove_controls: expected=control identity for "
                f"eventID {eid!r}, got=control-looking but unclassified "
                f"({reason}), fallback=kept as a biological sample in the GBIF "
                f"output -- verify it is not a control before submission",
                flush=True,
            )

        mask = event_ids.isin(control_ids)
        n_removed = int(mask.sum())
        if n_removed > 0:
            logger.info(
                f"Removed {n_removed} control sample row(s) "
                f"(controls: {sorted(control_ids)})"
            )
        return df[~mask].reset_index(drop=True)

    @staticmethod
    def _sum_reads(df: pd.DataFrame) -> pd.DataFrame:
        """Add nb_reads_total column (total reads per eventID)."""
        totals = df.groupby("eventID")["nb_reads"].transform("sum")
        df = df.copy()
        df["nb_reads_total"] = totals
        return df

    @staticmethod
    def _create_occurrence_id(df: pd.DataFrame, marker: str) -> pd.Series:
        """Build a stable occurrenceID: marker:eventID:sha256(sequence)[:8].

        The hash makes the ID deterministic across re-runs of the same data
        (G5 fix). The previous scheme used a per-run `taxon_seqindex` which
        meant resubmitting the same dataset to GBIF created new occurrences
        instead of replacing the old ones.
        """
        seq = df.get("sequence", pd.Series("", index=df.index))
        seq_hash = seq.fillna("").astype(str).apply(
            lambda s: hashlib.sha256(s.upper().encode("ascii")).hexdigest()[:8]
            if s
            else "NOSEQ"
        )
        return (
            f"{marker}:" + df["eventID"].astype(str) + ":" + seq_hash
        )

    @staticmethod
    def _map_env_medium(term: str) -> str:
        """Convert environment medium term to ENVO standard.

        Unrecognised values raise (G3). Previously they silently defaulted to
        'water', which silently mislabelled terrestrial samples in GBIF.
        """
        if pd.isna(term):
            return ""
        key = str(term).lower()
        result = _ENVO_TERMS.get(key)
        if result is None:
            raise ValueError(
                f"Unknown env_medium '{term}'. Recognised values: "
                f"{sorted(_ENVO_TERMS.keys())}. Add a mapping in "
                f"darwincore_builder._ENVO_TERMS if you need a new one."
            )
        return result

    # ------------------------------------------------------------------
    # G3: input metadata validation
    # ------------------------------------------------------------------

    @staticmethod
    def _validate_sample_metadata(
        sample_meta: pd.DataFrame, sample_metadata_path: Path
    ) -> None:
        """Sanity-check sample metadata before merging.

        Raises ValueError on missing required columns, out-of-range coordinates, or
        unknown env_medium.
        """
        # env_medium is required: build() maps it to an ENVO term by direct indexing,
        # so a missing column must fail here with a clear message, not a raw KeyError later.
        required = ("eventID", "eventDate", "env_medium")
        missing = [c for c in required if c not in sample_meta.columns]
        if missing:
            raise ValueError(
                f"Sample metadata CSV '{sample_metadata_path}' is missing required "
                f"columns: {missing}. GBIF needs 'eventID' (must match the per-sample "
                f"column names carried in the taxonomy table), 'eventDate' (format yyyy, "
                f"yyyy.mm, or yyyy.mm.dd), and 'env_medium' (one of the known ENVO terms: "
                f"water, soil, river, marine, sediment). Rename your headers to these exact "
                f"names. Recognized optional columns: decimalLatitude, decimalLongitude, "
                f"depth, size_frac, samp_size (legacy 'volume' is auto-renamed to samp_size)."
            )

        # Latitude / longitude range checks
        if "decimalLatitude" in sample_meta.columns:
            lat = pd.to_numeric(sample_meta["decimalLatitude"], errors="coerce")
            bad_lat = lat[(lat < -90) | (lat > 90)]
            if len(bad_lat) > 0:
                raise ValueError(
                    f"Invalid decimalLatitude (must be in [-90, 90]): "
                    f"{bad_lat.tolist()}"
                )
        if "decimalLongitude" in sample_meta.columns:
            lon = pd.to_numeric(sample_meta["decimalLongitude"], errors="coerce")
            bad_lon = lon[(lon < -180) | (lon > 180)]
            if len(bad_lon) > 0:
                raise ValueError(
                    f"Invalid decimalLongitude (must be in [-180, 180]): "
                    f"{bad_lon.tolist()}"
                )

        # env_medium values must all map to known ENVO terms (or be NaN)
        if "env_medium" in sample_meta.columns:
            non_null = sample_meta["env_medium"].dropna().astype(str).str.lower().unique()
            unknown = [v for v in non_null if v not in _ENVO_TERMS]
            if unknown:
                raise ValueError(
                    f"Unknown env_medium values in sample metadata: {unknown}. "
                    f"Recognised values: {sorted(_ENVO_TERMS.keys())}."
                )

    @staticmethod
    def _validate_project_metadata(
        project_meta: pd.DataFrame, project_metadata_path: Path
    ) -> None:
        """Sanity-check project metadata before building the GBIF output."""
        required = ("marker", "recordedby", "identificationRemarks", "identificationReferences")
        missing = [c for c in required if c not in project_meta.columns]
        if missing:
            raise ValueError(
                f"Project metadata CSV '{project_metadata_path}' is missing required "
                f"column(s): {missing}. The project-metadata CSV needs all of (exact "
                f"lowercase headers): marker, recordedby, identificationRemarks, "
                f"identificationReferences. Optional columns: seqmet, otu_seq_comp_appr, "
                f"otu_db, chimera_check. It is a single-row table: one header row, one data "
                f"row describing the run. Note headers are case-sensitive (e.g. "
                f"'identificationReferences', not 'IdentificationReferences')."
            )
        if len(project_meta) == 0:
            raise ValueError(
                f"Project metadata CSV '{project_metadata_path}' has the required column "
                f"headers (marker, recordedby, identificationRemarks, "
                f"identificationReferences) but no data rows. GBIF needs exactly one "
                f"project-metadata row. Add a data row describing this project beneath the "
                f"header line."
            )
        empty = [
            c for c in required
            if project_meta[c].iloc[0] in ("", None) or pd.isna(project_meta[c].iloc[0])
        ]
        if empty:
            raise ValueError(
                f"Project metadata fields are empty (required for GBIF): {empty}. Fill in a "
                f"value for each of these columns in the data row of the project metadata CSV "
                f"'{project_metadata_path}' (recordedby = the data contributor / recorder; "
                f"identificationReferences = the reference-DB or method citation), then re-run."
            )

    @staticmethod
    def _check_required_fields(out: pd.DataFrame) -> None:
        """Verify every DwC-required column has data before writing (G1)."""
        # Zero-row check FIRST: on an empty frame `(series == "").all()` is True
        # for every column, so the blank-field loop would otherwise report all
        # required fields as 'blank in every row' and misattribute the cause to a
        # missing otu_db. The real cause is that no occurrence rows survived.
        if len(out) == 0:
            raise ValueError(
                "GBIF output has no occurrence rows: every record was dropped by "
                "control removal and non-target filtering (or all reads were "
                "zero-filtered upstream). There is nothing to submit. Check that "
                "the taxonomy results table is non-empty and that controls/"
                "non-target taxa did not remove every row."
            )
        empty_fields = []
        for col in _DWC_REQUIRED_FIELDS:
            if col not in out.columns:
                empty_fields.append(col)
                continue
            series = out[col].astype(str)
            if (series == "").all() or series.isna().all():
                empty_fields.append(col)
        if empty_fields:
            raise ValueError(
                f"GBIF output has these required DarwinCore fields blank in every row: "
                f"{empty_fields}. GBIF rejects submissions with empty required fields. The "
                f"usual cause is 'otu_db' missing or left blank in your project-metadata CSV "
                f"(the third create-gbif argument) -- it is an optional input column but a "
                f"required output field, and it is not checked at load time. Add an 'otu_db' "
                f"value (e.g. the reference database name, like 'CRABS MitoFish 2025') to that "
                f"CSV. If a primer field is listed (target_gene/pcr_primer_forward/"
                f"pcr_primer_reverse), the marker's row in the bundled primers_list.csv has a "
                f"blank cell for it -- fix that template entry."
            )
