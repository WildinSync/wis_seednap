"""Build DarwinCore-compliant GBIF occurrence CSVs from eDNA pipeline outputs."""

import logging
import re
from importlib import resources
from pathlib import Path
from typing import Union

import pandas as pd

from seednap.steps.formatting.non_target_filter import NonTargetFilter
from seednap.steps.formatting.taxonomy_enricher import TaxonomyEnricher

logger = logging.getLogger(__name__)

# ENVO environment medium mapping
_ENVO_TERMS = {
    "water": "liquid water [ENVO_00002006]",
    "soil": "soil [ENVO:00001998]",
    "river": "liquid water [ENVO_00002006]",
}


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

        # Marker info
        marker = project_meta["marker"].iloc[0]
        info_marker = primers[primers["name"] == marker]
        if info_marker.empty:
            logger.warning(f"Marker '{marker}' not found in primers list")

        # Filter non-target taxa
        results = NonTargetFilter().filter(results, marker)

        # Sum reads per eventID and create taxon_seqindex
        results = self._sum_reads(results)
        results = self._create_taxon_seqindex(results)

        # Map env_medium to ENVO terms
        sample_meta["env_medium_envo"] = sample_meta["env_medium"].map(
            self._map_env_medium
        )

        # Join sample metadata
        merged = results.merge(sample_meta, on="eventID", how="left")

        # Add project-level fields
        merged["recordedby"] = project_meta["recordedby"].iloc[0]

        # Create occurrenceID
        merged["occurrenceID"] = self._create_occurrence_id(merged)

        # Protocol (optional column)
        protocol = ""
        if "protocol" in project_meta.columns:
            protocol = project_meta["protocol"].iloc[0]

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
        out["DNA_sequence"] = merged.get("sequence", "").str.upper()

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

        group_cols = [
            c
            for c in [
                "kingdom", "phylum", "class", "order", "family",
                "genus", "species", "taxon", "rank", "sequence", "sampleID",
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
        """Remove control samples (blank, CNEG, CMET, CEXT)."""
        mask = df["eventID"].str.contains(
            r"blank|CNEG|CMET|CEXT", case=False, na=False
        )
        n_removed = mask.sum()
        if n_removed > 0:
            logger.info(f"Removed {n_removed} control sample row(s)")
        return df[~mask].reset_index(drop=True)

    @staticmethod
    def _sum_reads(df: pd.DataFrame) -> pd.DataFrame:
        """Add nb_reads_total column (total reads per eventID)."""
        totals = df.groupby("eventID")["nb_reads"].transform("sum")
        df = df.copy()
        df["nb_reads_total"] = totals
        return df

    @staticmethod
    def _create_taxon_seqindex(df: pd.DataFrame) -> pd.DataFrame:
        """Add taxon_seqindex column — unique index per sequence within each taxon."""
        df = df.copy()
        df["taxon_seqindex"] = df.groupby("taxon")["sequence"].transform(
            lambda s: [f"{s.name}-{i + 1}" for i in pd.factorize(s)[0]]
        )
        return df

    @staticmethod
    def _create_occurrence_id(df: pd.DataFrame) -> pd.Series:
        """Build occurrenceID: eventID:class;order;family;genus;species;taxon_seqindex."""
        return (
            df["eventID"].astype(str)
            + ":"
            + df.get("class", pd.Series("", index=df.index)).fillna("").astype(str)
            + ";"
            + df.get("order", pd.Series("", index=df.index)).fillna("").astype(str)
            + ";"
            + df.get("family", pd.Series("", index=df.index)).fillna("").astype(str)
            + ";"
            + df.get("genus", pd.Series("", index=df.index)).fillna("").astype(str)
            + ";"
            + df.get("species", pd.Series("", index=df.index)).fillna("").astype(str)
            + ";"
            + df.get("taxon_seqindex", pd.Series("", index=df.index))
            .fillna("")
            .astype(str)
        )

    @staticmethod
    def _map_env_medium(term: str) -> str:
        """Convert environment medium term to ENVO standard."""
        if pd.isna(term):
            return ""
        result = _ENVO_TERMS.get(str(term).lower())
        if result is None:
            logger.warning(f"Unknown env_medium '{term}' — defaulting to water")
            return _ENVO_TERMS["water"]
        return result
