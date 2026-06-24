"""Outputs config: GBIF/DarwinCore export, metrics, run reporting."""

from pathlib import Path
from typing import Optional

from pydantic import Field, field_validator

from seednap.config.models.base import StrictModel


# ===========================================================================
# OUTPUTS: GBIF/DarwinCore export, metrics, run reporting
# ===========================================================================


class GbifExportConfig(StrictModel):
    """GBIF / DarwinCore export options (the 'export' step runs iff listed in pipeline.steps).

    GBIF (the Global Biodiversity Information Facility) is where these occurrence datasets are
    published; export reshapes the taxonomy-annotated table into a DarwinCore-aligned form.

    Attributes:
        add_rank: Add a column naming the lowest taxonomic rank that was resolved per feature.
        add_taxon: Add a column carrying the lowest available taxon name per feature.
    """

    add_rank: bool = Field(default=True, description="Add taxonomic rank column")
    add_taxon: bool = Field(default=True, description="Add lowest available taxon column")


class DarwinCoreExportConfig(StrictModel):
    """DarwinCore occurrence export options (the 'darwincore' step runs iff listed in
    pipeline.steps).

    This is the GBIF-ready occurrence file: it joins the long-format export to per-sample +
    project metadata, fills DarwinCore fields, and (unless skipped) enriches higher ranks from
    NCBI/WoRMS. The metadata it needs comes from ``report.sample_metadata`` (per-sample) and
    ``report.project_metadata`` (per-project); the 'darwincore' step requires both to be set.

    Attributes:
        summarise_pcr_replicates: Collapse PCR-replicate suffixes, summing reads per sample.
        skip_enrichment: Skip the NCBI/WoRMS kingdom/phylum enrichment (offline / faster).
    """

    summarise_pcr_replicates: bool = Field(
        default=False, description="Collapse PCR-replicate suffixes, summing reads per sample"
    )
    skip_enrichment: bool = Field(
        default=False, description="Skip the NCBI/WoRMS kingdom/phylum enrichment"
    )


class ExportConfig(StrictModel):
    """Export step configuration.

    Attributes:
        gbif: GBIF long-format export options (the 'export' step).
        darwincore: DarwinCore occurrence export options (the 'darwincore' step).
    """

    gbif: GbifExportConfig = Field(
        default_factory=GbifExportConfig, description="GBIF long-format export settings"
    )
    darwincore: DarwinCoreExportConfig = Field(
        default_factory=DarwinCoreExportConfig, description="DarwinCore occurrence export settings"
    )


class ReportConfig(StrictModel):
    """Run reporting parameters (the 'report' step runs iff listed in pipeline.steps).

    When the 'report' step runs it always writes the per-step read/sequence tracking table and
    step summary; ``html_report`` additionally toggles the self-contained HTML document.
    """

    html_report: bool = Field(
        default=True,
        description="Generate a self-contained HTML run report with charts (on by default)",
    )
    output_dir: Optional[Path] = Field(
        default=None,
        description="Base directory for report artifacts; a per-marker subdirectory is created "
                    "inside it. Defaults to '<paths.output>/04_report' when unset.",
    )
    warn_below_retention_pct: float = Field(
        default=30.0, ge=0, le=100,
        description="Warn for samples whose final non-chimeric reads fall below this % of raw reads",
    )
    warn_step_loss_pct: float = Field(
        default=70.0, ge=0, le=100,
        description="Warn when a single pipeline step drops more than this % of a sample's reads",
    )
    sample_metadata: Optional[Path] = Field(
        default=None,
        description="Per-sample (field) metadata CSV for the report's Dataset/provenance section "
                    "(location, dates, sites, institution); optional",
    )
    project_metadata: Optional[Path] = Field(
        default=None,
        description="Project metadata CSV for the report's Dataset section "
                    "(recorder, sequencing method, reference DB); optional",
    )

    @field_validator("output_dir", "sample_metadata", "project_metadata")
    @classmethod
    def expand_optional_path(cls, v: Optional[Path]) -> Optional[Path]:
        """Expand ``~`` and resolve the path when set; leave None unchanged.

        Args:
            v: A configured output_dir / sample_metadata / project_metadata path, or None.

        Returns:
            The path with ``~`` expanded and resolved to absolute, or None unchanged.
        """
        return v.expanduser().resolve() if v is not None else v
