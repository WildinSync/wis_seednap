"""Configuration data models using Pydantic for type-safe configuration."""

import logging
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator

logger = logging.getLogger(__name__)


class StrictModel(BaseModel):
    """Base model that rejects unknown fields to catch config typos."""

    model_config = ConfigDict(extra="forbid")


class PrimerConfig(StrictModel):
    """Primer pair configuration."""

    forward: str = Field(..., min_length=10, description="Forward primer sequence (5' to 3')")
    reverse: str = Field(..., min_length=10, description="Reverse primer sequence (5' to 3')")
    name: Optional[str] = Field(None, description="Primer set name (e.g., 'Teleo')")
    target: Optional[str] = Field(None, description="Target region (e.g., '12S rRNA')")
    amplicon_length: Optional[tuple[int, int]] = Field(
        None, description="Expected amplicon length range [min, max]"
    )

    @field_validator("forward", "reverse")
    @classmethod
    def validate_dna_sequence(cls, v: str) -> str:
        """Validate that the sequence contains only valid DNA bases including IUPAC ambiguity codes."""
        valid_bases = set("ACGTRYMKSWHBVDN")
        v_upper = v.upper()
        if not all(base in valid_bases for base in v_upper):
            invalid_bases = set(v_upper) - valid_bases
            raise ValueError(
                f"Invalid DNA sequence. Contains invalid bases: {invalid_bases}. "
                f"Valid bases are: {', '.join(sorted(valid_bases))}"
            )
        return v_upper

    def reverse_complement(self) -> tuple[str, str]:
        """Get reverse complement of primers (calculated, not stored)."""
        from seednap.utils.sequences import reverse_complement

        return reverse_complement(self.forward), reverse_complement(self.reverse)


class MarkerConfig(StrictModel):
    """Marker-specific configuration (e.g., teleo, 16S, COI)."""

    name: str = Field(..., description="Marker name (lowercase)")
    description: Optional[str] = Field(None, description="Marker description")
    primers: PrimerConfig = Field(..., description="Primer pair configuration")


class PathsConfig(StrictModel):
    """File paths configuration."""

    raw_data: Path = Field(default=Path("data/raw"), description="Raw FASTQ data directory")
    output: Path = Field(default=Path("outputs"), description="Output directory")
    logs: Path = Field(default=Path("logs"), description="Log files directory")
    references: Path = Field(default=Path("references"), description="Reference databases directory")

    @field_validator("raw_data", "output", "logs", "references")
    @classmethod
    def expand_path(cls, v: Path) -> Path:
        """Expand ~ and relative paths to absolute paths."""
        return v.expanduser().resolve()


class DemultiplexConfig(StrictModel):
    """Demultiplexing configuration."""

    enabled: bool = Field(default=False, description="Whether demultiplexing is enabled")
    protocol: Literal["ligation", "standard", "none"] = Field(
        default="none", description="Demultiplexing protocol type"
    )
    metadata: Optional[Path] = Field(None, description="Path to metadata CSV file")
    # When raw inputs are already demultiplexed (one FASTQ per sample),
    # set skip=true so the orchestrator records the step as skipped rather
    # than running the demultiplex protocol against pre-demultiplexed data.
    skip: bool = Field(
        default=False,
        description="Skip the demultiplex step (use when raw inputs are pre-demultiplexed)",
    )
    # If more than this fraction of samples fail during demultiplexing, abort.
    # Otherwise log the failures and continue. Default 0.5 = abort if >50% fail.
    max_sample_failure_rate: float = Field(
        default=0.5, ge=0.0, le=1.0,
        description="Abort demultiplex if more than this fraction of samples fail",
    )

    @field_validator("metadata")
    @classmethod
    def validate_metadata_path(cls, v: Optional[Path], info: Any) -> Optional[Path]:
        """Validate that metadata file exists if demultiplexing is enabled."""
        if v is not None:
            v = v.expanduser().resolve()
            # Note: we don't check file existence here since config might be created
            # before the file exists. Validation happens at runtime.
        return v


class TrimmingConfig(StrictModel):
    """Primer trimming configuration."""

    tool: Literal["cutadapt"] = Field(default="cutadapt", description="Trimming tool to use")
    min_length: int = Field(default=20, ge=1, description="Minimum read length after trimming")
    max_error_rate: float = Field(
        default=0.1, ge=0.0, le=1.0, description="Maximum error rate for primer matching"
    )
    cores: int = Field(default=1, ge=1, description="Number of CPU cores to use")
    discard_untrimmed: bool = Field(
        default=True, description="Discard reads without detected primers"
    )
    overlap: int = Field(default=3, ge=1, description="Minimum overlap for primer detection")


class Dada2FilterConfig(StrictModel):
    """DADA2 filtering parameters."""

    max_ee: float = Field(default=2.0, ge=0, description="Maximum expected errors")
    trunc_q: int = Field(default=11, ge=0, description="Truncate reads at first base with quality <= trunc_q")
    max_n: int = Field(default=0, ge=0, description="Maximum number of N bases allowed")
    rm_phix: bool = Field(default=True, description="Remove PhiX reads")
    min_len: Optional[int] = Field(None, ge=1, description="Minimum read length (optional)")
    max_len: Optional[int] = Field(None, ge=1, description="Maximum read length (optional)")


class Dada2MergeConfig(StrictModel):
    """DADA2 read merging parameters."""

    min_overlap: int = Field(default=20, ge=1, description="Minimum overlap for merging")
    max_mismatch: int = Field(default=0, ge=0, description="Maximum mismatches in overlap region")


class Dada2ChimeraConfig(StrictModel):
    """DADA2 chimera removal parameters."""

    method: Literal["consensus", "pooled", "none"] = Field(
        default="consensus", description="Chimera detection method"
    )


class Dada2Config(StrictModel):
    """DADA2 processing configuration."""

    filter: Dada2FilterConfig = Field(default_factory=Dada2FilterConfig)
    merge: Dada2MergeConfig = Field(default_factory=Dada2MergeConfig)
    chimera: Dada2ChimeraConfig = Field(default_factory=Dada2ChimeraConfig)
    pool: bool = Field(default=False, description="Pool samples for denoising")
    multithread: bool = Field(default=True, description="Use multithreading")


class SwarmMergeConfig(StrictModel):
    """vsearch read merging parameters for SWARM pipeline."""

    fastq_maxdiffs: int = Field(default=10, ge=0, description="Max differences in overlap region")
    fastq_minovlen: int = Field(default=10, ge=1, description="Min overlap length for merging")
    allow_stagger: bool = Field(default=False, description="Allow merging of staggered reads")


class SwarmClusteringConfig(StrictModel):
    """SWARM clustering algorithm parameters."""

    d: int = Field(default=1, ge=1, description="Clustering distance threshold")
    fastidious: bool = Field(default=True, description="Enable fastidious mode (refine singletons)")
    boundary: int = Field(default=3, ge=1, description="Min mass for large OTUs in fastidious mode")
    threads: int = Field(default=4, ge=1, description="Number of threads")


class SwarmChimeraConfig(StrictModel):
    """SWARM chimera detection parameters."""

    method: Literal["denovo", "none"] = Field(
        default="denovo", description="Chimera detection method"
    )


class SwarmConfig(StrictModel):
    """SWARM OTU clustering pipeline configuration."""

    merge: SwarmMergeConfig = Field(default_factory=SwarmMergeConfig)
    clustering: SwarmClusteringConfig = Field(default_factory=SwarmClusteringConfig)
    chimera: SwarmChimeraConfig = Field(default_factory=SwarmChimeraConfig)
    min_sequence_length: int = Field(
        default=20, ge=1, description="Min sequence length after merging"
    )


class Dada2DatabaseConfig(StrictModel):
    """DADA2 taxonomic database configuration."""

    all: Path = Field(..., description="Path to database with all taxonomic ranks")
    species: Optional[Path] = Field(None, description="Path to species-level database")
    # Naive Bayesian classifier bootstrap threshold (Wang 2007 RDP standard).
    # Below this confidence, the rank is nulled and every finer rank is cascaded.
    # 80 is the published recommendation for short rRNA reads (<= 250 bp);
    # eDNA convention generally uses 80 or higher. 50 is too permissive.
    bootstrap_threshold: int = Field(
        default=80, ge=0, le=100,
        description="Minimum RDP bootstrap (%) for a rank to be retained (Wang 2007)"
    )

    @field_validator("all", "species")
    @classmethod
    def expand_path(cls, v: Optional[Path]) -> Optional[Path]:
        """Expand paths."""
        if v is not None:
            return v.expanduser().resolve()
        return v


class BlastDatabaseConfig(StrictModel):
    """BLAST database configuration."""

    fasta: Path = Field(..., description="Path to reference FASTA database")
    perc_identity: float = Field(default=80.0, ge=0, le=100, description="Minimum percent identity")
    qcov_hsp_perc: float = Field(
        default=80.0, ge=0, le=100, description="Minimum query coverage per HSP"
    )
    evalue: float = Field(default=1e-25, gt=0, description="Maximum e-value")
    max_target_seqs: int = Field(default=5, ge=1, description="Maximum number of target sequences")
    # blastn task. 'megablast' (word_size 28) is fastest and the right call for short,
    # high-identity vertebrate amplicons against curated reference DBs. Switch to 'blastn'
    # (word_size 11) for divergent references where the family/order tier of hits matters.
    task: Literal["megablast", "blastn", "dc-megablast", "blastn-short"] = Field(
        default="megablast", description="blastn task type"
    )
    # Thresholds for filtering by taxonomic rank (cascade: below threshold for rank R nulls
    # R and all finer ranks). Defaults follow the field-standard from Pappalardo 2025
    # (Methods Ecol. Evol. 16:2380-2394) with rRNA-marker tweaks (family raised vs eDNAFlow).
    threshold_species: float = Field(default=99.0, ge=0, le=100, description="Species-level identity threshold")
    threshold_genus: float = Field(default=96.0, ge=0, le=100, description="Genus-level identity threshold")
    threshold_family: float = Field(default=90.0, ge=0, le=100, description="Family-level identity threshold")
    threshold_order: float = Field(default=80.0, ge=0, le=100, description="Order-level identity threshold")
    threshold_class: float = Field(default=70.0, ge=0, le=100, description="Class-level identity threshold")
    # LCA top-bitscore band (MEGAN-LR style): hits within this percent of the best
    # bitscore are considered together for LCA resolution. 0 = exact ties only.
    top_bitscore_pct: float = Field(
        default=10.0, ge=0, le=100,
        description="LCA bitscore band as percent of best hit (MEGAN-LR topPercent default: 10.0)"
    )

    @field_validator("fasta")
    @classmethod
    def expand_path(cls, v: Path) -> Path:
        """Expand path."""
        return v.expanduser().resolve()


class EcotagDatabaseConfig(StrictModel):
    """Ecotag database configuration."""

    tree: Path = Field(..., description="Path to NCBI taxonomy tree directory")
    fasta: Path = Field(..., description="Path to reference FASTA database")

    @field_validator("tree", "fasta")
    @classmethod
    def expand_path(cls, v: Path) -> Path:
        """Expand path."""
        return v.expanduser().resolve()


class DecipherDatabaseConfig(StrictModel):
    """DECIPHER database configuration."""

    trained: Path = Field(..., description="Path to trained DECIPHER RDS file")
    threshold: int = Field(
        default=60, ge=0, le=100, description="Confidence threshold for assignment"
    )
    processors: int = Field(default=8, ge=1, description="Number of CPU cores to use")

    @field_validator("trained")
    @classmethod
    def expand_path(cls, v: Path) -> Path:
        """Expand path."""
        return v.expanduser().resolve()


class TaxonomicAssignmentConfig(StrictModel):
    """Taxonomic assignment configuration."""

    method: Literal["dada2", "blast", "ecotag", "decipher"] = Field(
        ..., description="Taxonomic assignment method"
    )
    databases: Dict[str, Any] = Field(default_factory=dict, description="Database configurations")
    # Marker-level contaminant list applied to whichever method is selected.
    # Species names matched against the assigned `species` column get an
    # `is_contaminant_candidate=True` annotation in the output. Rows are
    # NEVER deleted; downstream decides. Use the underscore-separated CRABS
    # format (e.g. "Homo_sapiens"). See Whitmore et al. 2023, Nat. Ecol. Evol.
    contaminants: List[str] = Field(
        default_factory=list,
        description="Species to flag as candidate contaminants (CRABS underscore format)",
    )

    @field_validator("databases")
    @classmethod
    def validate_databases(cls, v: Dict[str, Any], info: Any) -> Dict[str, Any]:
        """Validate that appropriate database config exists for the selected method."""
        if "method" in info.data:
            method = info.data["method"]
            if method not in v:
                raise ValueError(
                    f"No database configuration found for method '{method}'. "
                    f"Please add '{method}' section to databases configuration."
                )
        return v

    def get_database_config(self) -> Any:
        """Get the database config for the selected method."""
        db_config = self.databases.get(self.method, {})

        # Parse into appropriate model based on method
        if self.method == "dada2":
            return Dada2DatabaseConfig(**db_config)
        elif self.method == "blast":
            return BlastDatabaseConfig(**db_config)
        elif self.method == "ecotag":
            return EcotagDatabaseConfig(**db_config)
        elif self.method == "decipher":
            return DecipherDatabaseConfig(**db_config)
        else:
            return db_config


class GbifExportConfig(StrictModel):
    """GBIF export configuration."""

    enabled: bool = Field(default=True, description="Whether to generate GBIF format output")
    add_rank: bool = Field(default=True, description="Add taxonomic rank column")
    add_taxon: bool = Field(default=True, description="Add lowest available taxon column")


class ExportConfig(StrictModel):
    """Output export configuration."""

    formats: List[str] = Field(default=["csv"], description="Output formats to generate")
    gbif: GbifExportConfig = Field(
        default_factory=GbifExportConfig, description="GBIF export settings"
    )


class MetricsConfig(StrictModel):
    """Quality control metrics configuration."""

    generate_plots: bool = Field(default=True, description="Generate QC plots")
    plot_format: Literal["png", "pdf", "svg"] = Field(
        default="png", description="Plot output format"
    )
    metrics: List[str] = Field(
        default=["read_counts", "quality_scores", "length_distribution"],
        description="Metrics to calculate",
    )


class LoggingConfig(StrictModel):
    """Logging configuration."""

    level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = Field(
        default="INFO", description="Logging level"
    )
    format: Literal["simple", "detailed", "json"] = Field(
        default="detailed", description="Log format"
    )
    file: bool = Field(default=True, description="Write logs to file")
    console: bool = Field(default=True, description="Write logs to console")


class PipelineStepsConfig(StrictModel):
    """Pipeline steps configuration."""

    steps: List[str] = Field(
        default=["trim", "dada2", "taxonomy", "export"],
        description="Pipeline steps to execute in order",
    )
    skip: List[str] = Field(default_factory=list, description="Steps to skip")


class PipelineConfig(StrictModel):
    """Complete pipeline configuration."""

    version: str = Field(default="0.1.0", description="Config format version")
    marker: MarkerConfig = Field(..., description="Marker configuration")
    paths: PathsConfig = Field(default_factory=PathsConfig, description="Path configuration")
    demultiplex: DemultiplexConfig = Field(
        default_factory=DemultiplexConfig, description="Demultiplexing configuration"
    )
    trimming: TrimmingConfig = Field(
        default_factory=TrimmingConfig, description="Primer trimming configuration"
    )
    dada2: Dada2Config = Field(default_factory=Dada2Config, description="DADA2 configuration")
    swarm: SwarmConfig = Field(default_factory=SwarmConfig, description="SWARM clustering configuration")
    taxonomy: TaxonomicAssignmentConfig = Field(
        ..., description="Taxonomic assignment configuration"
    )
    export: ExportConfig = Field(
        default_factory=ExportConfig, description="Export configuration"
    )
    metrics: MetricsConfig = Field(
        default_factory=MetricsConfig, description="Metrics configuration"
    )
    logging: LoggingConfig = Field(
        default_factory=LoggingConfig, description="Logging configuration"
    )
    pipeline: PipelineStepsConfig = Field(
        default_factory=PipelineStepsConfig, description="Pipeline steps configuration"
    )

    def model_post_init(self, __context: Any) -> None:
        """Post-initialization validation."""
        # Create output directories if they don't exist
        for path_name in ["output", "logs"]:
            path = getattr(self.paths, path_name)
            path.mkdir(parents=True, exist_ok=True)
