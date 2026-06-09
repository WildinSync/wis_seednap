"""Taxonomy config: per-method database models and the assignment config."""

from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

from pydantic import Field, ValidationError, field_validator

from seednap.config.models.base import StrictModel

# Each taxonomy method's required database path(s), surfaced in error text so a user
# whose block fails validation is told exactly which key(s) the method needs.
_REQUIRED_DB_PATHS: Dict[str, str] = {
    "blast": "fasta",
    "dada2": "all",
    "ecotag": "tree + fasta",
    "decipher": "trained",
}


def _flatten_db_errors(exc: ValidationError) -> List[str]:
    """Turn a per-method DB-block ValidationError into one readable bullet per problem.

    Replaces dumping Pydantic's raw nested repr (which carries internal ``[type=...]``
    noise) with a flat, declarative list keyed by the offending field.
    """
    bullets: List[str] = []
    for err in exc.errors():
        loc = err.get("loc", ())
        field = ".".join(str(p) for p in loc) if loc else "(block)"
        etype = err.get("type", "")
        ctx = err.get("ctx") or {}
        if etype == "missing":
            problem = "required path is missing"
        elif etype == "extra_forbidden":
            problem = "unknown key (typo? SeeDNAP rejects unrecognised keys)"
        elif etype in ("greater_than_equal", "less_than_equal", "greater_than", "less_than"):
            lo = ctx.get("ge", ctx.get("gt"))
            hi = ctx.get("le", ctx.get("lt"))
            if lo is not None and hi is not None:
                problem = f"out of range (must be between {lo} and {hi})"
            else:
                problem = err.get("msg", "value is out of range")
        else:
            problem = err.get("msg", "invalid value")
        bullets.append(f"  - {field}: {problem}")
    return bullets


# ===========================================================================
# TAXONOMY: one DB model per method, then the assignment config. A config fills
# only the SELECTED method's database block (taxonomy.databases.<method>); the
# others are ignored. _DATABASE_MODELS dispatches both validation and runtime.
# ===========================================================================


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
    # An in-band hit must also be within this many percent-identity points of the best
    # in-band hit to count toward the LCA. Guards against a single near-identity off-target
    # hit (e.g. a 98.6% worm beside 100% Bos hits on a short marker) collapsing the LCA to a
    # high rank. Default 1.0 (eDNAFlow "diff 1"); 0 disables (bitscore band only).
    lca_pident_delta: float = Field(
        default=1.0, ge=0, le=100,
        description="LCA pident floor: in-band hits must be within this many %id points of the best"
    )
    # LCA algorithm. 'cascade' (default) is the header-derived per-rank/MEGAN-LR resolver.
    # 'collapsed_taxonomy' is the eDNAFlow/OceanOmics %identity-window collapse-to-LCA, also
    # header-based (no taxids/taxdump) and fully offline, tuned by lca_pid/lca_diff below.
    # 'fishbase_tiered' is not implemented (fish-specific, needs a bundled WoRMS file + a staged
    # Fishbase table) and raises if selected.
    lca_algorithm: Literal["cascade", "collapsed_taxonomy", "fishbase_tiered"] = Field(
        default="cascade", description="BLAST LCA algorithm (default: cascade = current behavior)"
    )
    # Parameters for lca_algorithm="collapsed_taxonomy" (eDNAFlow/OceanOmics). lca_pid is the
    # hard %identity floor; lca_diff is the top-identity window width within which disagreeing
    # hits are collapsed to their LCA. eDNAFlow defaults: lca_pid=90, lca_diff=1. (Query
    # coverage is enforced separately at the blastn step via qcov_hsp_perc.)
    lca_pid: float = Field(default=90.0, ge=0, le=100, description="collapsed_taxonomy %identity floor")
    lca_diff: float = Field(default=1.0, ge=0, le=100, description="collapsed_taxonomy identity-window width")

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


# Maps each taxonomy method to the strict model that validates its database block. Single source
# of truth for both load-time validation (validate_databases) and runtime dispatch
# (get_database_config), so the two cannot drift.
_DATABASE_MODELS: Dict[str, type] = {
    "dada2": Dada2DatabaseConfig,
    "blast": BlastDatabaseConfig,
    "ecotag": EcotagDatabaseConfig,
    "decipher": DecipherDatabaseConfig,
}


class TaxonomicAssignmentConfig(StrictModel):
    """Taxonomic assignment configuration."""

    method: Literal["dada2", "blast", "ecotag", "decipher"] = Field(
        ..., description="Taxonomic assignment method"
    )
    # Open dict keyed by method name ("blast"/"dada2"/"ecotag"/"decipher"); each value is that
    # method's database block. Only the selected method's block is used at run time
    # (get_database_config), but validate_databases parses EVERY present block into its strict
    # model at load time so a typo or missing path errors during `seednap validate`, not mid-run.
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
        """Validate the database configurations at load time.

        Two checks: the selected ``method`` must have a database block, and every recognised
        block present is parsed into its strict model now (not lazily at the taxonomy step), so a
        typo or a missing required path surfaces during ``seednap validate`` instead of mid-run.
        ``databases`` is an open dict, so ``extra="forbid"`` does not otherwise reach inside it.
        """
        method = info.data.get("method")
        if method is not None and method not in v:
            raise ValueError(
                f"No database configuration found for method '{method}'. "
                f"Please add '{method}' section to databases configuration."
            )
        for name, block in v.items():
            model = _DATABASE_MODELS.get(name)
            if model is None:
                continue  # not a recognised method block; left untouched
            try:
                model(**block)
            except ValidationError as exc:
                bullets = "\n".join(_flatten_db_errors(exc))
                required = _REQUIRED_DB_PATHS.get(name, "its method-specific path(s)")
                raise ValueError(
                    f"Invalid taxonomy.databases.{name} block:\n{bullets}\n"
                    f"The '{name}' database block must list its required path(s) ({required}) "
                    f"and use only recognised keys. For a fully-annotated reference template run: "
                    f"seednap init --full. Note: SeeDNAP validates EVERY database block present, "
                    f"not just the one named by taxonomy.method, so a leftover block for an unused "
                    f"method ('{name}' here) must also be valid -- delete it if it is not needed."
                ) from exc
            except TypeError as exc:
                # model(**block) needs a mapping; a non-dict value (e.g. a bare path or list)
                # cannot be unpacked. Name the offending block and what it must be.
                raise ValueError(
                    f"Invalid taxonomy.databases.{name}: expected a block of key/value settings, "
                    f"got {type(block).__name__}. Make taxonomy.databases.{name} a mapping that "
                    f"lists its required path(s) ({_REQUIRED_DB_PATHS.get(name, 'its paths')}); "
                    f"run `seednap init --full` for a reference template."
                ) from exc
        return v

    def get_database_config(self) -> Any:
        """Return the parsed database config model for the selected method."""
        db_config = self.databases.get(self.method, {})
        model = _DATABASE_MODELS.get(self.method)
        return model(**db_config) if model is not None else db_config
