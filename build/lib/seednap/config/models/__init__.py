"""Pydantic models for the SeeDNAP marker configuration (one YAML per marker).

Every model is a ``StrictModel`` (``extra="forbid"``), so an unknown key errors at load time.
Config is merged over these defaults by ``loader.load_config``, so a YAML only needs to specify
what differs from the defaults.

Required keys (the only ``Field(...)`` without a default): ``marker.name``,
``marker.primers.forward``/``reverse``, ``taxonomy.method``, and the required path(s) inside the
SELECTED database block (``blast.fasta``; ``dada2.all``; ``ecotag.tree`` + ``fasta``;
``decipher.trained``). Everything else has a default and may be omitted.

Section -> consumer map (which pipeline step reads each top-level section):
    marker, paths        -> all steps
    demultiplex          -> run_demultiplex   [runs iff "demultiplex" in pipeline.steps]
    trimming             -> run_trim
    dada2                -> run_dada2          [ASV path; listed in pipeline.steps]
    swarm                -> run_swarm          [OTU path; listed in pipeline.steps]
    taxonomy             -> run_taxonomy / get_database_config
    cleaning             -> run_clean          [runs iff "clean" in pipeline.steps]
    export, report, logging -> their respective steps
    pipeline.steps       -> the single source of truth for what runs and in what order
                            (a stage runs iff listed; dada2 and swarm are mutually exclusive)
"""

from seednap.config.models.base import StrictModel
from seednap.config.models.dada2 import (
    Dada2ChimeraConfig,
    Dada2Config,
    Dada2FilterConfig,
    Dada2MergeConfig,
)
from seednap.config.models.input import (
    DemultiplexConfig,
    MarkerConfig,
    PathsConfig,
    PrimerConfig,
)
from seednap.config.models.operational import (
    CleaningConfig,
    LoggingConfig,
    PipelineStepsConfig,
)
from seednap.config.models.outputs import (
    ExportConfig,
    GbifExportConfig,
    ReportConfig,
)
from seednap.config.models.pipeline import PipelineConfig
from seednap.config.models.swarm import (
    SwarmChimeraConfig,
    SwarmClusteringConfig,
    SwarmConfig,
    SwarmMergeConfig,
)
from seednap.config.models.taxonomy import (
    _DATABASE_MODELS,
    BlastDatabaseConfig,
    Dada2DatabaseConfig,
    DecipherDatabaseConfig,
    EcotagDatabaseConfig,
    TaxonomicAssignmentConfig,
)
from seednap.config.models.trimming import TrimmingConfig

__all__ = [
    "StrictModel",
    "PrimerConfig",
    "MarkerConfig",
    "PathsConfig",
    "DemultiplexConfig",
    "TrimmingConfig",
    "Dada2FilterConfig",
    "Dada2MergeConfig",
    "Dada2ChimeraConfig",
    "Dada2Config",
    "SwarmMergeConfig",
    "SwarmClusteringConfig",
    "SwarmChimeraConfig",
    "SwarmConfig",
    "Dada2DatabaseConfig",
    "BlastDatabaseConfig",
    "EcotagDatabaseConfig",
    "DecipherDatabaseConfig",
    "TaxonomicAssignmentConfig",
    "GbifExportConfig",
    "ExportConfig",
    "ReportConfig",
    "CleaningConfig",
    "LoggingConfig",
    "PipelineStepsConfig",
    "PipelineConfig",
]
