"""Load, merge, and validate a marker's YAML configuration into a PipelineConfig.

This is the entry point that turns a per-marker YAML file (under ``config/markers/``) into the
validated :class:`~seednap.config.models.PipelineConfig` the orchestrator runs from. It reads
the YAML, optionally merges it over a defaults file (user keys win; see :func:`merge_configs`),
and validates the result against the strict Pydantic models, humanising any validation error
into an actionable message. It also provides the helpers behind the ``seednap validate`` and
``seednap init`` CLI commands.
"""
from pathlib import Path
from typing import Any, Dict, Optional

import yaml
from pydantic import ValidationError

from seednap.config.models import PipelineConfig


class ConfigError(Exception):
    """Raised for any unreadable, malformed, or invalid pipeline configuration.

    Carries a human-readable, actionable message (e.g. the offending file, the closest-match
    key suggestion, or the fix) so a configuration mistake stops the run at load time with a
    clear explanation rather than a mid-run crash.
    """

    pass


def load_yaml(file_path: Path) -> Dict[str, Any]:
    """
    Load YAML file and return as dictionary.

    Args:
        file_path: Path to YAML file

    Returns:
        Dictionary with configuration

    Raises:
        ConfigError: If file cannot be read or parsed
    """
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            config_dict = yaml.safe_load(f)

        if config_dict is None:
            raise ConfigError(
                f"Config file is empty: {file_path}. A marker config must at minimum "
                f"define marker.name, marker.primers.forward/reverse, taxonomy.method, "
                f"and the selected method's databases.<method> block (you will also "
                f"normally set paths.raw_data to your FASTQ directory). Generate a "
                f"minimal starting template with: seednap init -o {file_path} --force"
            )

        if not isinstance(config_dict, dict):
            raise ConfigError(
                f"Config file {file_path} is not a YAML mapping: its top level parsed "
                f"as a {type(config_dict).__name__}, not key/value pairs. A SeeDNAP "
                f"config must be a mapping with top-level keys like marker:, paths:, "
                f"taxonomy:. Make sure this is a marker config (not a CSV, metadata, or "
                f"sequence file) and that the first non-comment line is a key "
                f"(e.g. 'marker:'), not a '-' list item or a bare value."
            )

        return config_dict

    except ConfigError:
        # Already an actionable message (empty / not-a-mapping); do not re-wrap it
        # into the generic 'Error reading config file' string below.
        raise
    except FileNotFoundError as e:
        raise ConfigError(
            f"Config file not found: {file_path}. Check the path is correct (it is "
            f"resolved relative to the current directory unless absolute). Marker "
            f"configs live under config/markers/; create a new one with "
            f"`seednap init -o {file_path}`."
        ) from e
    except yaml.YAMLError as e:
        raise ConfigError(
            f"Invalid YAML in config file {file_path}: {e} "
            f"[SDN-CFG-009] (run `seednap explain SDN-CFG-009` for more)"
        ) from e
    except Exception as e:
        raise ConfigError(
            f"Could not read config file {file_path}: {e}. Expected a single YAML "
            f"file; if you passed a directory (a common mistake, e.g. "
            f"`config/markers/` instead of `config/markers/teleo.yaml`), point at the "
            f"specific .yaml file instead. Otherwise this is a filesystem error "
            f"reading the path; check `ls -l {file_path}`."
        ) from e


def merge_configs(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively merge two configuration dictionaries.

    The override dictionary takes precedence. Nested dictionaries are merged
    recursively; lists and scalars are REPLACED wholesale, not appended. This is
    why a marker YAML need only specify what differs from the model defaults (any
    field with a default may be omitted), but also means that to change one entry
    of a list (e.g. ``pipeline.steps`` or ``taxonomy.contaminants``) you must
    restate the whole list, not just the changed element.

    Args:
        base: Base configuration dictionary
        override: Override configuration dictionary

    Returns:
        Merged configuration dictionary
    """
    merged = base.copy()

    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            # Recursively merge nested dictionaries
            merged[key] = merge_configs(merged[key], value)
        else:
            # Override value (including lists and primitives)
            merged[key] = value

    return merged


def load_config(
    config_path: Path, defaults_path: Optional[Path] = None, validate: bool = True
) -> PipelineConfig:
    """
    Load and validate pipeline configuration.

    This function:
    1. Loads default configuration (if provided)
    2. Loads user configuration
    3. Merges them (user config takes precedence)
    4. Validates using Pydantic model (only when validate=True)
    5. Returns the validated config object, or the raw merged dict

    Side effect: when validate=True, constructing PipelineConfig runs its
    model_post_init, which calls mkdir on paths.output and paths.logs. Loading
    a config therefore creates those two directories (and can raise on an
    unwritable path); it is not a read-only operation.

    Args:
        config_path: Path to user configuration YAML file
        defaults_path: Optional path to default configuration file
        validate: Whether to validate the configuration (default: True)

    Returns:
        A validated PipelineConfig when validate=True. When validate=False, the
        raw merged config dict is returned instead (the return annotation still
        reads PipelineConfig; this path is for debugging/special cases only).

    Raises:
        ConfigError: If configuration is invalid or cannot be loaded
    """
    # Load user config
    config_dict = load_yaml(config_path)

    # Load and merge defaults if provided
    if defaults_path is not None and defaults_path.exists():
        defaults_dict = load_yaml(defaults_path)
        config_dict = merge_configs(defaults_dict, config_dict)

    # Validate using Pydantic model
    if validate:
        try:
            config = PipelineConfig(**config_dict)
            return config
        except ValidationError as e:
            # Humanize into what/why/fix messages (closest-match suggestions, migration hints
            # for removed keys, valid-key listings). See seednap.errors.config.
            from seednap.errors import humanize_validation_error

            raise ConfigError(humanize_validation_error(e, Path(config_path))) from e
    else:
        # Return unvalidated dict (for debugging or special cases)
        return config_dict  # type: ignore


def validate_config_file(config_path: Path) -> tuple[bool, Optional[str]]:
    """
    Validate a configuration file by fully loading it.

    This performs a complete load_config(validate=True): it merges defaults,
    runs the full Pydantic model, and (via PipelineConfig.model_post_init)
    creates paths.output and paths.logs. It is not a lightweight, side-effect-
    free check; it reports any failure as (False, message) instead of raising.
    Used by CLI validation commands.

    Args:
        config_path: Path to configuration file

    Returns:
        Tuple of (is_valid, error_message)
        - is_valid: True if valid, False otherwise
        - error_message: Error message if invalid, None if valid
    """
    try:
        load_config(config_path, validate=True)
        return True, None
    except ConfigError as e:
        return False, str(e)
    except OSError as e:
        # PipelineConfig.model_post_init calls mkdir on paths.output and paths.logs
        # while loading; a raw PermissionError/OSError here means those paths point
        # somewhere the user cannot create or write.
        return False, (
            f"Could not create the output/log directories declared in the config: "
            f"{e}. SeeDNAP creates paths.output and paths.logs when it loads a config "
            f"(it calls mkdir on both). Check that those two paths point at a location "
            f"you can create and write to: not a read-only mount, not another user's "
            f"directory, and not a path whose parent is missing or is a file rather "
            f"than a directory. On the eDNA server, set paths.output/paths.logs under "
            f"a directory you own (e.g. your home or a run directory you created)."
        )
    except Exception as e:
        return False, f"Unexpected error: {e}"


INIT_TEMPLATES = ("minimal", "full")

# Used when the marker is not in the bundled primers_list.csv.
_PLACEHOLDER_PRIMERS = ("ACACCGCCCGTCACTCT", "CTTCCGGTACACTTACCATG")


def _lookup_primers(marker: str) -> tuple[str, str, bool]:
    """Return the marker's primers from the bundled primers_list.csv.

    Args:
        marker: Marker name, matched case-insensitively against the ``name`` column.

    Returns:
        ``(forward, reverse, found)``. When the marker is not listed, the teleo primers are
        returned as placeholders with ``found=False``.
    """
    import csv
    from importlib import resources

    ref = resources.files("seednap.data.templates").joinpath("primers_list.csv")
    with ref.open(encoding="utf-8") as f:
        for row in csv.DictReader(f, skipinitialspace=True):
            if row["name"].strip().lower() == marker.lower():
                return row["pcr_primer_forward"].strip(), row["pcr_primer_reverse"].strip(), True
    return _PLACEHOLDER_PRIMERS[0], _PLACEHOLDER_PRIMERS[1], False


def create_example_config(
    output_path: Path, marker: str = "teleo", template: str = "minimal"
) -> None:
    """
    Create an example configuration file.

    Args:
        output_path: Where to write the example config
        marker: Marker name for the example (default: 'teleo'). Its primers are filled in
            from the bundled primers_list.csv when the marker is listed there.
        template: ``"minimal"`` (default) writes the fields you normally edit and leaves the
            rest on built-in defaults; ``"full"`` writes every parameter with its
            default value and a short comment.

    Raises:
        ConfigError: If ``template`` is unknown or the file cannot be written
    """
    if template not in INIT_TEMPLATES:
        raise ConfigError(
            f"Unknown init template '{template}'. Choose one of: {', '.join(INIT_TEMPLATES)}."
        )

    forward, reverse, found = _lookup_primers(marker)
    primer_note = "" if found else "   # placeholder: replace with your primer"

    if template == "minimal":
        example_config = f"""# SeeDNAP config for {marker} (minimal). Anything not listed uses its built-in default;
# run `seednap init --full` to see every parameter. Check with `seednap validate <this file>`.

# Steps to run, in order. Use "swarm" instead of "dada2" for OTUs.
# Available: demultiplex (before trim), trim, dada2 | swarm, taxonomy, clean, export, report
pipeline:
  steps:
    - "trim"
    - "dada2"
    - "taxonomy"
    - "export"
    - "report"

marker:
  name: {marker}
  primers:                                   # 5'->3', IUPAC codes allowed
    forward: "{forward}"{primer_note}
    reverse: "{reverse}"{primer_note}

paths:
  raw_data: "/path/to/raw_fastq"             # paired-end FASTQ (<sample>_R1/_R2.fastq.gz) or libraries to demultiplex
  output: "outputs"
  logs: "logs"

# demultiplex:                               # only if "demultiplex" is in pipeline.steps
#   protocol: "ligation"                     # ligation | standard | none
#   metadata: "/path/to/metadata_lab.csv"

trimming:
  cores: 4

taxonomy:
  method: "blast"                            # blast | dada2 | ecotag
  databases:
    blast:
      fasta: "/path/to/{marker}_reference.fasta"
"""
    else:
        example_config = f"""# SeeDNAP config for {marker} (full). Every parameter is listed with its default value.
# Only marker, primers, paths.raw_data, taxonomy.method and that method's database block are
# required: any other line can be deleted and its default is used.
# Check with `seednap validate <this file>`.

# Steps to run: a step runs only if it is listed, in this order. dada2 and swarm are
# mutually exclusive. Available: demultiplex (before trim), trim, dada2 | swarm, taxonomy,
# clean (after dada2/swarm), export (after taxonomy), report
pipeline:
  steps:
    - "trim"
    - "dada2"
    - "taxonomy"
    - "clean"
    - "export"
    - "report"

marker:
  name: {marker}
  description: "{marker} eDNA metabarcoding"
  primers:                                   # 5'->3', IUPAC codes allowed (I is read as N)
    forward: "{forward}"{primer_note}
    reverse: "{reverse}"{primer_note}

paths:
  raw_data: "/path/to/raw_fastq"             # paired-end FASTQ (<sample>_R1/_R2.fastq.gz) or libraries to demultiplex
  output: "outputs"
  logs: "logs"

# Only used if "demultiplex" is added to pipeline.steps (before "trim")
demultiplex:
  protocol: "ligation"                       # ligation | standard | none
  metadata: "/path/to/metadata_lab.csv"      # tag/sample sheet
  max_sample_failure_rate: 0.5               # abort if more than this fraction of samples fail

# Primer trimming (cutadapt)
trimming:
  min_length: 20                             # min read length after trimming
  max_error_rate: 0.1                        # max error rate in the primer match
  cores: 4
  discard_untrimmed: true                    # drop reads without primers
  overlap: 3                                 # min primer overlap

# ASV path: used if "dada2" is in pipeline.steps
dada2:
  filter:
    max_ee: 2.0                              # max expected errors
    trunc_q: 11                              # truncate at first base with quality <= trunc_q
    max_n: 0
    rm_phix: true
    # min_len: 100                           # optional
    # max_len: 200                           # optional
  merge:
    min_overlap: 20
    max_mismatch: 0
  chimera:
    method: "consensus"                      # consensus | pooled | none
  pool: false                                # pool samples for denoising
  multithread: true
  per_library: false                         # learn error models per library, then merge
  collect_metrics: true                      # ASV summary stats to metrics.json/csv

# OTU path: used if "swarm" is in pipeline.steps
swarm:
  merge:
    fastq_maxdiffs: 10                       # max differences in the overlap
    fastq_minovlen: 10                       # min overlap length
    allow_stagger: false
  clustering:
    d: 1                                     # distance threshold
    fastidious: true                         # refine singletons
    boundary: 3                              # min mass of large OTUs (fastidious)
    threads: 4
  chimera:
    method: "denovo"                         # denovo | none
  min_sequence_length: 20                    # min length after merging

taxonomy:
  method: "blast"                            # blast | dada2 | ecotag
  # Flagged in the is_contaminant_candidate column, never removed. CRABS format (Genus_species).
  contaminants:
    - "Homo_sapiens"
    - "Bos_taurus"
    - "Sus_scrofa"
    - "Canis_lupus_familiaris"
    - "Canis_familiaris"
    - "Felis_catus"
    - "Mus_musculus"
    - "Gallus_gallus"
    - "Ovis_aries"
  # Every block present here is validated: keep only the one for your method.
  databases:
    blast:
      fasta: "/path/to/{marker}_reference.fasta"
      task: "megablast"                      # megablast | blastn | dc-megablast | blastn-short
      perc_identity: 80.0
      qcov_hsp_perc: 80.0
      evalue: 1.0e-25
      max_target_seqs: 5
      # Per-rank identity floors: below the floor for a rank, that rank and finer ones are empty
      threshold_species: 99.0
      threshold_genus: 96.0
      threshold_family: 90.0
      threshold_order: 80.0
      threshold_class: 70.0
      lca_algorithm: "cascade"               # cascade | collapsed_taxonomy
      top_bitscore_pct: 10.0                 # cascade: hits within this % of the best bitscore enter the LCA
      lca_pident_delta: 1.0                  # cascade: and within this %id of the best hit
      lca_pid: 90.0                          # collapsed_taxonomy: %id floor
      lca_diff: 1.0                          # collapsed_taxonomy: %id window collapsed to the LCA
    # dada2:
    #   all: "/path/to/dada2_all_ranks.fasta"
    #   species: "/path/to/dada2_species.fasta"
    #   bootstrap_threshold: 80
    # ecotag:
    #   tree: "/path/to/taxonomy_tree/"
    #   fasta: "/path/to/ecotag_db.fasta"
    #   threads: 8                           # parallel ecotag processes (query split in chunks)
    #   timeout: 14400                       # per-process timeout, seconds

# Control decontamination: used if "clean" is in pipeline.steps (controls come from the manifest)
cleaning:
  mode: "flag"                               # flag (annotate only) | subtract (remove control reads)

# GBIF table: used if "export" is in pipeline.steps
export:
  gbif:
    add_rank: true
    add_taxon: true

# Read tracking + HTML report: used if "report" is in pipeline.steps
report:
  html_report: true
  warn_below_retention_pct: 30.0             # warn for samples keeping < this % of raw reads
  warn_step_loss_pct: 70.0                   # warn when one step drops > this % of a sample's reads
  # output_dir: "/path/to/reports"           # default: <output>/04_report/<marker>
  # sample_metadata: "/path/to/metadata_field.csv"
  # project_metadata: "/path/to/metadata_proj.csv"

logging:
  level: "INFO"                              # DEBUG | INFO | WARNING | ERROR
  format: "detailed"                         # simple | detailed | json
  file: true
  console: true
"""

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(example_config)
    except Exception as e:
        raise ConfigError(
            f"Failed to write example config to {output_path}: {e}. SeeDNAP could not "
            f"create or write that file; its parent directory is likely not writable "
            f"by you (owned by another user, or on a read-only mount). Re-run with an "
            f"output path under a directory you own, e.g. `seednap init -o "
            f"~/myconfig.yaml`, and make sure its parent directory exists and is "
            f"writable."
        ) from e
