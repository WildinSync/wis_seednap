"""Command-line interface for the seednap eDNA metabarcoding pipeline.

This is the single Click entrypoint a user types at the shell. It exposes the whole
pipeline as a set of subcommands: the end-to-end ``run-pipeline`` plus per-step
commands (``trim``, ``demultiplex``, ``dada2``, ``swarm``, ``assign-taxonomy``),
config helpers (``init``, ``validate``, ``explain``), and post-processing commands
(``format-gbif``, ``create-gbif``, ``report``, ``manifest``, ``clean``, ``monitor``).

In pipeline terms this module is the thin user-facing shell. Each command parses
options, sets up logging, then delegates the actual biology (primer trimming, ASV
denoising with DADA2, OTU clustering with SWARM, BLAST/ecotag/DECIPHER taxonomy,
DarwinCore/GBIF export) to the processors and runners under ``seednap.steps`` and the
orchestrator under ``seednap.pipeline``. No biological computation happens here; this
file only wires arguments, reports progress to the console, persists [WARN] safety
messages to a per-command log file, and turns exceptions into actionable error text.
"""

import sys
from pathlib import Path
from typing import Any, Dict, Optional

import click
from rich.console import Console
from rich.table import Table

from seednap.__version__ import __version__
from seednap.config import ConfigError, create_example_config, load_config, validate_config_file
from seednap.utils.logging import get_logger, setup_logging

console = Console()

# Module logger. Routing print_warning through this (rather than console-only) means every
# user-facing warning is also written to whatever log file the running command configured,
# so the no-silent-fallback [WARN]s are persisted, not just printed. The root logger always
# has a console handler (set up in main()), so logger.warning stays visible on the console too.
logger = get_logger(__name__)

# Set by main() from the -v/--verbose flag. Command error handlers use it to show a full
# Python traceback only when the user asked for verbose output; otherwise the actionable
# message (often the external tool's own stderr) is what they see, not a buried stack trace.
# Newer commands call _maybe_traceback() (which reads this global); older command handlers
# (format-gbif, create-gbif, blast, manifest, clean) inline the same traceback dump but read
# ctx.obj["verbose"] instead. Both carry the same flag, set together in main().
_VERBOSE = False

# Logging level/console settings chosen in main() from -v/-q. Stored so a per-command
# call to _add_command_log_file() can re-run setup_logging with the same level/console
# behavior while adding a file handler.
_LOG_LEVEL = "INFO"
_LOG_CONSOLE = True


def _add_command_log_file(output_dir: Optional[Path], name: str) -> None:
    """Persist this command's logs (including [WARN]s) to a file under the run tree.

    main() only configures a console handler, so without this every standalone command's
    safety warnings would be console-only. Each command that knows where its outputs go
    calls this to add a file handler at <output_dir>/logs/<name>.log (or logs/<name>.log
    when the command has no output directory). Level and console behavior are preserved
    from main()'s -v/-q choice.

    A failure to set up the file handler must not abort the command, but it also must not
    be silent (CLAUDE.md no-silent-fallback rule): on error we warn and keep logging to the
    console only.

    Args:
        output_dir: The command's run-output directory; the log goes to
            ``<output_dir>/logs/<name>.log``. Pass ``None`` for commands that have no
            output tree, in which case the log goes to ``logs/<name>.log`` under the
            current working directory.
        name: Base name for the log file (no extension), typically the command name or
            ``<command>_<marker>``.

    Returns:
        None. Reconfigures the root logger in place to add the file handler.

    Raises:
        Nothing. An ``OSError`` while creating the log file is caught and downgraded to a
        printed [WARN]; logging then continues to the console only.
    """
    base = Path(output_dir) if output_dir is not None else Path("logs")
    log_dir = base / "logs" if output_dir is not None else base
    log_file = log_dir / f"{name}.log"
    try:
        setup_logging(level=_LOG_LEVEL, log_file=log_file, console_output=_LOG_CONSOLE)
    except OSError as e:
        print(
            f"[WARN] command log setup: expected=writable log file at {log_file}, "
            f"got={type(e).__name__}: {e}, fallback=console-only logging (warnings not persisted)",
            flush=True,
        )


def _maybe_traceback() -> None:
    """Print the current exception's full Python traceback, but only in verbose (-v) mode.

    Lets command error handlers show the actionable message (often the external tool's own
    stderr) by default, and surface the buried stack trace only when the user asked for it.
    Must be called from within an ``except`` block, as it reads the active exception.

    Returns:
        None. Writes the formatted traceback to the console when ``_VERBOSE`` is set;
        otherwise does nothing.
    """
    if _VERBOSE:
        import traceback

        console.print(traceback.format_exc())


def print_error(message: str) -> None:
    """Print an error message to the console, prefixed and styled in red.

    Args:
        message: The human-readable error text to show. Console-only; this does not
            persist to the log file and does not exit the process.

    Returns:
        None.
    """
    console.print(f"[bold red]Error:[/bold red] {message}")


def print_success(message: str) -> None:
    """Print a success message to the console with a green check mark.

    Args:
        message: The human-readable success text to show. Console-only.

    Returns:
        None.
    """
    console.print(f"[bold green]✓[/bold green] {message}")


def print_warning(message: str) -> None:
    """Emit a warning that is both shown on the console and persisted to the log file.

    Routes through the configured logger instead of writing straight to the rich console, so
    that the warning lands in whatever per-command log file is active (the no-silent-fallback
    rule: warnings must be captured, not console-only). The root logger's console handler keeps
    it visible to the user; the "[WARN]" prefix matches the convention used elsewhere in the
    codebase (CLAUDE.md section 4) so all warning channels read consistently in the log file.

    Args:
        message: The warning text. A "[WARN] " prefix is prepended before logging, so do
            not include one in the caller's message.

    Returns:
        None.
    """
    logger.warning(f"[WARN] {message}")


def _assign_kwargs_from_config(config: Any, method: str) -> Dict[str, Any]:
    """Build TaxonomicAssigner.assign_taxonomy kwargs from a marker config's selected
    method block, so the standalone `assign-taxonomy --config` matches `run-pipeline`.

    This mirrors the per-method kwargs the orchestrator passes (see
    pipeline/orchestrator.py run_taxonomy): for BLAST that includes the BLAST search
    params (perc_identity, qcov_hsp_perc, evalue, max_target_seqs, task) that the
    standalone command otherwise leaves at the assigner defaults -- the divergence this
    fixes. `contaminants` is the marker-level list, applied to every method just as the
    orchestrator does. Caller layers explicit CLI overrides on top of this dict.

    Args:
        config: A loaded ``PipelineConfig`` for the marker. Only its ``taxonomy`` block is
            read (the selected database config and the marker-level ``contaminants`` list).
        method: Which taxonomy method's parameter block to extract: one of ``"blast"``,
            ``"dada2"``, ``"ecotag"``, or ``"decipher"``.

    Returns:
        A dict of keyword arguments for ``TaxonomicAssigner.assign_taxonomy`` matching the
        requested method (e.g. reference paths, per-rank percent-identity thresholds and
        LCA parameters for BLAST; RDP and species DB paths for DADA2). Returns an empty
        dict for any unrecognized method. ``contaminants`` is included for every known
        method. Names of taxa to treat as laboratory contaminants and drop from results.

    Raises:
        AttributeError: If the resolved database config has no attribute one of the
            method branches reads (e.g. the method's database block is missing or only a
            partial dict, so ``db.fasta`` and similar accesses fail).
        pydantic.ValidationError: Propagated from ``get_database_config()`` if the
            method's database block in the config is malformed.
    """
    db = config.taxonomy.get_database_config()
    contaminants = config.taxonomy.contaminants
    if method == "blast":
        return {
            "reference_fasta": db.fasta,
            "threshold_species": db.threshold_species,
            "threshold_genus": db.threshold_genus,
            "threshold_family": db.threshold_family,
            "threshold_order": db.threshold_order,
            "threshold_class": db.threshold_class,
            "top_bitscore_pct": db.top_bitscore_pct,
            "lca_pident_delta": db.lca_pident_delta,
            "lca_algorithm": db.lca_algorithm,
            "lca_pid": db.lca_pid,
            "lca_diff": db.lca_diff,
            "contaminants": contaminants,
            "perc_identity": db.perc_identity,
            "qcov_hsp_perc": db.qcov_hsp_perc,
            "evalue": db.evalue,
            "max_target_seqs": db.max_target_seqs,
            "task": db.task,
        }
    if method == "dada2":
        return {
            "rdp_db_path": db.all,
            "species_db_path": db.species,
            "bootstrap_threshold": db.bootstrap_threshold,
            "contaminants": contaminants,
        }
    if method == "ecotag":
        return {
            "taxonomy_db": db.tree,
            "reference_db": db.fasta,
            "contaminants": contaminants,
        }
    if method == "decipher":
        return {
            "trained_classifier_path": db.trained,
            "threshold": db.threshold,
            "processors": db.processors,
            "contaminants": contaminants,
        }
    return {}


@click.group()
@click.version_option(version=__version__, prog_name="seednap")
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Enable verbose logging (DEBUG level)",
)
@click.option(
    "--quiet",
    "-q",
    is_flag=True,
    help="Suppress console output (only errors)",
)
@click.pass_context
def main(ctx: click.Context, verbose: bool, quiet: bool) -> None:
    """
    seednap: Modern eDNA metabarcoding pipeline with DADA2.

    A pipeline for processing eDNA metabarcoding data with support
    for multiple taxonomic assignment methods.

    This is the top-level Click group: it runs before any subcommand to capture the global
    verbosity flags, stash them on the Click context and in module globals (so command error
    handlers and per-command log setup can read them), and configure logging once. Each
    subcommand may add its own file handler on top.

    Args:
        ctx: The Click context. Used to ensure and populate ``ctx.obj`` (a dict carrying
            ``verbose`` and ``quiet`` for subcommands).
        verbose: ``-v/--verbose``. Enable DEBUG-level logging and full tracebacks on error.
        quiet: ``-q/--quiet``. Suppress console output below WARNING level. ``verbose``
            takes precedence over ``quiet`` for the log level.

    Returns:
        None.
    """
    # Store options in context for subcommands
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet
    global _VERBOSE, _LOG_LEVEL, _LOG_CONSOLE
    _VERBOSE = verbose

    # Setup basic logging (subcommands may reconfigure)
    level = "DEBUG" if verbose else "WARNING" if quiet else "INFO"
    _LOG_LEVEL = level
    _LOG_CONSOLE = not quiet
    setup_logging(level=level, console_output=not quiet)


@main.command()
@click.argument("config_file", type=click.Path(exists=True, path_type=Path))
@click.pass_context
def validate(ctx: click.Context, config_file: Path) -> None:
    """
    Validate a configuration file.

    CONFIG_FILE: Path to the configuration YAML file to validate.

    This command checks:
    - YAML syntax is valid
    - All required fields are present
    - Field types and values are correct
    - Referenced paths and files exist (where applicable)

    A "valid" config that loads but points at missing inputs (raw FASTQ, a taxonomy
    database file) is still rejected here, via a read-only preflight check, so the user
    learns about a broken path before launching a long run rather than mid-run.

    Args:
        ctx: The Click context (carries the global verbose/quiet flags).
        config_file: Path to the marker configuration YAML to validate. Must already exist
            (enforced by Click).

    Returns:
        None. Exits the process with status 0 if the config is well-formed and all
        referenced inputs are usable, or status 1 otherwise.

    Raises:
        SystemExit: Always, via ``sys.exit``. Code 0 on success; code 1 if the schema is
            invalid or preflight finds unusable referenced inputs.
    """
    console.print(f"\n[bold]Validating configuration:[/bold] {config_file}\n")

    is_valid, error_message = validate_config_file(config_file)

    if is_valid:
        # Schema is well-formed; now show the summary and run preflight (referenced files /
        # database block) so a config that loads but points at missing inputs FAILS here rather
        # than mid-run. The success banner is deferred until preflight also passes.
        preflight_problems: list = []
        try:
            config = load_config(config_file)

            table = Table(title="Configuration Summary", show_header=True, header_style="bold cyan")
            table.add_column("Setting", style="cyan")
            table.add_column("Value", style="white")

            table.add_row("Marker", config.marker.name)
            table.add_row("Taxonomic Method", config.taxonomy.method)
            table.add_row("Output Directory", str(config.paths.output))
            table.add_row("Trimming Cores", str(config.trimming.cores))

            if "demultiplex" in config.pipeline.steps:
                table.add_row("Demultiplexing", f"Enabled ({config.demultiplex.protocol})")
            else:
                table.add_row("Demultiplexing", "Disabled")

            # Surface the database actually used for the selected method, and flag any referenced
            # path missing on disk (a config can be valid yet point at a file that is not there).
            # Read-only checks; nothing is created.
            def _exists(p: Path) -> str:
                """Return a colored found/MISSING status string for a path on disk.

                Args:
                    p: Filesystem path to check.

                Returns:
                    A rich-markup string: green "found" if the path exists, else red
                    "MISSING".
                """
                return "[green]found[/green]" if Path(p).exists() else "[red]MISSING[/red]"

            try:
                db_cfg = config.taxonomy.get_database_config()
                for field in type(db_cfg).model_fields:
                    val = getattr(db_cfg, field)
                    if isinstance(val, Path):
                        table.add_row(
                            f"DB ({config.taxonomy.method}.{field})", f"{val}  [{_exists(val)}]"
                        )
            except Exception as exc:  # never crash the summary; load validation already passed
                table.add_row(f"DB ({config.taxonomy.method})", f"[red]unresolved: {exc}[/red]")

            table.add_row("Raw data", f"{config.paths.raw_data}  [{_exists(config.paths.raw_data)}]")

            console.print(table)
            console.print()

            from seednap.errors import preflight_checks

            preflight_problems = preflight_checks(config)

        except Exception as e:
            print_warning(f"Could not load config for summary: {e}")

        if preflight_problems:
            print_error(
                "The config is well-formed, but referenced inputs are not usable "
                "(this would fail mid-run):\n"
            )
            for problem in preflight_problems:
                console.print(problem.render())
                console.print()
            sys.exit(1)

        print_success("Configuration is valid!")
        sys.exit(0)
    else:
        print_error("Configuration validation failed!\n")
        console.print(error_message)
        sys.exit(1)


@main.command()
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("config/markers/example.yaml"),
    help="Output path for example config",
)
@click.option(
    "--marker",
    "-m",
    default="teleo",
    help="Marker name for the example config",
)
@click.option(
    "--minimal/--full",
    default=True,
    help="Emit only the required fields (default) or the fully-annotated reference template",
)
@click.option(
    "--force",
    "-f",
    is_flag=True,
    help="Overwrite existing file",
)
def init(output: Path, marker: str, minimal: bool, force: bool) -> None:
    """
    Create an example configuration file.

    By default this writes a minimal config containing only the required fields (everything
    else uses built-in defaults); pass --full for the fully-annotated reference template.

    A marker config is the per-marker YAML that drives a whole run (primers, paths,
    trimming/DADA2/SWARM parameters, the taxonomy method and its reference databases). This
    command scaffolds one to edit rather than writing it by hand.

    Args:
        output: Path to write the example config to. Defaults to
            ``config/markers/example.yaml``.
        marker: Marker name to seed the example with (e.g. ``teleo``). Sets the marker
            block in the generated config.
        minimal: If True (``--minimal``, the default), emit only required fields; if False
            (``--full``), emit the fully-annotated reference template.
        force: If True (``--force``), overwrite an existing file at ``output``.

    Returns:
        None. Writes the config file and prints next-step hints on success.

    Raises:
        SystemExit: Code 1 if ``output`` already exists and ``force`` is not set, or if
            creating the example config raises ``ConfigError``.
    """
    if output.exists() and not force:
        print_error(f"File already exists: {output}")
        console.print("Use --force to overwrite.")
        sys.exit(1)

    try:
        create_example_config(output, marker=marker, minimal=minimal)
        print_success(f"Created example configuration: {output}")
        console.print("\nEdit this file to customize for your analysis.")
        console.print(f"Validate it with: [bold]seednap validate {output}[/bold]")
    except ConfigError as e:
        print_error(f"Failed to create config: {e}")
        sys.exit(1)


@main.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--format",
    "-f",
    "format_type",
    type=click.Choice(["dada2", "ecotag", "blast", "decipher"]),
    required=True,
    help="Input format type",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help="Output file path (default: input_file with _gbif_input suffix)",
)
@click.pass_context
def format_gbif(ctx: click.Context, input_file: Path, format_type: str, output: Optional[Path]) -> None:
    """
    Convert taxonomic assignment results to GBIF format.

    INPUT_FILE: Path to the taxonomic assignment CSV file.

    Transforms the wide-format table to long-format GBIF-compatible output.
    Adds 'rank' and 'taxon' columns, filters zero counts, and renames columns
    to match GBIF standards (eventID instead of filter_code).

    GBIF is the Global Biodiversity Information Facility, the public repository the lab's
    occurrence records are submitted to. This is the first of the two export steps: it
    reshapes one method's taxonomy+counts table into the long, per-observation form GBIF
    expects (one row per taxon per sample), keeping only non-zero read counts.

    Args:
        ctx: The Click context (carries the global verbose flag, read for traceback depth).
        input_file: Path to the taxonomic assignment CSV produced by an ``assign-taxonomy``
            run. Must already exist (enforced by Click).
        format_type: Which producing method wrote ``input_file``: one of ``"dada2"``,
            ``"ecotag"``, ``"blast"``, or ``"decipher"``. Selects the matching reshape
            logic and must match the method that generated the file.
        output: Path for the long-format GBIF CSV. If ``None``, defaults to the input
            file's directory with a ``_gbif_input`` suffix.

    Returns:
        None. Writes the reshaped CSV and prints record/eventID counts and a rank
        distribution on success.

    Raises:
        SystemExit: Code 1 if the input file is not found, if its columns do not match the
            chosen ``format_type``, or on any other conversion failure.
    """
    from seednap.steps.formatting.gbif_formatter import GBIFFormatter

    _add_command_log_file(output.parent if output is not None else input_file.parent, "format_gbif")

    console.print(f"\n[bold]Converting to GBIF format:[/bold] {input_file}")
    console.print(f"Input format: {format_type}\n")

    try:
        # Determine output path if not provided
        if output is None:
            output = input_file.parent / f"{input_file.stem}_gbif_input.csv"

        formatter = GBIFFormatter()
        df_out = formatter.from_method(format_type, input_file, output)

        # Print success message with stats
        print_success("Converted to GBIF format!")
        console.print(f"\nOutput file: [cyan]{output}[/cyan]")
        console.print(f"Total records: [green]{len(df_out)}[/green]")
        console.print(f"Unique eventIDs: [green]{df_out['eventID'].nunique()}[/green]")

        # Show rank distribution
        if "rank" in df_out.columns:
            rank_counts = df_out["rank"].value_counts()
            console.print("\n[bold]Rank distribution:[/bold]")
            for rank, count in rank_counts.items():
                console.print(f"  {rank}: {count}")

    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except ValueError as e:
        print_error(
            f"Invalid input for format '{format_type}': {e}. The CSV does not have the "
            f"columns this format expects. -f dada2/blast/decipher all expect a wide "
            f"taxonomy table with columns kingdom,phylum,class,order,family,genus,species,"
            f"sequence plus one numeric column per sample; -f ecotag expects an "
            f"ecotag-derived CSV with *_name columns. Re-check that the file came from the "
            f"matching `assign-taxonomy` run and that -f matches the producing method."
        )
        sys.exit(1)
    except Exception as e:
        print_error(f"Failed to convert file: {e}")
        if ctx.obj.get("verbose"):
            import traceback

            console.print(traceback.format_exc())
        sys.exit(1)


@main.command("create-gbif")
@click.argument("taxonomy_results", type=click.Path(exists=True, path_type=Path))
@click.argument("sample_metadata", type=click.Path(exists=True, path_type=Path))
@click.argument("project_metadata", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
@click.option(
    "--summarise-pcr/--no-summarise-pcr",
    default=False,
    help="Summarise PCR replicates by sample before building",
)
@click.option(
    "--skip-enrichment",
    is_flag=True,
    default=False,
    help="Skip NCBI/WORMS taxonomy enrichment (kingdom/phylum lookup)",
)
@click.pass_context
def create_gbif(
    ctx: click.Context,
    taxonomy_results: Path,
    sample_metadata: Path,
    project_metadata: Path,
    output: Path,
    summarise_pcr: bool,
    skip_enrichment: bool,
) -> None:
    """
    Build a DarwinCore-compliant GBIF occurrence CSV.

    Takes three input files and produces a single DarwinCore CSV with all
    required columns populated.

    \b
    TAXONOMY_RESULTS: Taxonomy CSV from format-gbif step (long format with
                      class, order, family, genus, species, taxon, rank,
                      sequence, nb_reads, eventID columns).
    SAMPLE_METADATA:  Per-sample metadata CSV (eventID, lat/lon, eventDate,
                      env_medium, samp_size, depth, size_frac).
    PROJECT_METADATA: Per-project metadata CSV (marker, recordedby, seqmet,
                      identificationRemarks, identificationReferences,
                      otu_seq_comp_appr, otu_db, chimera_check).
    OUTPUT:           Path for the output DarwinCore CSV file.

    \b
    Set the NCBI_API_KEY environment variable (or in a .env file) to enable
    automatic kingdom/phylum enrichment via NCBI Entrez and WORMS.

    DarwinCore is the standardized biodiversity-data vocabulary GBIF ingests. This is the
    second export step: it joins the long-format taxonomy table with per-sample field
    metadata (where, when, how each water sample was taken) and per-project metadata
    (marker, sequencing method, reference DB), optionally enriching higher taxonomy
    (kingdom/phylum) from NCBI/WORMS, into one occurrence CSV ready for submission.

    Args:
        ctx: The Click context (carries the global verbose flag, read for traceback depth).
        taxonomy_results: Long-format taxonomy CSV from the ``format-gbif`` step (columns
            class, order, family, genus, species, taxon, rank, sequence, nb_reads,
            eventID). Must exist.
        sample_metadata: Per-sample (field) metadata CSV keyed by eventID, with lat/lon,
            eventDate, env_medium, samp_size, depth and size_frac. Must exist.
        project_metadata: Per-project metadata CSV with marker, recordedby, seqmet and the
            identification/OTU/chimera provenance fields. Must exist.
        output: Path for the output DarwinCore occurrence CSV.
        summarise_pcr: If True (``--summarise-pcr``), collapse PCR replicates to one row per
            sample before building. PCR replicates are repeat amplifications of the same
            sample. Defaults to False.
        skip_enrichment: If True (``--skip-enrichment``), do not look up kingdom/phylum via
            NCBI Entrez / WORMS. Defaults to False.

    Returns:
        None. Writes the DarwinCore CSV and prints its path on success.

    Raises:
        SystemExit: Code 1 if an input file is missing, on a metadata validation error, or
            on any other build failure.
    """
    from seednap.steps.formatting.darwincore_builder import DarwinCoreBuilder

    _add_command_log_file(output.parent, "create_gbif")

    console.print("\n[bold]Building DarwinCore GBIF CSV[/bold]")
    console.print(f"  Taxonomy results:  {taxonomy_results}")
    console.print(f"  Sample metadata:   {sample_metadata}")
    console.print(f"  Project metadata:  {project_metadata}")
    console.print(f"  Output:            {output}")
    if summarise_pcr:
        console.print("  PCR replicates:    will be summarised")
    if skip_enrichment:
        console.print("  Enrichment:        skipped")
    console.print()

    try:
        builder = DarwinCoreBuilder(
            taxonomy_results_path=taxonomy_results,
            sample_metadata_path=sample_metadata,
            project_metadata_path=project_metadata,
            output_path=output,
            summarise_pcr_replicates=summarise_pcr,
            skip_enrichment=skip_enrichment,
        )
        result_path = builder.build()

        print_success(f"DarwinCore CSV written to {result_path}")

    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except ValueError as e:
        print_error(f"Validation error: {e}")
        sys.exit(1)
    except Exception as e:
        print_error(f"Failed to build DarwinCore CSV: {e}")
        if ctx.obj.get("verbose"):
            import traceback

            console.print(traceback.format_exc())
        sys.exit(1)


@main.command("wis-metadata")
@click.option(
    "--database-url",
    envvar="WIS_DATABASE_URL",
    required=True,
    help="SQLAlchemy URL for the WIS database, or set WIS_DATABASE_URL "
    "(e.g. postgresql://user:pass@host:5432/wis).",
)
@click.option(
    "--marker",
    required=True,
    help="Marker name (e.g. teleo); used for the project row and the output filenames.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    required=True,
    help="Directory for <marker>_sample_metadata.csv and <marker>_project_metadata.csv.",
)
@click.option(
    "--monitoring",
    default=None,
    help="Restrict to one WIS monitoring_id (the site / long-term project).",
)
@click.option(
    "--mission",
    default=None,
    help="Restrict to one WIS mission_id (the sampling campaign).",
)
@click.option(
    "--event-id-field",
    type=click.Choice(["sample_id", "material_sample_id"]),
    default="sample_id",
    show_default=True,
    help="Which WIS identifier becomes eventID; match your FASTQ/sample naming.",
)
@click.option(
    "--recorded-by",
    required=True,
    help="DwC recordedBy (data contributor) for the project row.",
)
@click.option(
    "--identification-remarks",
    required=True,
    help="Identification-method note for the project row.",
)
@click.option(
    "--identification-references",
    required=True,
    help="Reference-DB / method citation for the project row.",
)
@click.option(
    "--seq-meth", default="", help="Optional sequencing-method description (DwC seq_meth)."
)
@click.option(
    "--otu-seq-comp-appr", default="", help="Optional OTU/ASV sequence-comparison approach."
)
@click.pass_context
def wis_metadata(
    ctx: click.Context,
    database_url: str,
    marker: str,
    output_dir: Path,
    monitoring: Optional[str],
    mission: Optional[str],
    event_id_field: str,
    recorded_by: str,
    identification_remarks: str,
    identification_references: str,
    seq_meth: str,
    otu_seq_comp_appr: str,
) -> None:
    """Generate the GBIF export's metadata CSVs from the WIS database.

    Reads per-sample field metadata (eventID, eventDate, coordinates, env_medium, depth, size)
    from the WIS PostgreSQL/PostGIS database and writes the two CSVs the DarwinCore export
    consumes: ``<marker>_sample_metadata.csv`` (one row per sample) and
    ``<marker>_project_metadata.csv`` (one project row). Point ``report.sample_metadata`` /
    ``report.project_metadata`` (or the ``create-gbif`` arguments) at the generated files.

    \b
    Requires the optional database extra:  pip install 'seednap[wis]'
    (adds SQLAlchemy + psycopg2). The reference-database and chimera-removal provenance are
    filled by the 'darwincore' pipeline step from the run config, so they are not written here.

    Args:
        ctx: Click context (carries the global verbose flag).
        database_url: SQLAlchemy URL for the WIS database (or the WIS_DATABASE_URL env var).
        marker: Marker name for the project row and the output filenames.
        output_dir: Directory for the two output CSVs.
        monitoring: Optional ``monitoring_id`` (site/project) filter.
        mission: Optional ``mission_id`` (campaign) filter.
        event_id_field: WIS identifier used as eventID (``sample_id`` or ``material_sample_id``).
        recorded_by: DwC ``recordedBy`` for the project row.
        identification_remarks: Identification-method note for the project row.
        identification_references: Reference-DB / method citation for the project row.
        seq_meth: Optional sequencing-method description.
        otu_seq_comp_appr: Optional OTU/ASV sequence-comparison approach.

    Returns:
        None. Writes the two CSVs and prints their paths.

    Raises:
        SystemExit: Code 1 if the optional dependency is missing, no samples match the
            selector, or any other failure occurs.
    """
    from seednap.steps.formatting.wis_metadata import WisMetadataExporter

    _add_command_log_file(output_dir, "wis_metadata")

    console.print("\n[bold]Generating GBIF metadata from the WIS database[/bold]")
    console.print(f"  Marker:        {marker}")
    console.print(f"  Output dir:    {output_dir}")
    console.print(f"  Filter:        monitoring={monitoring or '—'}, mission={mission or '—'}")
    console.print(f"  eventID field: {event_id_field}")
    console.print()

    try:
        sample_csv, project_csv = WisMetadataExporter(database_url).export(
            output_dir=output_dir,
            marker=marker,
            recorded_by=recorded_by,
            identification_remarks=identification_remarks,
            identification_references=identification_references,
            monitoring=monitoring,
            mission=mission,
            event_id_field=event_id_field,
            seq_meth=seq_meth,
            otu_seq_comp_appr=otu_seq_comp_appr,
        )
        print_success(f"Wrote sample metadata to {sample_csv}")
        print_success(f"Wrote project metadata to {project_csv}")
        console.print(
            "\nNext: set report.sample_metadata / report.project_metadata to these files "
            "(or pass them to create-gbif)."
        )
    except ValueError as e:
        print_error(f"Validation error: {e}")
        sys.exit(1)
    except RuntimeError as e:
        print_error(str(e))
        sys.exit(1)
    except Exception as e:
        print_error(f"Failed to generate WIS metadata: {e}")
        if ctx.obj.get("verbose"):
            import traceback

            console.print(traceback.format_exc())
        sys.exit(1)


@main.command()
@click.argument("query_fasta", type=click.Path(exists=True, path_type=Path))
@click.argument("ref_fasta", type=click.Path(exists=True, path_type=Path))
@click.argument("asv_count", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help="Output CSV file path (default: query_fasta with _blast_taxonomy suffix)",
)
@click.option(
    "--perc-identity",
    default=80.0,
    type=float,
    help="Minimum percent identity for BLAST hits (default: 80.0)",
)
@click.option(
    "--qcov-hsp-perc",
    default=80.0,
    type=float,
    help="Minimum query coverage per HSP (default: 80.0)",
)
@click.option(
    "--evalue",
    default=1e-25,
    type=float,
    help="Maximum e-value for BLAST hits (default: 1e-25)",
)
@click.option(
    "--threshold-species",
    default=99.0,
    type=float,
    help="Minimum percent identity for species-level assignment (default: 99.0)",
)
@click.option(
    "--threshold-genus",
    default=96.0,
    type=float,
    help="Minimum percent identity for genus-level assignment (default: 96.0)",
)
@click.option(
    "--threshold-family",
    default=90.0,
    type=float,
    help="Minimum percent identity for family-level assignment (default: 90.0)",
)
@click.option(
    "--threshold-order",
    default=80.0,
    type=float,
    help="Minimum percent identity for order-level assignment (default: 80.0)",
)
@click.option(
    "--threshold-class",
    default=70.0,
    type=float,
    help="Minimum percent identity for class-level assignment (default: 70.0)",
)
@click.option(
    "--top-bitscore-pct",
    default=10.0,
    type=float,
    help="cascade LCA: include hits within this % of the best bitscore (MEGAN-LR; default: 10.0)",
)
@click.option(
    "--lca-pident-delta",
    default=1.0,
    type=float,
    help="cascade LCA: in-band hits must be within this %id of the best in-band hit (default: 1.0)",
)
@click.option(
    "--task",
    default="megablast",
    type=click.Choice(["megablast", "blastn", "dc-megablast", "blastn-short"]),
    help="blastn task (default: megablast)",
)
@click.option(
    "--lca-algorithm",
    default="cascade",
    type=click.Choice(["cascade", "collapsed_taxonomy"]),
    help="LCA algorithm: cascade (per-rank thresholds, default) or collapsed_taxonomy "
    "(eDNAFlow/OceanOmics %identity-window collapse-to-LCA)",
)
@click.option(
    "--lca-pid",
    default=90.0,
    type=float,
    help="collapsed_taxonomy only: hard %identity floor for hits (default: 90.0)",
)
@click.option(
    "--lca-diff",
    default=1.0,
    type=float,
    help="collapsed_taxonomy only: identity-window width collapsed to the LCA (default: 1.0)",
)
@click.pass_context
def blast(
    ctx: click.Context,
    query_fasta: Path,
    ref_fasta: Path,
    asv_count: Path,
    output: Optional[Path],
    perc_identity: float,
    qcov_hsp_perc: float,
    evalue: float,
    threshold_species: float,
    threshold_genus: float,
    threshold_family: float,
    threshold_order: float,
    threshold_class: float,
    top_bitscore_pct: float,
    lca_pident_delta: float,
    task: str,
    lca_algorithm: str,
    lca_pid: float,
    lca_diff: float,
) -> None:
    """
    Run BLAST taxonomic assignment with LCA resolution.

    QUERY_FASTA: Path to query sequences (ASVs from DADA2, e.g., query.fasta)
    REF_FASTA: Path to reference database FASTA file
    ASV_COUNT: Path to ASV count table CSV (seqtab_clean.csv from DADA2)

    Reported as three console stages:
    1. Run blastn search (builds the BLAST DB from REF_FASTA if needed).
    2. Process BLAST results: extract lineage from the reference headers,
       apply the per-rank percent-identity thresholds (species/genus/family/
       order/class), resolve ambiguous hits with the selected LCA algorithm,
       and left-merge taxonomy onto the ASV abundance table.
    3. Finalize: write the taxonomy+counts CSV and print a resolution summary.

    BLAST aligns each query sequence (an ASV/OTU representative) against a reference
    database; LCA (lowest common ancestor) resolution then assigns a sequence to the most
    specific rank its top hits agree on, so an ambiguous match is reported at, say, genus
    rather than guessing a single species. The per-rank percent-identity thresholds gate
    how high a sequence's similarity to references must be before that rank is accepted.

    Args:
        ctx: The Click context (carries the global verbose flag, read for traceback depth).
        query_fasta: Path to the query sequences (ASV representatives from DADA2). Must
            exist.
        ref_fasta: Path to the reference database FASTA. Must exist. The BLAST DB is built
            from it if not already present.
        asv_count: Path to the ASV count table CSV (e.g. ``seqtab_clean.csv`` from DADA2).
            Must exist. Taxonomy is left-merged onto these abundances.
        output: Path for the taxonomy+counts CSV. If ``None``, defaults to the query
            FASTA's directory with a ``_blast_taxonomy`` suffix.
        perc_identity: Minimum percent identity for a BLAST hit to be kept (blastn search
            filter). Percent.
        qcov_hsp_perc: Minimum query coverage per HSP for a hit to be kept. Percent.
        evalue: Maximum e-value (expected number of chance hits) for a BLAST hit.
        threshold_species: Minimum percent identity to assign at species rank.
        threshold_genus: Minimum percent identity to assign at genus rank.
        threshold_family: Minimum percent identity to assign at family rank.
        threshold_order: Minimum percent identity to assign at order rank.
        threshold_class: Minimum percent identity to assign at class rank.
        top_bitscore_pct: Cascade LCA: include hits whose bitscore is within this percent
            of the best hit's bitscore (MEGAN-LR-style band).
        lca_pident_delta: Cascade LCA: in-band hits must be within this percent identity of
            the best in-band hit to be considered.
        task: The blastn task/algorithm: one of ``megablast``, ``blastn``,
            ``dc-megablast``, or ``blastn-short``.
        lca_algorithm: Which LCA algorithm to use: ``cascade`` (per-rank thresholds) or
            ``collapsed_taxonomy`` (eDNAFlow/OceanOmics percent-identity-window collapse).
        lca_pid: collapsed_taxonomy only: hard percent-identity floor below which hits are
            discarded.
        lca_diff: collapsed_taxonomy only: width (percent identity) of the window of hits
            collapsed down to their LCA.

    Returns:
        None. Writes the taxonomy CSV and prints per-rank resolution percentages on
        success.

    Raises:
        SystemExit: Code 1 if a required file is missing or if the BLAST search or
            assignment fails.
    """
    from seednap.steps.taxonomic_assignment import BlastRunner, BlastTaxonomicAssigner

    _add_command_log_file(output.parent if output is not None else query_fasta.parent, "blast")

    console.print("\n[bold]Running BLAST taxonomic assignment[/bold]")
    console.print(f"Query: {query_fasta}")
    console.print(f"Reference: {ref_fasta}")
    console.print(f"ASV counts: {asv_count}\n")

    try:
        # Determine output path
        if output is None:
            output = query_fasta.parent / f"{query_fasta.stem}_blast_taxonomy.csv"

        # Create a unique temporary directory for BLAST output. Using
        # tempfile.mkdtemp (not a fixed "blast_temp" name with exist_ok=True)
        # guarantees this command never adopts -- and then rmtree's -- a
        # pre-existing user directory that happens to be named blast_temp.
        import tempfile

        blast_output_dir = Path(
            tempfile.mkdtemp(prefix="blast_temp_", dir=query_fasta.parent)
        )

        # Run BLAST search
        console.print("[cyan]Step 1/3:[/cyan] Running BLAST search...")
        runner = BlastRunner(
            perc_identity=perc_identity, qcov_hsp_perc=qcov_hsp_perc, evalue=evalue, task=task
        )

        blast_tsv = runner.run_blast_pipeline(
            query_fasta=query_fasta,
            db_fasta=ref_fasta,
            output_dir=blast_output_dir,
            marker="temp",
        )

        # Run taxonomic assignment
        console.print("[cyan]Step 2/3:[/cyan] Processing BLAST results...")
        assigner = BlastTaxonomicAssigner(
            reference_fasta=ref_fasta,
            threshold_species=threshold_species,
            threshold_genus=threshold_genus,
            threshold_family=threshold_family,
            threshold_order=threshold_order,
            threshold_class=threshold_class,
            top_bitscore_pct=top_bitscore_pct,
            lca_pident_delta=lca_pident_delta,
            lca_algorithm=lca_algorithm,
            lca_pid=lca_pid,
            lca_diff=lca_diff,
        )

        result = assigner.assign_taxonomy(
            blast_tsv=blast_tsv, asv_count_csv=asv_count, asv_fasta=query_fasta, output_path=output
        )

        # Print summary
        console.print("[cyan]Step 3/3:[/cyan] Finalizing results...")
        print_success("BLAST taxonomic assignment completed!")
        console.print(f"\nOutput file: [cyan]{output}[/cyan]")
        console.print(f"Total ASVs/OTUs: [green]{len(result)}[/green]")

        # Show taxonomic resolution summary, excluding unassigned ranks from the assigned
        # count (otherwise the summary would report 100% at every rank). No-hit OTUs get the
        # literal "Unassigned" (see BlastOutputFormatter), but the LCA cascade also NULLs
        # ranks below the resolved level; .astype(str) renders those None/NaN cells as
        # "None"/"nan", so the empty/NA-like strings here are real producers, not dead guards.
        unassigned = {"Unassigned", "nan", "", "NA", "None"}
        taxonomic_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        console.print("\n[bold]Taxonomic resolution:[/bold]")
        for rank in taxonomic_ranks:
            if rank in result.columns:
                n_assigned = int((~result[rank].astype(str).isin(unassigned)).sum())
                pct = (n_assigned / len(result) * 100) if len(result) else 0.0
                console.print(f"  {rank.capitalize()}: {n_assigned} ({pct:.1f}%)")

        # Clean up temporary directory
        import shutil

        shutil.rmtree(blast_output_dir, ignore_errors=True)

    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except Exception as e:
        print_error(f"BLAST assignment failed: {e}")
        if ctx.obj.get("verbose"):
            import traceback

            console.print(traceback.format_exc())
        sys.exit(1)


@main.command()
@click.argument("input_dir", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--forward-primer",
    required=True,
    help="Forward primer sequence (5' to 3')",
)
@click.option(
    "--reverse-primer",
    required=True,
    help="Reverse primer sequence (5' to 3')",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Output directory for trimmed reads",
)
@click.option(
    "--cores",
    "-c",
    default=1,
    type=int,
    help="Number of CPU cores to use",
)
def trim(
    input_dir: Path, forward_primer: str, reverse_primer: str, output_dir: Path, cores: int
) -> None:
    """
    Trim primers from FASTQ files (two-pass cutadapt).

    INPUT_DIR: Directory containing raw FASTQ files (R1/R2 pairs).

    Performs two-pass primer trimming:
    1. Remove 5' primers (anchored search)
    2. Remove 3' primers from pass 1 output

    Primers are the short fixed oligonucleotides that flank the amplified marker region;
    they are PCR artefacts, not biological sequence, so they must be stripped before
    denoising or clustering. This standalone command runs the same two-pass cutadapt
    trimming the pipeline uses, over every R1/R2 read pair in a directory.

    Args:
        input_dir: Directory of raw FASTQ files as R1/R2 pairs. Must exist.
        forward_primer: Forward primer sequence, 5' to 3'.
        reverse_primer: Reverse primer sequence, 5' to 3'.
        output_dir: Directory to write trimmed reads to.
        cores: Number of CPU cores cutadapt may use.

    Returns:
        None. Writes trimmed FASTQs and prints the number of samples processed on success.

    Raises:
        SystemExit: Code 1 on any trimming failure.
    """
    from seednap.steps.trimming import StandardTrimmer

    _add_command_log_file(output_dir, "trim")

    console.print(f"\n[bold]Trimming primers from:[/bold] {input_dir}")
    console.print(f"Forward primer: {forward_primer}")
    console.print(f"Reverse primer: {reverse_primer}")
    console.print(f"Output directory: {output_dir}")
    console.print(f"Cores: {cores}\n")

    try:
        trimmer = StandardTrimmer(cores=cores)

        results = trimmer.trim_directory(
            raw_reads_dir=input_dir,
            output_dir=output_dir,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            keep_untrimmed=False,
        )

        print_success(f"\nCompleted trimming {len(results)} samples!")
        console.print(f"Trimmed reads saved to: {output_dir}\n")

    except Exception as e:
        print_error(f"Trimming failed: {str(e)}")
        sys.exit(1)


@main.command()
@click.argument("raw_reads_dir", type=click.Path(exists=True, path_type=Path))
@click.argument("library_name", type=str)
@click.argument("metadata_csv", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--forward-primer",
    "-f",
    required=True,
    type=str,
    help="Forward primer sequence",
)
@click.option(
    "--reverse-primer",
    "-r",
    required=True,
    type=str,
    help="Reverse primer sequence",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Output base directory",
)
@click.option(
    "--cores",
    "-c",
    default=1,
    type=int,
    help="Number of CPU cores to use",
)
@click.option(
    "--no-gunzip",
    is_flag=True,
    help="Keep output files gzipped (default: gunzip)",
)
def demultiplex(
    raw_reads_dir: Path,
    library_name: str,
    metadata_csv: Path,
    forward_primer: str,
    reverse_primer: str,
    output_dir: Path,
    cores: int,
    no_gunzip: bool,
) -> None:
    """
    Demultiplex and trim ligation-based libraries.

    RAW_READS_DIR: Directory containing raw library FASTQ files.
    LIBRARY_NAME: Library identifier (matches filename prefix).
    METADATA_CSV: Metadata CSV with eventID, tag_demultiplex, and library columns.

    This command performs the complete ligation-based workflow:
    1. Generate tag files from metadata
    2. Demultiplex reads by tags
    3. Detect primers (both orientations)
    4. Merge and realign reads

    In a ligation-based library, many samples are pooled into one sequencing run, each
    sample marked by a short tag (barcode) ligated to its reads. Demultiplexing splits the
    pooled FASTQs back into per-sample reads by those tags before primer trimming, so each
    sample can be analysed separately downstream.

    Args:
        raw_reads_dir: Directory containing the pooled raw library FASTQ files. Must exist.
        library_name: Library identifier; matches the FASTQ filename prefix for this
            library.
        metadata_csv: Metadata CSV with eventID, tag_demultiplex and library columns,
            mapping tags to samples. Must exist.
        forward_primer: Forward primer sequence.
        reverse_primer: Reverse primer sequence.
        output_dir: Base output directory for the demultiplexed/realigned reads.
        cores: Number of CPU cores to use.
        no_gunzip: If True (``--no-gunzip``), leave the output files gzipped; otherwise the
            outputs are gunzipped.

    Returns:
        None. Writes the realigned per-sample reads and prints their directory on success.

    Raises:
        SystemExit: Code 1 on any demultiplexing failure.
    """
    from seednap.steps.trimming import LigationTrimmer

    _add_command_log_file(output_dir, f"demultiplex_{library_name}")

    console.print(f"\n[bold]Demultiplexing ligation library:[/bold] {library_name}")
    console.print(f"Raw reads: {raw_reads_dir}")
    console.print(f"Metadata: {metadata_csv}")
    console.print(f"Forward primer: {forward_primer}")
    console.print(f"Reverse primer: {reverse_primer}")
    console.print(f"Output directory: {output_dir}")
    console.print(f"Cores: {cores}\n")

    try:
        trimmer = LigationTrimmer(cores=cores)

        realigned_dir = trimmer.process_library(
            raw_reads_dir=raw_reads_dir,
            library_name=library_name,
            metadata_csv=metadata_csv,
            output_base_dir=output_dir,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            gunzip_output=not no_gunzip,
        )

        print_success("\nCompleted ligation library processing!")
        console.print(f"Realigned reads saved to: {realigned_dir}\n")

    except Exception as e:
        print_error(f"Demultiplexing failed: {str(e)}")
        sys.exit(1)


@main.command()
@click.argument("marker", type=str)
@click.argument("trimmed_reads_dir", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("outputs"),
    help="Base output directory (default: outputs/)",
)
@click.option(
    "--max-ee",
    default=2.0,
    type=float,
    help="Maximum expected errors for filtering (default: 2.0)",
)
@click.option(
    "--trunc-q",
    default=11,
    type=int,
    help="Truncate reads at first base with quality below this (default: 11)",
)
@click.option(
    "--min-overlap",
    default=20,
    type=int,
    help="Minimum overlap for merging paired reads (default: 20)",
)
@click.option(
    "--assign-taxonomy",
    is_flag=True,
    help="Run taxonomic assignment with DADA2 (requires --rdp-db and --species-db)",
)
@click.option(
    "--rdp-db",
    type=click.Path(exists=True, path_type=Path),
    help="Path to RDP-formatted taxonomy database (genus-level)",
)
@click.option(
    "--species-db",
    type=click.Path(exists=True, path_type=Path),
    help="Path to species-level taxonomy database",
)
@click.option(
    "--library-map",
    type=click.Path(exists=True, path_type=Path),
    help="CSV with 'sample,library' columns: learn DADA2 errors per library, then merge "
    "(DADA2-by-library). With 2+ libraries it denoises each separately and collapses identical "
    "ASVs; omit it (or a single library) for the standard single-batch path.",
)
def dada2(
    marker: str,
    trimmed_reads_dir: Path,
    output_dir: Path,
    max_ee: float,
    trunc_q: int,
    min_overlap: int,
    assign_taxonomy: bool,
    rdp_db: Optional[Path],
    species_db: Optional[Path],
    library_map: Optional[Path],
) -> None:
    """
    Run DADA2 processing on trimmed reads.

    MARKER: Marker name (e.g., teleo, amph).
    TRIMMED_READS_DIR: Directory containing primer-trimmed FASTQ files.

    This command performs:
    1. Quality control (before/after filtering)
    2. Filter and trim
    3. Error learning
    4. Denoising and sample inference
    5. Merge paired-end reads
    6. Chimera removal
    7. ASV table generation
    8. Metrics collection and reporting

    Optional: Taxonomic assignment using DADA2's naive Bayesian classifier

    DADA2 turns trimmed reads into ASVs (amplicon sequence variants): exact, denoised
    sequences inferred by modelling per-run sequencing error, as opposed to the fixed-%
    clusters that OTU methods produce. Along the way it filters by expected error,
    merges the paired forward/reverse reads into full amplicons, and removes chimeras
    (artefactual sequences fused from two real templates during PCR).

    Args:
        marker: Marker name (e.g. ``teleo``, ``amph``). Names the output subtree.
        trimmed_reads_dir: Directory of primer-trimmed FASTQ files (R1/R2 pairs). Must
            exist.
        output_dir: Base output directory. Defaults to ``outputs``.
        max_ee: Maximum expected errors allowed per read during filtering; reads above
            this are discarded.
        trunc_q: Truncate each read at the first base whose quality score falls below this
            value.
        min_overlap: Minimum number of overlapping bases required to merge a forward and
            reverse read.
        assign_taxonomy: If True (``--assign-taxonomy``), run DADA2's naive Bayesian
            taxonomy after denoising; requires ``rdp_db`` and ``species_db``.
        rdp_db: Path to the RDP-formatted (genus-level) taxonomy database. Required when
            ``assign_taxonomy`` is set. Must exist if given.
        species_db: Path to the species-level taxonomy database. Required when
            ``assign_taxonomy`` is set. Must exist if given.
        library_map: Optional CSV with ``sample,library`` columns. When 2+ libraries are
            present, errors are learned per library and identical ASVs are then collapsed
            (DADA2-by-library); omit it (or use a single library) for the standard
            single-batch path. Must exist if given.

    Returns:
        None. Writes the sequence table, query FASTA, ASV correspondence and metrics (and,
        if requested, the taxonomy table) and prints their paths on success.

    Raises:
        SystemExit: Code 1 if a required file is missing, if ``--assign-taxonomy`` is set
            without both DB paths, or on any other DADA2 failure.
    """
    from seednap.steps.dada2 import Dada2Processor

    _add_command_log_file(output_dir, marker)

    console.print("\n[bold]Running DADA2 processing:[/bold]")
    console.print(f"Marker: {marker}")
    console.print(f"Trimmed reads: {trimmed_reads_dir}")
    console.print(f"Output directory: {output_dir}")
    console.print(f"Parameters: maxEE={max_ee}, truncQ={trunc_q}, minOverlap={min_overlap}\n")

    try:
        # Initialize processor
        processor = Dada2Processor(
            marker=marker,
            trimmed_reads_dir=trimmed_reads_dir,
            output_base_dir=output_dir,
        )

        # Run DADA2 processing
        console.print("[bold]Running DADA2 workflow...[/bold]")
        outputs = processor.process(
            max_ee=max_ee,
            trunc_q=trunc_q,
            min_overlap=min_overlap,
            collect_metrics=True,
            library_map=library_map,
        )

        print_success("\nDADA2 processing completed successfully!")
        console.print("\nOutput files:")
        console.print(f"  Sequence table: {outputs['seqtab_clean']}")
        console.print(f"  Query FASTA: {outputs['query_fasta']}")
        console.print(f"  ASV correspondence: {outputs['corresp_seq']}")
        console.print(f"  Metrics: {outputs['metrics_dir']}")

        # Run taxonomy assignment if requested
        if assign_taxonomy:
            if not rdp_db or not species_db:
                print_error("--rdp-db and --species-db are required for taxonomy assignment")
                sys.exit(1)

            from seednap.steps.taxonomic_assignment.dada2_taxonomy_runner import Dada2TaxonomyRunner

            console.print("\n[bold]Running taxonomic assignment...[/bold]")
            taxonomy_runner = Dada2TaxonomyRunner()
            taxo_outputs = taxonomy_runner.run_dada2_taxonomy(
                marker=marker,
                output_dir=output_dir,
                rdp_db_path=rdp_db,
                species_db_path=species_db,
                query_fasta=outputs["query_fasta"],
                # "02_dada2" is the DADA2 step's output subdir convention (see the
                # orchestrator/processor and the `report` command); keep this literal in
                # sync with that layout if the output tree ever changes.
                log_file=output_dir / "02_dada2" / marker / "dada2_taxonomy.log",
            )

            print_success("\nTaxonomic assignment completed!")
            console.print(f"  Taxonomy table: {taxo_outputs['final_table']}")

        console.print()

    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except Exception as e:
        print_error(f"DADA2 processing failed: {e}")
        _maybe_traceback()
        sys.exit(1)


@main.command()
@click.argument("marker", type=str)
@click.argument("trimmed_reads_dir", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("outputs"),
    help="Base output directory (default: outputs/)",
)
@click.option(
    "--distance",
    "-d",
    default=1,
    type=int,
    help="SWARM clustering distance threshold (default: 1)",
)
@click.option(
    "--threads",
    "-t",
    default=4,
    type=int,
    help="Number of threads (default: 4)",
)
@click.option(
    "--no-fastidious",
    is_flag=True,
    help="Disable fastidious mode (singleton refinement)",
)
@click.option(
    "--no-chimera-filter",
    is_flag=True,
    help="Skip de novo chimera detection",
)
def swarm(
    marker: str,
    trimmed_reads_dir: Path,
    output_dir: Path,
    distance: int,
    threads: int,
    no_fastidious: bool,
    no_chimera_filter: bool,
) -> None:
    """
    Run SWARM OTU clustering on trimmed reads.

    MARKER: Marker name (e.g., teleo, amph).
    TRIMMED_READS_DIR: Directory containing primer-trimmed FASTQ files (R1/R2 pairs).

    This command performs:
    1. Merge paired-end reads (vsearch)
    2. Per-sample dereplication (vsearch)
    3. Global dereplication
    4. SWARM OTU clustering
    5. Sort representatives by abundance
    6. De novo chimera detection (vsearch UCHIME)
    7. OTU contingency table generation

    SWARM is an alternative to DADA2: instead of denoising to exact variants, it clusters
    near-identical sequences into OTUs (operational taxonomic units) by single-linkage
    growth at a small distance threshold ``d``, which avoids the arbitrary global %-identity
    cutoff of classic OTU pickers. Paired reads are merged first; the de novo chimera step
    (vsearch UCHIME) drops PCR-fused artefacts before the per-sample OTU count table is
    built. OTU counts legitimately differ from DADA2 ASV counts; that is not a bug.

    Args:
        marker: Marker name (e.g. ``teleo``, ``amph``). Names the output subtree.
        trimmed_reads_dir: Directory of primer-trimmed FASTQ files (R1/R2 pairs). Must
            exist.
        output_dir: Base output directory. Defaults to ``outputs``.
        distance: SWARM clustering distance threshold ``d`` (max nucleotide differences for
            single-linkage growth). Defaults to 1.
        threads: Number of threads to use. Defaults to 4.
        no_fastidious: If True (``--no-fastidious``), disable SWARM's fastidious singleton
            refinement; otherwise fastidious mode is on.
        no_chimera_filter: If True (``--no-chimera-filter``), skip de novo chimera
            detection; otherwise chimeras are removed.

    Returns:
        None. Writes the query FASTA, OTU tables and merged reads, and prints their paths
        on success.

    Raises:
        SystemExit: Code 1 if no R1/R2 FASTQ pairs are found in ``trimmed_reads_dir`` or on
            any other SWARM failure.
    """
    from seednap.steps.swarm import SwarmProcessor

    _add_command_log_file(output_dir, marker)

    console.print("\n[bold]Running SWARM OTU clustering:[/bold]")
    console.print(f"Marker: {marker}")
    console.print(f"Trimmed reads: {trimmed_reads_dir}")
    console.print(f"Output directory: {output_dir}")
    console.print(f"Parameters: d={distance}, threads={threads}, "
                  f"fastidious={not no_fastidious}, chimera={not no_chimera_filter}\n")

    try:
        processor = SwarmProcessor(
            marker=marker,
            trimmed_reads_dir=trimmed_reads_dir,
            output_base_dir=output_dir,
        )

        outputs = processor.process(
            d=distance,
            fastidious=not no_fastidious,
            threads=threads,
            chimera_detection=not no_chimera_filter,
        )

        print_success("\nSWARM OTU clustering completed successfully!")
        console.print("\nOutput files:")
        console.print(f"  Query FASTA: {outputs['query_fasta']}")
        console.print(f"  OTU table: {outputs['seqtab_clean_t']}")
        console.print(f"  Full OTU table: {outputs['otu_table_full']}")
        console.print(f"  Merged reads: {outputs['merged_dir']}")
        console.print()

    except FileNotFoundError as e:
        print_error(
            f"SWARM found no R1/R2 FASTQ pairs to process: {e}. Confirm that "
            f"{trimmed_reads_dir} is the trim step's output for marker '{marker}' (not the "
            f"raw/untrimmed reads), that the trim step completed, and that the files use a "
            f"recognized naming pattern: *_R1.fastq[.gz], *.R1.fastq[.gz], or "
            f"*_R1_001.fastq[.gz] with matching R2 files."
        )
        sys.exit(1)
    except Exception as e:
        print_error(f"SWARM processing failed: {e}")
        _maybe_traceback()
        sys.exit(1)


@main.command()
@click.argument("method", type=click.Choice(["blast", "dada2", "ecotag", "decipher"]))
@click.argument("marker", type=str)
@click.argument("query_fasta", type=click.Path(exists=True, path_type=Path))
@click.argument("asv_count_csv", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--config",
    "config_file",
    type=click.Path(exists=True, path_type=Path),
    help="Marker YAML config: use its taxonomy.databases.<method> block (db path + method "
    "params) so the result matches run-pipeline. Explicit options below override the config.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("outputs"),
    help="Base output directory (default: outputs/)",
)
@click.option(
    "--reference-fasta",
    type=click.Path(exists=True, path_type=Path),
    help="Reference database FASTA (for BLAST method)",
)
@click.option(
    "--rdp-db",
    type=click.Path(exists=True, path_type=Path),
    help="RDP-formatted taxonomy database (for DADA2 method)",
)
@click.option(
    "--species-db",
    type=click.Path(exists=True, path_type=Path),
    help="Species-level database (for DADA2 method)",
)
@click.option(
    "--taxonomy-db",
    type=click.Path(exists=True, path_type=Path),
    help="NCBI taxonomy database (for ecotag method)",
)
@click.option(
    "--reference-db",
    type=click.Path(exists=True, path_type=Path),
    help="Reference sequence database (for ecotag method)",
)
@click.option(
    "--trained-classifier",
    type=click.Path(exists=True, path_type=Path),
    help="Trained DECIPHER classifier .rds file (for DECIPHER method)",
)
@click.option(
    "--threshold-species",
    type=float,
    default=99.0,
    help="Percent identity threshold for species (BLAST, default: 99.0)",
)
@click.option(
    "--threshold-genus",
    type=float,
    default=96.0,
    help="Percent identity threshold for genus (BLAST, default: 96.0)",
)
@click.option(
    "--threshold-family",
    type=float,
    default=90.0,
    help="Percent identity threshold for family (BLAST, default: 90.0)",
)
@click.option(
    "--threshold-order",
    type=float,
    default=80.0,
    help="Percent identity threshold for order (BLAST, default: 80.0)",
)
@click.option(
    "--threshold-class",
    type=float,
    default=70.0,
    help="Percent identity threshold for class (BLAST, default: 70.0)",
)
@click.option(
    "--top-bitscore-pct",
    type=float,
    default=10.0,
    help="cascade LCA bitscore band as % of best hit (BLAST; default: 10.0)",
)
@click.option(
    "--lca-pident-delta",
    type=float,
    default=1.0,
    help="cascade LCA: in-band hits within this %id of the best in-band hit (BLAST; default: 1.0)",
)
@click.option(
    "--lca-algorithm",
    type=click.Choice(["cascade", "collapsed_taxonomy"]),
    default="cascade",
    help="BLAST LCA algorithm: cascade (default) or collapsed_taxonomy (eDNAFlow/OceanOmics)",
)
@click.option(
    "--lca-pid",
    type=float,
    default=90.0,
    help="collapsed_taxonomy: hard %identity floor (BLAST; default: 90.0)",
)
@click.option(
    "--lca-diff",
    type=float,
    default=1.0,
    help="collapsed_taxonomy: identity-window width collapsed to the LCA (BLAST; default: 1.0)",
)
@click.option(
    "--confidence-threshold",
    type=int,
    default=60,
    help="Confidence threshold for DECIPHER (0-100, default: 60)",
)
@click.option(
    "--processors",
    "-c",
    type=int,
    default=8,
    help="Number of CPU cores (default: 8)",
)
@click.pass_context
def assign_taxonomy(
    ctx: click.Context,
    method: str,
    marker: str,
    query_fasta: Path,
    asv_count_csv: Path,
    config_file: Optional[Path],
    output_dir: Path,
    reference_fasta: Optional[Path],
    rdp_db: Optional[Path],
    species_db: Optional[Path],
    taxonomy_db: Optional[Path],
    reference_db: Optional[Path],
    trained_classifier: Optional[Path],
    threshold_species: float,
    threshold_genus: float,
    threshold_family: float,
    threshold_order: float,
    threshold_class: float,
    top_bitscore_pct: float,
    lca_pident_delta: float,
    lca_algorithm: str,
    lca_pid: float,
    lca_diff: float,
    confidence_threshold: int,
    processors: int,
) -> None:
    """
    Assign taxonomy to ASVs using various methods.

    METHOD: Taxonomic assignment method (blast, dada2, ecotag, decipher).
    MARKER: Marker name (e.g., teleo, amph).
    QUERY_FASTA: Query FASTA file with ASV sequences.
    ASV_COUNT_CSV: ASV count table (seqtab_clean.csv or _t.csv).

    Each method requires specific database files:

    \b
    BLAST: --reference-fasta
    DADA2: --rdp-db and --species-db
    ecotag: --taxonomy-db and --reference-db
    DECIPHER: --trained-classifier

    Pass --config <marker.yaml> to use that marker's taxonomy.databases.<method>
    block (database path plus method parameters, e.g. BLAST evalue/task/thresholds)
    so the standalone result matches run-pipeline. Any explicit option above
    overrides the corresponding config value.

    Taxonomic assignment is the step that names each ASV/OTU sequence by comparing it to a
    reference database. The four methods differ in approach: BLAST (alignment + LCA),
    DADA2 (naive Bayesian classifier), ecotag (OBITools tree-based assignment), and
    DECIPHER (IdTaxa classifier). This command is the standalone counterpart to the
    taxonomy step of ``run-pipeline``.

    Args:
        ctx: The Click context, used both for the global verbose flag and to detect which
            options the user actually typed (so config values are only overridden by
            explicit CLI flags, not by click defaults).
        method: Assignment method: one of ``"blast"``, ``"dada2"``, ``"ecotag"``,
            ``"decipher"``.
        marker: Marker name (e.g. ``teleo``, ``amph``). Names the output subtree.
        query_fasta: Query FASTA of ASV/OTU sequences to assign. Must exist.
        asv_count_csv: ASV/OTU count table (``seqtab_clean.csv`` or ``_t.csv``). Must
            exist.
        config_file: Optional marker YAML; if given, its ``taxonomy.databases.<method>``
            block supplies DB paths and method parameters. Its ``taxonomy.method`` must
            equal ``method``. Must exist if given.
        output_dir: Base output directory. Defaults to ``outputs``.
        reference_fasta: Reference database FASTA (BLAST method). Must exist if given.
        rdp_db: RDP-formatted taxonomy database (DADA2 method). Must exist if given.
        species_db: Species-level taxonomy database (DADA2 method). Must exist if given.
        taxonomy_db: NCBI taxonomy database (ecotag method). Must exist if given.
        reference_db: Reference sequence database (ecotag method). Must exist if given.
        trained_classifier: Trained DECIPHER classifier ``.rds`` file (DECIPHER method).
            Must exist if given.
        threshold_species: BLAST: minimum percent identity to assign at species rank.
        threshold_genus: BLAST: minimum percent identity to assign at genus rank.
        threshold_family: BLAST: minimum percent identity to assign at family rank.
        threshold_order: BLAST: minimum percent identity to assign at order rank.
        threshold_class: BLAST: minimum percent identity to assign at class rank.
        top_bitscore_pct: BLAST cascade LCA: bitscore band as a percent of the best hit.
        lca_pident_delta: BLAST cascade LCA: in-band hits must be within this percent
            identity of the best in-band hit.
        lca_algorithm: BLAST LCA algorithm: ``cascade`` or ``collapsed_taxonomy``.
        lca_pid: BLAST collapsed_taxonomy: hard percent-identity floor.
        lca_diff: BLAST collapsed_taxonomy: identity-window width collapsed to the LCA.
        confidence_threshold: DECIPHER confidence threshold (0-100) for accepting a rank.
        processors: Number of CPU cores (used by DECIPHER).

    Returns:
        None. Writes the taxonomy outputs for the chosen method and prints their paths on
        success.

    Raises:
        SystemExit: Code 1 if ``--config`` is for a different method than ``method``, if a
            method's required DB option is neither passed nor present in the config, if an
            input file is missing, on a value error, or on any other assignment failure.
    """
    from seednap.steps.taxonomic_assignment import TaxonomicAssigner

    _add_command_log_file(output_dir, marker)

    console.print("\n[bold]Taxonomic Assignment:[/bold]")
    console.print(f"Method: {method}")
    console.print(f"Marker: {marker}")
    console.print(f"Query: {query_fasta}")
    console.print(f"ASV counts: {asv_count_csv}")
    console.print(f"Output directory: {output_dir}\n")

    try:
        # When --config is given, derive method parameters from that marker's
        # taxonomy.databases.<method> block (db path + blast evalue/task/thresholds/
        # lca params, etc.) so this standalone command matches run-pipeline. Without it,
        # the explicit options / model defaults are used as before.
        #
        # The config's taxonomy.method must equal the positional METHOD; using a config
        # written for one method while asking for another would silently apply the wrong
        # parameter set, so we raise instead of guessing (no-silent-fallback rule).
        config_kwargs: Dict[str, Any] = {}
        if config_file is not None:
            cfg = load_config(config_file)
            if cfg.taxonomy.method != method:
                print_error(
                    f"--config {config_file} is for taxonomy method "
                    f"'{cfg.taxonomy.method}', but you asked for '{method}'. "
                    f"Run `assign-taxonomy {cfg.taxonomy.method} ...` or pass a config "
                    f"whose taxonomy.method is '{method}'."
                )
                sys.exit(1)
            config_kwargs = _assign_kwargs_from_config(cfg, method)

        # Initialize assigner
        assigner = TaxonomicAssigner(
            method=method,
            marker=marker,
            output_dir=output_dir,
        )

        # Prepare method-specific arguments. Start from the config-derived values (if any)
        # then let any EXPLICITLY-passed CLI option override its config counterpart;
        # _given() distinguishes a user-typed flag from a value left at its click default.
        def _given(param_name: str) -> bool:
            """Return whether a CLI option was explicitly typed by the user.

            Distinguishes a user-typed flag from a value left at its click default, so a
            config-derived value is overridden only by an explicit CLI flag.

            Args:
                param_name: The Click parameter name (the function argument name) to check.

            Returns:
                True if the option's source is the command line; False if it came from a
                default, environment, or other source.
            """
            from click.core import ParameterSource

            return ctx.get_parameter_source(param_name) == ParameterSource.COMMANDLINE

        kwargs: Dict[str, Any] = dict(config_kwargs)

        if method == "blast":
            if _given("reference_fasta") or "reference_fasta" not in kwargs:
                if not reference_fasta:
                    print_error(
                        "--reference-fasta is required for BLAST method "
                        "(or pass --config with a taxonomy.databases.blast.fasta path)"
                    )
                    sys.exit(1)
                kwargs["reference_fasta"] = reference_fasta
            for opt in (
                "threshold_species", "threshold_genus", "threshold_family",
                "threshold_order", "threshold_class", "top_bitscore_pct",
                "lca_pident_delta", "lca_algorithm", "lca_pid", "lca_diff",
            ):
                if _given(opt) or opt not in kwargs:
                    kwargs[opt] = locals()[opt]

        elif method == "dada2":
            if _given("rdp_db") or "rdp_db_path" not in kwargs:
                if not rdp_db:
                    print_error(
                        "--rdp-db is required for DADA2 method "
                        "(or pass --config with a taxonomy.databases.dada2.all path)"
                    )
                    sys.exit(1)
                kwargs["rdp_db_path"] = rdp_db
            if _given("species_db") or "species_db_path" not in kwargs:
                if not species_db:
                    print_error(
                        "--species-db is required for DADA2 method "
                        "(or pass --config with a taxonomy.databases.dada2.species path)"
                    )
                    sys.exit(1)
                kwargs["species_db_path"] = species_db

        elif method == "ecotag":
            if _given("taxonomy_db") or "taxonomy_db" not in kwargs:
                if not taxonomy_db:
                    print_error(
                        "--taxonomy-db is required for ecotag method "
                        "(or pass --config with a taxonomy.databases.ecotag.tree path)"
                    )
                    sys.exit(1)
                kwargs["taxonomy_db"] = taxonomy_db
            if _given("reference_db") or "reference_db" not in kwargs:
                if not reference_db:
                    print_error(
                        "--reference-db is required for ecotag method "
                        "(or pass --config with a taxonomy.databases.ecotag.fasta path)"
                    )
                    sys.exit(1)
                kwargs["reference_db"] = reference_db

        elif method == "decipher":
            if _given("trained_classifier") or "trained_classifier_path" not in kwargs:
                if not trained_classifier:
                    print_error(
                        "--trained-classifier is required for DECIPHER method "
                        "(or pass --config with a taxonomy.databases.decipher.trained path)"
                    )
                    sys.exit(1)
                kwargs["trained_classifier_path"] = trained_classifier
            if _given("confidence_threshold") or "threshold" not in kwargs:
                kwargs["threshold"] = confidence_threshold
            if _given("processors") or "processors" not in kwargs:
                kwargs["processors"] = processors

        # Run taxonomic assignment
        console.print(f"[bold]Running {method.upper()} taxonomic assignment...[/bold]")
        outputs = assigner.assign_taxonomy(
            query_fasta=query_fasta,
            asv_count_csv=asv_count_csv,
            **kwargs,
        )

        print_success("\nTaxonomic assignment completed!")
        console.print("\nOutput files:")
        for key, path in outputs.items():
            if path and Path(path).exists():
                console.print(f"  {key}: {path}")

        console.print()

    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except ValueError as e:
        print_error(str(e))
        sys.exit(1)
    except Exception as e:
        from seednap.steps.taxonomic_assignment.blast_runner import BlastError
        from seednap.steps.taxonomic_assignment.dada2_taxonomy_runner import Dada2TaxonomyError
        from seednap.steps.taxonomic_assignment.decipher_runner import DecipherError
        from seednap.steps.taxonomic_assignment.ecotag_runner import EcotagError

        if isinstance(e, (EcotagError, BlastError, DecipherError, Dada2TaxonomyError)):
            # These carry self-contained, actionable messages (e.g. the OBITools-missing
            # what/why/fix block); print verbatim rather than mislabeling them as a crash.
            print_error(str(e))
        else:
            print_error(f"Taxonomic assignment failed unexpectedly: {e}")
        _maybe_traceback()
        sys.exit(1)


@main.command()
@click.argument("config", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--resume",
    is_flag=True,
    help="Resume pipeline from previous run",
)
@click.option(
    "--state-file",
    type=click.Path(path_type=Path),
    help="Path to state file (default: auto-generated in output dir)",
)
@click.option(
    "--stop-on-error/--continue-on-error",
    default=True,
    help="Whether to stop on first error or continue (default: stop)",
)
def run_pipeline(
    config: Path,
    resume: bool,
    state_file: Optional[Path],
    stop_on_error: bool,
) -> None:
    """
    Run complete SeeDNAP eDNA metabarcoding pipeline.

    This command orchestrates the full pipeline from raw reads to taxonomic assignments:
    1. Demultiplexing (optional)
    2. Primer trimming with cutadapt
    3. DADA2 processing (filtering, denoising, merging, chimera removal)
    4. Taxonomic assignment (DADA2/BLAST/ecotag/DECIPHER)
    5. Export to GBIF format

    CONFIG: Path to pipeline configuration YAML file

    Examples:

        # Run complete pipeline
        seednap run-pipeline config/markers/teleo.yaml

        # Resume from previous run
        seednap run-pipeline config/markers/teleo.yaml --resume

        # Continue on errors
        seednap run-pipeline config/markers/teleo.yaml --continue-on-error

        # Use custom state file
        seednap run-pipeline config/markers/teleo.yaml --state-file my_state.json

    This is the main entrypoint: it loads the marker config, runs a read-only preflight to
    fail fast on missing inputs, then hands control to the orchestrator, which executes each
    step and records its status, duration and outputs in a state JSON so the run can be
    resumed. A summary table of per-step status is printed at the end.

    Args:
        config: Path to the pipeline configuration YAML. Must exist (enforced by Click).
        resume: If True (``--resume``), continue from a previous run's state file, skipping
            already-completed steps.
        state_file: Path to the run-state JSON. If ``None``, the orchestrator
            auto-generates one in the output directory.
        stop_on_error: If True (``--stop-on-error``, the default), abort on the first failed
            step; if False (``--continue-on-error``), keep running subsequent steps.

    Returns:
        None. Prints a completion summary on success.

    Raises:
        SystemExit: Code 1 if preflight finds unusable inputs, on a config error, on a
            pipeline-step value error, on a missing required file, or on any other
            unexpected failure.
    """
    from seednap.pipeline.orchestrator import PipelineOrchestrator

    try:
        console.print("\n[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]")
        console.print("[bold cyan]     SeeDNAP eDNA Metabarcoding Pipeline[/bold cyan]")
        console.print("[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]\n")

        # Load config to show marker info
        from seednap.config.loader import load_config

        config_obj = load_config(config)

        # Preflight: fail before trimming/clustering if referenced inputs are missing, rather
        # than wasting compute and failing mid-run.
        from seednap.errors import preflight_checks

        problems = preflight_checks(config_obj)
        if problems:
            print_error("Cannot start the pipeline -- referenced inputs are not usable:\n")
            for problem in problems:
                console.print(problem.render())
                console.print()
            sys.exit(1)

        console.print(f"[bold]Marker:[/bold] {config_obj.marker.name}")
        console.print(f"[bold]Description:[/bold] {config_obj.marker.description}")
        console.print(f"[bold]Taxonomy method:[/bold] {config_obj.taxonomy.method}")

        if resume:
            console.print("\n[yellow]Resuming from previous run[/yellow]")

        console.print()

        # Create and run orchestrator
        orchestrator = PipelineOrchestrator(
            config=config, state_file=state_file, resume=resume
        )

        final_state = orchestrator.run(stop_on_error=stop_on_error)

        # Print summary
        summary = final_state.get_summary()

        console.print("\n[bold green]✓ Pipeline Completed Successfully![/bold green]\n")
        console.print("[bold]Summary:[/bold]")
        console.print(f"  Total duration: {summary['total_duration_seconds']:.1f}s")
        console.print(f"  Completed steps: {summary['completed']}/{summary['total_steps']}")
        console.print(f"  Failed steps: {summary['failed']}")
        console.print(f"  Skipped steps: {summary['skipped']}")

        if summary["completed"] > 0:
            console.print("\n[bold]Completed steps:[/bold]")
            for step_name, step_info in summary["steps"].items():
                if step_info["status"] == "completed":
                    duration = step_info["duration_seconds"]
                    console.print(
                        f"  [green]✓[/green] {step_name}: {duration:.1f}s"
                        if duration
                        else f"  [green]✓[/green] {step_name}"
                    )

        if summary["failed"] > 0:
            console.print("\n[bold yellow]Failed steps:[/bold yellow]")
            for step_name, step_info in summary["steps"].items():
                if step_info["status"] == "failed":
                    error = step_info.get("error", "Unknown error")
                    console.print(f"  [red]✗[/red] {step_name}: {error}")

        console.print(f"\n[bold]Output directory:[/bold] {config_obj.paths.output}")
        console.print(f"[bold]Log directory:[/bold] {config_obj.paths.logs}")
        console.print()

    except ConfigError as e:
        # Already a self-contained, actionable config message; print verbatim, not as a crash.
        print_error(str(e))
        sys.exit(1)
    except ValueError as e:
        # Config errors are caught above (ConfigError), so a ValueError here comes from a
        # pipeline step; config_obj is bound by now. Point the user at the run's state file
        # (default <paths.output>/.<marker>_state.json) and log directory.
        state_loc = config_obj.paths.output / f".{config_obj.marker.name}_state.json"
        print_error(
            f"Pipeline aborted: {e}\n"
            f"If a step failed, its status and error are recorded in the state file "
            f"({state_loc}) and the run log under {config_obj.paths.logs}; fix the cause and "
            f"re-run with --resume to continue from the failed step."
        )
        sys.exit(1)
    except FileNotFoundError as e:
        print_error(
            f"Pipeline could not find a required file: {e}. A path configured in your YAML "
            f"(raw_data, or a taxonomy database under taxonomy.databases.<method>) points at "
            f"something that is not on disk. Run `seednap validate {config}` to see each "
            f"configured path flagged found/MISSING, fix the offending entry, then re-run."
        )
        sys.exit(1)
    except Exception as e:
        print_error(f"Pipeline failed unexpectedly: {e}")
        _maybe_traceback()
        sys.exit(1)


@main.command()
@click.argument("marker", type=str)
@click.option(
    "--output-dir", "-o", type=click.Path(path_type=Path), default=Path("outputs"),
    help="Base output directory (default: outputs/)",
)
@click.option("--html", "html_report", is_flag=True, help="Also generate the self-contained HTML run report")
@click.option("--warn-retention", type=float, default=30.0, help="Warn below this overall retention % (default: 30)")
@click.option("--warn-step-loss", type=float, default=70.0, help="Warn when a step drops more than this % (default: 70)")
@click.option("--field-metadata", type=click.Path(exists=True, path_type=Path),
              help="Per-sample (field) metadata CSV for the Dataset/provenance section (location, dates, sites)")
@click.option("--project-metadata", type=click.Path(exists=True, path_type=Path),
              help="Project metadata CSV for the Dataset section (recorder, sequencing, reference DB)")
@click.option("--log-file", type=click.Path(exists=True, path_type=Path),
              help="Pipeline run log to embed in the HTML report (auto-located from logs/ if omitted)")
def report(
    marker: str, output_dir: Path, html_report: bool,
    warn_retention: float, warn_step_loss: float,
    field_metadata: Optional[Path], project_metadata: Optional[Path],
    log_file: Optional[Path],
) -> None:
    """
    Build the read/sequence tracking report from existing run outputs.

    MARKER: marker name (e.g. teleo, mam07).

    Rebuilds the per-step read-loss table (raw -> trimmed -> ... -> final) from
    the on-disk Cutadapt logs and the DADA2/SWARM outputs. Pass --html for the
    full self-contained visual report.

    Read tracking is the per-sample audit of how many reads survive each step (raw,
    trimmed, filtered, merged, non-chimeric, and so on). A sample that loses most of its
    reads at one step signals a quality or parameter problem, so this table is a key
    sanity check before trusting downstream abundances. The command only reads existing
    outputs; it never modifies the run and can be regenerated any time.

    Args:
        marker: Marker name (e.g. ``teleo``, ``mam07``). Selects the run's output subtree.
        output_dir: Base output directory of the run to report on. Defaults to ``outputs``.
        html_report: If True (``--html``), also write the self-contained HTML report.
        warn_retention: Warn when a sample's overall retention falls below this percent.
            Defaults to 30.
        warn_step_loss: Warn when any single step drops more than this percent of reads.
            Defaults to 70.
        field_metadata: Optional per-sample field metadata CSV for the HTML Dataset
            section. Auto-located near the output if omitted. Must exist if given.
        project_metadata: Optional project metadata CSV for the HTML Dataset section.
            Auto-located near the output if omitted. Must exist if given.
        log_file: Optional pipeline run log to embed in the HTML report. Auto-located from
            ``logs/`` if omitted. Must exist if given.

    Returns:
        None. Writes the read-tracking and step-summary CSVs (and, with ``--html``, the
        HTML report), prints the tables to the console, and emits [WARN]s for missing
        inputs or data-loss thresholds.

    Raises:
        SystemExit: Code 1 if no samples are found for the marker under ``output_dir``, or
            on any other reporting failure.
    """
    import pandas as pd

    from seednap.steps.report import HTMLReportBuilder, ReadTrackingBuilder

    out = output_dir
    _add_command_log_file(out, marker)
    dada2_dir = out / "02_dada2" / marker
    swarm_otu = out / "02_swarm" / marker / "otu_table.csv"

    kwargs: Dict[str, Any] = {
        # Cutadapt per-sample logs live under <output>/01_trim/<marker>/logs (written by
        # the trim step); read them there so raw/trimmed counts (and % retained) populate.
        "marker": marker, "logs_dir": out / "01_trim" / marker / "logs",
        "warn_below_retention_pct": warn_retention, "warn_step_loss_pct": warn_step_loss,
    }
    method = None
    if (dada2_dir / "track_reads.csv").exists():
        kwargs["dada2_dir"] = dada2_dir
        method = "DADA2"
    elif swarm_otu.exists():
        kwargs["swarm_otu_table"] = swarm_otu
        method = "SWARM"
    else:
        print_warning(
            f"No DADA2 track_reads.csv or SWARM otu_table.csv found under {out} "
            f"for marker '{marker}'; reporting raw/trimmed counts from logs only."
        )

    console.print(
        f"\n[bold]Read tracking report:[/bold] {marker}"
        + (f"  [cyan]({method})[/cyan]" if method else "")
    )
    try:
        builder = ReadTrackingBuilder(**kwargs)
        df = builder.build()
        if df.empty:
            print_error(
                f"No samples found for marker '{marker}' under {out}. This report is built "
                f"from the per-sample Cutadapt trim logs in {out}/01_trim/{marker}/logs/ "
                f"(<sample>_trim_pass1.txt), and none were found there. Most likely "
                f"--output-dir does not point at the directory used for `run-pipeline` "
                f"(default: outputs/), or trimming has not run yet. For the per-step ASV/OTU "
                f"columns the report also reads {out}/02_dada2/{marker}/track_reads.csv "
                f"(DADA2) or {out}/02_swarm/{marker}/otu_table.csv (SWARM); if those are "
                f"missing only raw/trimmed are reported. Check --output-dir, then that "
                f"'{marker}' matches the subdirectory name under 02_dada2/ or 02_swarm/."
            )
            sys.exit(1)

        report_dir = out / "04_report" / marker
        paths = builder.write(report_dir, df=df)
        step_summary_df = builder.step_summary(df)
        builder.write_step_summary(report_dir, summary_df=step_summary_df)
        warns = builder.warnings(df, log=False)  # keep the console clean; shown below + in HTML

        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("sample")
        for step in builder.steps:
            table.add_column(step, justify="right")
        table.add_column("% retained", justify="right")
        for _, row in df.iterrows():
            cells = [str(row["sample"])]
            for step in builder.steps:
                cells.append("NA" if pd.isna(row[step]) else f"{int(row[step]):,}")
            pr = row["pct_retained"]
            cells.append("NA" if pd.isna(pr) else f"{pr:.1f}%")
            table.add_row(*cells)
        console.print(table)

        # Run-level step summary: total reads + ASV/OTU count after each step.
        feat_label = "ASVs" if "nonchim" in builder.steps else "OTUs"
        ss_table = Table(title="Step summary (run totals)", show_header=True, header_style="bold cyan")
        ss_table.add_column("step")
        ss_table.add_column("total reads", justify="right")
        ss_table.add_column(feat_label, justify="right")
        for _, srow in step_summary_df.iterrows():
            tr, nf = srow["total_reads"], srow["n_features"]
            ss_table.add_row(
                str(srow["step"]),
                "NA" if pd.isna(tr) else f"{int(tr):,}",
                "-" if pd.isna(nf) else f"{int(nf):,}",
            )
        console.print(ss_table)

        print_success(f"Wrote {paths['read_tracking_csv']}")
        print_success(f"Wrote {report_dir / 'step_summary.csv'}")
        if warns:
            print_warning(f"{len(warns)} data-loss/measurement warning(s) — see the run log.")

        if html_report:
            import json as _json

            state = None
            state_file = out / f".{marker}_state.json"
            if state_file.exists():
                try:
                    state = _json.loads(state_file.read_text())
                except (ValueError, OSError):
                    print_warning(f"Could not read state file {state_file}; timeline omitted.")
            # Locate the final taxonomy table (for the taxonomy/contamination panels).
            taxo = None
            for suffix in ("blast", "dada2RDP", "dada2", "ecotag", "decipher"):
                cand = out / f"{marker}_{suffix}.csv"
                if cand.exists():
                    taxo = cand
                    break
            if taxo is None:
                globbed = sorted(out.glob(f"{marker}_*.csv"))
                taxo = globbed[0] if globbed else None
            otu_full = swarm_otu.parent / "otu_table_full.csv" if method == "SWARM" else None
            # Auto-locate dataset metadata near the output if not given explicitly.
            if field_metadata is None:
                for cand in (out / f"metadata_field_{marker}.csv", out / "metadata" / f"metadata_field_{marker}.csv"):
                    if cand.exists():
                        field_metadata = cand
                        break
            if project_metadata is None:
                for cand in (out / f"metadata_proj_{marker}.csv", out / "metadata" / f"metadata_proj_{marker}.csv"):
                    if cand.exists():
                        project_metadata = cand
                        break
            if field_metadata or project_metadata:
                console.print(f"  dataset metadata: field={field_metadata or '—'}, project={project_metadata or '—'}")
            else:
                print_warning("No dataset metadata found; the Dataset section will note it is absent. "
                              "Pass --field-metadata / --project-metadata to include provenance.")
            # Locate the pipeline run log to embed (the colorized transcript section).
            if log_file is None:
                for cand_dir in (Path("logs"), out / "logs", out.parent / "logs"):
                    exact = cand_dir / f"{marker}_pipeline_run.log"
                    if exact.exists():
                        log_file = exact
                        break
                    globbed = sorted(cand_dir.glob(f"{marker}_pipeline_*.log"),
                                     key=lambda p: p.stat().st_mtime, reverse=True)
                    if globbed:
                        log_file = globbed[0]
                        break
            if log_file:
                console.print(f"  run log: {log_file}")
            else:
                print_warning("No run log found under logs/; the HTML report's Run-log section "
                              "will note it is absent. Pass --log-file to embed a specific log.")
            html_path = HTMLReportBuilder(
                marker, df, warnings=warns, steps=builder.steps,
                state=state, taxonomy_csv=taxo,
                otu_table_full=otu_full if (otu_full and otu_full.exists()) else None,
                field_metadata_csv=field_metadata, project_metadata_csv=project_metadata,
                log_file=log_file,
                step_summary_df=step_summary_df,
                summary={
                    "warn_below_retention_pct": warn_retention,
                    "subtitle": f"{len(df)} samples · marker {marker}",
                    "provenance": {"dataset_name": marker, "marker": marker},
                },
            ).write(report_dir / "report.html")
            print_success(f"Wrote HTML report: {html_path}")
        console.print()

    except Exception as e:
        print_error(
            f"Report failed for marker '{marker}': {e}. This command only reads existing run "
            f"outputs under {out} and never modifies the pipeline, so it is safe to retry. "
            f"Check that --output-dir points at a completed run and that its inputs are intact: "
            f"the Cutadapt logs under {out}/logs/, the DADA2 "
            f"{out}/02_dada2/{marker}/track_reads.csv or SWARM "
            f"{out}/02_swarm/{marker}/otu_table.csv, and (with --html) the auto-located "
            f"taxonomy CSV {out}/{marker}_*.csv; the report writes to {out}/04_report/{marker}/, "
            f"which must be writable. Most failures here mean a prior step did not finish or its "
            f"output was moved or truncated; re-run the pipeline step that produces the missing "
            f"input, or confirm --output-dir points at a completed run."
        )
        sys.exit(1)


@main.command()
def version() -> None:
    """Show version information.

    Prints the installed seednap version string and the repository URL to the console.

    Returns:
        None.
    """
    console.print(f"\n[bold]seednap[/bold] version [cyan]{__version__}[/cyan]\n")
    console.print("eDNA metabarcoding pipeline with DADA2")
    console.print("Repository: https://github.com/WildinSync/wis_seednap\n")


@main.command()
@click.argument("code", required=False)
def explain(code: Optional[str]) -> None:
    """Explain a seednap error code in depth.

    CODE: a code shown in an error message, e.g. SDN-CFG-001. With no CODE, lists all codes.

    seednap tags many of its errors with a stable code so a user can look up a fuller
    what/why/fix explanation than fits in the original message. This command is that lookup.

    Args:
        code: An error code (e.g. ``SDN-CFG-001``) to explain in detail. If ``None`` or
            empty, the command lists every known code with its short title instead.

    Returns:
        None. Prints either the full explanation for ``code`` or the list of all codes.

    Raises:
        SystemExit: Code 1 if ``code`` is given but is not a known error code.
    """
    from seednap.errors import all_codes
    from seednap.errors.catalog import explain as explain_code

    if not code:
        console.print("\n[bold]seednap error codes[/bold] (run `seednap explain <CODE>`):\n")
        for c, title in all_codes().items():
            console.print(f"  [cyan]{c}[/cyan]  {title}")
        console.print()
        return

    detail = explain_code(code)
    if detail is None:
        print_error(
            f"Unknown error code '{code}'. Run `seednap explain` (no argument) to list all codes."
        )
        sys.exit(1)
    console.print(f"\n{detail}\n")


@main.command()
@click.argument("field_metadata", type=click.Path(exists=True, path_type=Path))
@click.option("--project-metadata", type=click.Path(exists=True, path_type=Path),
              help="Project metadata CSV (supplies the marker -> target_gene/assay_name)")
@click.option("--lab-metadata", type=click.Path(exists=True, path_type=Path),
              help="Legacy demux metadata CSV (supplies seq_run_id/library and tag barcodes)")
@click.option("--seq-run-id", type=str, default=None,
              help="Sequencing-run id for the whole dataset (overrides lab/derived value)")
@click.option("--target-gene", type=str, default=None,
              help="Marker / target_gene (overrides the project metadata)")
@click.option("--date-order", type=click.Choice(["ymd", "dmy", "mdy"], case_sensitive=False),
              default=None,
              help="Force the eventDate field order for genuinely-ambiguous dotted dates "
                   "(every date has day and month <=12); without it such files raise rather "
                   "than be guessed")
@click.option("-o", "--output", type=click.Path(path_type=Path), default=None,
              help="Write the canonical manifest CSV here")
@click.option("--abundance", type=click.Path(exists=True, path_type=Path), default=None,
              help="Validate the manifest's eventIDs against this abundance/OTU table")
@click.option("--strict", is_flag=True, default=False,
              help="Raise if the abundance table has sample columns absent from the manifest "
                   "(default: warn)")
@click.pass_context
def manifest(
    ctx: click.Context,
    field_metadata: Path,
    project_metadata: Optional[Path],
    lab_metadata: Optional[Path],
    seq_run_id: Optional[str],
    target_gene: Optional[str],
    date_order: Optional[str],
    output: Optional[Path],
    abundance: Optional[Path],
    strict: bool,
) -> None:
    """
    Build (and optionally validate) a canonical FAIRe sample manifest.

    Derives a single canonical, strictly-validated per-sample-library manifest from the
    lab's existing CSVs and, with --abundance, cross-checks its eventIDs against an
    abundance/OTU table (the up-front silent-ID-mismatch guard). Standalone: it reads
    on-disk inputs and never runs or alters the pipeline.

    \b
    FIELD_METADATA: per-sample field metadata CSV (metadata_field_*.csv), or a legacy
                    demux lab CSV (metadata_lab_*.csv) carrying library/tag columns.

    Every assumption (a synthesised seq_run_id, an ambiguous control, a dropped column,
    an orphan eventID) is logged as a [WARN]; ambiguous dates and a missing sample key
    raise. Pass -o to write the manifest CSV.

    A manifest is the single authoritative per-sample-library table that ties an eventID
    (one biological sample, or a negative control) to its sequencing run, tags, dates and
    location, in the FAIRe convention. With ``--abundance`` it doubles as an up-front guard
    against the silent ID-mismatch bug: it checks that every sample column in an
    abundance/OTU table has a matching manifest row before that table is trusted downstream.

    Args:
        ctx: The Click context (carries the global verbose flag, read for traceback depth).
        field_metadata: Per-sample field metadata CSV (``metadata_field_*.csv``), or a
            legacy demux lab CSV (``metadata_lab_*.csv``) carrying library/tag columns.
            Must exist.
        project_metadata: Optional project metadata CSV supplying the marker ->
            target_gene/assay_name mapping. Must exist if given.
        lab_metadata: Optional legacy demux metadata CSV supplying seq_run_id/library and
            tag barcodes. Must exist if given.
        seq_run_id: Optional sequencing-run id for the whole dataset; overrides any
            lab-derived value.
        target_gene: Optional marker / target_gene; overrides the project metadata.
        date_order: Optional forced eventDate field order (``ymd``, ``dmy`` or ``mdy``) for
            genuinely-ambiguous dotted dates; without it such files raise rather than be
            guessed.
        output: Optional path to write the canonical manifest CSV to.
        abundance: Optional abundance/OTU table to cross-check the manifest's eventIDs
            against. Must exist if given.
        strict: If True (``--strict``), raise when the abundance table has sample columns
            absent from the manifest; otherwise warn. Defaults to False.

    Returns:
        None. Prints row/sample/control/run counts, writes the manifest CSV if ``output``
        is given, and prints the abundance cross-check result if ``abundance`` is given.

    Raises:
        SystemExit: Code 1 if an input file is missing, on an ambiguous date or missing
            sample-key value error, on an abundance cross-check value error, or on any
            other build failure.
    """
    from seednap.config.manifest import validate_against_abundance
    from seednap.config.manifest_migrate import migrate_to_manifest

    _add_command_log_file(output.parent if output is not None else None, "manifest")

    console.print("\n[bold]Building sample manifest[/bold]")
    console.print(f"  Field metadata:   {field_metadata}")
    if project_metadata:
        console.print(f"  Project metadata: {project_metadata}")
    if lab_metadata:
        console.print(f"  Lab metadata:     {lab_metadata}")
    console.print()

    try:
        m = migrate_to_manifest(
            field_csv=field_metadata,
            project_csv=project_metadata,
            lab_csv=lab_metadata,
            seq_run_id=seq_run_id,
            target_gene=target_gene,
            date_order=date_order,
        )
    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except ValueError as e:
        print_error(
            f"Manifest build failed: {e}. If this names an ambiguous eventDate (day and month "
            f"both <=12), re-run with --date-order ymd|dmy|mdy to force the field order instead "
            f"of editing the CSV. If it reports no eventID/samp_name column, add one of those "
            f"two columns to the field-metadata CSV ({field_metadata}). This build only reads "
            f"on-disk inputs and never runs or alters the pipeline, so fix the file or flag and "
            f"re-run."
        )
        sys.exit(1)
    except Exception as e:
        print_error(f"Failed to build manifest: {e}")
        if ctx.obj.get("verbose"):
            import traceback

            console.print(traceback.format_exc())
        sys.exit(1)

    n_ctrl = len(m.controls())
    console.print(
        f"  {len(m)} rows: {len(m.biological_samples())} sample(s), {n_ctrl} control(s), "
        f"{len(m.seq_run_ids())} sequencing run(s)"
    )

    if output is not None:
        m.to_csv(output)
        print_success(f"Manifest written to {output}")

    if abundance is not None:
        try:
            result = validate_against_abundance(m, abundance, raise_on_orphan=strict)
        except ValueError as e:
            print_error(str(e))
            sys.exit(1)
        if result.ok:
            print_success(
                f"Cross-check OK: all {len(result.abundance_samples)} abundance sample(s) "
                f"have a manifest row"
            )
        else:
            print_warning(
                f"{len(result.orphan_abundance_columns)} abundance sample column(s) have no "
                f"manifest row: {result.orphan_abundance_columns}"
            )
        if result.manifest_extra_rows:
            console.print(
                f"  {len(result.manifest_extra_rows)} manifest eventID(s) absent from the "
                f"abundance table (likely dropped by the pipeline)"
            )


@main.command()
@click.argument("abundance_csv", type=click.Path(exists=True, path_type=Path))
@click.argument("field_metadata", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
@click.option("--mode", type=click.Choice(["flag", "subtract"], case_sensitive=False),
              default="flag", help="flag (annotate, default) or subtract (remove control reads)")
@click.option("--project-metadata", type=click.Path(exists=True, path_type=Path), default=None,
              help="Project metadata CSV (marker -> target_gene)")
@click.option("--id-col", type=str, default=None,
              help="OTU/ASV identifier column in the abundance table (default: first column)")
@click.option("--report", "report_path", type=click.Path(path_type=Path), default=None,
              help="Per-sample cleaning report CSV (default: <output stem>_report.csv)")
@click.pass_context
def clean(
    ctx: click.Context,
    abundance_csv: Path,
    field_metadata: Path,
    output: Path,
    mode: str,
    project_metadata: Optional[Path],
    id_col: Optional[str],
    report_path: Optional[Path],
) -> None:
    """
    Decontaminate an abundance table against its negative controls.

    Derives control identity (extraction vs PCR blanks, extraction batches) from a FAIRe
    manifest migrated from FIELD_METADATA, then flags (default) or subtracts control reads:
    extraction blanks clean their own extraction_ID batch, PCR blanks clean the whole
    dataset. Standalone and read-only on its inputs; every assumption is logged as a [WARN].

    \b
    ABUNDANCE_CSV:  OTU/ASV x sample table (e.g. 02_swarm/<marker>/otu_table.csv).
    FIELD_METADATA: per-sample field metadata CSV (metadata_field_*.csv).
    OUTPUT:         path for the cleaned abundance CSV.

    Negative controls (blanks) are samples carried through extraction or PCR with no
    biological template; any reads they accumulate are contamination, so they mark which
    OTUs/reads to distrust in real samples. Extraction blanks clean only their own
    extraction batch; PCR blanks clean the whole dataset. ``flag`` annotates suspect OTUs,
    ``subtract`` removes the control read counts.

    Args:
        ctx: The Click context (carries the global verbose flag, read for traceback depth).
        abundance_csv: OTU/ASV-by-sample abundance table to decontaminate. Must exist.
        field_metadata: Per-sample field metadata CSV (``metadata_field_*.csv``) from which
            control identity and extraction batches are derived. Must exist.
        output: Path for the cleaned abundance CSV.
        mode: ``flag`` (annotate suspect OTUs, the default) or ``subtract`` (remove control
            read counts). Case-insensitive.
        project_metadata: Optional project metadata CSV (marker -> target_gene). Must exist
            if given.
        id_col: OTU/ASV identifier column name in the abundance table. Defaults to the
            first column if omitted.
        report_path: Optional path for the per-sample cleaning report CSV. Defaults to the
            output stem with a ``_report.csv`` suffix.

    Returns:
        None. Writes the cleaned abundance CSV and the per-sample report, and prints
        control/sample/flagged-OTU counts (and reads removed in ``subtract`` mode).

    Raises:
        SystemExit: Code 1 if an input file is missing, on a value error (e.g. the id
            column or a required metadata column is absent), or on any other cleaning
            failure.
    """
    import pandas as pd

    from seednap.config.manifest_migrate import migrate_to_manifest
    from seednap.steps.cleaning import CleaningProcessor

    _add_command_log_file(Path(output).parent, "clean")

    console.print("\n[bold]Cleaning abundance table[/bold]")
    console.print(f"  Abundance: {abundance_csv}")
    console.print(f"  Metadata:  {field_metadata}")
    console.print(f"  Mode:      {mode}")
    try:
        manifest = migrate_to_manifest(field_metadata, project_csv=project_metadata)
        abundance = pd.read_csv(abundance_csv)
        col = id_col or str(abundance.columns[0])
        cleaned, report, result = CleaningProcessor(mode=mode.lower()).clean(
            abundance, manifest, id_col=col
        )
    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except ValueError as e:
        print_error(
            f"Cleaning failed: {e}. Most often the OTU/ASV id column was not found: by default "
            f"`clean` uses the first column of {abundance_csv} as the id, so pass "
            f"--id-col <name> if the identifier is elsewhere. If the error names the field "
            f"metadata, fix {field_metadata} (it needs an eventID/samp_name column and ISO-8601 "
            f"dates). This command is read-only on its inputs; correct the named file or flag "
            f"and re-run."
        )
        sys.exit(1)
    except Exception as e:
        print_error(f"Failed to clean: {e}")
        if ctx.obj.get("verbose"):
            import traceback

            console.print(traceback.format_exc())
        sys.exit(1)

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    cleaned.to_csv(output, index=False)
    rep_path = Path(report_path) if report_path else output.with_name(f"{output.stem}_report.csv")
    report.to_csv(rep_path, index=False)

    console.print(
        f"  {result.n_controls} control(s), {result.n_samples} sample(s), "
        f"{result.n_otus_flagged} OTU(s) flagged"
        + (f", {result.total_reads_removed} reads removed" if mode == "subtract" else "")
    )
    print_success(f"Cleaned table -> {output}")
    print_success(f"Per-sample report -> {rep_path}")


@main.command()
@click.argument("marker")
@click.option("-o", "--output-dir", type=click.Path(path_type=Path), default=Path("outputs"),
              help="Base output directory of the run (default: outputs)")
@click.option("--state-file", type=click.Path(path_type=Path), default=None,
              help="Run state JSON (default: <output-dir>/.<marker>_state.json)")
@click.pass_context
def monitor(ctx: click.Context, marker: str, output_dir: Path, state_file: Optional[Path]) -> None:
    """
    Summarise a finished or in-progress run from its state JSON.

    Reads <output-dir>/.<marker>_state.json (written after every step) and prints a
    per-step status/duration table plus the read-tracking headline (raw -> final reads,
    mean retention, warnings). When per-sample counts are present it also writes a
    monitoring_summary.csv. Standalone and read-only -- regenerable any time, no re-run.

    MARKER: marker name (e.g. teleo, mam07).

    The state JSON is the source of truth for "did this run finish?": the orchestrator
    writes it after every step with that step's status, duration and outputs. This command
    renders it as a status table and a read-tracking headline, so a user can check on a
    run (finished or still going) without re-executing anything.

    Args:
        ctx: The Click context (carries the global verbose/quiet flags).
        marker: Marker name (e.g. ``teleo``, ``mam07``). Used to locate the default state
            file and name the per-sample summary output.
        output_dir: Base output directory of the run. Defaults to ``outputs``.
        state_file: Path to the run-state JSON. Defaults to
            ``<output_dir>/.<marker>_state.json``. Must exist (at the resolved path) for
            the command to proceed.

    Returns:
        None. Prints the per-step status/duration table and read-tracking headline, and
        writes ``monitoring_summary.csv`` when per-sample counts are present in the state.

    Raises:
        SystemExit: Code 1 if the state file is not found at the resolved path, or if it
            cannot be loaded (truncated or incompatible schema).
    """
    from seednap.pipeline.state import PipelineState

    _add_command_log_file(output_dir, marker)

    sf = Path(state_file) if state_file else (output_dir / f".{marker}_state.json")
    if not sf.exists():
        print_error(
            f"State file not found: {sf}. `monitor` reads the run-state JSON written by "
            f"`run-pipeline`, expected at <output-dir>/.<marker>_state.json (output-dir "
            f"default: outputs). Note the run writes it under its config's paths.output, so "
            f"for the shipped configs that is outputs_test/<marker>/ (e.g. outputs_test/mam07). "
            f"Check the MARKER spelling ('{marker}'), pass --output-dir matching the run's "
            f"paths.output, or pass --state-file pointing straight at the JSON."
        )
        sys.exit(1)
    try:
        state = PipelineState.load(sf)
    except Exception as e:
        print_error(
            f"Could not load state file {sf}: {e}. It is likely truncated (the run was "
            f"interrupted mid-write) or from an incompatible older schema version. Re-run the "
            f"pipeline to regenerate it; deleting the file removes monitor's input and the "
            f"--resume checkpoint, so only delete it if you do not need to resume."
        )
        sys.exit(1)

    console.print(f"\n[bold]Run monitor:[/bold] [cyan]{state.marker}[/cyan]")
    console.print(f"  started:   {state.started_at}")
    console.print(f"  completed: {state.completed_at if state.completed_at else '(not finished)'}")

    table = Table(title="Pipeline steps")
    for col in ("step", "status", "duration", "read tracking"):
        table.add_column(col)
    per_sample_step = None
    for name, step in state.steps.items():
        status = getattr(step.status, "value", str(step.status))
        dur = f"{step.duration_seconds:.1f}s" if step.duration_seconds is not None else "-"
        rt = step.metadata.get("read_tracking") if isinstance(step.metadata, dict) else None
        if rt:
            summary = (
                f"{rt.get('n_samples', '?')} samples, "
                f"raw {rt.get('raw_reads_total', '?')} -> {rt.get('final_step', '')} "
                f"{rt.get('final_reads_total', '?')} "
                f"({rt.get('mean_retention_pct', '?')}% mean), {rt.get('n_warnings', 0)} warn"
            )
        else:
            summary = "-"
        table.add_row(name, status, dur, summary)
        if isinstance(step.metadata, dict) and step.metadata.get("read_tracking_per_sample"):
            per_sample_step = step
    console.print(table)

    # Write the per-sample monitoring summary (the E4 artifact) when counts are present.
    if per_sample_step is not None:
        import pandas as pd

        ps = per_sample_step.metadata["read_tracking_per_sample"]
        rows = [{"eventID": ev, **counts} for ev, counts in ps.items()]
        df = pd.DataFrame(rows)
        report_dir = output_dir / "04_report" / marker
        report_dir.mkdir(parents=True, exist_ok=True)
        out_csv = report_dir / "monitoring_summary.csv"
        df.to_csv(out_csv, index=False)
        print_success(f"Per-sample monitoring summary written to {out_csv} ({len(df)} samples)")
    else:
        print_warning(
            "No per-sample read-tracking counts in the state JSON "
            "(include the 'report' step in pipeline.steps to populate them)."
        )


if __name__ == "__main__":
    main()
