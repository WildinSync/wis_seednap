"""Command-line interface for seednap pipeline."""

import sys
from pathlib import Path
from typing import Any, Dict, Optional

import click
from rich.console import Console
from rich.table import Table

from seednap.__version__ import __version__
from seednap.config import ConfigError, create_example_config, load_config, validate_config_file
from seednap.utils.logging import setup_logging

console = Console()


def print_error(message: str) -> None:
    """Print error message in red."""
    console.print(f"[bold red]Error:[/bold red] {message}")


def print_success(message: str) -> None:
    """Print success message in green."""
    console.print(f"[bold green]✓[/bold green] {message}")


def print_warning(message: str) -> None:
    """Print warning message in yellow."""
    console.print(f"[bold yellow]Warning:[/bold yellow] {message}")


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
    """
    # Store options in context for subcommands
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet

    # Setup basic logging (subcommands may reconfigure)
    level = "DEBUG" if verbose else "WARNING" if quiet else "INFO"
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
    """
    console.print(f"\n[bold]Validating configuration:[/bold] {config_file}\n")

    is_valid, error_message = validate_config_file(config_file)

    if is_valid:
        print_success("Configuration is valid!")

        # Load and display config summary
        try:
            config = load_config(config_file)

            table = Table(title="Configuration Summary", show_header=True, header_style="bold cyan")
            table.add_column("Setting", style="cyan")
            table.add_column("Value", style="white")

            table.add_row("Marker", config.marker.name)
            table.add_row("Taxonomic Method", config.taxonomy.method)
            table.add_row("Output Directory", str(config.paths.output))
            table.add_row("Trimming Cores", str(config.trimming.cores))

            if config.demultiplex.enabled:
                table.add_row("Demultiplexing", f"Enabled ({config.demultiplex.protocol})")
            else:
                table.add_row("Demultiplexing", "Disabled")

            # Surface the database actually used for the selected method, and flag any referenced
            # path missing on disk (a config can be valid yet point at a file that is not there).
            # Read-only checks; nothing is created.
            def _exists(p: Path) -> str:
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

        except Exception as e:
            print_warning(f"Could not load config for summary: {e}")

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
    """
    if output.exists() and not force:
        print_error(f"File already exists: {output}")
        console.print("Use --force to overwrite.")
        sys.exit(1)

    try:
        create_example_config(output, marker=marker, minimal=minimal)
        print_success(f"Created example configuration: {output}")
        console.print(f"\nEdit this file to customize for your analysis.")
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
    """
    from seednap.steps.formatting.gbif_formatter import GBIFFormatter

    console.print(f"\n[bold]Converting to GBIF format:[/bold] {input_file}")
    console.print(f"Input format: {format_type}\n")

    try:
        # Determine output path if not provided
        if output is None:
            output = input_file.parent / f"{input_file.stem}_gbif_input.csv"

        formatter = GBIFFormatter()
        df_out = formatter.from_method(format_type, input_file, output)

        # Print success message with stats
        print_success(f"Converted to GBIF format!")
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
        print_error(f"Invalid input file: {e}")
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
    """
    from seednap.steps.formatting.darwincore_builder import DarwinCoreBuilder

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

    This command:
    1. Creates BLAST database (if needed) from reference FASTA
    2. Runs blastn search with configurable parameters
    3. Extracts phylogeny from reference database headers
    4. Filters hits by percent identity thresholds (species/genus/family)
    5. Resolves ambiguous hits using LCA (Lowest Common Ancestor)
    6. Merges taxonomy with ASV abundance table
    7. Outputs final table with taxonomy and counts
    """
    from seednap.steps.taxonomic_assignment import BlastRunner, BlastTaxonomicAssigner

    console.print(f"\n[bold]Running BLAST taxonomic assignment[/bold]")
    console.print(f"Query: {query_fasta}")
    console.print(f"Reference: {ref_fasta}")
    console.print(f"ASV counts: {asv_count}\n")

    try:
        # Determine output path
        if output is None:
            output = query_fasta.parent / f"{query_fasta.stem}_blast_taxonomy.csv"

        # Create temporary directory for BLAST output
        blast_output_dir = query_fasta.parent / "blast_temp"
        blast_output_dir.mkdir(exist_ok=True)

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
        console.print(f"Total ASVs with taxonomy: [green]{len(result)}[/green]")

        # Show taxonomic resolution summary
        taxonomic_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        console.print("\n[bold]Taxonomic resolution:[/bold]")
        for rank in taxonomic_ranks:
            if rank in result.columns:
                n_assigned = result[rank].notna().sum()
                pct = (n_assigned / len(result)) * 100
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
    """
    from seednap.steps.trimming import StandardTrimmer

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
    """
    from seednap.steps.trimming import LigationTrimmer

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
    """
    from seednap.steps.dada2 import Dada2Processor

    console.print(f"\n[bold]Running DADA2 processing:[/bold]")
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
        )

        print_success("\nDADA2 processing completed successfully!")
        console.print(f"\nOutput files:")
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
        import traceback

        console.print(traceback.format_exc())
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
    """
    from seednap.steps.swarm import SwarmProcessor

    console.print(f"\n[bold]Running SWARM OTU clustering:[/bold]")
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
        console.print(f"\nOutput files:")
        console.print(f"  Query FASTA: {outputs['query_fasta']}")
        console.print(f"  OTU table: {outputs['seqtab_clean_t']}")
        console.print(f"  Full OTU table: {outputs['otu_table_full']}")
        console.print(f"  Merged reads: {outputs['merged_dir']}")
        console.print()

    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except Exception as e:
        print_error(f"SWARM processing failed: {e}")
        import traceback

        console.print(traceback.format_exc())
        sys.exit(1)


@main.command()
@click.argument("method", type=click.Choice(["blast", "dada2", "ecotag", "decipher"]))
@click.argument("marker", type=str)
@click.argument("query_fasta", type=click.Path(exists=True, path_type=Path))
@click.argument("asv_count_csv", type=click.Path(exists=True, path_type=Path))
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
def assign_taxonomy(
    method: str,
    marker: str,
    query_fasta: Path,
    asv_count_csv: Path,
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
    """
    from seednap.steps.taxonomic_assignment import TaxonomicAssigner

    console.print(f"\n[bold]Taxonomic Assignment:[/bold]")
    console.print(f"Method: {method}")
    console.print(f"Marker: {marker}")
    console.print(f"Query: {query_fasta}")
    console.print(f"ASV counts: {asv_count_csv}")
    console.print(f"Output directory: {output_dir}\n")

    try:
        # Initialize assigner
        assigner = TaxonomicAssigner(
            method=method,
            marker=marker,
            output_dir=output_dir,
        )

        # Prepare method-specific arguments
        kwargs = {}

        if method == "blast":
            if not reference_fasta:
                print_error("--reference-fasta is required for BLAST method")
                sys.exit(1)
            kwargs.update({
                "reference_fasta": reference_fasta,
                "threshold_species": threshold_species,
                "threshold_genus": threshold_genus,
                "threshold_family": threshold_family,
                "threshold_order": threshold_order,
                "threshold_class": threshold_class,
                "top_bitscore_pct": top_bitscore_pct,
                "lca_pident_delta": lca_pident_delta,
                "lca_algorithm": lca_algorithm,
                "lca_pid": lca_pid,
                "lca_diff": lca_diff,
            })

        elif method == "dada2":
            if not rdp_db or not species_db:
                print_error("--rdp-db and --species-db are required for DADA2 method")
                sys.exit(1)
            kwargs.update({
                "rdp_db_path": rdp_db,
                "species_db_path": species_db,
            })

        elif method == "ecotag":
            if not taxonomy_db or not reference_db:
                print_error("--taxonomy-db and --reference-db are required for ecotag method")
                sys.exit(1)
            kwargs.update({
                "taxonomy_db": taxonomy_db,
                "reference_db": reference_db,
            })

        elif method == "decipher":
            if not trained_classifier:
                print_error("--trained-classifier is required for DECIPHER method")
                sys.exit(1)
            kwargs.update({
                "trained_classifier_path": trained_classifier,
                "threshold": confidence_threshold,
                "processors": processors,
            })

        # Run taxonomic assignment
        console.print(f"[bold]Running {method.upper()} taxonomic assignment...[/bold]")
        outputs = assigner.assign_taxonomy(
            query_fasta=query_fasta,
            asv_count_csv=asv_count_csv,
            **kwargs,
        )

        print_success(f"\nTaxonomic assignment completed!")
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
        print_error(f"Taxonomic assignment failed: {e}")
        import traceback

        console.print(traceback.format_exc())
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
    """
    from seednap.pipeline.orchestrator import PipelineOrchestrator

    try:
        console.print("\n[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]")
        console.print("[bold cyan]     SeeDNAP eDNA Metabarcoding Pipeline[/bold cyan]")
        console.print("[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]\n")

        # Load config to show marker info
        from seednap.config.loader import load_config

        config_obj = load_config(config)
        console.print(f"[bold]Marker:[/bold] {config_obj.marker.name}")
        console.print(f"[bold]Description:[/bold] {config_obj.marker.description}")
        console.print(f"[bold]Taxonomy method:[/bold] {config_obj.taxonomy.method}")

        if resume:
            console.print(f"\n[yellow]Resuming from previous run[/yellow]")

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
            console.print(f"\n[bold]Completed steps:[/bold]")
            for step_name, step_info in summary["steps"].items():
                if step_info["status"] == "completed":
                    duration = step_info["duration_seconds"]
                    console.print(
                        f"  [green]✓[/green] {step_name}: {duration:.1f}s"
                        if duration
                        else f"  [green]✓[/green] {step_name}"
                    )

        if summary["failed"] > 0:
            console.print(f"\n[bold yellow]Failed steps:[/bold yellow]")
            for step_name, step_info in summary["steps"].items():
                if step_info["status"] == "failed":
                    error = step_info.get("error", "Unknown error")
                    console.print(f"  [red]✗[/red] {step_name}: {error}")

        console.print(f"\n[bold]Output directory:[/bold] {config_obj.paths.output}")
        console.print(f"[bold]Log directory:[/bold] {config_obj.paths.logs}")
        console.print()

    except ValueError as e:
        print_error(str(e))
        sys.exit(1)
    except FileNotFoundError as e:
        print_error(str(e))
        sys.exit(1)
    except Exception as e:
        print_error(f"Pipeline failed: {e}")
        import traceback

        console.print(traceback.format_exc())
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
    """
    import pandas as pd

    from seednap.steps.report import HTMLReportBuilder, ReadTrackingBuilder

    out = output_dir
    dada2_dir = out / "02_dada2" / marker
    swarm_otu = out / "02_swarm" / marker / "otu_table.csv"

    kwargs: Dict[str, Any] = {
        "marker": marker, "logs_dir": out / "logs",
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
            print_error(f"No samples found under {out} (need logs/ and a cluster output).")
            sys.exit(1)

        report_dir = out / "04_report" / marker
        paths = builder.write(report_dir, df=df)
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

        print_success(f"Wrote {paths['read_tracking_csv']}")
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
                summary={
                    "warn_below_retention_pct": warn_retention,
                    "subtitle": f"{len(df)} samples · marker {marker}",
                    "provenance": {"dataset_name": marker, "marker": marker},
                },
            ).write(report_dir / "report.html")
            print_success(f"Wrote HTML report: {html_path}")
        console.print()

    except Exception as e:
        print_error(f"Report failed: {e}")
        sys.exit(1)


@main.command()
def version() -> None:
    """Show version information."""
    console.print(f"\n[bold]seednap[/bold] version [cyan]{__version__}[/cyan]\n")
    console.print("eDNA metabarcoding pipeline with DADA2")
    console.print("Repository: https://github.com/WildinSync/wis_seednap\n")


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
    """
    from seednap.config.manifest import validate_against_abundance
    from seednap.config.manifest_migrate import migrate_to_manifest

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
        print_error(f"Manifest build failed: {e}")
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
    """
    import pandas as pd

    from seednap.config.manifest_migrate import migrate_to_manifest
    from seednap.steps.cleaning import CleaningProcessor

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
        print_error(f"Cleaning failed: {e}")
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
    """
    from seednap.pipeline.state import PipelineState

    sf = Path(state_file) if state_file else (output_dir / f".{marker}_state.json")
    if not sf.exists():
        print_error(f"State file not found: {sf}")
        sys.exit(1)
    try:
        state = PipelineState.load(sf)
    except Exception as e:
        print_error(f"Could not load state file {sf}: {e}")
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
            "(run the pipeline with report.read_tracking enabled to populate them)."
        )


if __name__ == "__main__":
    main()
