"""Command-line interface for seednap pipeline."""

import sys
from pathlib import Path
from typing import Optional

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
            table.add_row("CPU Cores", str(config.resources.max_cores))

            if config.demultiplex.enabled:
                table.add_row("Demultiplexing", f"Enabled ({config.demultiplex.protocol})")
            else:
                table.add_row("Demultiplexing", "Disabled")

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
    "--force",
    "-f",
    is_flag=True,
    help="Overwrite existing file",
)
def init(output: Path, marker: str, force: bool) -> None:
    """
    Create an example configuration file.

    This generates a template configuration file with sensible defaults
    that you can customize for your analysis.
    """
    if output.exists() and not force:
        print_error(f"File already exists: {output}")
        console.print("Use --force to overwrite.")
        sys.exit(1)

    try:
        create_example_config(output, marker=marker)
        print_success(f"Created example configuration: {output}")
        console.print(f"\nEdit this file to customize for your analysis.")
        console.print(f"Validate it with: [bold]seednap validate {output}[/bold]")
    except ConfigError as e:
        print_error(f"Failed to create config: {e}")
        sys.exit(1)


@main.command()
@click.argument("config_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--resume-from",
    type=click.Choice(["trim", "dada2", "taxonomy", "export"]),
    help="Resume pipeline from a specific step",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show what would be run without executing",
)
@click.pass_context
def run(ctx: click.Context, config_file: Path, resume_from: Optional[str], dry_run: bool) -> None:
    """
    Run the full seednap pipeline.

    CONFIG_FILE: Path to the configuration YAML file.

    This will execute all configured pipeline steps:
    1. Primer trimming (cutadapt)
    2. DADA2 denoising and merging
    3. Taxonomic assignment
    4. Export to various formats
    """
    print_warning("Pipeline execution not yet implemented (Phase 6)")
    console.print("\nThis command will be available after Phase 6 implementation.")
    console.print("For now, use the legacy main.sh script.")
    sys.exit(1)


@main.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--format",
    "-f",
    "format_type",
    type=click.Choice(["dada2", "ecotag"]),
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
    from seednap.steps.format_gbif import format_dada2_to_gbif, format_ecotag_to_gbif

    console.print(f"\n[bold]Converting to GBIF format:[/bold] {input_file}")
    console.print(f"Input format: {format_type}\n")

    try:
        # Determine output path if not provided
        if output is None:
            output = input_file.parent / f"{input_file.stem}_gbif_input.csv"

        # Call appropriate formatter
        if format_type == "dada2":
            df_out = format_dada2_to_gbif(input_file, output)
        elif format_type == "ecotag":
            df_out = format_ecotag_to_gbif(input_file, output)
        else:
            print_error(f"Unknown format type: {format_type}")
            sys.exit(1)

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
    default=98.0,
    type=float,
    help="Minimum percent identity for species-level assignment (default: 98.0)",
)
@click.option(
    "--threshold-genus",
    default=96.0,
    type=float,
    help="Minimum percent identity for genus-level assignment (default: 96.0)",
)
@click.option(
    "--threshold-family",
    default=86.5,
    type=float,
    help="Minimum percent identity for family-level assignment (default: 86.5)",
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
            perc_identity=perc_identity, qcov_hsp_perc=qcov_hsp_perc, evalue=evalue
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
    default=2,
    type=int,
    help="Maximum expected errors for filtering (default: 2)",
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
    max_ee: int,
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

            console.print("\n[bold]Running taxonomic assignment...[/bold]")
            taxo_outputs = processor.assign_taxonomy(
                rdp_db_path=rdp_db,
                species_db_path=species_db,
            )

            print_success("\nTaxonomic assignment completed!")
            console.print(f"  Taxonomy table: {taxo_outputs['complete']}")

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
def version() -> None:
    """Show version information."""
    console.print(f"\n[bold]seednap[/bold] version [cyan]{__version__}[/cyan]\n")
    console.print("eDNA metabarcoding pipeline with DADA2")
    console.print("Repository: https://gitlab.ethz.ch/ele-projects/edna/edna-app/seednap\n")


if __name__ == "__main__":
    main()
