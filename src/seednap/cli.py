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
    type=click.Choice(["dada2", "ecotag"]),
    required=True,
    help="Input format type",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help="Output file path (default: input_file with _gbif suffix)",
)
def format_gbif(input_file: Path, format: str, output: Optional[Path]) -> None:
    """
    Convert taxonomic assignment results to GBIF format.

    INPUT_FILE: Path to the taxonomic assignment CSV file.

    Transforms the wide-format table to long-format GBIF-compatible output.
    """
    print_warning("GBIF formatting not yet implemented (Phase 1)")
    console.print("\nThis command will be available after Phase 1 implementation.")
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
    print_warning("Primer trimming not yet implemented (Phase 3)")
    console.print("\nThis command will be available after Phase 3 implementation.")
    sys.exit(1)


@main.command()
def version() -> None:
    """Show version information."""
    console.print(f"\n[bold]seednap[/bold] version [cyan]{__version__}[/cyan]\n")
    console.print("eDNA metabarcoding pipeline with DADA2")
    console.print("Repository: https://github.com/eth-edna/seednap\n")


if __name__ == "__main__":
    main()
