"""OBITools ecotag subprocess wrapper for taxonomic assignment.

This module provides a Python wrapper around OBITools (ecotag, obiannotate, obitab)
for taxonomic assignment of eDNA sequences.

OBITools v1 has Python 2 dependencies that conflict with Python 3 tools (Cutadapt
5.x, etc.) so it lives in a separate conda env. The runner auto-discovers the
env by probing common install locations; users can override via the
`SEEDNAP_OBITOOLS_BIN` environment variable.
"""

import logging
import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Union

from seednap.utils.subprocess import run_subprocess

logger = logging.getLogger(__name__)


# Standard OBITools installation locations to probe when ecotag is not on PATH.
# Order matters: we prefer the explicit override, then well-known shared envs,
# then the user's home conda envs.
_OBITOOLS_CANDIDATE_BINS = [
    "/opt/anaconda3/envs/obitools/bin",
    "/opt/conda/envs/obitools/bin",
    str(Path.home() / "miniconda3" / "envs" / "obitools" / "bin"),
    str(Path.home() / ".conda" / "envs" / "obitools" / "bin"),
    str(Path.home() / "anaconda3" / "envs" / "obitools" / "bin"),
]

_REQUIRED_OBITOOLS = ("ecotag", "obiannotate", "obitab")


class EcotagError(Exception):
    """Exception raised for ecotag command errors."""

    pass


def _find_obitools_bin() -> Optional[Path]:
    """Locate a directory containing the required OBITools binaries.

    Returns the first directory that has every required tool, or None if
    OBITools cannot be found anywhere we know to look.

    Resolution order:
        1. ``SEEDNAP_OBITOOLS_BIN`` env var (user override).
        2. Whatever's already on PATH (so an activated obitools env wins).
        3. The well-known install locations in :data:`_OBITOOLS_CANDIDATE_BINS`.
    """
    override = os.environ.get("SEEDNAP_OBITOOLS_BIN")
    if override:
        candidates: List[Path] = [Path(override)]
    else:
        candidates = []

    # If everything is already on PATH, use that.
    if all(shutil.which(t) for t in _REQUIRED_OBITOOLS):
        # Use the directory containing the first tool as the bin dir.
        first = shutil.which(_REQUIRED_OBITOOLS[0])
        if first is not None:
            candidates.insert(0, Path(first).parent)

    candidates.extend(Path(p) for p in _OBITOOLS_CANDIDATE_BINS)

    for d in candidates:
        if d.is_dir() and all((d / t).is_file() for t in _REQUIRED_OBITOOLS):
            return d
    return None


class EcotagRunner:
    """
    Run OBITools ecotag for taxonomic assignment.

    This class wraps the OBITools command-line tools:
    - ecotag: Taxonomic assignment against reference database
    - obiannotate: Clean/filter FASTA annotations
    - obitab: Convert FASTA to TSV table
    """

    def __init__(self, timeout: int = 3600, bin_dir: Optional[Union[str, Path]] = None):
        """
        Initialize ecotag runner.

        Args:
            timeout: Command timeout in seconds (default: 3600 = 1 hour)
            bin_dir: Optional path to the directory containing the OBITools
                binaries. If not provided, the runner auto-discovers from
                PATH / SEEDNAP_OBITOOLS_BIN / well-known conda env paths.
        """
        self.timeout = timeout
        self.bin_dir = self._resolve_bin_dir(bin_dir)

    @staticmethod
    def _resolve_bin_dir(bin_dir: Optional[Union[str, Path]]) -> Path:
        """Find a usable OBITools bin directory or raise with a clear message."""
        if bin_dir is not None:
            d = Path(bin_dir)
            missing = [t for t in _REQUIRED_OBITOOLS if not (d / t).is_file()]
            if missing:
                raise EcotagError(
                    f"OBITools bin_dir '{d}' does not contain the required "
                    f"tools: {missing}. Found: {[t for t in _REQUIRED_OBITOOLS if (d / t).is_file()]}"
                )
            return d

        discovered = _find_obitools_bin()
        if discovered is None:
            probed = [
                "$SEEDNAP_OBITOOLS_BIN (env var override)",
                "PATH",
                *_OBITOOLS_CANDIDATE_BINS,
            ]
            raise EcotagError(
                "OBITools not found. ecotag/obiannotate/obitab were not on "
                "PATH and none of the standard install locations contained "
                "all three binaries.\n\n"
                "OBITools v1 has Python 2 dependencies that conflict with the "
                "main seednap conda env, so it must be installed in a separate "
                "env. On the ETH ELE eDNA server it lives at "
                "/opt/anaconda3/envs/obitools.\n\n"
                "To use ecotag taxonomy:\n"
                "  conda activate obitools  # before running seednap\n"
                "or set:\n"
                "  export SEEDNAP_OBITOOLS_BIN=/opt/anaconda3/envs/obitools/bin\n\n"
                "Probed locations:\n  - " + "\n  - ".join(probed)
            )
        logger.info(f"Discovered OBITools at {discovered}")
        return discovered

    def _tool(self, name: str) -> str:
        """Return the absolute path to an OBITools binary."""
        return str(self.bin_dir / name)

    def _run_command(
        self,
        cmd: list,
        log_file: Optional[Union[str, Path]] = None,
    ) -> str:
        """
        Execute command.

        Args:
            cmd: Command list to execute
            log_file: Optional path to log file for stdout/stderr

        Returns:
            stdout from command

        Raises:
            EcotagError: If command fails
        """
        return run_subprocess(
            cmd,
            timeout=self.timeout,
            log_file=log_file,
            error_class=EcotagError,
        )

    def run_ecotag(
        self,
        query_fasta: Union[str, Path],
        taxonomy_db: Union[str, Path],
        reference_db: Union[str, Path],
        output_fasta: Union[str, Path],
        log_file: Optional[Union[str, Path]] = None,
    ) -> Path:
        """
        Run ecotag taxonomic assignment.

        Args:
            query_fasta: Path to query FASTA file (ASVs)
            taxonomy_db: Path to taxonomy database (NCBI format)
            reference_db: Path to reference sequence database
            output_fasta: Path to output FASTA file with taxonomy annotations
            log_file: Optional path to log file

        Returns:
            Path to output FASTA file

        Raises:
            FileNotFoundError: If input files don't exist
            EcotagError: If ecotag command fails
        """
        query_fasta = Path(query_fasta)
        taxonomy_db = Path(taxonomy_db)
        reference_db = Path(reference_db)
        output_fasta = Path(output_fasta)

        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")
        if not taxonomy_db.exists():
            raise FileNotFoundError(f"Taxonomy database not found: {taxonomy_db}")
        if not reference_db.exists():
            raise FileNotFoundError(f"Reference database not found: {reference_db}")

        output_fasta.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            self._tool("ecotag"),
            "-t",
            str(taxonomy_db),
            "-R",
            str(reference_db),
            str(query_fasta),
        ]

        logger.info(f"Running ecotag on {query_fasta}")
        stdout = self._run_command(cmd, log_file)

        # ecotag writes to stdout, redirect to file
        with open(output_fasta, "w") as f:
            f.write(stdout)

        logger.info(f"Ecotag completed: {output_fasta}")
        return output_fasta

    def clean_annotations(
        self,
        input_fasta: Union[str, Path],
        output_fasta: Union[str, Path],
        tags_to_delete: Optional[list] = None,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Path:
        """
        Clean FASTA annotations using obiannotate.

        Removes unnecessary tags from ecotag output to simplify downstream processing.

        Args:
            input_fasta: Path to input FASTA file (from ecotag)
            output_fasta: Path to output cleaned FASTA file
            tags_to_delete: List of tag names to delete (default: common ecotag tags)
            log_file: Optional path to log file

        Returns:
            Path to cleaned FASTA file

        Raises:
            EcotagError: If obiannotate command fails
        """
        input_fasta = Path(input_fasta)
        output_fasta = Path(output_fasta)

        if not input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

        output_fasta.parent.mkdir(parents=True, exist_ok=True)

        # Default OBITools annotation tags removed before tabulation.
        if tags_to_delete is None:
            tags_to_delete = [
                "scientific_name_by_db",
                "obiclean_samplecount",
                "obiclean_count",
                "obiclean_singletoncount",
                "obiclean_cluster",
                "obiclean_internalcount",
                "obiclean_head",
                "obiclean_headcount",
                "id_status",
                "rank_by_db",
                "obiclean_status",
                "seq_length_ori",
                "sminL",
                "sminR",
                "reverse_score",
                "reverse_primer",
                "reverse_match",
                "reverse_tag",
                "forward_tag",
                "forward_score",
                "forward_primer",
                "forward_match",
                "tail_quality",
            ]

        cmd = [self._tool("obiannotate")]
        for tag in tags_to_delete:
            cmd.extend(["--delete-tag", tag])
        cmd.append(str(input_fasta))

        logger.info(f"Cleaning annotations from {input_fasta}")
        stdout = self._run_command(cmd, log_file)

        # obiannotate writes to stdout
        with open(output_fasta, "w") as f:
            f.write(stdout)

        logger.info(f"Annotations cleaned: {output_fasta}")
        return output_fasta

    def convert_to_table(
        self,
        input_fasta: Union[str, Path],
        output_tsv: Union[str, Path],
        log_file: Optional[Union[str, Path]] = None,
    ) -> Path:
        """
        Convert FASTA to TSV table using obitab.

        Args:
            input_fasta: Path to input FASTA file
            output_tsv: Path to output TSV file
            log_file: Optional path to log file

        Returns:
            Path to output TSV file

        Raises:
            EcotagError: If obitab command fails
        """
        input_fasta = Path(input_fasta)
        output_tsv = Path(output_tsv)

        if not input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

        output_tsv.parent.mkdir(parents=True, exist_ok=True)

        cmd = [self._tool("obitab"), "-o", str(input_fasta)]

        logger.info(f"Converting {input_fasta} to table")
        stdout = self._run_command(cmd, log_file)

        # obitab writes to stdout
        with open(output_tsv, "w") as f:
            f.write(stdout)

        logger.info(f"Table created: {output_tsv}")
        return output_tsv

    def run_complete_workflow(
        self,
        query_fasta: Union[str, Path],
        taxonomy_db: Union[str, Path],
        reference_db: Union[str, Path],
        output_dir: Union[str, Path],
        marker: str,
        log_file: Optional[Union[str, Path]] = None,
    ) -> Dict[str, Path]:
        """
        Run complete ecotag workflow.

        Workflow:
        1. Run ecotag taxonomic assignment
        2. Clean annotations with obiannotate
        3. Convert to table with obitab

        Args:
            query_fasta: Path to query FASTA file (ASVs from DADA2)
            taxonomy_db: Path to taxonomy database
            reference_db: Path to reference sequence database
            output_dir: Output directory
            marker: Marker name (for output file naming)
            log_file: Optional path to log file

        Returns:
            Dictionary with paths to output files:
            - ecotag_fasta: FASTA with taxonomy annotations
            - cleaned_fasta: FASTA with cleaned annotations
            - taxonomy_tsv: TSV table with taxonomy

        Raises:
            EcotagError: If any step fails
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Starting ecotag workflow for {marker}")

        # Step 1: Run ecotag
        ecotag_fasta = output_dir / "query_ecotag.fasta"
        self.run_ecotag(
            query_fasta=query_fasta,
            taxonomy_db=taxonomy_db,
            reference_db=reference_db,
            output_fasta=ecotag_fasta,
            log_file=log_file,
        )

        # Step 2: Clean annotations
        cleaned_fasta = output_dir / "query_ecotag_temp.fasta"
        self.clean_annotations(
            input_fasta=ecotag_fasta,
            output_fasta=cleaned_fasta,
            log_file=log_file,
        )

        # Step 3: Convert to table
        taxonomy_tsv = output_dir / "query_ecotag.tsv"
        self.convert_to_table(
            input_fasta=cleaned_fasta,
            output_tsv=taxonomy_tsv,
            log_file=log_file,
        )

        logger.info(f"Ecotag workflow completed for {marker}")

        return {
            "ecotag_fasta": ecotag_fasta,
            "cleaned_fasta": cleaned_fasta,
            "taxonomy_tsv": taxonomy_tsv,
        }

    def link_with_abundance_table(
        self,
        taxonomy_tsv: Union[str, Path],
        abundance_csv: Union[str, Path],
        output_csv: Union[str, Path],
        sequence_col: str = "sequence",
        contaminants: Optional[list] = None,
    ) -> Path:
        """
        Link ecotag taxonomy with DADA2/SWARM abundance table.

        Delegates to the shared taxonomy post-processor so ecotag, DECIPHER,
        DADA2 RDP, and BLAST all share the same output schema and the same
        correctness guarantees.

        Args:
            taxonomy_tsv: Path to ecotag taxonomy TSV
            abundance_csv: Path to abundance table (seqtab_clean_t.csv)
            output_csv: Path to output CSV file
            sequence_col: Name of sequence column (default: 'sequence')
            contaminants: Optional list of species to flag as contaminants

        Returns:
            Path to output CSV file with merged taxonomy and abundances
        """
        from seednap.utils.taxonomy import link_taxonomy_with_abundance

        return link_taxonomy_with_abundance(
            taxonomy_path=taxonomy_tsv,
            abundance_path=abundance_csv,
            output_path=output_csv,
            sequence_col=sequence_col,
            taxonomy_sep="\t",
            contaminants=contaminants,
        )
