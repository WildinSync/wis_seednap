"""BLAST runner for executing makeblastdb and blastn commands."""

import logging
import subprocess
from pathlib import Path
from typing import List, Union

logger = logging.getLogger(__name__)


class BlastDatabaseError(Exception):
    """Exception raised for BLAST database errors."""

    pass


class BlastRunner:
    """Run BLAST commands (makeblastdb, blastn) via subprocess."""

    def __init__(
        self,
        perc_identity: float = 80.0,
        qcov_hsp_perc: float = 80.0,
        evalue: float = 1e-25,
        max_target_seqs: int = 5,
    ):
        """
        Initialize BLAST runner with search parameters.

        Args:
            perc_identity: Minimum percent identity for hits (default: 80.0)
            qcov_hsp_perc: Minimum query coverage per HSP (default: 80.0)
            evalue: Maximum e-value for hits (default: 1e-25)
            max_target_seqs: Maximum number of target sequences to keep (default: 5)
        """
        self.perc_identity = perc_identity
        self.qcov_hsp_perc = qcov_hsp_perc
        self.evalue = evalue
        self.max_target_seqs = max_target_seqs

    def check_blast_db_exists(self, fasta_path: Union[str, Path]) -> bool:
        """
        Check if BLAST database files exist for given FASTA.

        Args:
            fasta_path: Path to FASTA file

        Returns:
            True if database files exist, False otherwise
        """
        fasta_path = Path(fasta_path)

        # BLAST database files have extensions: .nhr, .nin, .nsq (and .njs for newer versions)
        required_extensions = [".nhr", ".nin", ".nsq"]

        return all((fasta_path.parent / f"{fasta_path.name}{ext}").exists() for ext in required_extensions)

    def create_blast_db(self, fasta_path: Union[str, Path]) -> None:
        """
        Create BLAST database from FASTA file using makeblastdb.

        Args:
            fasta_path: Path to input FASTA file

        Raises:
            FileNotFoundError: If FASTA file does not exist
            BlastDatabaseError: If makeblastdb command fails
        """
        fasta_path = Path(fasta_path)

        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        logger.info(f"Creating BLAST database for {fasta_path}")

        cmd = ["makeblastdb", "-dbtype", "nucl", "-in", str(fasta_path)]

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True, timeout=600  # 10 minute timeout
            )

            logger.debug(f"makeblastdb stdout: {result.stdout}")
            logger.info("BLAST database created successfully")

        except subprocess.CalledProcessError as e:
            error_msg = f"makeblastdb failed: {e.stderr}"
            logger.error(error_msg)
            raise BlastDatabaseError(error_msg) from e

        except subprocess.TimeoutExpired as e:
            error_msg = "makeblastdb timed out after 10 minutes"
            logger.error(error_msg)
            raise BlastDatabaseError(error_msg) from e

        except FileNotFoundError as e:
            error_msg = "makeblastdb command not found. Is BLAST installed?"
            logger.error(error_msg)
            raise BlastDatabaseError(error_msg) from e

    def run_blastn(
        self, query_fasta: Union[str, Path], db_fasta: Union[str, Path], output_tsv: Union[str, Path]
    ) -> None:
        """
        Run blastn search against database.

        Args:
            query_fasta: Path to query sequences FASTA
            db_fasta: Path to database FASTA (database files must exist)
            output_tsv: Path to output TSV file

        Raises:
            FileNotFoundError: If query or database files do not exist
            BlastDatabaseError: If blastn command fails
        """
        query_fasta = Path(query_fasta)
        db_fasta = Path(db_fasta)
        output_tsv = Path(output_tsv)

        # Validate inputs
        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")

        if not db_fasta.exists():
            raise FileNotFoundError(f"Database FASTA not found: {db_fasta}")

        # Ensure database exists
        if not self.check_blast_db_exists(db_fasta):
            logger.info("BLAST database files not found, creating...")
            self.create_blast_db(db_fasta)

        # Create output directory
        output_tsv.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Running BLAST search: {query_fasta} vs {db_fasta}")
        logger.info(
            f"Parameters: pident={self.perc_identity}, qcov={self.qcov_hsp_perc}, "
            f"evalue={self.evalue}, max_targets={self.max_target_seqs}"
        )

        # Build blastn command
        cmd = [
            "blastn",
            "-query",
            str(query_fasta),
            "-db",
            str(db_fasta),
            "-out",
            str(output_tsv),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq",
            "-perc_identity",
            str(self.perc_identity),
            "-qcov_hsp_perc",
            str(self.qcov_hsp_perc),
            "-evalue",
            str(self.evalue),
            "-max_target_seqs",
            str(self.max_target_seqs),
        ]

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True, timeout=3600  # 1 hour timeout
            )

            logger.debug(f"blastn stdout: {result.stdout}")
            logger.info(f"BLAST search completed, output saved to {output_tsv}")

        except subprocess.CalledProcessError as e:
            error_msg = f"blastn failed: {e.stderr}"
            logger.error(error_msg)
            raise BlastDatabaseError(error_msg) from e

        except subprocess.TimeoutExpired as e:
            error_msg = "blastn timed out after 1 hour"
            logger.error(error_msg)
            raise BlastDatabaseError(error_msg) from e

        except FileNotFoundError as e:
            error_msg = "blastn command not found. Is BLAST installed?"
            logger.error(error_msg)
            raise BlastDatabaseError(error_msg) from e

    def run_blast_pipeline(
        self,
        query_fasta: Union[str, Path],
        db_fasta: Union[str, Path],
        output_dir: Union[str, Path],
        marker: str,
    ) -> Path:
        """
        Run complete BLAST pipeline: makeblastdb (if needed) + blastn.

        Args:
            query_fasta: Path to query sequences (ASVs from DADA2)
            db_fasta: Path to reference database FASTA
            output_dir: Directory for BLAST outputs
            marker: Marker name (for output file naming)

        Returns:
            Path to BLAST output TSV file

        Raises:
            FileNotFoundError: If input files do not exist
            BlastDatabaseError: If BLAST commands fail
        """
        output_dir = Path(output_dir)
        output_tsv = output_dir / f"{marker}_blastn_output.tsv"

        # Run BLAST search
        self.run_blastn(query_fasta, db_fasta, output_tsv)

        return output_tsv
