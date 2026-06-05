"""Unified taxonomic assignment interface for all methods.

This module provides a high-level interface for taxonomic assignment that supports
multiple methods (BLAST, DADA2, ecotag, DECIPHER) with a consistent API.
"""

import logging
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, Union

from seednap.steps.taxonomic_assignment.blast_runner import BlastRunner, BlastTaxonomicAssigner
from seednap.steps.taxonomic_assignment.decipher_runner import DecipherRunner
from seednap.steps.taxonomic_assignment.ecotag_runner import EcotagRunner

logger = logging.getLogger(__name__)


class TaxonomyMethod(str, Enum):
    """Supported taxonomic assignment methods."""

    BLAST = "blast"
    DADA2 = "dada2"
    ECOTAG = "ecotag"
    DECIPHER = "decipher"


class TaxonomicAssigner:
    """
    Unified interface for taxonomic assignment.

    This class provides a consistent API for all taxonomic assignment methods,
    handling method selection, configuration, and output formatting.

    Supported methods:
    - BLAST: BLAST search with LCA resolution
    - DADA2: Naive Bayesian classifier (RDP)
    - ecotag: OBITools taxonomic assignment
    - DECIPHER: DECIPHER IdTaxa classifier
    """

    def __init__(
        self,
        method: Union[str, TaxonomyMethod],
        marker: str,
        output_dir: Union[str, Path],
    ):
        """
        Initialize taxonomic assigner.

        Args:
            method: Assignment method ('blast', 'dada2', 'ecotag', or 'decipher')
            marker: Marker name (e.g., 'teleo', 'amph')
            output_dir: Base output directory
        """
        if isinstance(method, str):
            method = TaxonomyMethod(method.lower())

        self.method = method
        self.marker = marker.lower()
        self.output_dir = Path(output_dir)
        self.taxo_dir = self.output_dir / "03_taxo" / self.marker
        self.taxo_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Initialized taxonomic assigner: {self.method.value} for {self.marker}")

    def assign_taxonomy(
        self,
        query_fasta: Union[str, Path],
        asv_count_csv: Union[str, Path],
        **method_specific_kwargs,
    ) -> Dict[str, Path]:
        """
        Assign taxonomy using the specified method.

        This is the main entry point. It routes to the appropriate method-specific
        implementation based on self.method.

        Args:
            query_fasta: Path to query FASTA file (ASVs from DADA2)
            asv_count_csv: Path to ASV count table (seqtab_clean.csv or _t.csv)
            **method_specific_kwargs: Method-specific parameters

        Returns:
            Dictionary with paths to output files (varies by method)

        Raises:
            ValueError: If required method-specific parameters are missing
            FileNotFoundError: If input files don't exist
        """
        query_fasta = Path(query_fasta)
        asv_count_csv = Path(asv_count_csv)

        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")
        if not asv_count_csv.exists():
            raise FileNotFoundError(f"ASV count table not found: {asv_count_csv}")

        logger.info(f"Running taxonomic assignment: {self.method.value}")

        if self.method == TaxonomyMethod.BLAST:
            return self._assign_blast(query_fasta, asv_count_csv, **method_specific_kwargs)
        elif self.method == TaxonomyMethod.DADA2:
            return self._assign_dada2(query_fasta, asv_count_csv, **method_specific_kwargs)
        elif self.method == TaxonomyMethod.ECOTAG:
            return self._assign_ecotag(query_fasta, asv_count_csv, **method_specific_kwargs)
        elif self.method == TaxonomyMethod.DECIPHER:
            return self._assign_decipher(query_fasta, asv_count_csv, **method_specific_kwargs)
        else:
            raise ValueError(f"Unsupported taxonomy method: {self.method}")

    def _assign_blast(
        self,
        query_fasta: Path,
        asv_count_csv: Path,
        reference_fasta: Optional[Union[str, Path]] = None,
        threshold_species: float = 98.0,
        threshold_genus: float = 96.0,
        threshold_family: float = 86.5,
        threshold_order: float = 80.0,
        threshold_class: float = 70.0,
        top_bitscore_pct: float = 10.0,
        lca_pident_delta: float = 1.0,
        lca_algorithm: str = "cascade",
        lca_pid: float = 90.0,
        lca_diff: float = 1.0,
        contaminants: Optional[list] = None,
        perc_identity: float = 80.0,
        qcov_hsp_perc: float = 80.0,
        evalue: float = 1e-25,
        max_target_seqs: int = 5,
        task: str = "megablast",
        **kwargs,
    ) -> Dict[str, Path]:
        """
        Assign taxonomy using BLAST.

        Args:
            query_fasta: Query FASTA file
            asv_count_csv: ASV count table
            reference_fasta: Reference database FASTA (required)
            threshold_species: Percent identity threshold for species (default: 98.0)
            threshold_genus: Percent identity threshold for genus (default: 96.0)
            threshold_family: Percent identity threshold for family (default: 86.5)
            threshold_order: Percent identity threshold for order (default: 80.0)
            threshold_class: Percent identity threshold for class (default: 70.0)
            top_bitscore_pct: LCA bitscore band as percent of best (default: 10.0)
            perc_identity: Minimum percent identity for BLAST hits (default: 80.0)
            qcov_hsp_perc: Minimum query coverage per HSP (default: 80.0)
            evalue: Maximum e-value for BLAST hits (default: 1e-25)
            max_target_seqs: Maximum target sequences to keep (default: 5)

        Returns:
            Dictionary with 'final_table' key pointing to output CSV
        """
        if reference_fasta is None:
            raise ValueError("reference_fasta is required for BLAST method")

        reference_fasta = Path(reference_fasta)

        logger.info(f"Running BLAST taxonomic assignment")

        # Create BLAST database if needed
        runner = BlastRunner(
            perc_identity=perc_identity,
            qcov_hsp_perc=qcov_hsp_perc,
            evalue=evalue,
            max_target_seqs=max_target_seqs,
            task=task,
        )
        if not runner.check_blast_db_exists(reference_fasta):
            logger.info("Creating BLAST database...")
            runner.create_blast_db(reference_fasta)

        # Run BLAST search
        blast_output = self.taxo_dir / "output_blastn.tsv"
        runner.run_blastn(
            query_fasta=query_fasta,
            db_fasta=reference_fasta,
            output_tsv=blast_output,
        )

        # Process BLAST results
        assigner = BlastTaxonomicAssigner(
            reference_fasta=reference_fasta,
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
            contaminants=contaminants,
        )

        final_output = self.output_dir / f"{self.marker}_blast.csv"
        result_df = assigner.assign_taxonomy(
            blast_tsv=blast_output,
            asv_count_csv=asv_count_csv,
            asv_fasta=query_fasta,
            output_path=final_output,
        )

        logger.info(f"BLAST assignment completed: {final_output}")

        return {
            "blast_output": blast_output,
            "final_table": final_output,
        }

    def _assign_dada2(
        self,
        query_fasta: Path,
        asv_count_csv: Path,
        rdp_db_path: Optional[Union[str, Path]] = None,
        species_db_path: Optional[Union[str, Path]] = None,
        multithread: bool = True,
        bootstrap_threshold: int = 80,
        contaminants: Optional[list] = None,
        **kwargs,
    ) -> Dict[str, Path]:
        """
        Assign taxonomy using DADA2 naive Bayesian classifier.

        Note: This method requires the DADA2 runner from the dada2 module.
        It's typically called as part of the DADA2 workflow.

        Args:
            query_fasta: Query FASTA file (not used directly - uses seqtab_clean.rds)
            asv_count_csv: ASV count table
            rdp_db_path: Path to RDP-formatted database (required)
            species_db_path: Path to species database (required)
            multithread: Use multithreading (default: True)

        Returns:
            Dictionary with 'taxonomy' and 'final_table' keys
        """
        if rdp_db_path is None or species_db_path is None:
            raise ValueError("rdp_db_path and species_db_path are required for DADA2 method")

        from seednap.steps.taxonomic_assignment.dada2_taxonomy_runner import Dada2TaxonomyRunner
        from seednap.utils.taxonomy import link_taxonomy_with_abundance

        runner = Dada2TaxonomyRunner()
        outputs = runner.run_dada2_taxonomy(
            marker=self.marker,
            output_dir=self.output_dir,
            rdp_db_path=rdp_db_path,
            species_db_path=species_db_path,
            query_fasta=query_fasta,
            multithread=multithread,
            bootstrap_threshold=bootstrap_threshold,
        )

        # Merge per-sequence taxonomy with the abundance table via the shared
        # post-processor (LEFT-merge, cascade null already done in R, contaminant
        # flag, BLAST-compatible schema). The R script writes `bootstrap_min`
        # which we expose as `pident` for schema parity with BLAST.
        link_taxonomy_with_abundance(
            taxonomy_path=outputs["taxonomy"],
            abundance_path=asv_count_csv,
            output_path=outputs["final_table"],
            sequence_col="sequence",
            contaminants=contaminants,
            pident_col="bootstrap_min",
        )

        logger.info(f"DADA2 assignment completed: {outputs['final_table']}")

        return outputs

    def _assign_ecotag(
        self,
        query_fasta: Path,
        asv_count_csv: Path,
        taxonomy_db: Optional[Union[str, Path]] = None,
        reference_db: Optional[Union[str, Path]] = None,
        contaminants: Optional[list] = None,
        **kwargs,
    ) -> Dict[str, Path]:
        """
        Assign taxonomy using ecotag (OBITools).

        Args:
            query_fasta: Query FASTA file
            asv_count_csv: ASV count table
            taxonomy_db: Path to taxonomy database (NCBI format, required)
            reference_db: Path to reference sequence database (required)
            contaminants: Optional list of species to flag as contaminants

        Returns:
            Dictionary with 'taxonomy_tsv' and 'final_table' keys
        """
        if taxonomy_db is None or reference_db is None:
            raise ValueError("taxonomy_db and reference_db are required for ecotag method")

        taxonomy_db = Path(taxonomy_db)
        reference_db = Path(reference_db)

        logger.info(f"Running ecotag taxonomic assignment")

        runner = EcotagRunner()
        outputs = runner.run_complete_workflow(
            query_fasta=query_fasta,
            taxonomy_db=taxonomy_db,
            reference_db=reference_db,
            output_dir=self.taxo_dir,
            marker=self.marker,
        )

        # Link with abundance table via the shared post-processor
        # (LEFT-merge, cascade null, contaminant flag, BLAST-compatible schema)
        complete_output = self.output_dir / f"{self.marker}_ecotag.csv"
        runner.link_with_abundance_table(
            taxonomy_tsv=outputs["taxonomy_tsv"],
            abundance_csv=asv_count_csv,
            output_csv=complete_output,
            contaminants=contaminants,
        )

        outputs["final_table"] = complete_output

        logger.info(f"Ecotag assignment completed: {complete_output}")

        return outputs

    def _assign_decipher(
        self,
        query_fasta: Path,
        asv_count_csv: Path,
        trained_classifier_path: Optional[Union[str, Path]] = None,
        threshold: int = 60,
        processors: int = 8,
        contaminants: Optional[list] = None,
        **kwargs,
    ) -> Dict[str, Path]:
        """
        Assign taxonomy using DECIPHER.

        Args:
            query_fasta: Query FASTA file (not used directly - uses seqtab_clean.rds)
            asv_count_csv: ASV count table
            trained_classifier_path: Path to trained DECIPHER classifier (.rds, required)
            threshold: Confidence threshold (0-100, default: 60)
            processors: Number of CPU cores (default: 8)
            contaminants: Optional list of species to flag as contaminants

        Returns:
            Dictionary with 'taxonomy' and 'final_table' keys
        """
        if trained_classifier_path is None:
            raise ValueError("trained_classifier_path is required for DECIPHER method")

        trained_classifier_path = Path(trained_classifier_path)

        logger.info(f"Running DECIPHER taxonomic assignment")

        runner = DecipherRunner()
        outputs = runner.run_decipher_assignment(
            marker=self.marker,
            output_dir=self.output_dir,
            trained_classifier_path=trained_classifier_path,
            query_fasta=query_fasta,
            threshold=threshold,
            processors=processors,
        )

        # The R script writes only the per-sequence taxonomy CSV; do the merge
        # with the abundance table in Python via the shared post-processor so
        # we get LEFT-merge + cascade-null + contaminant flag + BLAST schema.
        runner.link_with_abundance_table(
            taxonomy_csv=outputs["taxonomy"],
            abundance_csv=asv_count_csv,
            output_csv=outputs["final_table"],
            contaminants=contaminants,
        )

        logger.info(f"DECIPHER assignment completed: {outputs['final_table']}")

        return outputs

    @staticmethod
    def get_method_requirements(method: Union[str, TaxonomyMethod]) -> Dict[str, str]:
        """
        Get required parameters for a taxonomic assignment method.

        Args:
            method: Assignment method

        Returns:
            Dictionary mapping parameter names to descriptions
        """
        if isinstance(method, str):
            method = TaxonomyMethod(method.lower())

        requirements = {
            TaxonomyMethod.BLAST: {
                "reference_fasta": "Path to reference database FASTA file",
                "threshold_species": "Percent identity threshold for species (optional, default: 98.0)",
                "threshold_genus": "Percent identity threshold for genus (optional, default: 96.0)",
                "threshold_family": "Percent identity threshold for family (optional, default: 86.5)",
            },
            TaxonomyMethod.DADA2: {
                "rdp_db_path": "Path to RDP-formatted taxonomy database",
                "species_db_path": "Path to species-level database",
            },
            TaxonomyMethod.ECOTAG: {
                "taxonomy_db": "Path to NCBI taxonomy database",
                "reference_db": "Path to reference sequence database",
            },
            TaxonomyMethod.DECIPHER: {
                "trained_classifier_path": "Path to trained DECIPHER classifier (.rds file)",
                "threshold": "Confidence threshold 0-100 (optional, default: 60)",
                "processors": "Number of CPU cores (optional, default: 8)",
            },
        }

        return requirements[method]
