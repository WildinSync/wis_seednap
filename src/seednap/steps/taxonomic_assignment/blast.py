"""BLAST-based taxonomic assignment with LCA (Lowest Common Ancestor) resolution.

This module provides functionality for processing BLAST output, extracting phylogenetic
information from reference databases, and resolving ambiguous hits using LCA.

Author: Théophile Sanchez (original), refactored for seednap v0.1.0
"""

from pathlib import Path
from typing import Dict, List, Union

import pandas as pd

from seednap.utils.sequences import fasta_to_df


class BlastOutputFormatter:
    """Format BLAST output by adding phylogenetic information from reference database."""

    BLAST_COLUMNS = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qseq",
        "sseq",
    ]

    TAXONOMIC_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    def __init__(self, reference_fasta: Union[str, Path]):
        """
        Initialize BLAST output formatter.

        Args:
            reference_fasta: Path to reference database FASTA file with phylogeny in headers

        Raises:
            FileNotFoundError: If reference FASTA file does not exist
        """
        self.reference_fasta = Path(reference_fasta)
        if not self.reference_fasta.exists():
            raise FileNotFoundError(f"Reference FASTA not found: {self.reference_fasta}")

        self._phylo_dict = self._load_phylogeny()

    def _load_phylogeny(self) -> Dict[str, str]:
        """
        Load phylogenetic information from reference FASTA headers.

        Returns:
            Dictionary mapping sequence IDs to full header lines
        """
        phylo_dict = {}
        with open(self.reference_fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Extract sequence ID (first part of header)
                    seq_id = line.split()[0].strip(">")
                    phylo_dict[seq_id] = line.strip()
        return phylo_dict

    def format_blast_output(
        self, blast_tsv: Union[str, Path], output_path: Union[str, Path, None] = None
    ) -> pd.DataFrame:
        """
        Format BLAST output by adding phylogenetic information.

        Reads BLAST TSV output and extracts phylogeny from reference database headers.
        Expected header format: >seq_id<TAB>kingdom;phylum;class;order;family;genus;species

        Args:
            blast_tsv: Path to BLAST output TSV file (format 6)
            output_path: Optional path to save formatted output

        Returns:
            DataFrame with BLAST results and phylogenetic columns

        Raises:
            FileNotFoundError: If BLAST TSV file does not exist
            KeyError: If sequence ID from BLAST output not found in reference database
        """
        blast_tsv = Path(blast_tsv)
        if not blast_tsv.exists():
            raise FileNotFoundError(f"BLAST output not found: {blast_tsv}")

        # Read BLAST TSV
        df = pd.read_csv(blast_tsv, sep="\t", header=None, names=self.BLAST_COLUMNS)

        # Add taxonomic columns
        for rank in self.TAXONOMIC_RANKS:
            df[rank] = None

        # Extract phylogeny for each hit
        for i, row in df.iterrows():
            seq_id = row["sseqid"]

            if seq_id not in self._phylo_dict:
                raise KeyError(f"Sequence ID '{seq_id}' not found in reference database")

            # Parse header: >seq_id<TAB>kingdom;phylum;class;order;family;genus;species
            header = self._phylo_dict[seq_id]
            phylo_string = header.replace("\n", "").split("\t")[1]
            phylo_values = phylo_string.split(";")

            # Assign to columns
            phylo = dict(zip(self.TAXONOMIC_RANKS, phylo_values))
            for rank in self.TAXONOMIC_RANKS:
                df.at[i, rank] = phylo[rank]

        # Add blast rank (1 = best hit, 2 = second best, etc.)
        df["blast_rank"] = df.groupby("qseqid").cumcount() + 1

        # Save if output path provided
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, sep="\t", header=True, index=False)

        return df


class BlastPhyloFilter:
    """Filter BLAST hits by percent identity thresholds for each taxonomic rank."""

    def __init__(
        self,
        threshold_species: float = 98.0,
        threshold_genus: float = 96.0,
        threshold_family: float = 86.5,
    ):
        """
        Initialize phylogenetic filter with thresholds.

        Args:
            threshold_species: Minimum percent identity for species-level assignment
            threshold_genus: Minimum percent identity for genus-level assignment
            threshold_family: Minimum percent identity for family-level assignment
        """
        self.thresholds = {
            "species": threshold_species,
            "genus": threshold_genus,
            "family": threshold_family,
        }

    def filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter phylogenetic assignments by percent identity thresholds.

        Sets taxonomic rank to None if percent identity is below threshold.

        Args:
            df: DataFrame with 'pident' and taxonomic rank columns

        Returns:
            Filtered DataFrame
        """
        df = df.copy()

        for phylo_level, threshold in self.thresholds.items():
            # Set to None if below threshold
            df.loc[pd.to_numeric(df["pident"]) < float(threshold), phylo_level] = None

        return df


class BlastLCAResolver:
    """Resolve ambiguous BLAST hits using LCA (Lowest Common Ancestor)."""

    TAXONOMIC_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    @staticmethod
    def resolve_ambiguous_hits(group: pd.DataFrame) -> pd.DataFrame:
        """
        Resolve ambiguous hits with identical bitscores using LCA.

        When multiple hits have the same bitscore, this function finds the most recent
        common ancestor (LCA) by identifying the lowest taxonomic rank where all hits agree.

        Args:
            group: DataFrame of BLAST hits for a single query sequence

        Returns:
            DataFrame with ambiguous hits resolved, marked with 'keep_for_analysis' column
        """
        group = group.reset_index(drop=True)

        if len(group) <= 1:
            # No ambiguity, keep the single hit
            group["keep_for_analysis"] = True
            return group

        # Get best bitscore
        best_bitscore = group["bitscore"].iloc[0]

        # Find all hits with best bitscore
        ambiguous_hits = group[group["bitscore"] == best_bitscore].copy()

        if len(ambiguous_hits) <= 1:
            # Only one hit with best score, no ambiguity
            group["keep_for_analysis"] = group["bitscore"] == best_bitscore
            return group

        # Check if all ambiguous hits have the same phylogeny
        same_phylo = all(
            ambiguous_hits[col][ambiguous_hits[col].notna()].nunique() < 2
            for col in BlastLCAResolver.TAXONOMIC_RANKS
        )

        if same_phylo:
            # All hits agree, keep first one
            group["keep_for_analysis"] = group["bitscore"] == best_bitscore
            return group

        # Phylogeny differs - perform LCA
        # Mark all original rows as not to keep
        group["keep_for_analysis"] = False

        # Create combined row with LCA
        combined_row_data = []
        for col in ambiguous_hits.columns:
            if ambiguous_hits[col].nunique() == 1:
                # All hits agree on this column
                combined_row_data.append(ambiguous_hits[col].iloc[0])
            else:
                # Hits disagree - set to None (this is the LCA point)
                combined_row_data.append(None)

        combined_row = pd.DataFrame([combined_row_data], columns=ambiguous_hits.columns)
        combined_row["keep_for_analysis"] = True

        # Append combined row
        result = pd.concat([group, combined_row], ignore_index=True)

        return result


class BlastTaxonomicAssigner:
    """Complete BLAST-based taxonomic assignment pipeline."""

    def __init__(
        self,
        reference_fasta: Union[str, Path],
        threshold_species: float = 98.0,
        threshold_genus: float = 96.0,
        threshold_family: float = 86.5,
    ):
        """
        Initialize BLAST taxonomic assigner.

        Args:
            reference_fasta: Path to reference database FASTA file
            threshold_species: Minimum percent identity for species-level assignment
            threshold_genus: Minimum percent identity for genus-level assignment
            threshold_family: Minimum percent identity for family-level assignment
        """
        self.formatter = BlastOutputFormatter(reference_fasta)
        self.filter = BlastPhyloFilter(threshold_species, threshold_genus, threshold_family)
        self.lca_resolver = BlastLCAResolver()

    def assign_taxonomy(
        self,
        blast_tsv: Union[str, Path],
        asv_count_csv: Union[str, Path],
        asv_fasta: Union[str, Path],
        output_path: Union[str, Path, None] = None,
    ) -> pd.DataFrame:
        """
        Complete BLAST taxonomic assignment workflow.

        This function:
        1. Formats BLAST output with phylogeny from reference DB
        2. Filters hits by percent identity thresholds
        3. Resolves ambiguous hits using LCA
        4. Merges with ASV count table
        5. Outputs final table with taxonomy and abundances

        Args:
            blast_tsv: Path to BLAST output TSV file
            asv_count_csv: Path to ASV count table CSV (samples x ASVs)
            asv_fasta: Path to ASV sequences FASTA file
            output_path: Optional path to save final table

        Returns:
            DataFrame with taxonomic assignments and ASV counts

        Raises:
            FileNotFoundError: If input files do not exist
        """
        # Format BLAST output
        formatted = self.formatter.format_blast_output(blast_tsv)

        # Replace 'None' strings with actual None
        formatted = formatted.replace("None", None)

        # Filter by thresholds
        filtered = self.filter.filter(formatted)

        # Keep only best hit initially
        filtered["keep_for_analysis"] = filtered["blast_rank"] == 1

        # Resolve ambiguous hits with LCA
        filtered = (
            filtered.groupby("qseqid", group_keys=False)
            .apply(self.lca_resolver.resolve_ambiguous_hits)
        )

        # Keep only resolved hits
        filtered = filtered[filtered["keep_for_analysis"] == True]  # noqa: E712

        # Select columns for output
        phylo_cols = BlastLCAResolver.TAXONOMIC_RANKS
        result = filtered[["qseqid", "pident"] + phylo_cols].rename(columns={"qseqid": "ASV_ID"})

        # Load ASV count table
        asv_count = pd.read_csv(asv_count_csv, sep=",", index_col=0).T

        # Load ASV sequences
        asv_sequences = fasta_to_df(asv_fasta)
        asv_sequences = asv_sequences.rename(columns={"id": "ASV_ID", "sequence": "Sequence"})

        # Merge count table with sequences
        asv_count = pd.merge(
            asv_count, asv_sequences, how="inner", left_index=True, right_on="Sequence"
        ).drop(columns="Sequence")

        # Merge taxonomy with counts
        final_table = pd.merge(result, asv_count, how="inner", on="ASV_ID")

        # Sort by ASV number
        final_table["asv_num"] = final_table["ASV_ID"].str.extract(r"(\d+)").astype(int)
        final_table = final_table.sort_values("asv_num").drop(columns="asv_num")

        # Add sequence back
        final_table = pd.merge(final_table, asv_sequences, how="inner", on="ASV_ID")

        # Save if output path provided
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            final_table.to_csv(output_path, header=True, index=False)

        return final_table