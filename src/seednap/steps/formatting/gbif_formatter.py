"""GBIF formatter for converting taxonomic assignments to GBIF-compatible format."""

import logging
from pathlib import Path
from typing import Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


class GBIFFormatter:
    """
    Format taxonomic assignment outputs to GBIF-compatible format.

    This class converts outputs from different taxonomic assignment methods
    (DADA2, ecotag, BLAST, DECIPHER) into a standardized GBIF format suitable
    for biodiversity databases.
    """

    def __init__(self):
        """Initialize GBIF formatter."""
        self.taxonomic_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    def _add_rank(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add taxonomic rank column based on lowest available taxonomic level.

        For DADA2 outputs, "/" in species column indicates ambiguous hits at genus level.

        Args:
            df: DataFrame with taxonomic columns

        Returns:
            DataFrame with added 'rank' column and cleaned 'species' column

        Logic:
            - species: If species column has no "/" and is not NA
            - genus: If genus is not NA but species has "/" or is NA
            - family: If family is not NA but genus is NA
            - higher: For any higher taxonomic level
        """
        df = df.copy()

        def determine_rank(row: pd.Series) -> str:
            # Check if species is valid (no "/" and not NA)
            species_valid = pd.notna(row.get("species")) and "/" not in str(
                row.get("species")
            )

            if species_valid:
                return "species"
            elif pd.notna(row.get("genus")):
                return "genus"
            elif pd.notna(row.get("family")):
                return "family"
            else:
                return "higher"

        df["rank"] = df.apply(determine_rank, axis=1)

        # Clean species column: if rank is not species and species contains "/", set to NA
        if "species" in df.columns:
            df.loc[
                (df["rank"] != "species") & (df["species"].str.contains("/", na=False)),
                "species",
            ] = pd.NA

        return df

    def _add_taxon(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add taxon column containing the lowest available taxonomic assignment.

        Args:
            df: DataFrame with 'rank' and taxonomic columns

        Returns:
            DataFrame with added 'taxon' column

        Logic:
            - If rank is 'species': use species column
            - If rank is 'genus': use genus column
            - If rank is 'family': use family column
            - If rank is 'higher': use first available of order, class, phylum, kingdom
        """
        df = df.copy()

        def get_taxon(row: pd.Series) -> Union[str, None]:
            rank = row.get("rank")

            if rank == "species":
                return row.get("species")
            elif rank == "genus":
                return row.get("genus")
            elif rank == "family":
                return row.get("family")
            elif rank == "higher":
                # Return first non-NA value from higher ranks
                for col in ["order", "class", "phylum", "kingdom"]:
                    if pd.notna(row.get(col)):
                        return row.get(col)
            return None

        df["taxon"] = df.apply(get_taxon, axis=1)

        return df

    def _transform_to_long_format(
        self, df: pd.DataFrame, taxonomic_cols: list
    ) -> pd.DataFrame:
        """
        Transform wide format (samples as columns) to long format (samples as rows).

        Args:
            df: DataFrame in wide format
            taxonomic_cols: List of taxonomic column names

        Returns:
            DataFrame in long format with 'eventID' and 'nb_reads' columns
        """
        # Get sample columns (all columns except taxonomic ones)
        sample_cols = [col for col in df.columns if col not in taxonomic_cols]

        if len(sample_cols) == 0:
            raise ValueError("No sample columns found in input file")

        # Transform to long format
        df_long = df.melt(
            id_vars=taxonomic_cols,
            value_vars=sample_cols,
            var_name="eventID",
            value_name="nb_reads",
        )

        # Filter out zero counts
        df_long = df_long[df_long["nb_reads"] > 0]

        return df_long

    def from_dada2_rdp(
        self,
        input_path: Union[str, Path],
        output_path: Optional[Union[str, Path]] = None,
        add_rank: bool = True,
        add_taxon: bool = True,
    ) -> pd.DataFrame:
        """
        Convert DADA2 RDP taxonomic assignment output to GBIF format.

        This method:
        1. Reads DADA2 output CSV (wide format with sample columns)
        2. Keeps taxonomic columns: kingdom, phylum, class, order, family, genus, species, sequence
        3. Transforms to long format (samples become rows)
        4. Filters out zero read counts
        5. Optionally adds 'rank' column (determined from taxonomic assignment)
        6. Optionally adds 'taxon' column (lowest available taxonomic level)
        7. Exports to CSV (if output_path specified)

        Args:
            input_path: Path to DADA2 output CSV file
            output_path: Optional path to output GBIF CSV file
            add_rank: Whether to add 'rank' column (default: True)
            add_taxon: Whether to add 'taxon' column (default: True)

        Returns:
            DataFrame in GBIF-compatible long format

        Raises:
            FileNotFoundError: If input file does not exist
            ValueError: If required columns are missing
        """
        input_path = Path(input_path)

        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")

        logger.info(f"Converting DADA2 output to GBIF format: {input_path}")

        # Read CSV
        df = pd.read_csv(input_path)

        # Remove X column if present (R index column)
        if "X" in df.columns:
            df = df.drop(columns=["X"])

        # Define taxonomic columns to keep
        taxonomic_cols = [
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "sequence",
        ]

        # Check required columns exist
        missing_cols = [col for col in taxonomic_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Transform to long format
        df_long = self._transform_to_long_format(df, taxonomic_cols)

        # Add rank and taxon if requested
        if add_rank:
            df_long = self._add_rank(df_long)

        if add_taxon:
            if "rank" not in df_long.columns:
                df_long = self._add_rank(df_long)
            df_long = self._add_taxon(df_long)

        # Reorder columns
        final_cols = [
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
        if "taxon" in df_long.columns:
            final_cols.append("taxon")
        if "rank" in df_long.columns:
            final_cols.append("rank")
        final_cols.extend(["sequence", "nb_reads", "eventID"])

        df_out = df_long[final_cols]

        # Write output if path specified
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df_out.to_csv(output_path, index=False)
            logger.info(f"Wrote GBIF output to {output_path}")

        return df_out

    def from_ecotag(
        self,
        input_path: Union[str, Path],
        output_path: Optional[Union[str, Path]] = None,
        add_rank: bool = True,
        add_taxon: bool = True,
    ) -> pd.DataFrame:
        """
        Convert ecotag (OBITools) output to GBIF format.

        Ecotag output has different column names:
        - family_name → family
        - genus_name → genus
        - species_name → species
        - order_name → order

        This method:
        1. Reads ecotag output CSV
        2. Renames ecotag-specific columns to standard names
        3. Adds placeholder kingdom, phylum, class columns (as NA)
        4. Drops ecotag-specific metadata columns
        5. Performs same transformation as DADA2 (long format, add rank/taxon)
        6. Exports to CSV (if output_path specified)

        Args:
            input_path: Path to ecotag output CSV file
            output_path: Optional path to output GBIF CSV file
            add_rank: Whether to add 'rank' column (default: True)
            add_taxon: Whether to add 'taxon' column (default: True)

        Returns:
            DataFrame in GBIF-compatible long format

        Raises:
            FileNotFoundError: If input file does not exist
            ValueError: If required columns are missing
        """
        input_path = Path(input_path)

        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")

        logger.info(f"Converting ecotag output to GBIF format: {input_path}")

        # Read CSV
        df = pd.read_csv(input_path)

        # Remove X column if present
        if "X" in df.columns:
            df = df.drop(columns=["X"])

        # Rename ecotag columns to standard names
        rename_map = {
            "family_name": "family",
            "genus_name": "genus",
            "species_name": "species",
            "order_name": "order",
        }
        df = df.rename(columns=rename_map)

        # Add placeholder columns for kingdom, phylum, class
        df["kingdom"] = pd.NA
        df["phylum"] = pd.NA
        df["class"] = pd.NA

        # Define columns to keep
        taxonomic_cols = [
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "sequence",
        ]

        # Check required columns exist
        missing_cols = [col for col in taxonomic_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns after renaming: {missing_cols}")

        # Drop ecotag-specific metadata columns
        cols_to_drop = ["id", "definition", "count", "scientific_name"]

        # Also drop columns matching certain patterns
        pattern_cols_to_drop = [
            col
            for col in df.columns
            if any(
                pattern in col
                for pattern in [
                    "best_identity",
                    "best_match",
                    "match_count",
                    "species_list",
                    "taxid",
                ]
            )
        ]
        cols_to_drop.extend(pattern_cols_to_drop)

        # Drop columns that exist
        cols_to_drop = [col for col in cols_to_drop if col in df.columns]
        if cols_to_drop:
            df = df.drop(columns=cols_to_drop)

        # Transform to long format
        df_long = self._transform_to_long_format(df, taxonomic_cols)

        # Add rank and taxon if requested
        if add_rank:
            df_long = self._add_rank(df_long)

        if add_taxon:
            if "rank" not in df_long.columns:
                df_long = self._add_rank(df_long)
            df_long = self._add_taxon(df_long)

        # Reorder columns
        final_cols = [
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
        if "taxon" in df_long.columns:
            final_cols.append("taxon")
        if "rank" in df_long.columns:
            final_cols.append("rank")
        final_cols.extend(["sequence", "nb_reads", "eventID"])

        # Only include columns that exist
        final_cols = [col for col in final_cols if col in df_long.columns]
        df_out = df_long[final_cols]

        # Write output if path specified
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df_out.to_csv(output_path, index=False)
            logger.info(f"Wrote GBIF output to {output_path}")

        return df_out

    def from_blast(
        self,
        input_path: Union[str, Path],
        output_path: Optional[Union[str, Path]] = None,
        add_rank: bool = True,
        add_taxon: bool = True,
    ) -> pd.DataFrame:
        """
        Convert BLAST output to GBIF format.

        Args:
            input_path: Path to BLAST output CSV file
            output_path: Optional path to output GBIF CSV file
            add_rank: Whether to add 'rank' column (default: True)
            add_taxon: Whether to add 'taxon' column (default: True)

        Returns:
            DataFrame in GBIF-compatible long format
        """
        # BLAST output should already be in a similar format to DADA2
        return self.from_dada2_rdp(input_path, output_path, add_rank, add_taxon)

    def from_decipher(
        self,
        input_path: Union[str, Path],
        output_path: Optional[Union[str, Path]] = None,
        add_rank: bool = True,
        add_taxon: bool = True,
    ) -> pd.DataFrame:
        """
        Convert DECIPHER output to GBIF format.

        Args:
            input_path: Path to DECIPHER output CSV file
            output_path: Optional path to output GBIF CSV file
            add_rank: Whether to add 'rank' column (default: True)
            add_taxon: Whether to add 'taxon' column (default: True)

        Returns:
            DataFrame in GBIF-compatible long format
        """
        # DECIPHER output should be similar to DADA2
        return self.from_dada2_rdp(input_path, output_path, add_rank, add_taxon)
