"""Convert taxonomic assignment outputs to GBIF-compatible format."""

from pathlib import Path
from typing import Union

import pandas as pd


def add_rank_dada(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add taxonomic rank column based on lowest available taxonomic level.

    DADA2 uses "/" in species column to indicate ambiguous hits at genus level.
    This function determines the rank and cleans up the species column.

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

    # Determine rank
    def determine_rank(row: pd.Series) -> str:
        # Check if species is valid (no "/" and not NA)
        species_valid = pd.notna(row.get("species")) and "/" not in str(row.get("species"))

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
    df.loc[(df["rank"] != "species") & (df["species"].str.contains("/", na=False)), "species"] = pd.NA

    return df


def add_taxon_dada(df: pd.DataFrame) -> pd.DataFrame:
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

    def get_taxon(row: pd.Series) -> str | None:
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


def format_dada2_to_gbif(input_path: Union[str, Path], output_path: Union[str, Path] | None = None) -> pd.DataFrame:
    """
    Convert DADA2 taxonomic assignment output to GBIF-compatible format.

    This function:
    1. Reads DADA2 output CSV (wide format with sample columns)
    2. Keeps taxonomic columns: kingdom, phylum, class, order, family, genus, species, sequence
    3. Transforms to long format (samples become rows)
    4. Filters out zero read counts
    5. Adds 'rank' column (determined from taxonomic assignment)
    6. Adds 'taxon' column (lowest available taxonomic level)
    7. Renames 'filter_code' to 'eventID' (GBIF standard)
    8. Exports to CSV (if output_path specified)

    Args:
        input_path: Path to DADA2 output CSV file
        output_path: Optional path to output GBIF CSV file.
                     If None, returns DataFrame without writing file.
                     If provided but is just a directory, generates filename from input.

    Returns:
        DataFrame in GBIF-compatible long format

    Raises:
        FileNotFoundError: If input file does not exist
        ValueError: If required columns are missing

    Examples:
        >>> df = format_dada2_to_gbif('outputs/teleo_dada2.csv', 'outputs/teleo_gbif.csv')
    """
    # Convert to Path
    input_path = Path(input_path)

    # Check input exists
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Read CSV
    df = pd.read_csv(input_path)

    # Define taxonomic columns to keep
    taxonomic_cols = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "sequence"]

    # Check required columns exist
    missing_cols = [col for col in taxonomic_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Remove X column if present (R index column)
    if "X" in df.columns:
        df = df.drop(columns=["X"])

    # Get sample columns (all columns except taxonomic ones)
    sample_cols = [col for col in df.columns if col not in taxonomic_cols]

    if len(sample_cols) == 0:
        raise ValueError("No sample columns found in input file")

    # Transform to long format
    df_long = df.melt(
        id_vars=taxonomic_cols,
        value_vars=sample_cols,
        var_name="filter_code",
        value_name="nb_reads",
    )

    # Filter out zero counts
    df_long = df_long[df_long["nb_reads"] > 0]

    # Add rank and taxon
    df_long = add_rank_dada(df_long)
    df_long = add_taxon_dada(df_long)

    # Reorder columns and rename filter_code to eventID
    final_cols = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "taxon", "rank", "sequence", "nb_reads", "filter_code"]
    df_out = df_long[final_cols].rename(columns={"filter_code": "eventID"})

    # Write output if path specified
    if output_path:
        output_path = Path(output_path)

        # If output_path is a directory, generate filename
        if output_path.is_dir():
            output_filename = input_path.stem + "_gbif_input.csv"
            output_path = output_path / output_filename

        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Write CSV
        df_out.to_csv(output_path, index=False)

    return df_out


def format_ecotag_to_gbif(input_path: Union[str, Path], output_path: Union[str, Path] | None = None) -> pd.DataFrame:
    """
    Convert Ecotag taxonomic assignment output to GBIF-compatible format.

    Ecotag output has different column names than DADA2:
    - family_name → family
    - genus_name → genus
    - species_name → species
    - order_name → order

    This function:
    1. Reads Ecotag output CSV
    2. Adds placeholder kingdom, phylum, class columns (as NA)
    3. Renames Ecotag-specific columns to standard names
    4. Drops Ecotag-specific metadata columns
    5. Performs same transformation as DADA2 (long format, add rank/taxon)
    6. Exports to CSV (if output_path specified)

    Args:
        input_path: Path to Ecotag output CSV file
        output_path: Optional path to output GBIF CSV file

    Returns:
        DataFrame in GBIF-compatible long format

    Raises:
        FileNotFoundError: If input file does not exist
        ValueError: If required columns are missing

    Examples:
        >>> df = format_ecotag_to_gbif('outputs/teleo_ecotag.csv', 'outputs/teleo_gbif.csv')
    """
    # Convert to Path
    input_path = Path(input_path)

    # Check input exists
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Read CSV
    df = pd.read_csv(input_path)

    # Remove X column if present
    if "X" in df.columns:
        df = df.drop(columns=["X"])

    # Rename Ecotag columns to standard names
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
    taxonomic_cols = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "sequence", "rank"]

    # Check required columns exist
    missing_cols = [col for col in taxonomic_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns after renaming: {missing_cols}")

    # Drop Ecotag-specific metadata columns
    cols_to_drop = [
        "id",
        "definition",
        "count",
        "scientific_name",
    ]
    # Also drop any columns matching these patterns
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
    df = df.drop(columns=cols_to_drop)

    # Get sample columns (all columns except taxonomic ones)
    sample_cols = [col for col in df.columns if col not in taxonomic_cols]

    if len(sample_cols) == 0:
        raise ValueError("No sample columns found in input file")

    # Transform to long format
    df_long = df.melt(
        id_vars=taxonomic_cols,
        value_vars=sample_cols,
        var_name="filter_code",
        value_name="nb_reads",
    )

    # Filter out zero counts
    df_long = df_long[df_long["nb_reads"] > 0]

    # Add rank and taxon (rank already exists from ecotag, but we'll recalculate for consistency)
    df_long = add_rank_dada(df_long)
    df_long = add_taxon_dada(df_long)

    # Reorder columns and rename filter_code to eventID
    final_cols = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "taxon", "rank", "sequence", "nb_reads", "filter_code"]
    # Only include columns that exist
    final_cols = [col for col in final_cols if col in df_long.columns]
    df_out = df_long[final_cols]

    # Rename filter_code to eventID if it exists
    if "filter_code" in df_out.columns:
        df_out = df_out.rename(columns={"filter_code": "eventID"})

    # Write output if path specified
    if output_path:
        output_path = Path(output_path)

        # If output_path is a directory, generate filename
        if output_path.is_dir():
            output_filename = input_path.stem + "_gbif_input.csv"
            output_path = output_path / output_filename

        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Write CSV
        df_out.to_csv(output_path, index=False)

    return df_out
