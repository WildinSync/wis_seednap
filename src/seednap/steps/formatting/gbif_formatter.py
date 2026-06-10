"""GBIF formatter for converting taxonomic assignments to GBIF-compatible format.

First step of the formatting stage, run by the ``format-gbif`` command. Each
taxonomic-assignment method (DADA2 RDP, ecotag, BLAST, DECIPHER) writes its
result as a wide table: one row per OTU/ASV (the clustered or denoised sequence
variant standing in for a taxon) with the taxonomy columns plus one numeric
column per sample holding that OTU's read count. This module normalises the
differing column names into one schema and reshapes the table to GBIF's long
"occurrence" layout: one row per (OTU, sample) pair with a single ``nb_reads``
count, dropping the zero counts that a wide table necessarily contains. It also
derives the lowest confidently assigned rank and the matching taxon name for
each row. The long table it produces is the input to ``DarwinCoreBuilder``.
"""

import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional, Union, cast

import pandas as pd

from seednap.utils.taxonomy import TAXONOMIC_RANKS

logger = logging.getLogger(__name__)


class GBIFFormatter:
    """
    Format taxonomic assignment outputs to GBIF-compatible format.

    This class converts outputs from different taxonomic assignment methods
    (DADA2, ecotag, BLAST, DECIPHER) into a standardized GBIF format suitable
    for biodiversity databases.
    """

    def __init__(self) -> None:
        """Initialize GBIF formatter with the standard taxonomic rank list."""
        self.taxonomic_ranks = list(TAXONOMIC_RANKS)

    def _add_rank(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add taxonomic rank column based on lowest available taxonomic level.

        Records, for each OTU, the finest rank at which the assignment is
        considered reliable. A species name containing "/" is a tie between
        several species (an ambiguous hit, common with short markers) and is
        treated as resolved only to genus; the "/" species string is then blanked
        on rows whose rank is coarser than species.

        Args:
            df: DataFrame with taxonomic columns (``species``, ``genus``,
                ``family`` are read; others may be present).

        Returns:
            A copy of ``df`` with an added ``rank`` column (one of ``species``,
            ``genus``, ``family``, ``higher``) and the ``species`` column cleaned
            of "/"-ambiguous values on non-species rows.

        Logic:
            - species: If species column has no "/" and is not NA
            - genus: If genus is not NA but species has "/" or is NA
            - family: If family is not NA but genus is NA
            - higher: For any higher taxonomic level
        """
        df = df.copy()

        def determine_rank(row: pd.Series) -> str:
            """Return the lowest taxonomic rank with a valid assignment for one row.

            Args:
                row: One OTU row; read keys ``species``, ``genus``, ``family``.

            Returns:
                The rank label ``species``, ``genus``, ``family``, or ``higher``.
            """
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

        Picks the single name that best labels each OTU: the value at the row's
        ``rank``, or for ``higher`` rows the first populated rank walking down
        from order. This becomes ``scientificName`` in the DarwinCore output.

        Args:
            df: DataFrame with a ``rank`` column (as added by ``_add_rank``) and
                taxonomic columns.

        Returns:
            A copy of ``df`` with an added ``taxon`` column (the chosen name, or
            None if no rank is populated).

        Logic:
            - If rank is 'species': use species column
            - If rank is 'genus': use genus column
            - If rank is 'family': use family column
            - If rank is 'higher': use first available of order, class, phylum, kingdom
        """
        df = df.copy()

        def get_taxon(row: pd.Series) -> Optional[str]:
            """Return the taxon name at the row's rank, or the lowest filled higher rank.

            Args:
                row: One OTU row; read keys ``rank`` plus the taxonomic columns
                    ``species``/``genus``/``family``/``order``/``class``/
                    ``phylum``/``kingdom``.

            Returns:
                The taxon name string, or None if the relevant rank columns are
                all empty.
            """
            rank = row.get("rank")

            if rank == "species":
                return cast(Optional[str], row.get("species"))
            elif rank == "genus":
                return cast(Optional[str], row.get("genus"))
            elif rank == "family":
                return cast(Optional[str], row.get("family"))
            elif rank == "higher":
                # Return first non-NA value from higher ranks
                for col in ["order", "class", "phylum", "kingdom"]:
                    if pd.notna(row.get(col)):
                        return cast(Optional[str], row.get(col))
            return None

        df["taxon"] = df.apply(get_taxon, axis=1)

        return df

    def _transform_to_long_format(
        self, df: pd.DataFrame, taxonomic_cols: List[str]
    ) -> pd.DataFrame:
        """
        Transform wide format (samples as columns) to long format (samples as rows).

        Reshapes the OTU-by-sample read-count matrix into one row per
        (OTU, sample) pair: each sample's numeric column becomes an ``eventID``
        value and its count a ``nb_reads`` value, which is GBIF's expected
        occurrence layout. Sample columns are detected as the numeric columns
        that are neither taxonomy nor per-OTU annotation columns. Per-OTU
        annotations (``ASV_ID``, ``pident``, ``is_contaminant_candidate``) are
        carried onto every resulting sample row, and rows with zero reads are
        dropped (a wide matrix is mostly zeros: an OTU is absent from most
        samples).

        Args:
            df: DataFrame in wide format (taxonomy columns plus one numeric
                read-count column per sample).
            taxonomic_cols: Taxonomic column names to keep as identifier columns
                during the reshape.

        Returns:
            DataFrame in long format: the ``taxonomic_cols``, any present
            annotation columns, plus ``eventID`` (sample name) and ``nb_reads``
            (integer read count > 0).

        Raises:
            ValueError: If no per-sample numeric read-count columns can be found
                (i.e. every non-taxonomy, non-annotation column is non-numeric),
                which usually means the wrong table was passed.
        """
        # Sample columns are everything that's NOT a known non-sample column.
        # The post-Commit-F BLAST schema includes per-OTU annotations (ASV_ID,
        # pident, is_contaminant_candidate) that must NOT be treated as samples.
        # We identify samples by being numeric and not in the taxonomic / annotation set.
        non_sample_known = set(taxonomic_cols) | {
            "ASV_ID", "pident", "is_contaminant_candidate",
        }
        sample_cols = [
            col for col in df.columns
            if col not in non_sample_known
            and pd.api.types.is_numeric_dtype(df[col])
        ]

        if len(sample_cols) == 0:
            non_numeric_cols = [
                col for col in df.columns
                if col not in non_sample_known
                and not pd.api.types.is_numeric_dtype(df[col])
            ]
            raise ValueError(
                "No per-sample read-count columns found in the input table. GBIF "
                "export expects a wide-format CSV where each sample is its own "
                "numeric column (one column per eventID, holding integer read "
                "counts) alongside the taxonomy columns (kingdom, phylum, class, "
                "order, family, genus, species, sequence). After excluding the "
                "taxonomy and per-OTU annotation columns (ASV_ID, pident, "
                "is_contaminant_candidate), none of the remaining columns were "
                f"numeric. Remaining columns were: {non_numeric_cols}. Check that "
                "you passed the wide abundance/taxonomy table from the taxonomy "
                "step (e.g. <marker>_<method>.csv), and that sample columns contain "
                "numeric read counts rather than text."
            )

        # Carry annotation columns through the melt as id_vars so they survive
        # to the long-format output (per-OTU info should appear on every sample row).
        annotation_cols = [
            c for c in ("ASV_ID", "pident", "is_contaminant_candidate") if c in df.columns
        ]

        df_long = df.melt(
            id_vars=taxonomic_cols + annotation_cols,
            value_vars=sample_cols,
            var_name="eventID",
            value_name="nb_reads",
        )

        # Filter out zero counts
        df_long = df_long[df_long["nb_reads"] > 0]

        return df_long

    def from_method(
        self,
        method: str,
        input_path: Union[str, Path],
        output_path: Optional[Union[str, Path]] = None,
        add_rank: bool = True,
        add_taxon: bool = True,
    ) -> pd.DataFrame:
        """
        Dispatch to the correct formatter based on taxonomy method name.

        Args:
            method: Taxonomy method ('dada2', 'ecotag', 'blast', 'decipher').
            input_path: Path to taxonomy CSV file.
            output_path: Optional path to output GBIF CSV file.
            add_rank: Whether to add 'rank' column (default: True).
            add_taxon: Whether to add 'taxon' column (default: True).

        Returns:
            DataFrame in GBIF-compatible long format.

        Raises:
            ValueError: If method is not recognised.
        """
        dispatch: Dict[str, Callable[..., pd.DataFrame]] = {
            "dada2": self.from_dada2_rdp,
            "ecotag": self.from_ecotag,
            "blast": self.from_blast,
            "decipher": self.from_decipher,
        }
        formatter_fn = dispatch.get(method)
        if formatter_fn is None:
            raise ValueError(
                f"Unknown taxonomy method '{method}'. "
                f"Supported: {', '.join(dispatch)}"
            )
        return formatter_fn(input_path, output_path, add_rank, add_taxon)

    def from_dada2_rdp(
        self,
        input_path: Union[str, Path],
        output_path: Optional[Union[str, Path]] = None,
        add_rank: bool = True,
        add_taxon: bool = True,
        _source_format: str = "dada2",
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
            _source_format: Internal label naming the actual source format
                (dada2/blast/decipher) for the INFO log line. The BLAST and
                DECIPHER outputs share this DADA2-shaped path; this parameter
                keeps the log truthful about which method produced the input.

        Returns:
            DataFrame in GBIF-compatible long format

        Raises:
            FileNotFoundError: If input file does not exist
            ValueError: If required columns are missing
        """
        input_path = Path(input_path)

        if not input_path.exists():
            raise FileNotFoundError(
                f"Taxonomy results file not found: {input_path}. This is the "
                f"per-marker taxonomy CSV produced by the taxonomy step (e.g. "
                f"outputs/<marker>_<method>.csv). Check the "
                f"path is correct; if you are resuming a pipeline, the file may "
                f"have been moved or deleted -- re-run the taxonomy step, or point "
                f"--input at the existing taxonomy CSV."
            )

        logger.info(f"Converting {_source_format} output to GBIF format: {input_path}")

        # Read CSV
        try:
            df = pd.read_csv(input_path)
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            raise ValueError(
                f"Could not read the taxonomy CSV '{input_path}': the file is "
                f"empty or is not valid CSV ({e}). Confirm the path points at the "
                f"taxonomy-step output CSV (not a FASTA, log, or binary), and that "
                f"the file is not truncated or zero-length -- a failed or "
                f"interrupted upstream taxonomy step can leave an empty file behind."
            ) from e

        # Remove X column if present (R index column)
        if "X" in df.columns:
            df = df.drop(columns=["X"])

        # Normalize BLAST/post-processor schema. Commits A-G produce a
        # `Sequence` column (capital S) and use the literal string "Unassigned"
        # for missing taxonomy. Map them back to the lowercase / NaN form the
        # rank-determination logic expects.
        if "Sequence" in df.columns and "sequence" not in df.columns:
            df = df.rename(columns={"Sequence": "sequence"})
        for col in TAXONOMIC_RANKS:
            if col in df.columns:
                df[col] = df[col].replace("Unassigned", pd.NA)

        # Define taxonomic columns to keep
        taxonomic_cols = [*TAXONOMIC_RANKS, "sequence"]

        # Check required columns exist
        missing_cols = [col for col in taxonomic_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(
                f"Taxonomy CSV '{input_path}' is missing required columns: "
                f"{missing_cols}. GBIF formatting needs all of: kingdom, phylum, "
                f"class, order, family, genus, species, sequence (a capital-S "
                f"'Sequence' is auto-mapped to 'sequence'). The usual cause is "
                f"pointing at the wrong file (e.g. the raw ASV count table instead "
                f"of the taxonomy-merged table) or a taxonomy method whose column "
                f"names were not normalised. Pass the taxonomy-step output for this "
                f"marker, or rename the columns to match the standard schema."
            )

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
        final_cols = [*TAXONOMIC_RANKS]
        if "taxon" in df_long.columns:
            final_cols.append("taxon")
        if "rank" in df_long.columns:
            final_cols.append("rank")
        final_cols.extend(["sequence", "nb_reads", "eventID"])
        # Preserve the contaminant annotation so create-gbif can surface contamination_flag.
        # _transform_to_long_format carries it through the melt; the column is only present
        # when taxonomy.contaminants was set (else create-gbif defaults the flag to False).
        if "is_contaminant_candidate" in df_long.columns:
            final_cols.append("is_contaminant_candidate")

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
        4. Drops ecotag-specific metadata columns: the named columns
           id, definition, count, scientific_name, plus any column whose
           name contains best_identity, best_match, match_count,
           species_list, or taxid (pattern-matched)
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
            raise FileNotFoundError(
                f"Ecotag results file not found: {input_path}. This is the "
                f"per-marker ecotag taxonomy CSV produced by the taxonomy step "
                f"(e.g. outputs/<marker>_ecotag.csv). Check "
                f"the path is correct; if you are resuming a pipeline, the file may "
                f"have been moved or deleted -- re-run the ecotag taxonomy step, or "
                f"point --input at the existing ecotag CSV."
            )

        logger.info(f"Converting ecotag output to GBIF format: {input_path}")

        # Read CSV
        try:
            df = pd.read_csv(input_path)
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            raise ValueError(
                f"Could not read the ecotag taxonomy CSV '{input_path}': the file "
                f"is empty or is not valid CSV ({e}). Confirm the path points at "
                f"the ecotag taxonomy-step output CSV (not a FASTA, log, or "
                f"binary), and that the file is not truncated or zero-length -- a "
                f"failed or interrupted upstream taxonomy step can leave an empty "
                f"file behind."
            ) from e

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
        taxonomic_cols = [*TAXONOMIC_RANKS, "sequence"]

        # Check required columns exist
        missing_cols = [col for col in taxonomic_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(
                f"Ecotag CSV '{input_path}' is missing required GBIF columns after "
                f"renaming: {missing_cols}. from_ecotag maps "
                f"family_name/genus_name/species_name/order_name to "
                f"family/genus/species/order and auto-adds kingdom/phylum/class, so "
                f"the real gap is among: order, family, genus, species, sequence. A "
                f"missing 'sequence' usually means this is the raw obitab/ecotag "
                f"table that was not yet linked with the abundance table; run the "
                f"ecotag link-with-abundance step first, then format with "
                f"--format ecotag."
            )

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
        final_cols = [*TAXONOMIC_RANKS]
        if "taxon" in df_long.columns:
            final_cols.append("taxon")
        if "rank" in df_long.columns:
            final_cols.append("rank")
        final_cols.extend(["sequence", "nb_reads", "eventID"])
        # Preserve the contaminant annotation so create-gbif can surface contamination_flag
        # (same as from_dada2_rdp); the filter below drops it harmlessly when absent.
        final_cols.append("is_contaminant_candidate")

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

        BLAST taxonomy output already shares the DADA2 wide-table schema (the
        post-processor emits the same columns), so this simply delegates to
        ``from_dada2_rdp`` with a ``blast`` source label for the log.

        Args:
            input_path: Path to BLAST output CSV file
            output_path: Optional path to output GBIF CSV file
            add_rank: Whether to add 'rank' column (default: True)
            add_taxon: Whether to add 'taxon' column (default: True)

        Returns:
            DataFrame in GBIF-compatible long format

        Raises:
            FileNotFoundError: If the input file does not exist.
            ValueError: If the CSV is empty/invalid or required columns are
                missing.
        """
        # BLAST output should already be in a similar format to DADA2
        return self.from_dada2_rdp(
            input_path, output_path, add_rank, add_taxon, _source_format="blast"
        )

    def from_decipher(
        self,
        input_path: Union[str, Path],
        output_path: Optional[Union[str, Path]] = None,
        add_rank: bool = True,
        add_taxon: bool = True,
    ) -> pd.DataFrame:
        """
        Convert DECIPHER output to GBIF format.

        DECIPHER (IdTaxa) taxonomy output shares the DADA2 wide-table schema, so
        this delegates to ``from_dada2_rdp`` with a ``decipher`` source label for
        the log.

        Args:
            input_path: Path to DECIPHER output CSV file
            output_path: Optional path to output GBIF CSV file
            add_rank: Whether to add 'rank' column (default: True)
            add_taxon: Whether to add 'taxon' column (default: True)

        Returns:
            DataFrame in GBIF-compatible long format

        Raises:
            FileNotFoundError: If the input file does not exist.
            ValueError: If the CSV is empty/invalid or required columns are
                missing.
        """
        # DECIPHER output should be similar to DADA2
        return self.from_dada2_rdp(
            input_path, output_path, add_rank, add_taxon, _source_format="decipher"
        )
