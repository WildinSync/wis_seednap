"""BLAST-based taxonomic assignment: command execution, LCA resolution, and output formatting.

This module provides:
- BlastRunner: subprocess wrapper for makeblastdb/blastn
- BlastOutputFormatter: adds phylogenetic info from reference DB headers
- BlastPhyloFilter: filters hits by percent identity thresholds
- BlastLCAResolver: resolves ambiguous hits using Lowest Common Ancestor
- BlastTaxonomicAssigner: end-to-end BLAST taxonomy pipeline

Author: Théophile Sanchez (original), refactored for seednap v0.1.0
"""

import logging
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd

from seednap.utils.sequences import fasta_to_df
from seednap.utils.subprocess import run_subprocess

logger = logging.getLogger(__name__)

# CRABS reference FASTA headers write the literal string "NA" where a taxonomic rank is unknown
# (verified in the 2025 DBs: e.g. ;Actinopteri;NA;Centropomidae; on Lates records -- teleo 941
# headers, mamm07 3056). These placeholders are NOT taxa. They are normalized to None when the
# header lineage is parsed (BlastOutputFormatter), so neither LCA resolver treats "NA" as a real
# value -- which would (a) over-collapse a call when one in-band hit carries "NA" at a rank where
# the others agree, and (b) leak a taxon literally named "NA" into the GBIF export. Matched
# case-insensitively after strip(); no real taxon is named "na"/"nan".
MISSING_RANK_SENTINELS = ("", "na", "nan")


class BlastError(Exception):
    """Exception raised for BLAST errors."""

    pass


# Backwards-compatible alias
BlastDatabaseError = BlastError


# ---------------------------------------------------------------------------
# Command runner
# ---------------------------------------------------------------------------


class BlastRunner:
    """Run BLAST commands (makeblastdb, blastn) via subprocess."""

    VALID_TASKS = ("megablast", "blastn", "dc-megablast", "blastn-short")

    def __init__(
        self,
        perc_identity: float = 80.0,
        qcov_hsp_perc: float = 80.0,
        evalue: float = 1e-25,
        max_target_seqs: int = 5,
        task: str = "megablast",
    ):
        """
        Initialize BLAST runner with search parameters.

        Args:
            perc_identity: Minimum percent identity for hits (default: 80.0)
            qcov_hsp_perc: Minimum query coverage per HSP (default: 80.0)
            evalue: Maximum e-value for hits (default: 1e-25)
            max_target_seqs: Maximum number of target sequences to keep (default: 5)
            task: blastn task. 'megablast' (word_size 28) is the right default for
                short, high-identity vertebrate amplicons; 'blastn' (word_size 11)
                is more sensitive for divergent matches. Default: 'megablast'.
        """
        if task not in self.VALID_TASKS:
            raise ValueError(
                f"Invalid blastn task '{task}'. Must be one of: {self.VALID_TASKS}"
            )
        self.perc_identity = perc_identity
        self.qcov_hsp_perc = qcov_hsp_perc
        self.evalue = evalue
        self.max_target_seqs = max_target_seqs
        self.task = task

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
            BlastError: If makeblastdb command fails
        """
        fasta_path = Path(fasta_path)

        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        logger.info(f"Creating BLAST database for {fasta_path}")

        cmd = ["makeblastdb", "-dbtype", "nucl", "-in", str(fasta_path)]
        run_subprocess(cmd, timeout=600, error_class=BlastError)
        logger.info("BLAST database created successfully")

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
            BlastError: If blastn command fails
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
            f"Parameters: task={self.task}, pident={self.perc_identity}, "
            f"qcov={self.qcov_hsp_perc}, evalue={self.evalue}, "
            f"max_targets={self.max_target_seqs}"
        )

        # Build blastn command
        cmd = [
            "blastn",
            "-task",
            self.task,
            "-query",
            str(query_fasta),
            "-db",
            str(db_fasta),
            "-out",
            str(output_tsv),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-perc_identity",
            str(self.perc_identity),
            "-qcov_hsp_perc",
            str(self.qcov_hsp_perc),
            "-evalue",
            str(self.evalue),
            "-max_target_seqs",
            str(self.max_target_seqs),
        ]

        run_subprocess(cmd, timeout=3600, error_class=BlastError)
        logger.info(f"BLAST search completed, output saved to {output_tsv}")

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
            BlastError: If BLAST commands fail
        """
        output_dir = Path(output_dir)
        output_tsv = output_dir / f"{marker}_blastn_output.tsv"

        # Run BLAST search
        self.run_blastn(query_fasta, db_fasta, output_tsv)

        return output_tsv


# ---------------------------------------------------------------------------
# Output formatting and phylogenetic extraction
# ---------------------------------------------------------------------------


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

        # Empty BLAST output -> no hits at all. Return empty DF with full schema so the
        # downstream left-merge produces an all-Unassigned final table.
        if blast_tsv.stat().st_size == 0:
            logger.warning(
                f"BLAST output {blast_tsv} is empty: no hits at all. All OTUs will be "
                f"marked 'Unassigned'."
            )
            return pd.DataFrame(
                columns=self.BLAST_COLUMNS + self.TAXONOMIC_RANKS + ["blast_rank"]
            )

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
            header_parts = header.replace("\n", "").split("\t")
            if len(header_parts) < 2:
                raise ValueError(
                    f"Malformed reference header for '{seq_id}': expected tab-separated "
                    f"ID and taxonomy, got: {header[:80]}"
                )
            phylo_string = header_parts[1]
            phylo_values = phylo_string.split(";")
            if len(phylo_values) < len(self.TAXONOMIC_RANKS):
                raise ValueError(
                    f"Incomplete taxonomy for '{seq_id}': expected {len(self.TAXONOMIC_RANKS)} "
                    f"ranks (;-separated), got {len(phylo_values)}: {phylo_string[:80]}"
                )

            # Assign to columns. The CRABS missing-rank sentinel ("NA") is normalized to None
            # here, at the single point where the DB-format lineage is parsed, so every
            # downstream consumer (both LCA resolvers, the threshold cascade, the export) sees a
            # genuine missing rank rather than a taxon literally named "NA" (CLAUDE.md sec.4).
            phylo = dict(zip(self.TAXONOMIC_RANKS, phylo_values))
            for rank in self.TAXONOMIC_RANKS:
                value = phylo[rank]
                df.at[i, rank] = (
                    None
                    if str(value).strip().lower() in MISSING_RANK_SENTINELS
                    else value
                )

        # Add blast rank (1 = best hit, 2 = second best, etc.)
        df["blast_rank"] = df.groupby("qseqid").cumcount() + 1

        # Save if output path provided
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, sep="\t", header=True, index=False)

        return df


# ---------------------------------------------------------------------------
# Phylogenetic filtering
# ---------------------------------------------------------------------------


class BlastPhyloFilter:
    """Filter BLAST hits by percent identity thresholds for each taxonomic rank.

    Applies a *cascade* nulling rule: when a hit's percent identity falls below the
    threshold for rank R, R is nulled and so are all finer ranks below R. This avoids
    the orphan-rank problem (e.g. order populated but family/genus/species nulled
    individually with no consistency check).

    Cascade ranks (coarse -> fine): class -> order -> family -> genus -> species.
    Kingdom and phylum are never auto-nulled by threshold; they pass through if the
    hit cleared the absolute BLAST `perc_identity` cutoff.
    """

    # Cascade ranks, ordered coarse to fine
    CASCADE_RANKS = ["class", "order", "family", "genus", "species"]

    def __init__(
        self,
        threshold_species: float = 98.0,
        threshold_genus: float = 96.0,
        threshold_family: float = 86.5,
        threshold_order: float = 80.0,
        threshold_class: float = 70.0,
    ):
        """
        Initialize phylogenetic filter with thresholds.

        Args:
            threshold_species: Minimum percent identity for species-level assignment
            threshold_genus: Minimum percent identity for genus-level assignment
            threshold_family: Minimum percent identity for family-level assignment
            threshold_order: Minimum percent identity for order-level assignment
            threshold_class: Minimum percent identity for class-level assignment
        """
        self.thresholds = {
            "class": threshold_class,
            "order": threshold_order,
            "family": threshold_family,
            "genus": threshold_genus,
            "species": threshold_species,
        }

    def filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Cascade-null phylogenetic assignments below each rank's identity threshold.

        For rank R with threshold T, rows where pident < T have R and every finer
        rank set to None.

        Args:
            df: DataFrame with 'pident' and taxonomic rank columns

        Returns:
            DataFrame with cascade-nulled taxonomy
        """
        df = df.copy()
        if len(df) == 0:
            return df

        pident = pd.to_numeric(df["pident"], errors="coerce")
        for rank, threshold in self.thresholds.items():
            if rank not in self.CASCADE_RANKS:
                continue
            mask = pident < float(threshold)
            if not mask.any():
                continue
            rank_idx = self.CASCADE_RANKS.index(rank)
            for finer in self.CASCADE_RANKS[rank_idx:]:
                df.loc[mask, finer] = None
        return df


# ---------------------------------------------------------------------------
# LCA resolution
# ---------------------------------------------------------------------------


class BlastLCAResolver:
    """Resolve ambiguous BLAST hits using LCA (Lowest Common Ancestor).

    Implements a MEGAN-LR-style top-bitscore band: any hit whose bitscore is within
    `top_bitscore_pct` percent of the best bitscore for a query is considered when
    deciding the LCA. This is more robust than collapsing only over exact bitscore
    ties, which is brittle to near-duplicate references.

    When in-band hits all agree on every taxonomic rank, the best hit is kept as-is.
    When they disagree at any rank, a synthetic combined row is produced with the
    disagreed ranks (and every finer rank, by cascade) set to None. The combined
    row's pident reports the best (max) pident in the band.
    """

    TAXONOMIC_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    def __init__(self, top_bitscore_pct: float = 10.0, lca_pident_delta: float = 1.0):
        """
        Initialize LCA resolver.

        Args:
            top_bitscore_pct: Hits with bitscore within this percent of the best
                are included in LCA resolution. 0 = exact ties only. Default 10.0
                follows MEGAN-LR's `topPercent` default.
            lca_pident_delta: An in-band hit is additionally required to be within
                this many percent-identity points of the best in-band hit. Default
                1.0 (the eDNAFlow "diff 1" convention). On a short marker a single
                mismatch can leave a phylogenetically distant hit (e.g. a 98.6%
                worm next to 100% Bos hits) inside the bitscore band; without this
                floor that one hit collapses the LCA all the way to kingdom. 0
                disables the floor (bitscore band only).
        """
        if top_bitscore_pct < 0 or top_bitscore_pct > 100:
            raise ValueError(
                f"top_bitscore_pct must be in [0, 100]; got {top_bitscore_pct}"
            )
        if lca_pident_delta < 0 or lca_pident_delta > 100:
            raise ValueError(
                f"lca_pident_delta must be in [0, 100]; got {lca_pident_delta}"
            )
        self.top_bitscore_pct = float(top_bitscore_pct)
        self.lca_pident_delta = float(lca_pident_delta)

    def resolve_ambiguous_hits(self, group: pd.DataFrame) -> pd.DataFrame:
        """
        Resolve ambiguous hits within the top-bitscore band using LCA.

        Args:
            group: DataFrame of BLAST hits for a single query sequence

        Returns:
            DataFrame with `keep_for_analysis` set: at most one row per query is
            True (either a single in-band hit, the best of an agreeing in-band
            cohort, or a synthetic LCA-combined row).
        """
        group = group.reset_index(drop=True).copy()

        if len(group) == 0:
            group["keep_for_analysis"] = False
            return group
        if len(group) == 1:
            group["keep_for_analysis"] = True
            return group

        # Bitscore band
        bitscores = pd.to_numeric(group["bitscore"], errors="coerce")
        best_bitscore = bitscores.max()
        if pd.isna(best_bitscore) or best_bitscore <= 0:
            group["keep_for_analysis"] = False
            return group
        threshold = best_bitscore * (1.0 - self.top_bitscore_pct / 100.0)
        in_bitscore = bitscores >= threshold
        # pident floor: exclude in-bitscore hits whose identity is more than
        # lca_pident_delta below the best in-band identity, so one near-identity
        # off-target hit cannot collapse the LCA of an otherwise-agreeing cohort.
        pidents = pd.to_numeric(group["pident"], errors="coerce")
        best_pident = pidents[in_bitscore].max()
        if pd.isna(best_pident) or self.lca_pident_delta <= 0:
            in_band_mask = in_bitscore
        else:
            in_band_mask = in_bitscore & (pidents >= best_pident - self.lca_pident_delta)
        ambiguous_hits = group[in_band_mask]

        # Helper: keep exactly one row by label
        def _keep_one(group_df: pd.DataFrame, keep_label: int) -> pd.DataFrame:
            keep_mask = pd.Series(False, index=group_df.index)
            keep_mask.loc[keep_label] = True
            group_df["keep_for_analysis"] = keep_mask
            return group_df

        if len(ambiguous_hits) <= 1:
            # No ambiguity - just keep the single in-band hit (the best one)
            return _keep_one(group, ambiguous_hits.index[0])

        # Multiple in-band hits: do they all agree on every taxonomic rank?
        same_phylo = all(
            ambiguous_hits[col].dropna().nunique() < 2
            for col in self.TAXONOMIC_RANKS
        )
        if same_phylo:
            # All agree - keep the best (first/max-bitscore) ambiguous hit
            best_idx = bitscores[in_band_mask].idxmax()
            return _keep_one(group, best_idx)

        # Disagreement: build LCA combined row from the agreed-on columns
        combined_row = ambiguous_hits.iloc[[0]].copy()
        for col in ambiguous_hits.columns:
            non_null = ambiguous_hits[col].dropna()
            if non_null.nunique() > 1:
                combined_row[col] = None
        # Best pident in the band best represents the LCA's identity score
        combined_row["pident"] = pd.to_numeric(
            ambiguous_hits["pident"], errors="coerce"
        ).max()
        combined_row["bitscore"] = best_bitscore

        # Cascade: if a coarse rank is None, every finer rank is also None
        for i, rank in enumerate(self.TAXONOMIC_RANKS):
            if combined_row[rank].isna().all():
                for finer_rank in self.TAXONOMIC_RANKS[i + 1:]:
                    combined_row[finer_rank] = None
                break

        group["keep_for_analysis"] = False
        combined_row["keep_for_analysis"] = True
        # The combined row has columns set to None where in-band hits disagreed.
        # pandas 2.x emits a FutureWarning about all-NA column concat behavior;
        # the current behavior is what we want, so suppress the noise.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FutureWarning)
            return pd.concat([group, combined_row], ignore_index=True)


class CollapsedTaxonomyLCAResolver:
    """eDNAFlow / OceanOmics "collapsed taxonomy" LCA, reimplemented on header lineages.

    Per OTU: keep hits with ``pident >= lca_pid`` (the hard identity floor; query coverage is
    already enforced at the blastn step via ``-qcov_hsp_perc``). Among those, keep the
    top-identity window ``[best_pident - lca_diff, best_pident]``. Collapse the window's
    lineages to their lowest common ancestor: walk kingdom->species and keep a rank while all
    windowed hits share a single non-empty value; at the first rank where they disagree, null
    that rank and every finer rank.

    This is the algorithm from OceanOmics-amplicon-nf's ``LCA`` process (eDNAFlow, Mahsa
    Mousavi-Derazmahalleh), adapted to SeeDNAP: the lineage comes from the reference FASTA
    headers (SeeDNAP DBs are CRABS-formatted), so it needs **no NCBI taxids and no taxdump**
    and runs fully offline. Upstream defects are fixed: the identity-difference comparison is
    numeric (upstream compares formatted strings), taxa absent from the reference raise at the
    formatter (upstream silently drops taxids missing from the dump), and no hit is dropped
    without the count being reported.
    """

    TAXONOMIC_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    def __init__(self, lca_pid: float = 90.0, lca_diff: float = 1.0):
        if not (0 <= lca_pid <= 100):
            raise ValueError(f"lca_pid must be in [0, 100]; got {lca_pid}")
        if not (0 <= lca_diff <= 100):
            raise ValueError(f"lca_diff must be in [0, 100]; got {lca_diff}")
        self.lca_pid = float(lca_pid)
        self.lca_diff = float(lca_diff)

    @staticmethod
    def _distinct(series: pd.Series) -> List[str]:
        """Distinct real-taxon values of a lineage column. The formatter already normalizes the
        CRABS "NA" sentinel to None; this also drops it defensively (frames built directly in
        tests bypass the formatter), using the module-level MISSING_RANK_SENTINELS."""
        return [
            v
            for v in series.dropna().unique()
            if str(v).strip().lower() not in MISSING_RANK_SENTINELS
        ]

    def resolve_ambiguous_hits(self, group: pd.DataFrame) -> pd.DataFrame:
        """Collapse one OTU's hits to their LCA over the top-identity window. Returns the
        group plus a synthetic combined row with ``keep_for_analysis=True`` (the LCA)."""
        group = group.reset_index(drop=True).copy()
        if len(group) == 0:
            group["keep_for_analysis"] = False
            return group

        pident = pd.to_numeric(group["pident"], errors="coerce")
        passed = pident >= self.lca_pid
        combined = group.iloc[[0]].copy()

        if not passed.any():
            # Nothing clears the identity floor -> Unassigned (null every rank).
            for rank in self.TAXONOMIC_RANKS:
                combined[rank] = None
            combined["pident"] = float(pident.max()) if pident.notna().any() else None
        else:
            best = float(pident[passed].max())
            window = group[passed & (pident >= best - self.lca_diff)]
            combined = window.iloc[[0]].copy()
            disagreed = False
            for i, rank in enumerate(self.TAXONOMIC_RANKS):
                if disagreed:
                    combined[rank] = None
                    continue
                vals = self._distinct(window[rank])
                if len(vals) > 1:
                    disagreed = True
                    combined[rank] = None
                else:
                    combined[rank] = vals[0] if vals else None
            combined["pident"] = best

        group["keep_for_analysis"] = False
        combined["keep_for_analysis"] = True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FutureWarning)
            return pd.concat([group, combined], ignore_index=True)


# ---------------------------------------------------------------------------
# End-to-end taxonomic assigner
# ---------------------------------------------------------------------------


class BlastTaxonomicAssigner:
    """Complete BLAST-based taxonomic assignment pipeline."""

    UNASSIGNED_LABEL = "Unassigned"
    CONTAMINANT_FLAG_COL = "is_contaminant_candidate"

    def __init__(
        self,
        reference_fasta: Union[str, Path],
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
        contaminants: Optional[List[str]] = None,
    ):
        """
        Initialize BLAST taxonomic assigner.

        Args:
            reference_fasta: Path to reference database FASTA file
            threshold_species: Minimum percent identity for species-level assignment
            threshold_genus: Minimum percent identity for genus-level assignment
            threshold_family: Minimum percent identity for family-level assignment
            threshold_order: Minimum percent identity for order-level assignment
            threshold_class: Minimum percent identity for class-level assignment
            top_bitscore_pct: LCA bitscore band as percent of best hit
                (MEGAN-LR style; default 10.0)
            contaminants: List of species names (CRABS underscore format) to flag as
                candidate contaminants. Rows are NEVER deleted; only flagged.
        """
        self.formatter = BlastOutputFormatter(reference_fasta)
        self.filter = BlastPhyloFilter(
            threshold_species=threshold_species,
            threshold_genus=threshold_genus,
            threshold_family=threshold_family,
            threshold_order=threshold_order,
            threshold_class=threshold_class,
        )
        self.lca_algorithm = lca_algorithm
        self.lca_pid = lca_pid
        self.lca_resolver = self._make_lca_resolver(
            lca_algorithm,
            top_bitscore_pct=top_bitscore_pct, lca_pident_delta=lca_pident_delta,
            lca_pid=lca_pid, lca_diff=lca_diff,
        )
        self.contaminants = list(contaminants) if contaminants else []

    @staticmethod
    def _make_lca_resolver(
        algorithm: str, *, top_bitscore_pct: float, lca_pident_delta: float,
        lca_pid: float = 90.0, lca_diff: float = 1.0,
    ) -> Union["BlastLCAResolver", "CollapsedTaxonomyLCAResolver"]:
        """Resolver factory keyed on lca_algorithm.

        'cascade' (default) = the header-derived per-rank/MEGAN-LR resolver. 'collapsed_taxonomy'
        = the eDNAFlow/OceanOmics identity-window collapse-to-LCA (also header-based, offline).
        'fishbase_tiered' (Fishbase->WoRMS->NCBI) is not implemented (fish-specific, needs the
        bundled WoRMS file + staged Fishbase parquet) and raises until provisioned.
        """
        if algorithm == "cascade":
            return BlastLCAResolver(
                top_bitscore_pct=top_bitscore_pct, lca_pident_delta=lca_pident_delta
            )
        if algorithm == "collapsed_taxonomy":
            return CollapsedTaxonomyLCAResolver(lca_pid=lca_pid, lca_diff=lca_diff)
        if algorithm == "fishbase_tiered":
            raise NotImplementedError(
                "lca_algorithm='fishbase_tiered' is not implemented yet: it needs the bundled "
                "WoRMS file + a staged Fishbase parquet and is fish-specific. Use 'cascade' or "
                "'collapsed_taxonomy'."
            )
        raise ValueError(f"unknown lca_algorithm {algorithm!r}")

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
        2. Cascade-filters hits by percent identity thresholds
        3. Resolves ambiguous hits using LCA (top-bitscore band)
        4. LEFT-merges taxonomy onto the ASV count table so every OTU survives;
           OTUs with no BLAST hits are explicitly marked 'Unassigned'.
        5. Outputs final table with taxonomy, sequences, and per-sample abundances

        Args:
            blast_tsv: Path to BLAST output TSV file
            asv_count_csv: Path to ASV count table CSV (sequences as rows, samples as columns)
            asv_fasta: Path to ASV sequences FASTA file
            output_path: Optional path to save final table

        Returns:
            DataFrame with taxonomic assignments and ASV counts. The row count
            equals the number of OTUs in the abundance table (no silent drops).

        Raises:
            FileNotFoundError: If input files do not exist
        """
        phylo_cols = BlastLCAResolver.TAXONOMIC_RANKS

        # 1. Format BLAST output (handles empty BLAST gracefully)
        formatted = self.formatter.format_blast_output(blast_tsv)
        formatted = formatted.replace("None", None)

        # 2. Cascade-filter and 3. LCA-resolve, only if we have any hits
        if len(formatted) == 0:
            result = pd.DataFrame(columns=["ASV_ID", "pident"] + phylo_cols)
        else:
            # cascade applies per-rank identity thresholds before the LCA; collapsed_taxonomy
            # does its own identity-floor + window inside the resolver, so it skips the cascade.
            if self.lca_algorithm == "collapsed_taxonomy":
                filtered = formatted.copy()
            else:
                filtered = self.filter.filter(formatted)
            filtered["keep_for_analysis"] = filtered["blast_rank"] == 1
            filtered = (
                filtered.groupby("qseqid", group_keys=True)
                .apply(self.lca_resolver.resolve_ambiguous_hits, include_groups=False)
                .reset_index(level="qseqid")
                .reset_index(drop=True)
            )
            filtered = filtered[filtered["keep_for_analysis"] == True]  # noqa: E712
            result = (
                filtered[["qseqid", "pident"] + phylo_cols]
                .rename(columns={"qseqid": "ASV_ID"})
                .reset_index(drop=True)
            )
            # D1 guard: an OTU whose best hit is well above the genus identity threshold
            # but which still collapsed below genus (an LCA across disagreeing in-band
            # hits) is surfaced loudly, so a confident call is never silently lost to an
            # over-broad band (CLAUDE.md sec.4).
            genus_thr = getattr(self.filter, "threshold_genus", 96.0)
            res_pident = pd.to_numeric(result["pident"], errors="coerce")
            collapsed = result[(res_pident >= genus_thr) & (result["genus"].isna())]
            if len(collapsed) > 0:
                ids = list(collapsed["ASV_ID"].astype(str))
                shown = ids[:15]
                more = "" if len(ids) <= 15 else f", +{len(ids) - 15} more"
                logger.warning(
                    f"[WARN] blast_taxonomy: {len(collapsed)} OTU(s) with best pident >= the genus "
                    f"threshold ({genus_thr}) collapsed below genus via LCA -- their in-band hits "
                    f"span taxa (likely short-fragment cross-taxa matches in the reference DB). "
                    f"Reported, not silently assigned: {shown}{more}"
                )

            # collapsed_taxonomy floor guard: an OTU that HAD BLAST hits but none cleared
            # lca_pid is marked Unassigned (all ranks null) yet retains its best below-floor
            # pident -- it would otherwise be indistinguishable from a genuine no-hit OTU,
            # which IS warned below (n_unassigned). Surface the dropped-by-floor count so the
            # two cases are distinguishable and no hit is discarded silently (CLAUDE.md sec.4).
            if self.lca_algorithm == "collapsed_taxonomy":
                all_null = result[phylo_cols].isna().all(axis=1)
                below_floor = result[all_null & res_pident.notna()]
                if len(below_floor) > 0:
                    lo, hi = float(res_pident[below_floor.index].min()), float(
                        res_pident[below_floor.index].max()
                    )
                    logger.warning(
                        f"[WARN] collapsed_taxonomy LCA: {len(below_floor)} OTU(s) had BLAST "
                        f"hit(s) but none cleared lca_pid={self.lca_pid} "
                        f"(best pident {lo:.1f}-{hi:.1f}); marked Unassigned. "
                        f"expected=>= floor, got=all hits below floor, fallback=Unassigned"
                    )

        # 4. Load ASV count table (sequences as rows, samples as columns)
        asv_count = pd.read_csv(asv_count_csv, sep=",", index_col=0)
        sample_cols = list(asv_count.columns)

        # Attach ASV_ID via the FASTA
        asv_sequences = fasta_to_df(asv_fasta)
        asv_sequences = asv_sequences.rename(columns={"id": "ASV_ID", "sequence": "Sequence"})
        asv_count = pd.merge(
            asv_count, asv_sequences, how="inner", left_index=True, right_on="Sequence"
        )

        # Diagnostic: how many OTUs got no BLAST hits at all?
        otus_with_hits = set(result["ASV_ID"].dropna().unique()) if len(result) > 0 else set()
        n_unassigned = (~asv_count["ASV_ID"].isin(otus_with_hits)).sum()
        if n_unassigned > 0:
            logger.warning(
                f"{n_unassigned} of {len(asv_count)} OTUs had no BLAST hits and "
                f"will be marked '{self.UNASSIGNED_LABEL}' in the output."
            )

        # 5. LEFT-merge so every OTU in the abundance table reaches the output.
        # Without this, OTUs without BLAST hits were silently dropped.
        final_table = pd.merge(asv_count, result, how="left", on="ASV_ID")

        # Fill missing taxonomy with 'Unassigned' (rank columns only). pident stays NaN
        # for genuine no-hit rows so downstream can distinguish them.
        for rank in phylo_cols:
            final_table[rank] = final_table[rank].fillna(self.UNASSIGNED_LABEL)

        # Flag candidate contaminants by species match. Never delete rows.
        if self.contaminants:
            contam_set = set(self.contaminants)
            is_contam = final_table["species"].astype(str).isin(contam_set)
            final_table[self.CONTAMINANT_FLAG_COL] = is_contam
            n_flagged = int(is_contam.sum())
            if n_flagged > 0:
                breakdown = (
                    final_table.loc[is_contam, "species"].value_counts().to_dict()
                )
                logger.warning(
                    f"Flagged {n_flagged} OTUs as candidate contaminants "
                    f"(by species match): {breakdown}"
                )
        else:
            final_table[self.CONTAMINANT_FLAG_COL] = False

        # Sort by ASV number for deterministic output
        final_table["asv_num"] = (
            final_table["ASV_ID"].astype(str).str.extract(r"(\d+)").astype("Int64")
        )
        final_table = final_table.sort_values("asv_num").drop(columns="asv_num").reset_index(
            drop=True
        )

        # Stable column order: ASV_ID, pident, taxonomy, contaminant flag, samples, Sequence
        ordered = (
            ["ASV_ID", "pident"] + phylo_cols + [self.CONTAMINANT_FLAG_COL]
            + sample_cols + ["Sequence"]
        )
        final_table = final_table[[c for c in ordered if c in final_table.columns]]

        # Save if requested
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            final_table.to_csv(output_path, header=True, index=False)

        return final_table
