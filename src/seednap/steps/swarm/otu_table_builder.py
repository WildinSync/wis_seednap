"""Build the OTU contingency table from SWARM clustering outputs.

Final step of the SWARM OTU path (``steps/swarm/``), run after vsearch
dereplication, SWARM clustering, and chimera detection. Parses the SWARM
output files (representatives, stats, swarm membership, chimera detection)
and the per-sample FASTA files, then assembles a full OTU contingency table
(one row per OTU, one read-count column per sample) and the DADA2-compatible
normalized outputs (query FASTA + abundance CSV) that feed downstream
taxonomy assignment.

An OTU (Operational Taxonomic Unit) is a cluster of near-identical amplicon
sequences standing in for one putative taxon; the contingency table records
how many reads of each OTU were seen in each sample. The construction follows
the OTU contingency table approach of the published SWARM amplicon clustering
pipeline.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd

logger = logging.getLogger(__name__)

# Column names reserved for OTU metadata in the full contingency table.
# A sample whose name collides with one of these would overwrite the metadata
# value in build() and then be dropped from the abundance matrix by
# to_taxonomy_input() (sample_cols is a set-difference against this list plus
# "sequence"). That is a silent per-sample data loss, so build() refuses it.
_RESERVED_METADATA_COLS = (
    "OTU", "total", "cloud", "amplicon", "length",
    "abundance", "chimera", "spread", "sequence",
)


class OtuTableBuilder:
    """
    Build the OTU contingency table from SWARM clustering results.

    An OTU (Operational Taxonomic Unit) is a cluster of near-identical
    amplicon sequences treated as one putative taxon. After SWARM groups
    the dereplicated reads into OTUs, this class stitches the various SWARM
    output files back together with the per-sample abundance data to produce
    one wide table: one row per OTU, OTU metadata columns plus one read-count
    column per sample. That table is the SWARM-path equivalent of a DADA2 ASV
    table and feeds downstream taxonomy assignment.
    """

    def build(
        self,
        representatives_fasta: Union[str, Path],
        stats_file: Union[str, Path],
        swarm_file: Union[str, Path],
        uchime_file: Optional[Union[str, Path]],
        sample_fastas: Sequence[Union[str, Path]],
    ) -> pd.DataFrame:
        """
        Build the full OTU contingency table from the SWARM outputs.

        Joins the four SWARM/vsearch artefacts (representative sequences,
        per-OTU statistics, cluster membership, chimera calls) with the
        per-sample dereplicated FASTAs to assemble one row per OTU, ranked
        by total abundance, with a read count for every sample.

        Args:
            representatives_fasta: Path to the abundance-sorted seed
                (representative) sequences FASTA, one sequence per OTU.
            stats_file: Path to the SWARM ``-s`` statistics file (per-OTU
                mass, cloud size, seed ID, seed abundance).
            swarm_file: Path to the SWARM ``-o`` cluster membership file
                (seed amplicon followed by its member amplicons).
            uchime_file: Path to the vsearch UCHIME chimera-detection table,
                or None to skip chimera annotation (all OTUs get status "NA").
            sample_fastas: Per-sample dereplicated FASTA files carrying the
                ``;size=N`` abundance of each amplicon in each sample.

        Returns:
            DataFrame with one row per OTU, sorted by decreasing total
            abundance. Columns: OTU (1-based rank), total (summed read mass),
            cloud (count of unique amplicons in the cluster), amplicon (seed
            ID), length (seed sequence length in bp), abundance (seed read
            count), chimera (Y/N/?/NA), spread (number of samples the OTU
            occurs in), sequence (seed nucleotide sequence), and one integer
            column per sample holding that sample's read count for the OTU.

        Raises:
            ValueError: If any sample name collides with a reserved OTU-table
                metadata column (see ``_RESERVED_METADATA_COLS``); such a
                collision would overwrite metadata and silently drop the
                sample from the abundance matrix, so it fails loudly instead.
            FileNotFoundError: If any input file path does not exist.
        """
        representatives = self._parse_representatives(representatives_fasta)
        stats, sorted_seeds, seeds = self._parse_stats(stats_file)
        swarms = self._parse_swarms(swarm_file)
        uchime = self._parse_uchime(uchime_file) if uchime_file else {}
        amplicons2samples, samples = self._parse_sample_fastas(sample_fastas)

        # Guard against a sample name colliding with a reserved metadata column.
        # Such a sample would overwrite the metadata value here and then be
        # silently excluded from the abundance matrix in to_taxonomy_input()
        # (which derives sample_cols by set-difference against the metadata
        # column names). Fail loudly with the offending name instead.
        reserved = set(_RESERVED_METADATA_COLS)
        colliding = [s for s in samples if s in reserved]
        if colliding:
            raise ValueError(
                f"Sample name(s) collide with reserved OTU-table metadata "
                f"columns: {sorted(colliding)}. These names "
                f"({sorted(reserved)}) are used for OTU metadata; a sample "
                f"named like one of them would overwrite the metadata value "
                f"and be dropped from the abundance matrix passed to taxonomy. "
                f"Rename the offending sample(s) upstream."
            )

        rows = []
        for i, (seed, mass) in enumerate(sorted_seeds, start=1):
            sequence = representatives.get(seed)
            if not sequence:
                logger.warning(f"Seed {seed} not found in representatives FASTA, skipping OTU {i}")
                continue
            seed_abundance, cloud = seeds.get(seed, (0, 0))
            chimera_status = uchime.get(seed, "NA")

            # Sum per-sample abundances for all amplicons in this OTU
            occurrences = {sample: 0 for sample in samples}
            for amplicon in swarms.get(seed, []):
                if amplicon in amplicons2samples:
                    for sample, abundance in amplicons2samples[amplicon].items():
                        occurrences[sample] = occurrences.get(sample, 0) + abundance

            spread = sum(1 for v in occurrences.values() if v > 0)

            row = {
                "OTU": i,
                "total": mass,
                "cloud": cloud,
                "amplicon": seed,
                "length": len(sequence),
                "abundance": seed_abundance,
                "chimera": chimera_status,
                "spread": spread,
                "sequence": sequence,
            }
            for sample in samples:
                row[sample] = occurrences.get(sample, 0)

            rows.append(row)

        df = pd.DataFrame(rows)
        logger.info(f"Built OTU table: {len(df)} OTUs across {len(samples)} samples")
        return df

    def to_taxonomy_input(
        self,
        otu_table: pd.DataFrame,
        query_fasta_path: Union[str, Path],
        abundance_csv_path: Union[str, Path],
    ) -> Tuple[Path, Path]:
        """
        Write normalized outputs compatible with the taxonomy assignment step.

        Chimeras are artefactual sequences formed when two distinct templates
        fuse during PCR; keeping them would inflate the OTU count with false
        taxa. This drops only OTUs flagged as definite chimeras (chimera ==
        "Y"). Borderline ("?") and unscored ("NA") OTUs are intentionally
        KEPT, so tightening this filter would change the clean OTU count.
        Produces:
        - query.fasta: representative sequences for the taxonomy search.
        - otu_table.csv: abundance matrix (same format as DADA2
          seqtab_clean_t.csv) so downstream steps treat SWARM and DADA2
          outputs identically.

        Args:
            otu_table: Full OTU table DataFrame from build(), including the
                "chimera", "sequence", "OTU", and per-sample columns.
            query_fasta_path: Output path for the representative sequences
                FASTA (one record per kept OTU, header ``>OTU_<n>``).
            abundance_csv_path: Output path for the abundance CSV (first
                column "sequence", remaining columns are per-sample counts).

        Returns:
            Tuple of (query_fasta_path, abundance_csv_path) as Path objects,
            pointing at the two files just written.

        Raises:
            OSError: If either output path cannot be created or written.
        """
        query_fasta_path = Path(query_fasta_path)
        abundance_csv_path = Path(abundance_csv_path)
        query_fasta_path.parent.mkdir(parents=True, exist_ok=True)
        abundance_csv_path.parent.mkdir(parents=True, exist_ok=True)

        # Drop only definite chimeras (Y); keep borderline (?) and unscored (NA).
        non_chimeric = otu_table[otu_table["chimera"] != "Y"].copy()
        logger.info(
            f"Filtered chimeras: {len(otu_table)} → {len(non_chimeric)} OTUs"
        )

        # Write query FASTA
        with open(query_fasta_path, "w") as f:
            for _, row in non_chimeric.iterrows():
                f.write(f">OTU_{row['OTU']}\n{row['sequence']}\n")

        logger.info(f"Wrote representative sequences → {query_fasta_path}")

        # Write abundance CSV (DADA2 seqtab_clean_t format):
        # First column = sequence, other columns = sample abundances.
        # Reuse the single reserved-metadata list build() guards against, so the
        # exclusion set here cannot drift from the collision check there.
        sample_cols = [
            c for c in non_chimeric.columns if c not in _RESERVED_METADATA_COLS
        ]

        abundance_df = non_chimeric[["sequence"] + sample_cols].copy()
        abundance_df = abundance_df.set_index("sequence")
        abundance_df.to_csv(abundance_csv_path)

        logger.info(f"Wrote abundance table → {abundance_csv_path}")

        return query_fasta_path, abundance_csv_path

    # ------------------------------------------------------------------
    # Parsers for SWARM output and per-sample FASTA files
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_representatives(fasta_path: Union[str, Path]) -> Dict[str, str]:
        """Parse representative sequences FASTA: amplicon ID → sequence.

        Requires unwrapped (single-line-per-sequence) FASTA: the one line
        after each header is taken as the full sequence. This holds because
        the upstream vsearch steps run with ``--fasta_width 0`` (see
        vsearch_runner.py) and swarm ``--seeds`` output is unwrapped. If that
        flag is ever dropped, wrapped sequences would be silently truncated
        to their first line.

        Args:
            fasta_path: Path to the swarm representative-seeds FASTA (one OTU seed
                per record, header ``>amplicon;size=N``).

        Returns:
            Mapping of amplicon ID (the ``;size=`` suffix stripped) to its sequence.
        """
        fasta_path = Path(fasta_path)
        separator = ";size="
        representatives = {}

        with open(fasta_path) as f:
            for line in f:
                if line.startswith(">"):
                    amplicon = line.strip(">;\n").split(separator)[0]
                else:
                    representatives[amplicon] = line.strip()

        logger.debug(f"Parsed {len(representatives)} representative sequences")
        return representatives

    @staticmethod
    def _parse_stats(
        stats_path: Union[str, Path],
    ) -> Tuple[Dict[str, int], List[Tuple[str, int]], Dict[str, Tuple[int, int]]]:
        """
        Parse the SWARM ``-s`` statistics file.

        Reads positional, tab-separated columns from the SWARM ``-s`` output.
        Only the first four columns are used (column order is SWARM
        version-sensitive across releases, so verify against the installed swarm):
            - parts[0]: number of unique amplicons in the OTU (stored as
              ``cloud_size``; note this is the count of unique amplicons in
              the cluster, not an OTU-table "cloud" of reads).
            - parts[1]: total mass (sum of all amplicon abundances), the OTU
              total used to rank OTUs.
            - parts[2]: seed amplicon ID (the OTU representative).
            - parts[3]: seed amplicon abundance.
        The seed ID carries a ``;size=N`` suffix because swarm is invoked
        with ``--usearch-abundance`` (see swarm_runner.py); it is stripped so
        the seed matches IDs from the other parsers.

        Args:
            stats_path: Path to the SWARM ``-s`` statistics file (one row per OTU).

        Returns:
            Tuple of:
            - stats: seed → total mass
            - sorted_seeds: list of (seed, mass) sorted by decreasing mass
            - seeds: seed → (seed_abundance, cloud_size)
        """
        stats_path = Path(stats_path)
        stats = {}
        seeds = {}

        with open(stats_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                cloud, mass, seed_raw, seed_abundance = parts[0], parts[1], parts[2], parts[3]
                # Strip the ;size=N suffix (from --usearch-abundance) so the
                # seed ID matches the other parsers.
                seed = seed_raw.split(";size=")[0]
                stats[seed] = int(mass)
                seeds[seed] = (int(seed_abundance), int(cloud))

        sorted_seeds = sorted(stats.items(), key=lambda x: (x[1], x[0]), reverse=True)

        logger.debug(f"Parsed stats for {len(stats)} OTUs")
        return stats, sorted_seeds, seeds

    @staticmethod
    def _parse_swarms(swarm_path: Union[str, Path]) -> Dict[str, List[str]]:
        """Parse the SWARM membership file: seed amplicon -> its member amplicons.

        Each line of the swarm output lists the amplicons clustered into one OTU; the
        first token is the seed (OTU representative), the rest are members.

        Args:
            swarm_path: Path to the SWARM cluster-membership file (one OTU per line).

        Returns:
            Mapping of seed amplicon ID to the list of member amplicon IDs in that OTU
            (each ID has its ``;size=`` suffix stripped).
        """
        swarm_path = Path(swarm_path)
        swarms = {}

        with open(swarm_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Split on spaces to get individual amplicon entries,
                # then strip abundance annotations (;size=N; or ;size=N or _N)
                amplicons = []
                for entry in line.split(" "):
                    if not entry:
                        continue
                    if ";size=" in entry:
                        amp = entry.split(";size=")[0]
                    else:
                        amp = re.sub(r"_\d+$", "", entry)
                    if amp:
                        amplicons.append(amp)
                if amplicons:
                    seed = amplicons[0]
                    swarms[seed] = amplicons

        logger.debug(f"Parsed {len(swarms)} SWARM clusters")
        return swarms

    @staticmethod
    def _parse_uchime(uchime_path: Union[str, Path]) -> Dict[str, str]:
        """Parse UCHIME output: amplicon ID → chimera status (Y/N/?).

        Reads the tab-separated ``--uchimeout`` table emitted by vsearch
        ``--uchime_denovo`` (see vsearch_runner.py). Column positions are
        version-sensitive and load-bearing, so verify against the installed vsearch:
            - parts[1]: the query/sequence label (the OTU seed amplicon);
              its ``;size=N`` annotation is stripped to match other parsers.
            - parts[17]: the chimera classification flag, one of Y (chimera),
              N (not a chimera), or ? (borderline).
        Blank lines are skipped quietly; a non-blank line with no query-label
        column is logged with a [WARN] and skipped (never dropped silently); a
        line that has a label but no classification column defaults to "NA".

        Args:
            uchime_path: Path to the vsearch ``--uchimeout`` table.

        Returns:
            Mapping of amplicon ID (``;size=`` stripped) to its chimera status:
            ``"Y"`` (chimera), ``"N"`` (not a chimera), ``"?"`` (borderline), or
            ``"NA"`` (label present but unscored).
        """
        uchime_path = Path(uchime_path)
        uchime = {}

        with open(uchime_path) as f:
            for line in f:
                if not line.strip():
                    continue  # blank line: skip quietly
                parts = line.strip().split("\t")
                try:
                    seed = parts[1].split(";")[0]
                except IndexError:
                    # No silent drop: a non-blank line we cannot resolve to a seed label is
                    # logged, never dropped quietly (the SWARM/BLAST paths shipped a silent
                    # ID-mismatch zero-fill once; this catches the next one in the log).
                    logger.warning(
                        f"[WARN] _parse_uchime: skipping unparseable line in {uchime_path} "
                        f"(no query-label column): {line.strip()[:80]!r}"
                    )
                    continue
                try:
                    status = parts[17]
                except IndexError:
                    status = "NA"
                uchime[seed] = status

        logger.debug(f"Parsed chimera status for {len(uchime)} sequences")
        return uchime

    @staticmethod
    def _parse_sample_fastas(
        fasta_paths: Sequence[Union[str, Path]],
    ) -> Tuple[Dict[str, Dict[str, int]], List[str]]:
        """
        Parse per-sample FASTA files to get amplicon → sample → abundance mapping.

        Reads only header lines, taking the amplicon ID and its ``;size=N``
        abundance (written by vsearch ``--sizeout``); sequence body lines are
        ignored, so this parser is unaffected by FASTA line wrapping.

        Args:
            fasta_paths: List of per-sample dereplicated FASTA files

        Returns:
            Tuple of:
            - amplicons2samples: amplicon → {sample: abundance}
            - samples: sorted list of sample names
        """
        separator = ";size="
        samples_set: Dict[str, int] = {}
        amplicons2samples: Dict[str, Dict[str, int]] = {}

        for fasta_path in fasta_paths:
            fasta_path = Path(fasta_path)
            sample = fasta_path.stem
            samples_set[sample] = samples_set.get(sample, 0) + 1

            with open(fasta_path) as f:
                for line in f:
                    if line.startswith(">"):
                        parts = line.strip(">;\n").split(separator)
                        amplicon = parts[0]
                        abundance = int(parts[1]) if len(parts) > 1 else 1

                        if amplicon not in amplicons2samples:
                            amplicons2samples[amplicon] = {}
                        amplicons2samples[amplicon][sample] = (
                            amplicons2samples[amplicon].get(sample, 0) + abundance
                        )

        samples = sorted(samples_set.keys())
        logger.debug(
            f"Parsed {len(amplicons2samples)} amplicons across {len(samples)} samples"
        )
        return amplicons2samples, samples
