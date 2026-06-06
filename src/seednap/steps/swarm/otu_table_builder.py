"""Build OTU contingency table from SWARM clustering outputs.

Parses SWARM output files (representatives, stats, swarm membership, chimera
detection) and per-sample FASTA files to produce a full OTU contingency table
and DADA2-compatible normalized outputs for downstream taxonomy assignment.

Implements the OTU contingency table construction from the published SWARM
amplicon clustering pipeline.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd

logger = logging.getLogger(__name__)


class OtuTableBuilder:
    """
    Build OTU contingency table from SWARM clustering results.

    Parses all SWARM output files and per-sample abundance data to create
    a comprehensive OTU table with per-sample read counts.
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
        Build the full OTU contingency table.

        Args:
            representatives_fasta: Sorted representative sequences FASTA
            stats_file: SWARM statistics file
            swarm_file: SWARM cluster membership file
            uchime_file: UCHIME chimera detection results (None to skip)
            sample_fastas: Per-sample dereplicated FASTA files

        Returns:
            DataFrame with OTU table (sorted by decreasing total abundance)
        """
        representatives = self._parse_representatives(representatives_fasta)
        stats, sorted_seeds, seeds = self._parse_stats(stats_file)
        swarms = self._parse_swarms(swarm_file)
        uchime = self._parse_uchime(uchime_file) if uchime_file else {}
        amplicons2samples, samples = self._parse_sample_fastas(sample_fastas)

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

        Filters chimeric OTUs and produces:
        - query.fasta: representative sequences for taxonomy search
        - otu_table.csv: abundance matrix (same format as DADA2 seqtab_clean_t.csv)

        Args:
            otu_table: Full OTU table DataFrame from build()
            query_fasta_path: Output path for representative sequences FASTA
            abundance_csv_path: Output path for abundance CSV

        Returns:
            Tuple of (query_fasta_path, abundance_csv_path)
        """
        query_fasta_path = Path(query_fasta_path)
        abundance_csv_path = Path(abundance_csv_path)
        query_fasta_path.parent.mkdir(parents=True, exist_ok=True)
        abundance_csv_path.parent.mkdir(parents=True, exist_ok=True)

        # Filter chimeric OTUs
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
        # First column = sequence, other columns = sample abundances
        metadata_cols = [
            "OTU", "total", "cloud", "amplicon", "length",
            "abundance", "chimera", "spread",
        ]
        sample_cols = [c for c in non_chimeric.columns if c not in metadata_cols + ["sequence"]]

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
        """Parse representative sequences FASTA: amplicon ID → sequence."""
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
        Parse SWARM stats file.

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
                # Strip ;size=N; annotation so seed matches other parsers
                seed = seed_raw.split(";size=")[0]
                stats[seed] = int(mass)
                seeds[seed] = (int(seed_abundance), int(cloud))

        sorted_seeds = sorted(stats.items(), key=lambda x: (x[1], x[0]), reverse=True)

        logger.debug(f"Parsed stats for {len(stats)} OTUs")
        return stats, sorted_seeds, seeds

    @staticmethod
    def _parse_swarms(swarm_path: Union[str, Path]) -> Dict[str, List[str]]:
        """Parse SWARM membership file: seed → list of member amplicon IDs."""
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
        """Parse UCHIME output: amplicon ID → chimera status (Y/N/?)."""
        uchime_path = Path(uchime_path)
        uchime = {}

        with open(uchime_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                try:
                    seed = parts[1].split(";")[0]
                except IndexError:
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
