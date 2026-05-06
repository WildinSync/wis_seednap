"""Integration test for OtuTableBuilder: catches the silent-zero OTU table bug.

The pre-Mar-2026 SWARM pipeline silently produced an OTU contingency table
with all-zero counts when SWARM seed IDs (with `;size=N;` annotations) didn't
match the per-sample FASTA amplicon IDs. This test reproduces a small
realistic input and verifies the table contains non-zero per-sample counts.
"""

from pathlib import Path

import pandas as pd
import pytest

from seednap.steps.swarm.otu_table_builder import OtuTableBuilder


@pytest.fixture
def swarm_outputs(tmp_path: Path) -> Path:
    """Create realistic SWARM-style outputs with semicolon-annotated seed IDs."""
    # Sorted representatives FASTA: SWARM emits ;size=N; annotations on the headers
    # but the OTU table builder must strip them before joining.
    reps = tmp_path / "cluster_representatives.sorted.fasta"
    reps.write_text(
        ">SHA1A;size=20;\nACGTACGTACGTACGTACGT\n"
        ">SHA1B;size=10;\nTTTTAAAACCCCGGGGAAAA\n"
        ">SHA1C;size=5;\nGGGGCCCCAAAATTTTGGGG\n"
    )
    # Stats: cloud_size, total_mass, seed_id_with_size, seed_abundance, ...
    stats = tmp_path / "all.stats"
    stats.write_text(
        "1\t20\tSHA1A;size=20\t20\t\n"
        "1\t10\tSHA1B;size=10\t10\t\n"
        "1\t5\tSHA1C;size=5\t5\t\n"
    )
    # Swarm membership: seed [member member ...]
    swarm = tmp_path / "all.swarm"
    swarm.write_text(
        "SHA1A;size=10; SHA1AA;size=5; SHA1AB;size=5;\n"
        "SHA1B;size=10;\n"
        "SHA1C;size=5;\n"
    )
    # Per-sample dereplicated FASTAs - amplicon IDs match seeds (no size suffix here)
    sample_dir = tmp_path / "dereplicated"
    sample_dir.mkdir()
    (sample_dir / "S1.fasta").write_text(
        ">SHA1A;size=8;\nACGTACGTACGTACGTACGT\n"
        ">SHA1AA;size=3;\nACGTACGTACGTACGTACGT\n"
        ">SHA1B;size=6;\nTTTTAAAACCCCGGGGAAAA\n"
    )
    (sample_dir / "S2.fasta").write_text(
        ">SHA1A;size=2;\nACGTACGTACGTACGTACGT\n"
        ">SHA1AB;size=2;\nACGTACGTACGTACGTACGT\n"
        ">SHA1B;size=4;\nTTTTAAAACCCCGGGGAAAA\n"
        ">SHA1C;size=5;\nGGGGCCCCAAAATTTTGGGG\n"
    )
    return tmp_path


def test_otu_table_has_nonzero_counts(swarm_outputs: Path) -> None:
    """The OTU contingency table must have non-zero per-sample counts.

    Specifically guards against the pre-fix bug where seed-ID/amplicon-ID
    mismatches produced an all-zeros table even with valid input.
    """
    builder = OtuTableBuilder()
    sample_fastas = sorted((swarm_outputs / "dereplicated").glob("*.fasta"))
    df = builder.build(
        representatives_fasta=swarm_outputs / "cluster_representatives.sorted.fasta",
        stats_file=swarm_outputs / "all.stats",
        swarm_file=swarm_outputs / "all.swarm",
        uchime_file=None,
        sample_fastas=sample_fastas,
    )

    # Built 3 OTUs across 2 samples
    assert len(df) == 3
    sample_cols = ["S1", "S2"]
    for col in sample_cols:
        assert col in df.columns, f"missing sample column {col}"

    # The bug class: all per-sample counts were zero. Assert at least one cell
    # is non-zero in the abundance block.
    abundance_block = df[sample_cols].to_numpy()
    assert (abundance_block > 0).any(), (
        f"OTU table has no positive counts; this is the silent-zero bug class:\n{df}"
    )

    # Stronger: every OTU we built should be present in at least one sample
    per_otu_total = abundance_block.sum(axis=1)
    assert (per_otu_total > 0).all(), (
        f"At least one OTU has zero count in every sample:\n{df}"
    )


def test_seed_size_annotations_are_stripped(swarm_outputs: Path) -> None:
    """The seed parser must strip ;size=N; annotations so amplicon IDs match across files."""
    builder = OtuTableBuilder()
    df = builder.build(
        representatives_fasta=swarm_outputs / "cluster_representatives.sorted.fasta",
        stats_file=swarm_outputs / "all.stats",
        swarm_file=swarm_outputs / "all.swarm",
        uchime_file=None,
        sample_fastas=sorted((swarm_outputs / "dereplicated").glob("*.fasta")),
    )
    # The amplicon column holds the canonical seed IDs (no ;size=...; suffix)
    assert ";size=" not in str(df["amplicon"].tolist())
