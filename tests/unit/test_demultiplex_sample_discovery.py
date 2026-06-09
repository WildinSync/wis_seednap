"""Regression test for ligation demultiplex sample discovery.

LigationTrimmer.process_library lists demultiplexed samples by globbing the
demux directory and excluding cutadapt's catch-all 'unknown' bucket. The bucket
is named exactly 'unknown.R*.fastq.gz'; a substring filter ('unknown' in name)
silently drops legitimate samples whose eventID contains 'unknown' (e.g. the
lab's BB_Unknown_Svalbard_unknown dataset). This test reproduces the discovery
logic over a temp directory and asserts the exact-match behavior: such a sample
survives while the real catch-all bucket is excluded.
"""

from __future__ import annotations

from pathlib import Path

# Read the production source from the repo-under-test directly. Importing the
# module can resolve to a different installed copy on the eDNA server, so the
# regression guard below reads the file next to this test's package root.
_REPO_TRIMMING = (
    Path(__file__).resolve().parents[2]
    / "src"
    / "seednap"
    / "steps"
    / "trimming"
    / "trimming_pipeline.py"
)


def _discover_samples(demux_dir: Path) -> list[str]:
    """Mirror the sample-discovery selection in LigationTrimmer.process_library.

    Kept in lockstep with the production code: extract the {name} token, then
    exclude only the exact 'unknown' catch-all (case-insensitive), never a
    substring match.
    """
    all_names = sorted(
        {
            f.name.replace(".R1.fastq.gz", "").replace(".R2.fastq.gz", "")
            for f in demux_dir.glob("*.R*.fastq.gz")
        }
    )
    return [name for name in all_names if name.lower() != "unknown"]


def test_sample_with_unknown_in_name_is_not_dropped(tmp_path: Path) -> None:
    demux_dir = tmp_path / "demultiplex"
    demux_dir.mkdir()

    # A legitimate sample whose eventID contains 'unknown'/'Unknown'.
    legit = "BB_Unknown_Svalbard_unknown"
    for read in ("R1", "R2"):
        (demux_dir / f"{legit}.{read}.fastq.gz").touch()
    # A plain sample.
    for read in ("R1", "R2"):
        (demux_dir / f"sampleA.{read}.fastq.gz").touch()
    # cutadapt's exact catch-all bucket (should be excluded).
    for read in ("R1", "R2"):
        (demux_dir / f"unknown.{read}.fastq.gz").touch()

    samples = _discover_samples(demux_dir)

    # The substring filter would have dropped the BB_Unknown_Svalbard_unknown
    # sample; the exact-match filter keeps it.
    assert legit in samples
    assert "sampleA" in samples
    # Only the exact 'unknown' catch-all is excluded.
    assert "unknown" not in samples
    assert len(samples) == 2


def test_production_uses_exact_unknown_match_not_substring() -> None:
    """Guard against a regression to the substring filter.

    The substring form `"unknown" not in f.name.lower()` silently drops any
    sample whose name contains 'unknown'. process_library must compare the
    extracted sample name to the exact value 'unknown' instead.
    """
    src = _REPO_TRIMMING.read_text()
    assert 'name.lower() != "unknown"' in src
    assert '"unknown" not in f.name.lower()' not in src
