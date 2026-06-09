"""Unit tests for SwarmProcessor._find_sample_pairs R1/R2 pairing.

The pairing derives the R2 filename from the R1 filename. A naive
``r1.name.replace("R1", "R2", 1)`` rewrites the FIRST "R1" substring anywhere
in the name, so a sample-name prefix that itself contains "R1" (e.g.
``MR12_R1.fastq``, ``R1B-site_R1.fastq``, ``SAMPLE-R1A_R1_001.fastq.gz``) gets
corrupted: the derived R2 path does not exist and the whole biological sample
is silently dropped (or, worse, paired with another sample's R2). These tests
pin the read-token-only rewrite so only the matched R1 token is changed.
"""

from __future__ import annotations

from pathlib import Path

from seednap.steps.swarm.processor import SwarmProcessor


def _touch(path: Path) -> None:
    path.write_text("")


def test_simple_names_pair_correctly(tmp_path: Path) -> None:
    for name in ("sampleA.R1.fastq", "sampleB_R1.fastq", "sampleC_R1_001.fastq.gz"):
        _touch(tmp_path / name)
    for name in ("sampleA.R2.fastq", "sampleB_R2.fastq", "sampleC_R2_001.fastq.gz"):
        _touch(tmp_path / name)

    pairs = SwarmProcessor._find_sample_pairs(tmp_path)
    found = {sample: (r1.name, r2.name) for sample, r1, r2 in pairs}

    assert found["sampleA"] == ("sampleA.R1.fastq", "sampleA.R2.fastq")
    assert found["sampleB"] == ("sampleB_R1.fastq", "sampleB_R2.fastq")
    assert found["sampleC"] == ("sampleC_R1_001.fastq.gz", "sampleC_R2_001.fastq.gz")


def test_prefix_contains_r1_substring_not_corrupted(tmp_path: Path) -> None:
    """Sample names whose prefix contains 'R1' must still pair with the real R2."""
    # Each of these would mis-derive its R2 with a blind first-substring replace.
    cases = {
        "MR12_R1.fastq": "MR12_R2.fastq",
        "R1B-site_R1.fastq": "R1B-site_R2.fastq",
        "SAMPLE-R1A_R1_001.fastq.gz": "SAMPLE-R1A_R2_001.fastq.gz",
    }
    for r1_name, r2_name in cases.items():
        _touch(tmp_path / r1_name)
        _touch(tmp_path / r2_name)

    pairs = SwarmProcessor._find_sample_pairs(tmp_path)
    paired_r1_to_r2 = {r1.name: r2.name for _, r1, r2 in pairs}

    # Every R1 must be paired (no sample silently dropped) and to the CORRECT R2.
    assert len(pairs) == len(cases)
    for r1_name, r2_name in cases.items():
        assert paired_r1_to_r2.get(r1_name) == r2_name


def test_orphan_r1_with_corruptible_prefix_is_dropped_not_misderived(
    tmp_path: Path,
) -> None:
    """A genuinely orphaned R1 (no matching R2) is dropped, and the buggy
    first-substring rewrite must not accidentally pair it with an unrelated R2.

    ``MR12_R1.fastq`` would mis-derive to ``MR22_R1.fastq`` under the old code;
    if a file by that name existed it would be silently mis-paired. Here the
    correct R2 (``MR12_R2.fastq``) is absent, so the sample must simply not pair.
    """
    _touch(tmp_path / "MR12_R1.fastq")
    # A decoy that the buggy MR12->MR22 derivation could latch onto.
    _touch(tmp_path / "MR22_R1.fastq")
    _touch(tmp_path / "MR22_R2.fastq")

    pairs = SwarmProcessor._find_sample_pairs(tmp_path)
    paired = {sample: (r1.name, r2.name) for sample, r1, r2 in pairs}

    # MR12 has no real R2 -> not paired (and never mis-paired to MR22's R2).
    assert "MR12" not in paired
    # MR22 pairs correctly with its own R2.
    assert paired["MR22"] == ("MR22_R1.fastq", "MR22_R2.fastq")
