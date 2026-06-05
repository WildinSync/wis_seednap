"""Unit tests for the control-decontamination CleaningProcessor.

Synthetic fixture pins the validated standard: an extraction blank cleans only its own
extraction batch, a PCR blank cleans the whole dataset, flag mode never changes counts,
and every edge (orphan control, zero controls, unmatched batch) is warned, not silent.
"""

import logging

import pandas as pd
import pytest

from seednap.config.manifest import SampleManifest, SampleManifestRow
from seednap.steps.cleaning import CleaningProcessor


def _sample(ev, ext):
    return SampleManifestRow(eventID=ev, seq_run_id="run1", samp_category="sample",
                             eventDate="2025-01-01", extraction_ID=ext)


def _neg(ev, neg_type, ext=None):
    return SampleManifestRow(eventID=ev, seq_run_id="run1", samp_category="negative control",
                             neg_cont_type=neg_type, extraction_ID=ext)


@pytest.fixture
def manifest():
    return SampleManifest(rows=[
        _sample("S1", "EXP1"), _sample("S2", "EXP1"), _sample("S3", "EXP2"),
        _neg("Bext1", "extraction negative", "EXP1"),
        _neg("Bpcr", "PCR negative", None),
    ])


@pytest.fixture
def abundance():
    # O1 in the PCR blank (whole-dataset); O2 in the EXP1 extraction blank; O3 clean.
    return pd.DataFrame({
        "sequence": ["O1", "O2", "O3"],
        "S1": [100, 50, 20], "S2": [100, 50, 0], "S3": [100, 50, 20],
        "Bext1": [0, 5, 0], "Bpcr": [10, 0, 0],
    })


def test_invalid_mode():
    with pytest.raises(ValueError):
        CleaningProcessor(mode="bogus")


def test_flag_mode_never_changes_counts(manifest, abundance):
    df, rep, res = CleaningProcessor(mode="flag").clean(abundance, manifest, id_col="sequence")
    assert res.mode == "flag" and res.total_reads_removed == 0
    # counts identical to input
    for c in ("S1", "S2", "S3", "Bext1", "Bpcr"):
        assert list(df[c]) == list(abundance[c])
    assert (rep["reads_before"] == rep["reads_after"]).all()
    # O1 (PCR) and O2 (ext) are flagged; O3 is not
    flagged = dict(zip(df["sequence"], df["in_negative_control"]))
    assert flagged["O1"] and flagged["O2"] and not flagged["O3"]


def test_subtract_extraction_scoping(manifest, abundance):
    df, rep, res = CleaningProcessor(mode="subtract").clean(abundance, manifest, id_col="sequence")
    r = rep.set_index("eventID")
    # S1/S2 (EXP1): lose O1 (PCR, whole-ds) + O2 (Bext1, EXP1) = 2 OTUs, 150 reads
    assert r.loc["S1", "n_otus_removed"] == 2 and r.loc["S1", "n_reads_removed"] == 150
    assert r.loc["S2", "n_reads_removed"] == 150
    # S3 (EXP2): loses only O1 (PCR); O2's extraction blank is EXP1, so it is NOT removed
    assert r.loc["S3", "n_otus_removed"] == 1 and r.loc["S3", "n_reads_removed"] == 100
    # the actual matrix: O2 still present in S3, gone from S1/S2
    by = df.set_index("sequence")
    assert by.loc["O2", "S3"] == 50 and by.loc["O2", "S1"] == 0 and by.loc["O2", "S2"] == 0
    assert by.loc["O3", "S1"] == 20  # untouched


def test_orphan_control_classified_by_name_and_warned(manifest, abundance, caplog):
    """A control column absent from the manifest (e.g. Blank-PCR-3) is classified by name."""
    ab = abundance.copy()
    ab["Blank-PCR-3"] = [7, 0, 0]  # an orphan PCR blank, not in the manifest
    with caplog.at_level(logging.WARNING):
        df, rep, res = CleaningProcessor(mode="subtract").clean(ab, manifest, id_col="sequence")
    assert res.n_controls == 3  # Bext1, Bpcr, Blank-PCR-3
    assert any("Blank-PCR-3" in r.message and "absent" in r.message for r in caplog.records)


def test_zero_controls_warns(caplog):
    m = SampleManifest(rows=[_sample("S1", "EXP1"), _sample("S2", "EXP1")])
    ab = pd.DataFrame({"sequence": ["O1"], "S1": [10], "S2": [10]})
    with caplog.at_level(logging.WARNING):
        df, rep, res = CleaningProcessor(mode="subtract").clean(ab, m, id_col="sequence")
    assert res.n_controls == 0 and res.total_reads_removed == 0
    assert any("at least one negative control" in r.message for r in caplog.records)


def test_extraction_blank_matching_no_sample_warns(caplog):
    m = SampleManifest(rows=[_sample("S1", "EXP1"), _neg("Bext9", "extraction negative", "EXP9")])
    ab = pd.DataFrame({"sequence": ["O1"], "S1": [10], "Bext9": [5]})
    with caplog.at_level(logging.WARNING):
        CleaningProcessor(mode="subtract").clean(ab, m, id_col="sequence")
    assert any("matches no biological sample" in r.message for r in caplog.records)
