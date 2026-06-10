"""Unit test for implicit sample-column detection in CleaningProcessor.

Standalone ``seednap clean`` does not pass an explicit ``sample_cols`` list, so the
processor must auto-detect samples on a taxonomy/BLAST-shaped table: a sample column is
numeric AND not one of the known per-OTU annotation columns (taxonomic ranks,
sequence/Sequence, ASV_ID, OTU/OTU_ID, taxon, rank, pident, is_contaminant_candidate).
Without this, the numeric annotation columns ``pident`` and the boolean
``is_contaminant_candidate`` are mistaken for biological samples and corrupt the result.

This pins the auto-detection (the implicit path); the explicit-``sample_cols`` path is
covered separately in test_cleaning.py.
"""

import pandas as pd

from seednap.config.manifest import SampleManifest, SampleManifestRow
from seednap.steps.cleaning import CleaningProcessor


def _sample(ev, ext):
    return SampleManifestRow(eventID=ev, seq_run_id="run1", samp_category="sample",
                             eventDate="2025-01-01", extraction_ID=ext)


def _neg(ev, neg_type, ext=None):
    return SampleManifestRow(eventID=ev, seq_run_id="run1", samp_category="negative control",
                             neg_cont_type=neg_type, extraction_ID=ext)


def test_taxonomy_table_auto_detects_only_sample_columns():
    """A taxonomy-shaped table (ranks + pident + is_contaminant_candidate + 2 negative
    controls + 2 real samples) cleaned WITHOUT explicit sample_cols: only the 4 sample
    columns are treated as samples, pident is untouched, and the per-sample report lists
    only the 2 biological samples (not pident / is_contaminant_candidate / controls)."""
    manifest = SampleManifest(rows=[
        _sample("S1", "EXP1"), _sample("S2", "EXP1"),
        _neg("Bext1", "extraction negative", "EXP1"),
        _neg("Bpcr", "PCR negative", None),
    ])
    # O1 lives in the PCR blank (whole-dataset); O2 is clean.
    tax = pd.DataFrame({
        "ASV_ID": ["O1", "O2"],
        "pident": [100.0, 98.0],                  # numeric, NOT a sample
        "kingdom": ["Animalia", "Animalia"],
        "phylum": ["Chordata", "Chordata"],
        "class": ["Actinopteri", "Actinopteri"],
        "order": ["Cypriniformes", "Cypriniformes"],
        "family": ["Cyprinidae", "Cyprinidae"],
        "genus": ["Squalius", "Barbus"],
        "species": ["Squalius_cephalus", "Barbus_barbus"],
        "is_contaminant_candidate": [False, True],  # boolean -> numeric, NOT a sample
        "Bext1": [0, 0],                           # negative control
        "Bpcr": [10, 0],                           # negative control
        "S1": [100, 50],                           # real sample
        "S2": [100, 50],                           # real sample
        "Sequence": ["ACGT", "TGCA"],
    })
    pident_before = list(tax["pident"])
    contam_before = list(tax["is_contaminant_candidate"])

    cleaned, report, result = CleaningProcessor(mode="subtract").clean(
        tax, manifest, id_col="ASV_ID"
    )

    # Exactly the 2 biological sample columns are treated as samples; the 2 negative
    # controls are counted as controls, not samples.
    assert result.n_samples == 2
    assert result.n_controls == 2

    # pident and is_contaminant_candidate are NOT treated as samples (left untouched).
    assert list(cleaned["pident"]) == pident_before
    assert list(cleaned["is_contaminant_candidate"]) == contam_before

    # The per-sample report lists only the real samples -- never the numeric annotation
    # columns and never the controls.
    assert set(report["eventID"]) == {"S1", "S2"}

    # Sanity: the actual decontamination still happened on the real samples (O1 was in
    # the PCR blank, so it is zeroed in S1/S2; O2 survives).
    by = cleaned.set_index("ASV_ID")
    assert by.loc["O1", "S1"] == 0 and by.loc["O1", "S2"] == 0
    assert by.loc["O2", "S1"] == 50 and by.loc["O2", "S2"] == 50
