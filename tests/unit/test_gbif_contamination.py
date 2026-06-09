"""GBIFFormatter must carry is_contaminant_candidate through to the GBIF output.

Regression: the column was previously dropped at the final column selection, so create-gbif's
contamination_flag was always False end-to-end. The taxonomy step sets is_contaminant_candidate
from taxonomy.contaminants; it must survive format-gbif so create-gbif can surface it.
"""

import pandas as pd

from seednap.steps.formatting.gbif_formatter import GBIFFormatter


def _wide_table() -> pd.DataFrame:
    """A minimal wide-format taxonomy table with one contaminant row and one sample column."""
    return pd.DataFrame(
        {
            "kingdom": ["Animalia", "Animalia"],
            "phylum": ["Chordata", "Chordata"],
            "class": ["Actinopteri", "Mammalia"],
            "order": ["Cypriniformes", "Primates"],
            "family": ["Cyprinidae", "Hominidae"],
            "genus": ["Cyprinus", "Homo"],
            "species": ["Cyprinus carpio", "Homo sapiens"],
            "sequence": ["ACGTACGT", "TTGCTTGC"],
            "is_contaminant_candidate": [False, True],
            "sampleA": [12, 5],
        }
    )


def test_from_dada2_rdp_preserves_contaminant_flag(tmp_path):
    inp = tmp_path / "tax.csv"
    _wide_table().to_csv(inp, index=False)
    out = GBIFFormatter().from_dada2_rdp(inp)
    assert "is_contaminant_candidate" in out.columns
    homo = out[out["species"] == "Homo sapiens"]
    assert bool(homo["is_contaminant_candidate"].iloc[0]) is True


def test_from_dada2_rdp_without_contaminant_column_does_not_crash(tmp_path):
    """When taxonomy.contaminants is unset the column is absent; output omits it, no crash."""
    inp = tmp_path / "tax.csv"
    _wide_table().drop(columns=["is_contaminant_candidate"]).to_csv(inp, index=False)
    out = GBIFFormatter().from_dada2_rdp(inp)
    assert "is_contaminant_candidate" not in out.columns
    assert len(out) > 0


def _ecotag_table() -> pd.DataFrame:
    """Minimal ecotag-shaped table (the *_name columns from_ecotag renames)."""
    return pd.DataFrame(
        {
            "order_name": ["Cypriniformes", "Primates"],
            "family_name": ["Cyprinidae", "Hominidae"],
            "genus_name": ["Cyprinus", "Homo"],
            "species_name": ["Cyprinus carpio", "Homo sapiens"],
            "sequence": ["ACGTACGT", "TTGCTTGC"],
            "is_contaminant_candidate": [False, True],
            "sampleA": [12, 5],
        }
    )


def test_from_ecotag_preserves_contaminant_flag(tmp_path):
    inp = tmp_path / "ecotag.csv"
    _ecotag_table().to_csv(inp, index=False)
    out = GBIFFormatter().from_ecotag(inp)
    assert "is_contaminant_candidate" in out.columns
    homo = out[out["species"] == "Homo sapiens"]
    assert bool(homo["is_contaminant_candidate"].iloc[0]) is True


def test_from_ecotag_without_contaminant_column_does_not_crash(tmp_path):
    inp = tmp_path / "ecotag.csv"
    _ecotag_table().drop(columns=["is_contaminant_candidate"]).to_csv(inp, index=False)
    out = GBIFFormatter().from_ecotag(inp)
    assert "is_contaminant_candidate" not in out.columns
    assert len(out) > 0
