"""GBIFFormatter INFO log must name the actual source format, not always "DADA2".

Regression: from_blast and from_decipher delegate to from_dada2_rdp, which used to
hardcode "Converting DADA2 output to GBIF format". A blast/ecotag/decipher run therefore
logged a misleading DADA2 message. The conversion log line must reflect the real method.
"""

import pandas as pd

from seednap.steps.formatting.gbif_formatter import GBIFFormatter


def _wide_table() -> pd.DataFrame:
    """A minimal DADA2/BLAST-shaped wide table with one sample column."""
    return pd.DataFrame(
        {
            "kingdom": ["Animalia"],
            "phylum": ["Chordata"],
            "class": ["Actinopteri"],
            "order": ["Cypriniformes"],
            "family": ["Cyprinidae"],
            "genus": ["Cyprinus"],
            "species": ["Cyprinus carpio"],
            "sequence": ["ACGTACGT"],
            "sampleA": [12],
        }
    )


def _ecotag_table() -> pd.DataFrame:
    """Minimal ecotag-shaped table (the *_name columns from_ecotag renames)."""
    return pd.DataFrame(
        {
            "order_name": ["Cypriniformes"],
            "family_name": ["Cyprinidae"],
            "genus_name": ["Cyprinus"],
            "species_name": ["Cyprinus carpio"],
            "sequence": ["ACGTACGT"],
            "sampleA": [12],
        }
    )


def test_from_blast_log_names_blast_not_dada2(tmp_path, caplog):
    inp = tmp_path / "tax.csv"
    _wide_table().to_csv(inp, index=False)
    with caplog.at_level("INFO", logger="seednap.steps.formatting.gbif_formatter"):
        GBIFFormatter().from_blast(inp)
    convert_lines = [r.message for r in caplog.records if "to GBIF format" in r.message]
    assert convert_lines, "expected a 'Converting ... to GBIF format' INFO line"
    assert any("blast" in m for m in convert_lines)
    assert all("DADA2" not in m for m in convert_lines)


def test_from_decipher_log_names_decipher_not_dada2(tmp_path, caplog):
    inp = tmp_path / "tax.csv"
    _wide_table().to_csv(inp, index=False)
    with caplog.at_level("INFO", logger="seednap.steps.formatting.gbif_formatter"):
        GBIFFormatter().from_decipher(inp)
    convert_lines = [r.message for r in caplog.records if "to GBIF format" in r.message]
    assert convert_lines, "expected a 'Converting ... to GBIF format' INFO line"
    assert any("decipher" in m for m in convert_lines)
    assert all("DADA2" not in m for m in convert_lines)


def test_from_ecotag_log_names_ecotag_not_dada2(tmp_path, caplog):
    inp = tmp_path / "ecotag.csv"
    _ecotag_table().to_csv(inp, index=False)
    with caplog.at_level("INFO", logger="seednap.steps.formatting.gbif_formatter"):
        GBIFFormatter().from_ecotag(inp)
    convert_lines = [r.message for r in caplog.records if "to GBIF format" in r.message]
    assert convert_lines, "expected a 'Converting ... to GBIF format' INFO line"
    assert any("ecotag" in m for m in convert_lines)
    assert all("DADA2" not in m for m in convert_lines)


def test_from_dada2_rdp_log_still_names_dada2(tmp_path, caplog):
    inp = tmp_path / "tax.csv"
    _wide_table().to_csv(inp, index=False)
    with caplog.at_level("INFO", logger="seednap.steps.formatting.gbif_formatter"):
        GBIFFormatter().from_dada2_rdp(inp)
    convert_lines = [r.message for r in caplog.records if "to GBIF format" in r.message]
    assert convert_lines, "expected a 'Converting ... to GBIF format' INFO line"
    assert any("dada2" in m for m in convert_lines)
