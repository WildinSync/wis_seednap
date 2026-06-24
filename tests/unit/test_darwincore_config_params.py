"""Increment 2: the DarwinCore builder prefers config-supplied provenance over the CSV.

When the 'darwincore' pipeline step supplies otu_db / chimera_check derived from the run
config, those win over the project-metadata CSV (so pipeline parameters are not re-entered);
a differing CSV value is surfaced with a [WARN]. With no config value the CSV value is used,
so the standalone create-gbif command is unaffected.
"""

import logging

from seednap.steps.formatting.darwincore_builder import DarwinCoreBuilder


def test_config_value_wins():
    assert DarwinCoreBuilder._prefer_config("otu_db", "CONFIG_DB", "csv_db") == "CONFIG_DB"


def test_falls_back_to_csv_when_no_config_value():
    assert DarwinCoreBuilder._prefer_config("otu_db", None, "csv_db") == "csv_db"


def test_warns_on_disagreement(caplog):
    with caplog.at_level(logging.WARNING):
        value = DarwinCoreBuilder._prefer_config("otu_db", "CONFIG_DB", "csv_db")
    assert value == "CONFIG_DB"
    assert any("otu_db" in r.message and "differs" in r.message for r in caplog.records)


def test_no_warning_when_csv_matches_or_empty(caplog):
    with caplog.at_level(logging.WARNING):
        DarwinCoreBuilder._prefer_config("chimera_check", "X", "X")  # identical
        DarwinCoreBuilder._prefer_config("chimera_check", "X", "")   # empty CSV
        DarwinCoreBuilder._prefer_config("chimera_check", "X", None)  # missing CSV
    assert not any("differs" in r.message for r in caplog.records)
