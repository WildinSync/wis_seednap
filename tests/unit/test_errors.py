"""Unit tests for the error-explainability module (humanizer, preflight, explain, base)."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from seednap.config.models import PipelineConfig
from seednap.errors import (
    SeednapError,
    explain,
    humanize_validation_error,
    preflight_checks,
)


def _base(tmp_path) -> dict:
    """A schema-valid config dict (output/logs under tmp so model_post_init can mkdir)."""
    return {
        "marker": {
            "name": "t",
            "primers": {"forward": "ACACCGCCCGTCACTCT", "reverse": "CTTCCGGTACACTTACCATG"},
        },
        "paths": {
            "raw_data": str(tmp_path),
            "output": str(tmp_path / "o"),
            "logs": str(tmp_path / "o" / "l"),
        },
        "taxonomy": {"method": "blast", "databases": {"blast": {"fasta": str(tmp_path / "ref.fasta")}}},
        "pipeline": {"steps": ["trim", "swarm", "taxonomy", "report"]},
    }


def _humanize(d) -> str:
    with pytest.raises(ValidationError) as ei:
        PipelineConfig(**d)
    return humanize_validation_error(ei.value, Path("config.yaml"))


def test_humanizer_typo_suggests_closest_key(tmp_path):
    d = _base(tmp_path)
    d["paths"]["raw_dat"] = "/tmp"
    msg = _humanize(d)
    assert "paths.raw_dat" in msg
    assert "closest valid key is 'raw_data'" in msg
    assert "SDN-CFG-001" in msg


def test_humanizer_removed_key_gives_migration_hint(tmp_path):
    d = _base(tmp_path)
    d["cleaning"] = {"enabled": False, "mode": "flag"}
    msg = _humanize(d)
    assert "cleaning.enabled" in msg
    assert "removed" in msg and "pipeline.steps" in msg


def test_humanizer_bad_literal_lists_allowed(tmp_path):
    d = _base(tmp_path)
    d["taxonomy"]["method"] = "blastn"
    msg = _humanize(d)
    assert "Allowed" in msg and "blast" in msg
    assert "closest allowed value is 'blast'" in msg
    assert "SDN-CFG-005" in msg


def test_humanizer_missing_required(tmp_path):
    d = _base(tmp_path)
    del d["marker"]
    msg = _humanize(d)
    assert "marker" in msg and "missing" in msg.lower() and "SDN-CFG-002" in msg


def test_humanizer_out_of_range(tmp_path):
    d = _base(tmp_path)
    d["dada2"] = {"filter": {"max_ee": -1}}
    msg = _humanize(d)
    assert "max_ee" in msg and "out of range" in msg and "SDN-CFG-004" in msg


def test_humanizer_string_too_short_is_not_unknown_key(tmp_path):
    # A too-short primer is a valid key with a bad value: it must NOT be tagged
    # SDN-CFG-001 ('unknown key'), which would send the user down the typo path.
    d = _base(tmp_path)
    d["marker"]["primers"]["forward"] = "ACGT"  # min_length=10
    msg = _humanize(d)
    assert "marker.primers.forward" in msg
    assert "wrong length" in msg
    assert "SDN-CFG-004" in msg
    assert "SDN-CFG-001" not in msg


def test_humanizer_demux_cross_field_uses_pipeline_topic(tmp_path):
    # The demultiplex-protocol cross-field validator raises with an empty loc;
    # it must map to the pipeline.steps topic (006), not the generic 'invalid
    # choice' code (005).
    d = _base(tmp_path)
    d["pipeline"]["steps"] = ["demultiplex", "trim", "swarm", "taxonomy", "report"]
    msg = _humanize(d)
    assert "demultiplex" in msg and "pipeline.steps" in msg
    assert "SDN-CFG-006" in msg
    assert "SDN-CFG-005" not in msg


def test_preflight_flags_missing_inputs(tmp_path):
    d = _base(tmp_path)
    d["paths"]["raw_data"] = str(tmp_path / "does_not_exist")  # missing dir; fasta also missing
    cfg = PipelineConfig(**d)
    problems = preflight_checks(cfg)
    assert any("Raw-data directory does not exist" in p.summary for p in problems)
    assert any("fasta does not exist" in p.summary for p in problems)
    assert all(p.code == "SDN-CFG-007" for p in problems)


def test_preflight_passes_when_inputs_exist(tmp_path):
    (tmp_path / "ref.fasta").write_text(">x\nACGT\n")
    cfg = PipelineConfig(**_base(tmp_path))  # raw_data=tmp_path exists, fasta exists
    assert preflight_checks(cfg) == []


def test_explain_known_and_unknown():
    assert "Unknown configuration key" in (explain("SDN-CFG-001") or "")
    assert explain("sdn-cfg-001") is not None  # case-insensitive
    assert explain("SDN-NOPE-999") is None


def test_seednaperror_render():
    e = SeednapError("X failed", why="because Y", fix="do Z", code="SDN-CFG-001")
    s = e.render()
    assert "X failed" in s and "Why: because Y" in s and "Fix: do Z" in s
    assert "SDN-CFG-001" in s and "seednap explain SDN-CFG-001" in s
