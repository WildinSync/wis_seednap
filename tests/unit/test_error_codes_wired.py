"""Verify the SDN- error codes are actually emitted at runtime.

These guard against the orphaned-code regression: a code in the catalog
(used by `seednap explain`) must correspond to a message a user can really
hit. We exercise the three runtime paths and assert their messages carry the
matching code, and that each code resolves via `seednap explain`.

- SDN-TOOL-001: subprocess wrapper, tool not on PATH (FileNotFoundError)
- SDN-TOOL-002: subprocess wrapper, tool exits non-zero (CalledProcessError)
- SDN-CFG-009 : config loader, malformed YAML and empty config
"""

import sys
from pathlib import Path

import pytest

from seednap.config.loader import ConfigError, load_yaml
from seednap.errors import explain
from seednap.utils.subprocess import run_subprocess


# --- subprocess wrapper -------------------------------------------------------

def test_subprocess_tool_not_found_carries_tool_001() -> None:
    """A missing executable raises with the [SDN-TOOL-001] tag."""
    with pytest.raises(RuntimeError) as ei:
        run_subprocess(
            ["seednap_no_such_tool_xyz", "--version"],
            error_class=RuntimeError,
        )
    msg = str(ei.value)
    assert "[SDN-TOOL-001]" in msg
    assert "not installed or not on PATH" in msg


def test_subprocess_nonzero_exit_carries_tool_002() -> None:
    """A tool that runs but exits non-zero raises with the [SDN-TOOL-002] tag."""
    with pytest.raises(RuntimeError) as ei:
        run_subprocess(
            [sys.executable, "-c", "import sys; sys.stderr.write('boom\\n'); sys.exit(3)"],
            error_class=RuntimeError,
        )
    msg = str(ei.value)
    assert "[SDN-TOOL-002]" in msg
    assert "exited with status 3" in msg
    assert "boom" in msg  # the tool's own stderr is surfaced


# --- config loader ------------------------------------------------------------

def test_loader_malformed_yaml_carries_cfg_009(tmp_path: Path) -> None:
    """Unparseable YAML raises ConfigError tagged [SDN-CFG-009]."""
    bad = tmp_path / "bad.yaml"
    bad.write_text("marker:\n\tname: t\n")  # tab indentation is invalid YAML
    with pytest.raises(ConfigError) as ei:
        load_yaml(bad)
    msg = str(ei.value)
    assert "[SDN-CFG-009]" in msg
    assert "Invalid YAML" in msg


def test_loader_empty_config_is_not_tagged_malformed(tmp_path: Path) -> None:
    """An empty/comment-only config is VALID YAML (parses to None), not malformed, so it gets a
    clear 'empty' message but NOT the SDN-CFG-009 (malformed YAML) code."""
    empty = tmp_path / "empty.yaml"
    empty.write_text("# only a comment, no keys\n")
    with pytest.raises(ConfigError) as ei:
        load_yaml(empty)
    msg = str(ei.value)
    assert "empty" in msg
    assert "SDN-CFG-009" not in msg


# --- catalog round-trip -------------------------------------------------------

def test_wired_codes_resolve_via_explain() -> None:
    """Every code emitted above must be a real catalog entry, not an orphan."""
    for code in ("SDN-TOOL-001", "SDN-TOOL-002", "SDN-CFG-009"):
        assert explain(code) is not None, f"{code} has no extended explanation"
