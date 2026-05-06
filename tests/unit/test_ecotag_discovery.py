"""Unit tests for OBITools discovery in EcotagRunner (Issue #3 fix).

Validates that:
  - explicit `bin_dir` works (and rejects directories that don't have all 3 tools)
  - SEEDNAP_OBITOOLS_BIN env var is honoured
  - well-known install locations are probed
  - missing OBITools raises a clear, actionable error pointing at the docs

No real OBITools binaries are needed; we synthesise stub executables on disk.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path
from typing import List

import pytest

from seednap.steps.taxonomic_assignment.ecotag_runner import (
    EcotagError,
    EcotagRunner,
    _find_obitools_bin,
)


def _make_stub_bin(directory: Path, names: List[str]) -> None:
    """Create no-op shell scripts named like OBITools binaries."""
    directory.mkdir(parents=True, exist_ok=True)
    for n in names:
        f = directory / n
        f.write_text("#!/bin/sh\necho 'stub OBITools'\n")
        f.chmod(f.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)


def test_explicit_bin_dir_succeeds(tmp_path: Path) -> None:
    bin_dir = tmp_path / "obibin"
    _make_stub_bin(bin_dir, ["ecotag", "obiannotate", "obitab"])
    runner = EcotagRunner(bin_dir=bin_dir)
    assert runner.bin_dir == bin_dir
    assert runner._tool("ecotag") == str(bin_dir / "ecotag")


def test_explicit_bin_dir_missing_tool_raises(tmp_path: Path) -> None:
    bin_dir = tmp_path / "obibin"
    _make_stub_bin(bin_dir, ["ecotag", "obiannotate"])  # missing obitab
    with pytest.raises(EcotagError, match="missing.*obitab|does not contain"):
        EcotagRunner(bin_dir=bin_dir)


def test_env_var_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bin_dir = tmp_path / "obibin"
    _make_stub_bin(bin_dir, ["ecotag", "obiannotate", "obitab"])
    monkeypatch.setenv("SEEDNAP_OBITOOLS_BIN", str(bin_dir))
    # Force PATH lookup to fail by clearing PATH so only the env override matches.
    monkeypatch.setenv("PATH", "/nonexistent")
    found = _find_obitools_bin()
    assert found == bin_dir


def test_missing_obitools_raises_with_install_instructions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No OBITools anywhere -> error must mention conda + env var + bin paths."""
    monkeypatch.setenv("PATH", str(tmp_path))  # empty PATH
    monkeypatch.delenv("SEEDNAP_OBITOOLS_BIN", raising=False)
    # Patch the candidate list to a single nonexistent dir so we don't accidentally
    # discover the real OBITools install on the host.
    import seednap.steps.taxonomic_assignment.ecotag_runner as er

    monkeypatch.setattr(er, "_OBITOOLS_CANDIDATE_BINS", [str(tmp_path / "absent")])
    with pytest.raises(EcotagError) as exc_info:
        EcotagRunner()
    msg = str(exc_info.value)
    # The error message must give the user a path forward
    assert "OBITools not found" in msg
    assert "SEEDNAP_OBITOOLS_BIN" in msg
    assert "conda activate obitools" in msg
    assert "Probed locations" in msg
