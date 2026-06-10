"""R scripts now ship inside the importable ``seednap`` package.

The three R scripts moved from a repo-root ``scripts/`` directory (a sibling
of ``src/``, never bundled in a wheel) into ``seednap/scripts/`` so they ride
along in both editable and wheel installs. ``r_runner`` now resolves them via
``importlib.resources`` instead of walking ``__file__`` up to the repo root.

These tests pin that packaging contract:

  - all three bundled R scripts resolve to existing files,
  - the scripts directory is anchored inside the installed ``seednap`` package
    (not the old repo-root ``scripts/``),
  - a missing script raises a clear FileNotFoundError whose wording reflects
    "broken/incomplete install", not the dropped "not shipped in a wheel /
    only present when editable" explanation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from seednap.utils.r_runner import (
    SCRIPTS_DIR,
    RScriptError,
    RScriptRunner,
    r_script_path,
)

EXPECTED_SCRIPTS = (
    "dada2_process.R",
    "taxo_dada2_marker.R",
    "taxo_decipher_marker.R",
)


def test_all_bundled_scripts_resolve_and_exist() -> None:
    """Every R script the runners reference resolves to a real file."""
    for name in EXPECTED_SCRIPTS:
        path = r_script_path(name)
        assert path.exists(), f"bundled R script missing: {path}"
        assert path.name == name


def test_scripts_dir_is_inside_seednap_package() -> None:
    """The scripts directory ships inside the package, not at the repo root.

    Anchoring is via importlib.resources, so the resolved directory must be
    ``.../seednap/scripts`` rather than a sibling-of-src ``scripts/``.
    """
    assert SCRIPTS_DIR.name == "scripts"
    assert SCRIPTS_DIR.parent.name == "seednap"
    assert SCRIPTS_DIR.is_dir()
    resolved_names = sorted(p.name for p in SCRIPTS_DIR.iterdir() if p.suffix == ".R")
    assert resolved_names == sorted(EXPECTED_SCRIPTS)


def test_missing_script_error_wording_reflects_broken_install(tmp_path: Path) -> None:
    """A missing packaged script now means a broken install, not "no wheel".

    The old message told users the scripts were not shipped in a wheel and
    only present in editable installs. Now that they are packaged, that
    wording is gone; the message must point at a broken/incomplete install.
    """
    runner = RScriptRunner.__new__(RScriptRunner)
    runner.timeout = 60
    runner._error_class = RScriptError
    bogus = tmp_path / "not_a_real_script.R"
    with pytest.raises(FileNotFoundError) as exc_info:
        runner._run_r_script(bogus, args=[])
    msg = str(exc_info.value)
    lower = msg.lower()
    # The dropped wording claimed the scripts were NOT shipped in a wheel and
    # were only present in editable installs. Both stale phrasings must be gone.
    assert "not shipped in a wheel" not in lower, f"stale wheel wording leaked: {msg}"
    assert "only present when" not in lower, f"stale editable wording leaked: {msg}"
    # The message now frames a missing script as a broken/incomplete install.
    assert "broken or incomplete" in lower
