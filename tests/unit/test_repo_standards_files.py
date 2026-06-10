"""Smoke tests for the repository-standards meta files.

Guards CONTRIBUTING.md, CHANGELOG.md, CODEOWNERS, and LICENSE: they must exist
at the repo root and carry the key, load-bearing content (commit-prefix
convention, Keep-a-Changelog [Unreleased] section, a default CODEOWNERS owner,
and an MIT LICENSE). These are cheap structural checks, not prose review.
"""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(name: str) -> str:
    path = REPO_ROOT / name
    assert path.is_file(), f"{name} must exist at the repo root"
    return path.read_text(encoding="utf-8")


def test_license_exists_and_is_mit() -> None:
    """pyproject declares MIT, so a matching LICENSE must be present."""
    text = _read("LICENSE")
    assert "MIT License" in text


def test_contributing_covers_key_topics() -> None:
    """CONTRIBUTING must mention the supported dev workflow commands."""
    text = _read("CONTRIBUTING.md")
    for token in ("conda", "pytest", "ruff check", "mypy", "pip install -e ."):
        assert token in text, f"CONTRIBUTING.md should mention '{token}'"
    # The shared server env path is the documented way in on the eDNA server.
    assert "/home/shared/edna/envs/seednap" in text


@pytest.mark.parametrize(
    "prefix", ["[FIX]", "[FEAT]", "[REFACTOR]", "[DOCS]", "[CONFIG]", "[TEST]"]
)
def test_contributing_documents_commit_prefixes(prefix: str) -> None:
    """Every commit-title prefix the repo uses must be documented."""
    text = _read("CONTRIBUTING.md")
    assert prefix in text


def test_changelog_keep_a_changelog_shape() -> None:
    """CHANGELOG must follow Keep-a-Changelog and have an [Unreleased] section."""
    text = _read("CHANGELOG.md")
    assert "Keep a Changelog" in text
    assert "[Unreleased]" in text
    # At least one of the standard section headers should be present.
    assert any(h in text for h in ("### Added", "### Changed", "### Fixed"))


def test_codeowners_has_default_org_owner() -> None:
    """CODEOWNERS must map a catch-all to the GitHub org handle."""
    text = _read("CODEOWNERS")
    lines = [
        ln.strip()
        for ln in text.splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    assert lines, "CODEOWNERS must have at least one non-comment rule"
    assert any(ln.startswith("*") and "@WildinSync" in ln for ln in lines)
