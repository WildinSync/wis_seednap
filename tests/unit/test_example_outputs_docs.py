"""Structural checks for the committed worked-example outputs under docs/.

The full pipeline outputs live under outputs/ (gitignored), so a fresh clone
cannot see what a finished run produces. docs/example-outputs/ ships a README
plus trimmed copies of the four headline tables. These tests guard that the
trimmed copies stay well-formed (a real header plus a few rows, consistent
column counts, small file sizes) and that the README links to them resolve.
"""

import csv
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = REPO_ROOT / "docs" / "example-outputs"

EXPECTED_CSVS = [
    "1_read_tracking.csv",
    "2_otu_table.csv",
    "3_taxonomy_blast.csv",
    "4_sample_manifest_FAIRe.csv",
]

# Trimmed copies must stay small so they survive a clone cheaply.
MAX_CSV_BYTES = 16 * 1024


def _read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return [row for row in csv.reader(handle) if row]


def test_example_dir_and_readme_present() -> None:
    """The example directory ships with a README so a clone can read it."""
    assert EXAMPLE_DIR.is_dir(), f"missing {EXAMPLE_DIR}"
    readme = EXAMPLE_DIR / "README.md"
    assert readme.is_file(), "docs/example-outputs/README.md must be committed"
    text = readme.read_text(encoding="utf-8")
    for name in EXPECTED_CSVS:
        assert name in text, f"README does not mention {name}"


@pytest.mark.parametrize("name", EXPECTED_CSVS, ids=lambda n: n)
def test_trimmed_csv_is_well_formed(name: str) -> None:
    """Each trimmed CSV is a real header plus a few rows, all the same width."""
    path = EXAMPLE_DIR / name
    assert path.is_file(), f"missing {path}"
    assert path.stat().st_size <= MAX_CSV_BYTES, f"{name} is too large to be a trimmed sample"

    rows = _read_rows(path)
    assert len(rows) >= 2, f"{name} must have a header and at least one data row"
    # Trimmed to roughly five representative rows; keep it small.
    assert len(rows) <= 7, f"{name} should be a trimmed excerpt, got {len(rows)} rows"

    header_width = len(rows[0])
    assert header_width >= 2, f"{name} header looks empty"
    for i, row in enumerate(rows[1:], start=1):
        assert len(row) == header_width, (
            f"{name} row {i} has {len(row)} columns, header has {header_width}"
        )


def test_readme_links_to_example_outputs() -> None:
    """The top-level README points readers at the worked example."""
    text = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    assert "docs/example-outputs/" in text


def test_pipeline_steps_links_to_example_outputs() -> None:
    """docs/pipeline-steps.md points readers at the worked example."""
    text = (REPO_ROOT / "docs" / "pipeline-steps.md").read_text(encoding="utf-8")
    assert "example-outputs/" in text
