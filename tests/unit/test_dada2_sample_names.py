"""Unit test for sample-name extraction in scripts/dada2_process.R.

The trim step writes "<sample>.R1.fastq" / "<sample>.R2.fastq". DADA2 derives
the per-sample name from those filenames and uses it as the abundance-table
column header (seqtab_clean_t.csv) and the track_reads.csv label.

The old code split on the first "." (strsplit(basename(x), "\\.")[1]), which
silently truncated any sample whose name itself contains a "." and, worse,
collided two distinct samples sharing a prefix before the first dot
("Site.A" and "Site.B" both -> "Site"), fusing/mislabelling per-sample counts
in a GBIF-bound dataset. The fix strips the read-suffix token instead.

This test runs the EXACT extraction expression as it appears in the R script
through Rscript, so it stays coupled to the real fix and fails on a regression
back to the dot-split (or any change that re-introduces the collision).
"""

from __future__ import annotations

import re
import shutil
import subprocess

import pytest

# The R scripts ship inside the package; resolve via the same helper the runners use
# (robust to the package layout, not a hardcoded repo-root path).
from seednap.utils.r_runner import r_script_path

R_SCRIPT = r_script_path("dada2_process.R")

# The bug cases: distinct sample names that the old dot-split collapsed/collided.
SAMPLES_TO_FILES = {
    "Site.A": "Site.A.R1.fastq",
    "Site.B": "Site.B.R1.fastq",
    "Sample01": "Sample01.R1.fastq",
    "Blank-ext": "Blank-ext.R1.fastq.gz",
    "BB_Unknown_Svalbard": "BB_Unknown_Svalbard.R1.fastq",
}


def _extraction_regex_from_script() -> str:
    """Pull the read-suffix sub() pattern out of the R script verbatim."""
    text = R_SCRIPT.read_text(encoding="utf-8")
    m = re.search(r'sub\("(\\\\\.\[Rr\]\[12\][^"]*)",\s*""', text)
    assert m, "could not find the read-suffix sub() pattern in dada2_process.R"
    return m.group(1)


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not installed")
def test_sample_names_strip_read_suffix_not_first_dot():
    pattern = _extraction_regex_from_script()
    files = list(SAMPLES_TO_FILES.values())
    expected = list(SAMPLES_TO_FILES.keys())

    # Run the same extraction R uses, on the bug-case filenames.
    r_files = "c(" + ", ".join(f'"{f}"' for f in files) + ")"
    r_expr = f'cat(sub("{pattern}", "", basename({r_files})), sep="\\n")'
    proc = subprocess.run(
        ["Rscript", "-e", r_expr],
        capture_output=True,
        text=True,
        check=True,
    )
    got = [line for line in proc.stdout.splitlines() if line != ""]

    assert got == expected, f"expected {expected}, got {got}"
    # The load-bearing property: the two dotted names stay distinct (no collision).
    assert "Site.A" in got and "Site.B" in got
    assert len(set(got)) == len(got), "sample names collided after extraction"
