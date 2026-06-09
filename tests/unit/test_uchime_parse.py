"""OtuTableBuilder._parse_uchime: parse Y/N/? status, skip blanks quietly, never drop silently.

UCHIME --uchimeout columns: parts[1] = query label, parts[17] = chimera flag. A non-blank line
with no query-label column must be logged with a [WARN] (no-silent-fallback rule), not dropped
silently; blank lines are skipped quietly; a labelled line missing the flag column defaults "NA".
"""

import logging

from seednap.steps.swarm.otu_table_builder import OtuTableBuilder


def _line(label: str, status: str) -> str:
    """An 18-column UCHIME row: label at index 1, status at index 17."""
    fields = ["0.0", label] + ["x"] * 15 + [status]
    return "\t".join(fields)


def test_parse_uchime_status_blank_and_malformed(tmp_path, caplog):
    path = tmp_path / "uchime.txt"
    path.write_text(
        _line("seqA;size=5", "Y") + "\n"
        + _line("seqB;size=3", "N") + "\n"
        + "\n"                       # blank: skipped quietly
        + "0.0\tseqC;size=2\n"      # labelled but no status column -> NA
        + "garbage\n"                # malformed: no label column -> [WARN] + skip
    )
    with caplog.at_level(logging.WARNING):
        result = OtuTableBuilder._parse_uchime(path)

    assert result == {"seqA": "Y", "seqB": "N", "seqC": "NA"}
    warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert any("garbage" in w and "_parse_uchime" in w for w in warnings)
    # the blank line produced no warning (skipped quietly)
    assert sum("_parse_uchime" in w for w in warnings) == 1
