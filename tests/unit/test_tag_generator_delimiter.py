"""Unit tests for tag-file metadata loading delimiter handling.

Covers the lab-data usability fix in
``seednap.steps.trimming.tag_generator.TagFileGenerator``:

- semicolon-delimited metadata that already has the required header
  (eventID/tag_demultiplex/library) loads and produces tag files;
- a lab-style headerless Corr_tags file
  (well;library;sample;project;marker;tagseq) raises a clear error that
  names the required columns and the Corr_tags conversion mapping.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from seednap.steps.trimming.tag_generator import TagFileGenerator


def test_semicolon_metadata_with_header_loads(tmp_path: Path) -> None:
    """A semicolon-delimited metadata with the right header generates tag files."""
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "eventID;tag_demultiplex;library\n"
        "sampleA;ACGTACGT;LIB1\n"
        "sampleB;TGCATGCA;LIB1\n"
    )
    out_dir = tmp_path / "tags"

    gen = TagFileGenerator()
    tag_files = gen.generate_ligation_tag_files(
        metadata_csv=metadata, output_dir=out_dir
    )

    assert "LIB1" in tag_files
    fasta_text = Path(tag_files["LIB1"]).read_text()
    # Both samples written; tag sequence carried through into the adapter spec.
    assert ">sampleA" in fasta_text
    assert ">sampleB" in fasta_text
    assert "ACGTACGT" in fasta_text


def test_lab_corr_tags_headerless_raises_clear_guidance(tmp_path: Path) -> None:
    """A headerless lab Corr_tags file raises an error naming the required columns."""
    # Lab Corr_tags layout: well;library;sample;project;marker;tagseq, no header.
    corr_tags = tmp_path / "Corr_tags.csv"
    corr_tags.write_text(
        "A1;LIB1;sampleA;PROJ;teleo;ACGTACGT\n"
        "A2;LIB1;sampleB;PROJ;teleo;TGCATGCA\n"
    )
    out_dir = tmp_path / "tags"

    gen = TagFileGenerator()
    with pytest.raises(ValueError) as excinfo:
        gen.generate_ligation_tag_files(metadata_csv=corr_tags, output_dir=out_dir)

    msg = str(excinfo.value)
    # Names the required columns.
    assert "eventID" in msg
    assert "tag_demultiplex" in msg
    assert "library" in msg
    # Points at the Corr_tags conversion explicitly (no silent guessing).
    assert "Corr_tags" in msg
    assert "tagseq" in msg
