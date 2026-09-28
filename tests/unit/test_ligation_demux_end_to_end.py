"""Ligation demultiplex -> trim -> read tracking on a synthetic library, with real cutadapt.

Checks that the per-sample outputs are gzipped and correctly oriented, that no
intermediate file survives, that the cutadapt reports are kept, and that read
tracking follows each sample from its demultiplexed reads through primer
detection to trimming.
"""

import gzip
import shutil

import pytest

from seednap.steps.report import ReadTrackingBuilder
from seednap.steps.trimming import LigationTrimmer, StandardTrimmer
from seednap.utils.sequences import reverse_complement

pytestmark = pytest.mark.skipif(shutil.which("cutadapt") is None, reason="cutadapt not on PATH")

FWD = "ACACCGCCCGTCACTCT"
REV = "CTTCCGGTACACTTACCATG"
INSERT = "GATTAGATACCCCACTATGCTTAGCCCTAAACATAGATAATTTTACAACAAAATAATTCGCCAGAG"
TAGS = {"S1": "ACACACAC", "S2": "AGCTAGCT"}
UNKNOWN_TAG = "TTTTGGGG"
# (sample tag, pairs in the expected orientation, pairs in the reverse orientation)
LIBRARY = [("S1", 3, 2), ("S2", 4, 0)]
N_UNKNOWN = 5


def _pair(tag, reverse_orientation):
    fragment = tag + FWD + INSERT + reverse_complement(REV) + reverse_complement(tag)
    r1, r2 = fragment, reverse_complement(fragment)
    return (r2, r1) if reverse_orientation else (r1, r2)


def _write_fastq_gz(path, seqs):
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(seqs):
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def _read_seqs(path):
    with gzip.open(path, "rt") as fh:
        return fh.read().splitlines()[1::4]


@pytest.fixture
def library(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    pairs = []
    for sample, n_fwd, n_rev in LIBRARY:
        pairs += [_pair(TAGS[sample], False)] * n_fwd + [_pair(TAGS[sample], True)] * n_rev
    pairs += [_pair(UNKNOWN_TAG, False)] * N_UNKNOWN
    _write_fastq_gz(raw / "LIB1_S1_L001_R1.fastq.gz", [p[0] for p in pairs])
    _write_fastq_gz(raw / "LIB1_S1_L001_R2.fastq.gz", [p[1] for p in pairs])
    meta = tmp_path / "meta.csv"
    meta.write_text(
        "eventID,tag_demultiplex,library\n"
        + "".join(f"{s},{t},LIB1\n" for s, t in TAGS.items())
        + "X1,GGGGCCCC,OTHER_LIB\n"
    )
    return raw, meta


def test_demux_trim_and_read_tracking(tmp_path, library):
    raw, meta = library
    demux = tmp_path / "demux"
    trimmer = LigationTrimmer()
    tag_files = trimmer.generate_tag_files(meta, demux / "cutadapt_tags", libraries=["LIB1"])
    assert list(tag_files) == ["LIB1"]

    written = trimmer.process_library(
        raw_reads_dir=raw,
        library_name="LIB1",
        tag_file=trimmer.library_tag_file(tag_files, "LIB1", meta),
        output_dir=demux / "samples",
        work_dir=demux / "work" / "LIB1",
        log_dir=demux / "logs",
        forward_primer=FWD,
        reverse_primer=REV,
    )

    assert sorted(written) == ["S1", "S2"]
    assert not (demux / "work" / "LIB1").exists()
    assert sorted(p.name for p in (demux / "samples").iterdir()) == [
        "S1.R1.fastq.gz", "S1.R2.fastq.gz", "S2.R1.fastq.gz", "S2.R2.fastq.gz",
    ]
    # Both orientations end up forward-strand in R1, reverse-strand in R2.
    s1_r1 = _read_seqs(demux / "samples" / "S1.R1.fastq.gz")
    s1_r2 = _read_seqs(demux / "samples" / "S1.R2.fastq.gz")
    assert len(s1_r1) == len(s1_r2) == 5
    assert all(seq.startswith(FWD) for seq in s1_r1)
    assert all(seq.startswith(REV) for seq in s1_r2)
    assert sorted(p.name for p in (demux / "logs").iterdir()) == [
        "LIB1_demultiplex.txt",
        "S1_primer_round1.txt", "S1_primer_round2.txt",
        "S2_primer_round1.txt", "S2_primer_round2.txt",
    ]

    trim_dir = tmp_path / "trim"
    for sample in written:
        r1, r2 = StandardTrimmer().trim_sample(
            demux / "samples" / f"{sample}.R1.fastq.gz",
            demux / "samples" / f"{sample}.R2.fastq.gz",
            trim_dir, sample, FWD, REV,
        )
        assert r1.name.endswith(".fastq.gz") and r2.name.endswith(".fastq.gz")
        assert all(seq == INSERT for seq in _read_seqs(r1))
    assert not list(trim_dir.glob("*TEMPORARY*"))

    builder = ReadTrackingBuilder(
        marker="teleo", logs_dir=trim_dir / "logs", demux_logs_dir=demux / "logs"
    )
    assert builder.steps == ["raw", "primer_found", "trimmed"]
    df = builder.build().set_index("sample")
    assert df.loc["S1", ["raw", "primer_found", "trimmed"]].tolist() == [5, 5, 5]
    assert df.loc["S2", ["raw", "primer_found", "trimmed"]].tolist() == [4, 4, 4]

    summary = builder.demux_summary()
    assert summary.to_dict("records") == [
        {"library": "LIB1", "read_pairs": 14, "assigned": 9, "pct_assigned": 64.29}
    ]
    assert builder.write_demux_summary(tmp_path / "report").name == "demux_summary.csv"


def test_read_tracking_without_demux_is_unchanged(tmp_path):
    builder = ReadTrackingBuilder(
        marker="teleo", logs_dir=tmp_path, demux_logs_dir=tmp_path / "missing"
    )
    assert builder.steps == ["raw", "trimmed"]
    assert builder.demux_summary().empty
    assert builder.write_demux_summary(tmp_path) is None
