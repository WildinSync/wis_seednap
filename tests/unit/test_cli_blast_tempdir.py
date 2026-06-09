"""Unit test for the `blast` CLI command's scratch-directory handling.

Regression guard: the command used to derive its scratch dir as
`query_fasta.parent / "blast_temp"`, create it with mkdir(exist_ok=True) (thus
ADOPTING any pre-existing directory of that name), and on success unconditionally
rmtree it. A user directory that happened to be named `blast_temp` next to the
query FASTA was then deleted with its contents.

The fix uses tempfile.mkdtemp, so the command only ever removes a directory it
created in this invocation. This test stubs the heavy BLAST collaborators (no
blastn binary needed) and asserts a pre-existing `blast_temp` with user data
survives the run.
"""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pandas as pd
from click.testing import CliRunner

from seednap.cli import main


def test_blast_command_does_not_delete_preexisting_blast_temp(tmp_path: Path) -> None:
    query = tmp_path / "query.fasta"
    query.write_text(">ASV_1\nACGTACGTACGT\n")
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref1\nACGTACGTACGT\n")
    counts = tmp_path / "counts.csv"
    pd.DataFrame({"S1": [10]}, index=pd.Index(["ACGTACGTACGT"], name="sequence")).to_csv(
        counts
    )

    # Pre-existing user directory that collides with the old hardcoded name.
    user_dir = tmp_path / "blast_temp"
    user_dir.mkdir()
    sentinel = user_dir / "important_user_data.txt"
    sentinel.write_text("do not delete me")

    # Stub the BLAST runner + assigner so no blastn binary is needed. The runner
    # just writes a dummy TSV into whatever scratch dir the command made; the
    # assigner returns an empty-but-valid result frame.
    fake_runner = mock.MagicMock()

    def _fake_pipeline(query_fasta, db_fasta, output_dir, marker):
        tsv = Path(output_dir) / "blast.tsv"
        tsv.write_text("")
        return tsv

    fake_runner.run_blast_pipeline.side_effect = _fake_pipeline

    fake_assigner = mock.MagicMock()
    result_df = pd.DataFrame(
        columns=["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    )
    fake_assigner.assign_taxonomy.return_value = result_df

    with mock.patch(
        "seednap.steps.taxonomic_assignment.BlastRunner", return_value=fake_runner
    ), mock.patch(
        "seednap.steps.taxonomic_assignment.BlastTaxonomicAssigner",
        return_value=fake_assigner,
    ):
        runner = CliRunner()
        res = runner.invoke(
            main,
            ["blast", str(query), str(ref), str(counts)],
            catch_exceptions=False,
        )

    assert res.exit_code == 0, res.output
    # The user's pre-existing blast_temp directory and its file must survive.
    assert user_dir.is_dir()
    assert sentinel.exists()
    assert sentinel.read_text() == "do not delete me"
