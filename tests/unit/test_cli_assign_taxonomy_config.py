"""Unit tests for `assign-taxonomy --config`.

Regression guard for the bug where the standalone `assign-taxonomy` command ignored
the marker config and ran BLAST with the assigner's hardcoded defaults (e.g.
evalue=1e-25), diverging from what `run-pipeline` does for the same data (e.g. a
config with evalue=1e-10, task=megablast, max_target_seqs=10).

The fix adds `--config <marker.yaml>`: the command loads that marker's
taxonomy.databases.<method> block and threads its parameters through to
TaxonomicAssigner.assign_taxonomy, while keeping the explicit CLI options as
overrides. These tests stub TaxonomicAssigner (no blastn binary needed) and assert
the kwargs the CLI passes carry the config values, not the hardcoded defaults.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict
from unittest import mock

import pandas as pd
from click.testing import CliRunner

from seednap.cli import main

# Config blast params chosen to differ from the assigner's hardcoded defaults
# (evalue 1e-25, perc_identity 80, task megablast, max_target_seqs 5,
# threshold_species 98/99) so a passthrough is provably from the config.
_CONFIG_EVALUE = 1e-10
_CONFIG_TASK = "blastn"
_CONFIG_MAX_TARGET_SEQS = 10
_CONFIG_PERC_IDENTITY = 75.0
_CONFIG_THRESHOLD_SPECIES = 97.5


def _write_inputs(tmp_path: Path) -> Dict[str, Path]:
    """Write a query FASTA, an ASV count CSV, and a BLAST marker config."""
    query = tmp_path / "query.fasta"
    query.write_text(">ASV_1\nACGTACGTACGT\n")

    counts = tmp_path / "counts.csv"
    pd.DataFrame(
        {"S1": [10]}, index=pd.Index(["ACGTACGTACGT"], name="sequence")
    ).to_csv(counts)

    config = tmp_path / "teleo.yaml"
    config.write_text(
        "marker:\n"
        "  name: teleo\n"
        "  primers:\n"
        '    forward: "ACACCGCCCGTCACTCT"\n'
        '    reverse: "CTTCCGGTACACTTACCATG"\n'
        "paths:\n"
        '  raw_data: "data/raw"\n'
        "taxonomy:\n"
        '  method: "blast"\n'
        "  databases:\n"
        "    blast:\n"
        '      fasta: "references/teleo/blast_db.fasta"\n'
        f'      task: "{_CONFIG_TASK}"\n'
        f"      evalue: {_CONFIG_EVALUE}\n"
        f"      perc_identity: {_CONFIG_PERC_IDENTITY}\n"
        f"      max_target_seqs: {_CONFIG_MAX_TARGET_SEQS}\n"
        f"      threshold_species: {_CONFIG_THRESHOLD_SPECIES}\n"
    )
    return {"query": query, "counts": counts, "config": config}


def _run_with_captured_assigner(tmp_path: Path, extra_args: list) -> Dict[str, Any]:
    """Invoke assign-taxonomy with TaxonomicAssigner stubbed; return captured kwargs."""
    inputs = _write_inputs(tmp_path)

    captured: Dict[str, Any] = {}
    fake_assigner = mock.MagicMock()

    def _capture(query_fasta, asv_count_csv, **kwargs):
        captured.update(kwargs)
        return {"final_table": tmp_path / "teleo_blast.csv"}

    fake_assigner.assign_taxonomy.side_effect = _capture

    with mock.patch(
        "seednap.steps.taxonomic_assignment.TaxonomicAssigner",
        return_value=fake_assigner,
    ):
        runner = CliRunner()
        res = runner.invoke(
            main,
            [
                "assign-taxonomy",
                "blast",
                "teleo",
                str(inputs["query"]),
                str(inputs["counts"]),
                "--output-dir",
                str(tmp_path / "out"),
                *extra_args,
            ],
            catch_exceptions=False,
        )
    assert res.exit_code == 0, res.output
    assert captured, "TaxonomicAssigner.assign_taxonomy was never called"
    return captured


def test_config_blast_params_passed_through(tmp_path: Path) -> None:
    """With --config, the config's BLAST params reach the assigner, not the defaults."""
    inputs = _write_inputs(tmp_path)
    captured = _run_with_captured_assigner(
        tmp_path, ["--config", str(inputs["config"])]
    )

    # The bug: evalue used to be the hardcoded assigner default (1e-25). It must now
    # be the config value.
    assert captured["evalue"] == _CONFIG_EVALUE
    assert captured["evalue"] != 1e-25

    assert captured["task"] == _CONFIG_TASK
    assert captured["max_target_seqs"] == _CONFIG_MAX_TARGET_SEQS
    assert captured["perc_identity"] == _CONFIG_PERC_IDENTITY
    assert captured["threshold_species"] == _CONFIG_THRESHOLD_SPECIES
    # reference_fasta comes from the config's blast.fasta (resolved absolute path).
    assert Path(captured["reference_fasta"]).name == "blast_db.fasta"


def test_explicit_option_overrides_config(tmp_path: Path) -> None:
    """An explicitly-passed CLI option wins over the config value for the same param."""
    inputs = _write_inputs(tmp_path)
    captured = _run_with_captured_assigner(
        tmp_path,
        ["--config", str(inputs["config"]), "--threshold-species", "91.0"],
    )

    # CLI override beats config for threshold_species ...
    assert captured["threshold_species"] == 91.0
    # ... while params NOT given on the CLI still come from the config.
    assert captured["evalue"] == _CONFIG_EVALUE
    assert captured["task"] == _CONFIG_TASK


def test_without_config_evalue_stays_default(tmp_path: Path) -> None:
    """Without --config the command keeps its prior behavior: no config-derived blast
    search params (evalue/task/etc.) are injected, so the assigner default applies.

    --reference-fasta is supplied because BLAST still requires a reference; the point is
    that the search-level knobs (evalue/task/perc_identity) are NOT passed without a
    config (the standalone command never exposed them as options)."""
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref1\nACGTACGTACGT\n")
    captured = _run_with_captured_assigner(
        tmp_path, ["--reference-fasta", str(ref)]
    )

    # The standalone command has no --evalue/--task options, and without --config it
    # must not pass them, leaving the assigner's own default (the prior behavior).
    assert "evalue" not in captured
    assert "task" not in captured
    assert "perc_identity" not in captured


def test_config_method_mismatch_errors(tmp_path: Path) -> None:
    """A config whose taxonomy.method differs from the requested METHOD must error,
    not silently apply the wrong method's parameters (no-silent-fallback rule)."""
    inputs = _write_inputs(tmp_path)  # config taxonomy.method == "blast"

    with mock.patch(
        "seednap.steps.taxonomic_assignment.TaxonomicAssigner"
    ) as fake_cls:
        runner = CliRunner()
        res = runner.invoke(
            main,
            [
                "assign-taxonomy",
                "decipher",  # mismatched against the blast config
                "teleo",
                str(inputs["query"]),
                str(inputs["counts"]),
                "--output-dir",
                str(tmp_path / "out"),
                "--config",
                str(inputs["config"]),
            ],
            catch_exceptions=False,
        )

    assert res.exit_code == 1, res.output
    assert "blast" in res.output and "decipher" in res.output
    fake_cls.return_value.assign_taxonomy.assert_not_called()
