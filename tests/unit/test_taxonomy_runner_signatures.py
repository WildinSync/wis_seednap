"""Issue #2 fix: the DADA2 taxonomy runner now accepts query_fasta and does not
require seqtab_clean.rds.

These tests don't actually run R -- they verify the Python wrapper
signatures and the input-validation path so that:

  - calling without query_fasta is impossible (Python TypeError)
  - missing query_fasta raises a clear FileNotFoundError, not a cryptic R
    error 5 minutes into execution
  - the runner no longer enforces an seqtab_clean.rds presence check

Real R-side e2e validation lives in the validation report (Phase B for
DADA2 ASV->RDP, plus the new SWARM->RDP cross-mode test).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from seednap.steps.taxonomic_assignment.dada2_taxonomy_runner import (
    Dada2TaxonomyError,
    Dada2TaxonomyRunner,
)


@pytest.fixture
def fake_inputs(tmp_path: Path):
    rdp = tmp_path / "rdp.fasta"
    rdp.write_text(">REF1\tk;p;c;o;f;g;s\nACGT\n")
    species = tmp_path / "species.fasta"
    species.write_text(">REF1 k g s\nACGT\n")
    query = tmp_path / "query.fasta"
    query.write_text(">OTU_1\nACGTACGT\n")
    return {"rdp": rdp, "species": species, "query": query, "out": tmp_path / "out"}


def test_dada2_runner_requires_query_fasta() -> None:
    """The new signature makes query_fasta mandatory."""
    runner = Dada2TaxonomyRunner.__new__(Dada2TaxonomyRunner)
    with pytest.raises(TypeError):
        runner.run_dada2_taxonomy(
            marker="teleo",
            output_dir="/tmp",
            rdp_db_path="/tmp/rdp.fa",
            species_db_path="/tmp/species.fa",
        )


def test_dada2_runner_missing_query_raises_filenotfound(fake_inputs, tmp_path: Path) -> None:
    """Missing query.fasta is caught up-front, not 5 minutes into R execution."""
    runner = Dada2TaxonomyRunner.__new__(Dada2TaxonomyRunner)
    runner.timeout = 60
    with pytest.raises(FileNotFoundError, match="Query FASTA not found"):
        runner.run_dada2_taxonomy(
            marker="teleo",
            output_dir=fake_inputs["out"],
            rdp_db_path=fake_inputs["rdp"],
            species_db_path=fake_inputs["species"],
            query_fasta=tmp_path / "does_not_exist.fasta",
        )


def test_dada2_runner_no_longer_checks_seqtab_rds(fake_inputs, tmp_path: Path) -> None:
    """Issue #2 fix: pre-fix runner blocked SWARM->DADA2-RDP because it
    insisted on 02_dada2/{marker}/seqtab_clean.rds. Now only query.fasta
    matters; the runner gets past validation and only fails because we
    point it at a non-existent R script.
    """
    runner = Dada2TaxonomyRunner.__new__(Dada2TaxonomyRunner)
    runner.timeout = 60
    bogus_script = tmp_path / "absolutely_not_an_R_script.R"
    with pytest.raises((FileNotFoundError, Dada2TaxonomyError)) as exc_info:
        runner.run_dada2_taxonomy(
            marker="teleo",
            output_dir=fake_inputs["out"],
            rdp_db_path=fake_inputs["rdp"],
            species_db_path=fake_inputs["species"],
            query_fasta=fake_inputs["query"],
            script_path=bogus_script,
        )
    msg = str(exc_info.value)
    assert "seqtab_clean.rds" not in msg, (
        f"Issue #2 not fixed: error still mentions seqtab_clean.rds: {msg}"
    )

