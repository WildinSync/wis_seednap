"""Unit tests: the DADA2 taxonomy species DB is required, not optional.

The config schema (Dada2DatabaseConfig.species), the deployed docs, the
runner's error text, the assigner's runtime guard, and the R script
(addSpecies is always called) must all agree that the species-level DB is
required for the dada2 method. Previously the field was declared
Optional[Path]=None and the runner's FileNotFoundError told users they could
"remove the species key to skip species-level assignment" -- but no layer
actually supported skipping, so a config omitting it passed `seednap validate`
and then crashed mid-run with a misleading "required" error.

These tests pin the resolved behavior: omitting `species` from a dada2 DB
block is rejected at config-load time (fails before the fix, passes after),
and the runner's error text no longer advertises a skip path that does not
exist.
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from seednap.config.models.taxonomy import Dada2DatabaseConfig


def test_dada2_db_requires_species(tmp_path) -> None:
    """A dada2 DB block without `species` must fail validation.

    Before the fix `species` was Optional[Path]=None, so this constructed
    cleanly; after the fix it is a required field.
    """
    all_db = tmp_path / "all.fasta"
    all_db.write_text(">REF1\tk;p;c;o;f;g;s\nACGT\n")
    with pytest.raises(ValidationError):
        Dada2DatabaseConfig(all=all_db)


def test_dada2_db_with_species_ok(tmp_path) -> None:
    """A dada2 DB block that supplies both `all` and `species` validates."""
    all_db = tmp_path / "all.fasta"
    all_db.write_text(">REF1\tk;p;c;o;f;g;s\nACGT\n")
    species_db = tmp_path / "species.fasta"
    species_db.write_text(">REF1 k g s\nACGT\n")
    cfg = Dada2DatabaseConfig(all=all_db, species=species_db)
    assert cfg.species is not None


def test_runner_error_text_does_not_advertise_skip() -> None:
    """The species-DB FileNotFoundError no longer tells users they may skip
    species-level assignment by removing the key (no layer supports it)."""
    import inspect

    from seednap.steps.taxonomic_assignment import dada2_taxonomy_runner

    src = inspect.getsource(dada2_taxonomy_runner.Dada2TaxonomyRunner.run_dada2_taxonomy)
    assert "skip species-level assignment" not in src
    assert "remove the species key" not in src
