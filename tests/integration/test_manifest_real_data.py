"""Real-data integration test for the FAIRe manifest migrator + cross-CSV validator.

Synthetic unit tests prove the logic; this proves the *approach works on the actual lab
data* by migrating a real dataset and asserting known-true facts about it (dotted dates,
Blank-ext/Blank-PCR controls, the extraction batch, and the documented Blank-PCR-3 orphan
in the OTU table). It is skipped where the server data tree is not present, so it never
breaks CI elsewhere.

All paths under /home/shared/edna/workflows are READ-ONLY reference data.
"""

from pathlib import Path

import pytest

from seednap.config.manifest import validate_against_abundance
from seednap.config.manifest_migrate import migrate_to_manifest

# Real Rhône 2025 teleo dataset (lab metadata + the seednap OTU table from the run).
_FIELD = Path(
    "/home/shared/edna/workflows/fw_ch_rhone_2025/teleo/create_df_gbif/data/"
    "metadata_field_fw_ch_rhone_2025.csv"
)
_PROJ = Path(
    "/home/shared/edna/workflows/fw_ch_rhone_2025/teleo/create_df_gbif/data/"
    "metadata_proj_fw_ch_rhone_2025.csv"
)
_OTU = Path(
    "/home/shared/edna/seednap/outputs_test/teleo_rhone/02_swarm/teleo_rhone/otu_table.csv"
)

pytestmark = pytest.mark.skipif(
    not (_FIELD.exists() and _PROJ.exists()),
    reason="real Rhône lab metadata not present (server-only data)",
)


def test_migrate_real_rhone_teleo():
    m = migrate_to_manifest(_FIELD, project_csv=_PROJ)

    # 18 biological samples + Blank-ext-1 + Blank-PCR-1 = 20 rows, 2 controls.
    assert len(m) == 20
    assert len(m.biological_samples()) == 18
    assert len(m.controls()) == 2

    rows = {r.eventID: r for r in m.rows}

    # dotted lab date 2025.08.19 -> ISO-8601
    assert rows["DAR-2025-1103"].eventDate == "2025-08-19"
    # volume -> samp_size; marker -> target_gene
    assert rows["DAR-2025-1103"].samp_size == pytest.approx(15.0)
    assert rows["DAR-2025-1103"].target_gene == "teleo"

    # control identity (two-field FAIRe model) + the extraction-batch key
    assert rows["Blank-ext-1"].samp_category == "negative control"
    assert rows["Blank-ext-1"].neg_cont_type == "extraction negative"
    assert rows["Blank-ext-1"].extraction_ID == "EXP524"
    assert rows["Blank-PCR-1"].neg_cont_type == "PCR negative"
    assert rows["Blank-PCR-1"].extraction_ID is None  # PCR blanks legitimately have none


@pytest.mark.skipif(not _OTU.exists(), reason="seednap teleo_rhone OTU table not present")
def test_cross_check_real_rhone_reproduces_known_orphan():
    """The OTU table carries Blank-PCR-3 (no metadata row) while the metadata lists
    Blank-PCR-1: the validator must flag exactly that orphan, never silently."""
    m = migrate_to_manifest(_FIELD, project_csv=_PROJ)
    res = validate_against_abundance(m, _OTU)
    assert res.orphan_abundance_columns == ["Blank-PCR-3"]
    assert res.manifest_extra_rows == ["Blank-PCR-1"]
    assert not res.ok
