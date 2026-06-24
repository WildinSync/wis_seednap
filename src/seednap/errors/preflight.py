"""Pre-run checks for referenced files the schema validates but does not exist-check.

Pydantic validates the *shape* of a config and expands paths, but never checks that referenced
files exist or that the selected taxonomy method's database block resolves. Those gaps mean a
bad path passes ``seednap validate`` and only fails mid-run, after trimming and clustering have
already burned compute. ``preflight_checks`` closes that: it returns the problems as structured
``SeednapError``s so ``validate`` (and the start of ``run-pipeline``) can fail early with a
clear what/why/fix.
"""

import typing
from pathlib import Path
from typing import Any, List

from pydantic import BaseModel

from seednap.config.models import PipelineConfig
from seednap.errors.base import SeednapError


def _is_path_field(annotation: Any) -> bool:
    """Return True if the field annotation is ``Path`` or ``Optional[Path]``.

    Used to pick out the filesystem-path fields of a taxonomy database block (reference
    FASTA, RDP/species databases, and the like) so only those get an existence check.

    Args:
        annotation: A type annotation taken from a Pydantic field.

    Returns:
        True if the annotation is ``Path`` or a Union that includes ``Path`` (e.g.
        ``Optional[Path]``); False otherwise.
    """
    if annotation is Path:
        return True
    if typing.get_origin(annotation) is typing.Union:
        return Path in typing.get_args(annotation)
    return False


def preflight_checks(config: PipelineConfig) -> List[SeednapError]:
    """Check that the inputs a config references actually exist on disk.

    Pydantic confirms the config's shape but never opens the filesystem, so a wrong
    path passes ``seednap validate`` and only fails partway through a run. This catches
    those before any trimming or clustering: it verifies the raw-data directory exists,
    confirms the selected ``taxonomy.method``'s database block resolves, and checks that
    every reference-database file that block points at is present (the reference FASTA
    or RDP/species databases the classifier reads to assign taxonomy). Each ``SeednapError``
    is collected and returned (not raised) so the caller can report all problems at once.

    Args:
        config: The fully validated ``PipelineConfig`` whose referenced paths and
            taxonomy database block are to be existence-checked.

    Returns:
        A list of ``SeednapError`` objects, one per problem found, each carrying a
        what/why/fix triad and an ``SDN-CFG-007``/``SDN-CFG-008`` code. An empty list
        means every referenced input was found. If the taxonomy database block fails to
        resolve, the list is returned immediately with only that single problem.
    """
    problems: List[SeednapError] = []

    raw = Path(config.paths.raw_data)
    if not raw.is_dir():
        problems.append(
            SeednapError(
                f"Raw-data directory does not exist: {raw}",
                why="this is paths.raw_data; the pipeline reads input FASTQ pairs from here",
                fix="point paths.raw_data at the directory holding your <sample>_R1/_R2 FASTQ "
                "files, and confirm it exists and is readable",
                code="SDN-CFG-007",
            )
        )

    method = config.taxonomy.method
    try:
        block = config.taxonomy.get_database_config()
    except Exception as exc:  # the method's DB block is missing/incomplete
        first = str(exc).splitlines()[0]
        problems.append(
            SeednapError(
                f"taxonomy.method is '{method}' but its database block does not resolve",
                why=first,
                fix=f"add or complete the `{method}:` block under taxonomy.databases with the "
                f"reference path(s) that method needs",
                code="SDN-CFG-008",
            )
        )
        return problems

    if isinstance(block, BaseModel):
        for name, field in block.model_fields.items():
            if not _is_path_field(field.annotation):
                continue
            val = getattr(block, name, None)
            if val is not None and not Path(val).exists():
                problems.append(
                    SeednapError(
                        f"taxonomy.databases.{method}.{name} does not exist: {val}",
                        why=f"this reference database file is required for the '{method}' "
                        "taxonomy method",
                        fix=f"set taxonomy.databases.{method}.{name} to an existing file path",
                        code="SDN-CFG-007",
                    )
                )

    # The DarwinCore occurrence step joins the taxonomy table to per-sample + per-project
    # metadata; require both up front rather than after trim/cluster/taxonomy have already run.
    if "darwincore" in config.pipeline.steps:
        for label, val in (
            ("report.sample_metadata", config.report.sample_metadata),
            ("report.project_metadata", config.report.project_metadata),
        ):
            if val is None:
                problems.append(
                    SeednapError(
                        f"pipeline.steps includes 'darwincore' but {label} is not set",
                        why="the DarwinCore occurrence file is built by joining the taxonomy "
                        "table to per-sample and per-project metadata; without it every "
                        "occurrence would have blank date / coordinates / provenance",
                        fix=f"set {label} to this dataset's metadata CSV, or remove 'darwincore' "
                        "from pipeline.steps",
                        code="SDN-CFG-007",
                    )
                )
            elif not Path(val).exists():
                problems.append(
                    SeednapError(
                        f"{label} does not exist: {val}",
                        why="the 'darwincore' step reads this metadata CSV to fill the "
                        "DarwinCore occurrence fields",
                        fix=f"set {label} to an existing metadata CSV",
                        code="SDN-CFG-007",
                    )
                )

    return problems
