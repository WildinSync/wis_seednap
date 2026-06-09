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
    """True if the field annotation is Path or Optional[Path]."""
    if annotation is Path:
        return True
    if typing.get_origin(annotation) is typing.Union:
        return Path in typing.get_args(annotation)
    return False


def preflight_checks(config: PipelineConfig) -> List[SeednapError]:
    """Return structured problems for referenced inputs that do not exist (empty == all good)."""
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

    return problems
