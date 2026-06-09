"""SeeDNAP error-explainability layer.

A small, dependency-free module that turns failures into actionable what/why/fix messages:

- ``SeednapError``: a structured user-facing error (summary + why + fix + optional code/docs).
- ``humanize_validation_error``: Pydantic ValidationError -> friendly config messages with
  closest-match suggestions and migration hints for removed keys.
- ``preflight_checks``: catch referenced-file / database-block problems at load time, not mid-run.
- ``catalog``: stable error codes + extended explanations (used by ``seednap explain``).
"""

from seednap.errors.base import SeednapError
from seednap.errors.catalog import REMOVED_KEYS, all_codes, explain
from seednap.errors.config import humanize_validation_error
from seednap.errors.preflight import preflight_checks

__all__ = [
    "SeednapError",
    "humanize_validation_error",
    "preflight_checks",
    "explain",
    "all_codes",
    "REMOVED_KEYS",
]
