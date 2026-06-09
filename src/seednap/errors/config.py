"""Humanize Pydantic ValidationError into actionable, what/why/fix config messages.

Pydantic v2 exposes a structured error contract (``ValidationError.errors()`` -> a list of
dicts with stable ``type``/``loc``/``msg``/``input``/``ctx`` keys). We map each error ``type``
to a friendly template, introspect the model to list the valid keys at the offending location,
add a declarative closest-match suggestion (difflib; phrased as a statement, not a question, per
rustc style), and consult a registry of keys removed in past migrations so a stale config gets a
migration hint rather than a bare "unknown key". Each line carries a stable code for
``seednap explain``.
"""

import difflib
import typing
from pathlib import Path
from typing import Any, List, Optional, Tuple, Type, Union

from pydantic import BaseModel, ValidationError

from seednap.config.models import PipelineConfig
from seednap.errors.catalog import REMOVED_KEYS

# Per-field guidance for the most important required fields (loc dotted -> (why, fix)).
_MISSING_HINTS = {
    "marker": (
        "every run needs a marker identity and its primer pair to trim and label outputs",
        "add a `marker:` block with `name:` and `primers: {forward:, reverse:}` "
        "(see `seednap init`)",
    ),
    "marker.primers.forward": (
        "two-pass cutadapt trimming needs both 5'->3' primer sequences",
        "add `forward: <5'->3' sequence>` under `marker.primers`",
    ),
    "marker.primers.reverse": (
        "two-pass cutadapt trimming needs both 5'->3' primer sequences",
        "add `reverse: <5'->3' sequence>` under `marker.primers`",
    ),
    "taxonomy.method": (
        "method selects which classifier runs and which database block is required",
        "set `taxonomy.method:` to one of blast, dada2, ecotag, decipher",
    ),
}

_KIND = {
    "int_parsing": "a whole number", "int_type": "a whole number",
    "bool_parsing": "a true/false value", "bool_type": "a true/false value",
    "float_parsing": "a number", "float_type": "a number",
    "string_type": "a text string",
}


def _unwrap(annotation: Any) -> Any:
    """Strip Optional/Union to the first BaseModel arg, if any."""
    if typing.get_origin(annotation) is Union:
        for arg in typing.get_args(annotation):
            if isinstance(arg, type) and issubclass(arg, BaseModel):
                return arg
    return annotation


def _model_at(loc: Tuple[Any, ...]) -> Optional[Type[BaseModel]]:
    """Resolve the Pydantic model reached by walking ``loc`` from PipelineConfig."""
    model: Any = PipelineConfig
    for part in loc:
        if not (isinstance(model, type) and issubclass(model, BaseModel)):
            return None
        field = model.model_fields.get(str(part))
        if field is None:
            return None
        model = _unwrap(field.annotation)
    return model if isinstance(model, type) and issubclass(model, BaseModel) else None


def _closest(word: str, options: List[str]) -> Optional[str]:
    """Return the single closest match to ``word`` among ``options``, or None."""
    matches = difflib.get_close_matches(word, options, n=1, cutoff=0.6)
    return matches[0] if matches else None


def _one(err: Any) -> Tuple[str, str]:
    """Map one Pydantic error dict (ErrorDetails) to (message, code)."""
    loc: Tuple[Any, ...] = tuple(err.get("loc", ()))
    dotted = ".".join(str(p) for p in loc) or "(top level)"
    etype = err.get("type", "")
    ctx = err.get("ctx") or {}
    given = err.get("input")
    msg = err.get("msg", "")

    if etype == "extra_forbidden":
        # Removed-in-migration key?
        if dotted in REMOVED_KEYS:
            return f"Unknown config key '{dotted}': {REMOVED_KEYS[dotted]}", "SDN-CFG-001"
        parent = _model_at(loc[:-1])
        valid = list(parent.model_fields.keys()) if parent else []
        bad = str(loc[-1]) if loc else dotted
        suggestion = ""
        close = _closest(bad, valid)
        if close:
            suggestion = f" The closest valid key is '{close}'."
        valid_txt = (
            f" Valid keys under '{'.'.join(str(p) for p in loc[:-1]) or 'the top level'}': "
            f"{', '.join(valid)}." if valid else ""
        )
        return (
            f"Unknown config key '{dotted}' (seednap rejects unknown keys to catch typos)."
            f"{suggestion}{valid_txt}",
            "SDN-CFG-001",
        )

    if etype == "missing":
        if dotted in _MISSING_HINTS:
            why, fix = _MISSING_HINTS[dotted]
            return f"Required field '{dotted}' is missing. Why: {why}. Fix: {fix}.", "SDN-CFG-002"
        return (
            f"Required field '{dotted}' is missing; add it to the config "
            f"(`seednap init` shows a template).",
            "SDN-CFG-002",
        )

    if etype in _KIND:
        return (
            f"'{dotted}': expected {_KIND[etype]}, got {given!r}. "
            f"Fix: provide a valid {_KIND[etype]} (unquoted numbers; true/false for booleans).",
            "SDN-CFG-003",
        )

    if etype in ("greater_than_equal", "greater_than", "less_than_equal", "less_than"):
        op = {"greater_than_equal": ">=", "greater_than": ">",
              "less_than_equal": "<=", "less_than": "<"}[etype]
        bound = ctx.get("ge", ctx.get("gt", ctx.get("le", ctx.get("lt"))))
        return (
            f"'{dotted}': value {given!r} is out of range (must be {op} {bound}).",
            "SDN-CFG-004",
        )

    if etype == "literal_error":
        expected = str(ctx.get("expected", "")).replace("'", "")
        options = [o.strip() for o in expected.replace(" or ", ", ").split(",") if o.strip()]
        close = _closest(str(given), options)
        suggestion = f" The closest allowed value is '{close}'." if close else ""
        return (
            f"'{dotted}': {given!r} is not a valid value. Allowed: {', '.join(options)}."
            f"{suggestion}",
            "SDN-CFG-005",
        )

    if etype == "value_error":
        # Custom validators (pipeline.steps DAG, taxonomy.databases) already produce
        # self-contained messages; surface them verbatim minus pydantic's prefix.
        text = msg[len("Value error, "):] if msg.startswith("Value error, ") else msg
        # Pick the code from where in the loc path the failing validator sits, so
        # `seednap explain <CODE>` lands on the right topic: a validator rooted at
        # `pipeline.*` (the steps DAG) -> 006; any path containing `databases`
        # (a taxonomy database block) -> 008; anything else -> 005 (generic invalid value).
        code = "SDN-CFG-006" if loc and str(loc[0]) == "pipeline" else (
            "SDN-CFG-008" if "databases" in (str(p) for p in loc) else "SDN-CFG-005")
        return text, code

    # Fallback: keep pydantic's message but name the location clearly.
    return f"'{dotted}': {msg}", "SDN-CFG-001"


def humanize_validation_error(exc: ValidationError, config_path: Path) -> str:
    """Build a full, actionable config-error message from a Pydantic ValidationError."""
    blocks: List[str] = []
    seen = set()
    for err in exc.errors():
        text, code = _one(err)
        key = (text, code)
        if key in seen:
            continue
        seen.add(key)
        blocks.append(f"  - {text}  [{code}]")
    body = "\n".join(blocks)
    return (
        f"Configuration is not valid: {config_path}\n{body}\n"
        f"  (run `seednap explain <CODE>` for any code above; see docs/configuration.md.)"
    )
