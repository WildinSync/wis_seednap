"""Stable error-code catalog + removed-config-key registry.

Following rustc's model: each error class has a stable code (``SDN-XXX-NNN``) and a long-form
*extended explanation* that builds understanding of WHY the error occurs (not a copy-paste
quick fix). The inline message names what/why/fix for the specific case; ``seednap explain
<code>`` shows the deeper explanation. Codes are assigned only where the explanation adds more
than the inline message.
"""

from typing import Dict, Optional, Tuple

# ---------------------------------------------------------------------------
# Config keys removed in past migrations -> what to do instead. Used by the
# config humanizer so a stale config gets a migration hint, not just "unknown key".
# ---------------------------------------------------------------------------
REMOVED_KEYS: Dict[str, str] = {
    "version": "the config is no longer versioned; delete this key.",
    "metrics": "removed in the steps-model migration; ASV summary stats are now "
    "`dada2.collect_metrics`, and run reporting/plots are the `report` step.",
    "metrics.collect_asv_metrics": "moved to `dada2.collect_metrics`.",
    "metrics.generate_plots": "renamed and moved to `dada2.collect_metrics`.",
    "metrics.plot_format": "removed; the HTML report figures are produced by the `report` step.",
    "cleaning.enabled": "removed in the steps-model migration; cleaning runs only if `clean` is "
    "listed in `pipeline.steps`. Remove this key and add/remove `clean` from steps.",
    "export.gbif.enabled": "removed in the steps-model migration; GBIF export runs only if "
    "`export` is listed in `pipeline.steps`. Remove this key and add/remove `export` from steps.",
    "demultiplex.enabled": "removed in the steps-model migration; demultiplexing runs only if "
    "`demultiplex` is listed in `pipeline.steps` (before `trim`).",
    "demultiplex.skip": "removed; for pre-demultiplexed inputs just omit `demultiplex` from "
    "`pipeline.steps`.",
    "report.read_tracking": "removed; the read-tracking table is written whenever the `report` "
    "step runs (list/omit `report` in `pipeline.steps`).",
    "pipeline.skip": "removed; `pipeline.steps` is now the single source of truth -- just don't "
    "list a stage to skip it.",
    "trimming.tool": "removed; trimming always uses cutadapt.",
    "paths.references": "removed; reference databases are configured under "
    "`taxonomy.databases.<method>`, not `paths`.",
    "marker.primers.name": "removed; the primer block only needs `forward` and `reverse`.",
    "marker.primers.target": "removed; the primer block only needs `forward` and `reverse`.",
    "marker.primers.amplicon_length": "removed; the primer block only needs `forward` and "
    "`reverse`.",
}

# ---------------------------------------------------------------------------
# code -> (title, extended explanation). Kept terse; the explanation is the WHY.
# ---------------------------------------------------------------------------
CODES: Dict[str, Tuple[str, str]] = {
    "SDN-CFG-001": (
        "Unknown configuration key",
        "A key in your YAML is not part of the SeeDNAP config schema. SeeDNAP validates configs "
        "strictly (extra keys are rejected) so a typo or a key left over from an older version is "
        "caught at load time rather than silently ignored -- a silently-ignored setting would "
        "produce a dataset that looks valid but was processed with the wrong options. Compare the "
        "flagged key against the valid keys listed in the error (and `docs/configuration.md`); if "
        "it is a leftover from an older config, the error names what replaced it.",
    ),
    "SDN-CFG-002": (
        "Missing required configuration field",
        "A field SeeDNAP cannot run without is absent. Required fields (the marker identity and "
        "its primers, and the taxonomy method + its database block) have no safe default -- "
        "guessing them would mislabel or mis-trim your data. `seednap init` emits a template with "
        "every required field present.",
    ),
    "SDN-CFG-003": (
        "Wrong type for a configuration value",
        "A value is the wrong kind (e.g. text where a number or true/false is expected). YAML "
        "infers types from how a value is written; quote-wrapping a number, or a typo like "
        "`cores: lots`, yields a string. Use an unquoted number for numeric fields and "
        "true/false for booleans.",
    ),
    "SDN-CFG-004": (
        "Configuration value out of range",
        "A numeric value violates the allowed range for that field (e.g. a percentage above 100, "
        "or a count below 1). The bounds encode what is physically meaningful for the parameter; "
        "the error states the bound and what the field controls.",
    ),
    "SDN-CFG-005": (
        "Invalid choice for a configuration field",
        "A field that accepts only a fixed set of values was given something else. The error "
        "lists the allowed values; if the given value is a near-miss the closest valid one is "
        "named.",
    ),
    "SDN-CFG-006": (
        "Invalid pipeline.steps order or set",
        "`pipeline.steps` is the ordered list of stages to run; a stage runs iff listed. Stages "
        "have dependencies (trim before dada2/swarm; a feature step before taxonomy/clean; "
        "taxonomy before export) and dada2/swarm are mutually exclusive. The order is validated "
        "against this dependency graph so an impossible pipeline is rejected before it starts.",
    ),
    "SDN-CFG-007": (
        "Referenced file or directory not found",
        "A path in the config (raw-data directory, reference database) does not exist on disk. "
        "SeeDNAP checks referenced paths at `validate` / load time so a missing input is caught "
        "before trimming and clustering burn compute, rather than failing mid-run. A common cause "
        "is a config copied from another dataset whose paths point elsewhere.",
    ),
    "SDN-CFG-008": (
        "Taxonomy database block does not resolve",
        "`taxonomy.method` selects a classifier, and each method needs its own block under "
        "`taxonomy.databases` (blast -> a reference FASTA; dada2 -> RDP + species DBs; etc.). The "
        "selected method's block is missing or incomplete, so taxonomy could not run. This is "
        "checked at load time so it does not surface only after the feature step completes.",
    ),
    "SDN-CFG-009": (
        "Malformed configuration YAML",
        "The file is not valid YAML (often a tab used for indentation -- YAML requires spaces -- "
        "or an unclosed bracket/quote). The error gives the location; fix the syntax and re-run "
        "`seednap validate`.",
    ),
    "SDN-TOOL-001": (
        "External tool not installed or not on PATH",
        "SeeDNAP shells out to external bioinformatics tools (cutadapt, vsearch, swarm, "
        "blastn/makeblastdb, Rscript). One could not be launched, almost always because the wrong "
        "(or no) conda environment is active. Activate the environment that ships these tools and "
        "verify with `<tool> --version`.",
    ),
    "SDN-TOOL-002": (
        "External tool exited with an error",
        "An external tool ran but returned a non-zero status. The tool's own stderr (shown in the "
        "error) is the primary clue -- it is the tool's diagnostic, not a SeeDNAP bug. Common "
        "causes are a malformed/empty input, a wrong database path, or a tool-version mismatch.",
    ),
}


def explain(code: str) -> Optional[str]:
    """Return the extended explanation for an error code, or None if it is unknown.

    Backs the ``seednap explain <code>`` command: given a stable code shown in an error
    message, it returns the long-form title plus the WHY explanation from ``CODES``.

    Args:
        code: An error code such as ``SDN-CFG-001``. Surrounding whitespace is stripped
            and the code is upper-cased before lookup, so casing/padding does not matter.

    Returns:
        A formatted string ``"<CODE>: <title>\\n\\n<detail>"`` when the code is known,
        otherwise None.
    """
    entry = CODES.get(code.strip().upper())
    if entry is None:
        return None
    title, detail = entry
    return f"{code.strip().upper()}: {title}\n\n{detail}"


def all_codes() -> Dict[str, str]:
    """Return every known error code mapped to its short title, for listing.

    Used to render a directory of codes (e.g. ``seednap explain`` with no argument)
    without exposing the long-form explanations.

    Returns:
        A dict mapping each ``SDN-XXX-NNN`` code to its one-line title; the extended
        explanation is omitted.
    """
    return {c: t for c, (t, _) in CODES.items()}
