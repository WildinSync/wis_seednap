#!/usr/bin/env python3
"""Migrate a pre-"steps model" seednap config to the unified pipeline.steps schema.

The redesign made ``pipeline.steps`` the single source of truth: a stage runs iff listed,
and the per-section enable gates were removed. This migrator rewrites an old config so it
loads under the new schema, **preserving its existing behavior and comments**:

  - new ``pipeline.steps`` = the old steps (minus ``pipeline.skip``), then:
      * drop ``demultiplex`` unless it was both listed AND ``demultiplex.enabled: true``
      * drop ``export``      unless it was both listed AND ``export.gbif.enabled: true``
      * insert ``clean`` after ``taxonomy`` if ``cleaning.enabled: true``
      * append ``report``    if ``report.read_tracking`` was true (the old default)
  - strip removed keys: ``demultiplex.enabled``/``skip``, ``export.gbif.enabled``,
    ``cleaning.enabled``, ``report.read_tracking``, ``pipeline.skip``
  - fold ``metrics.collect_asv_metrics`` -> ``dada2.collect_metrics`` (remove ``metrics:``)

Text-based: every other line and comment is preserved byte-for-byte. A ``.bak`` is written.
Usage: python migrate_config_to_steps_model.py FILE [FILE ...]
"""
import re
import sys
from pathlib import Path

import yaml


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _new_steps(cfg: dict) -> list:
    pipe = cfg.get("pipeline") or {}
    skip = set(pipe.get("skip", []) or [])
    steps = [s for s in (pipe.get("steps", []) or []) if s not in skip]
    demux_on = (cfg.get("demultiplex") or {}).get("enabled", False)
    gbif_on = ((cfg.get("export") or {}).get("gbif") or {}).get("enabled", True)
    clean_on = (cfg.get("cleaning") or {}).get("enabled", False)
    read_tracking = (cfg.get("report") or {}).get("read_tracking", True)
    if "demultiplex" in steps and not demux_on:
        steps.remove("demultiplex")
    if "export" in steps and not gbif_on:
        steps.remove("export")
    if clean_on and "clean" not in steps:
        steps.insert(steps.index("taxonomy") + 1 if "taxonomy" in steps else len(steps), "clean")
    if read_tracking and "report" not in steps:
        steps.append("report")
    return steps


# (top-level-section, indent, key) lines to strip outright
_STRIP = [
    ("demultiplex", 2, "enabled"), ("demultiplex", 2, "skip"),
    ("export", 4, "enabled"),  # export.gbif.enabled
    ("cleaning", 2, "enabled"),
    ("report", 2, "read_tracking"),
    ("pipeline", 2, "skip"),
]


def migrate(text: str) -> str:
    cfg = yaml.safe_load(text) or {}
    new_steps = _new_steps(cfg)
    metrics_val = (cfg.get("metrics") or {}).get("collect_asv_metrics", True)

    lines = text.splitlines(keepends=True)
    out: list = []
    top = None
    i = 0
    while i < len(lines):
        line = lines[i]
        ind = _indent(line)
        if ind == 0 and re.match(r"[A-Za-z_]+\s*:", line):
            top = line.split(":")[0].strip()
        body = line[ind:]

        # strip removed keys (indentation- and section-anchored)
        if any(top == sec and ind == d and re.match(rf"{k}\s*:", body) for sec, d, k in _STRIP):
            i += 1
            continue
        # remove the whole top-level metrics: block (folded into dada2.collect_metrics)
        if ind == 0 and re.match(r"metrics\s*:", line):
            i += 1
            while i < len(lines) and (lines[i].strip() == "" or _indent(lines[i]) > 0):
                if lines[i].strip() != "" and _indent(lines[i]) == 0:
                    break
                i += 1
            continue
        # fold collect_metrics into dada2: emit it right after the dada2 section's
        # last top-level (indent-2) scalar, before the next section.
        if top == "dada2" and ind == 0 and re.match(r"dada2\s*:", line):
            out.append(line)
            i += 1
            while i < len(lines) and not (_indent(lines[i]) == 0 and lines[i].strip()):
                out.append(lines[i])
                i += 1
            out.append(f"  collect_metrics: {str(metrics_val).lower()}"
                       f"   # ASV summary stats to metrics.json/csv + console (DADA2 path only)\n")
            continue
        # rewrite pipeline.steps list
        if top == "pipeline" and ind == 2 and re.match(r"steps\s*:", body):
            out.append(line)  # keep the "  steps:" line (+ any trailing comment)
            i += 1
            while i < len(lines) and lines[i].lstrip().startswith("- "):
                i += 1  # drop old items
            for s in new_steps:
                out.append(f'    - "{s}"\n')
            continue

        out.append(line)
        i += 1
    return "".join(out)


def main(argv) -> None:
    for fp in argv:
        p = Path(fp)
        new = migrate(p.read_text())
        if new != p.read_text():
            p.with_suffix(p.suffix + ".bak").write_text(p.read_text())
            p.write_text(new)
            print(f"{p}: migrated to the pipeline.steps model (.bak written)")
        else:
            print(f"{p}: already current")


if __name__ == "__main__":
    main(sys.argv[1:])
