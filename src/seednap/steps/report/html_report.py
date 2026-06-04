"""Self-contained HTML run report (opt-in, ``report.html_report``).

Renders a single portable ``.html`` file styled like a typeset scientific
paper: a warm-paper page, serif (Computer Modern) typography, justified text,
``Figure N`` / ``Table N`` captions, and restrained monochrome publication
figures with a single SeeDNAP-green accent. Every section (Summary, Dataset,
Read tracking, Per-sample detail, Taxonomy, Feature QC, Controls, Provenance,
Run log, Notes) is a self-contained selectable panel behind a sticky top
navigation bar, implemented with pure CSS radio-tabs -- one panel visible at a
time on screen, all panels expanded when printed. Nothing is repeated across
tabs.

Charts are matplotlib PNGs embedded as base64, so there are no external
assets, no CDN, and no JavaScript (the tab switching is pure CSS). It is dataset-agnostic: every number and
label is derived from the data passed in (read-tracking df, optional taxonomy
CSV, optional ``otu_table_full``, optional state JSON, optional run-log file).
The optional run-log section embeds the pipeline's console transcript,
colorized by level via rich's own HTML export so the palette matches the live
console exactly. Optional sources are ``[WARN]``-guarded -- a missing one
yields an explanatory sentence rather than vanishing silently (CLAUDE.md
section 4) -- and matplotlib is imported lazily so the report still renders
(text + tables) if it is absent.
"""

import base64
import html as _html
import io
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
from jinja2 import Template

from seednap.steps.report.read_tracking import DADA2_STEPS
from seednap.utils.logging import get_logger

logger = get_logger(__name__)

_DEFAULT_STEPS = DADA2_STEPS
_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
_TAX_META = {"ASV_ID", "OTU_ID", "pident", "is_contaminant_candidate", "Sequence", "sequence"}
_UNASSIGNED = {"Unassigned", "unassigned", "", "NA", "nan", "None"}

# Log levels recognised in the run-log transcript. The file logger writes
# "TIME | LEVEL | name:lineno | message"; EVENT levels are the ones a scientist
# must always see, so they survive truncation of long logs (see _select_log_lines).
_LOG_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
_EVENT_LEVELS = {"WARNING", "ERROR", "CRITICAL"}

# --- publication figure palette (mostly ink/grey, one sea-green accent) -------
INK = "#222222"
GREY = "#8a8a8a"
ACCENT = "#2e8b57"   # seednap sea-green -- the single emphasis color
HAIR = "#cccccc"

PAPER_RC = {
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["cmr10", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "mathtext.rm": "serif",
    "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.linewidth": 0.6, "axes.edgecolor": INK, "axes.labelcolor": INK,
    "text.color": INK, "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlelocation": "left",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.color": INK, "ytick.color": INK,
    "axes.grid": False, "grid.color": HAIR, "grid.linewidth": 0.5, "grid.alpha": 0.7,
    "legend.frameon": False,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "figure.dpi": 150, "savefig.dpi": 190, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
    "lines.linewidth": 1.1,
}

_TEMPLATE = Template(
    """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SeeDNAP run report: {{ marker }}</title>
<style>
  :root {
    --serif:"Latin Modern Roman","CMU Serif","STIX Two Text","Nimbus Roman","Times New Roman",Georgia,Cambria,"DejaVu Serif",Times,serif;
    --mono:"DejaVu Sans Mono",ui-monospace,"SFMono-Regular",Consolas,"Liberation Mono",Menlo,monospace;
    --measure:70ch; --leading:1.55; --ink:#1a1a1a; --muted:#555; --paper:#fdfdfb;
    --rule:#111; --hair:#d8d8d8; --accent:#2e8b57; --link:#1a4f8b;
  }
  *{box-sizing:border-box;}
  html{font-size:17px;-webkit-text-size-adjust:100%;}
  body{font-family:var(--serif); font-size:1rem; line-height:var(--leading); color:var(--ink);
       background:var(--paper); margin:0; padding:0 0 5rem;
       text-rendering:optimizeLegibility; font-kerning:normal;
       font-feature-settings:"kern" 1,"liga" 1;}
  /* Centered reading column for prose/tables; the top bar and the terminal
     deliberately break out of it. */
  .panel{max-width:var(--measure); margin-left:auto; margin-right:auto;
       padding-left:1.25rem; padding-right:1.25rem;}
  .descriptor{color:var(--muted); font-size:.95rem; margin:0 0 1rem;}
  .summary-lead{font-size:1.02rem;}
  p{margin:0 0 .85em; text-align:justify; hyphens:auto; -webkit-hyphens:auto; hyphenate-limit-chars:6 3 3;}
  h2,h3{font-family:var(--serif); font-weight:700; line-height:1.2; text-align:left; hyphens:manual;}
  h2{font-size:1.28rem; margin:0 0 .7em; border-bottom:1px solid var(--hair); padding-bottom:.2em;}
  figure{margin:1.6rem 0; text-align:center;}
  figure img{max-width:100%; height:auto; display:block; margin:0 auto;}
  figcaption{font-size:.87rem; text-align:left; margin-top:.5rem; line-height:1.4;}
  figcaption b{font-weight:700;}
  table{border-collapse:collapse; margin:1.2rem auto; font-size:.9rem; width:100%;
        border-top:1.4px solid var(--rule); font-variant-numeric:tabular-nums lining-nums;}
  caption{caption-side:top; text-align:left; font-size:.87rem; margin-bottom:.4rem;}
  caption b{font-weight:700;}
  thead th{border-bottom:1px solid var(--rule); font-weight:700; text-align:right; padding:.3rem .7rem;}
  thead th:first-child, tbody td:first-child{text-align:left;}
  tbody td{padding:.26rem .7rem; text-align:right; border-bottom:1px solid var(--hair);}
  tbody tr:last-child td{border-bottom:1.4px solid var(--rule);}
  .scroll{max-height:30rem; overflow:auto;}
  .flag-low{font-weight:700;} .flag-low::after{content:" *";}
  .na{color:var(--muted); font-style:italic;}
  .warn-head{font-variant-caps:small-caps; letter-spacing:.06em; font-weight:700; font-size:.8rem;
             color:var(--accent); margin:.9rem 0 .25rem;}
  .warn-none{font-style:italic; color:var(--muted);}
  code{font-family:var(--mono); font-size:.85em;}
  /* Tabbed panels (pure CSS, no JS): hidden radios drive which panel shows.
     The labels live in a sticky top bar; one panel is visible at a time. */
  input.tab-radio{position:absolute; width:1px; height:1px; opacity:0; pointer-events:none;}
  .topbar{position:sticky; top:0; z-index:10; background:var(--paper);
          border-bottom:1px solid var(--hair); box-shadow:0 1px 6px rgba(0,0,0,.04);}
  .topbar-inner{display:flex; align-items:center; gap:.9rem; flex-wrap:wrap;
          max-width:1180px; margin:0 auto; padding:.5rem 1.25rem;}
  .brand{font-variant-caps:small-caps; letter-spacing:.12em; font-weight:700; color:var(--accent);
          white-space:nowrap; font-size:.95rem;}
  .brand .brand-sub{color:var(--muted); font-weight:400; font-variant-caps:normal; letter-spacing:0;
          font-size:.8rem; margin-left:.4rem;}
  .tabs{display:flex; flex-wrap:wrap; gap:.35rem;}
  .tabs label{font-family:var(--serif); font-size:.85rem; cursor:pointer; user-select:none; white-space:nowrap;
        padding:.3rem .7rem; border:1px solid var(--hair); border-radius:5px; color:var(--muted); background:#fff;}
  .tabs label:hover{color:var(--ink); border-color:var(--accent);}
  .panel{display:none;} .panel h2{margin-top:.2rem;} .panel{padding-top:1.4rem; padding-bottom:1rem;}
  {% for i in range(12) %}#tab-{{ i }}:checked ~ #panel-{{ i }}{% if not loop.last %},
  {% endif %}{% endfor %}{display:block;}
  {% for i in range(12) %}#tab-{{ i }}:checked ~ .topbar label[for="tab-{{ i }}"]{% if not loop.last %},
  {% endif %}{% endfor %}{color:#fff; background:var(--accent); border-color:var(--accent); font-weight:700;}
  {% for i in range(12) %}#tab-{{ i }}:focus-visible ~ .topbar label[for="tab-{{ i }}"]{% if not loop.last %},
  {% endif %}{% endfor %}{outline:2px solid var(--accent); outline-offset:2px;}
  /* Run-log transcript rendered as a real terminal window: dark chrome with
     traffic-light dots, a large dark body, and the bright ANSI palette that
     reads on black (info blue, warning yellow, error red). The window breaks
     out of the reading column so it is genuinely terminal-sized. */
  .runlog-meta{font-family:var(--mono); font-size:.78rem; color:var(--muted); margin:.2rem 0 .6rem;}
  .term-legend{font-family:var(--mono); font-size:.78rem; color:var(--muted);}
  .lvl-chip{font-family:var(--mono); padding:.04rem .42rem; border-radius:3px; background:#1b1d23; margin:0 .12rem;}
  .lvl-info{color:#6cb6ff;} .lvl-warning{color:#e3b341; font-weight:600;} .lvl-error{color:#ff6b6b; font-weight:700;}
  .terminal{width:min(95vw,1200px); margin:1.1rem 0; margin-left:50%; transform:translateX(-50%);
            border-radius:9px; overflow:hidden; border:1px solid #000;
            box-shadow:0 10px 34px rgba(0,0,0,.30); background:#15171c;}
  .term-bar{display:flex; align-items:center; gap:.5rem; padding:.55rem .8rem;
            background:linear-gradient(#3a3d44,#2b2e34); border-bottom:1px solid #000;}
  .term-dot{width:12px; height:12px; border-radius:50%; box-shadow:inset 0 0 0 .5px rgba(0,0,0,.25);}
  .term-dot.r{background:#ff5f56;} .term-dot.y{background:#ffbd2e;} .term-dot.g{background:#27c93f;}
  .term-title{margin-left:.5rem; color:#c8ccd2; font-family:var(--mono); font-size:.8rem; letter-spacing:.01em;}
  .term-body{height:72vh; min-height:24rem; overflow:auto; background:#15171c;}
  .term-body.compact{height:auto; max-height:26rem; min-height:0;}
  .runlog{font-family:var(--mono); font-size:.82rem; line-height:1.55; color:#d6d9df;
          background:transparent; margin:0; padding:1rem 1.15rem; white-space:pre; tab-size:2;}
  /* Pure-CSS Fullscreen toggle for the run-log terminal (no JS). */
  .term-max-toggle{position:absolute; width:1px; height:1px; opacity:0; pointer-events:none;}
  .term-max-btn{margin-left:auto; cursor:pointer; user-select:none; color:#c8ccd2;
            font-family:var(--mono); font-size:.74rem; border:1px solid #555; border-radius:4px; padding:.08rem .55rem;}
  .term-max-btn:hover{border-color:#9aa0a6; color:#fff;}
  .lbl-close{display:none;}
  #termmax:checked ~ .terminal{position:fixed; inset:0; width:100vw; height:100vh; max-width:none;
            margin:0; transform:none; border-radius:0; z-index:9999;}
  #termmax:checked ~ .terminal .term-body{height:calc(100vh - 2.7rem); max-height:none;}
  #termmax:checked ~ .terminal .runlog{font-size:.95rem; line-height:1.6;}
  #termmax:checked ~ .terminal .lbl-open{display:none;} #termmax:checked ~ .terminal .lbl-close{display:inline;}
  @media print{ body{font-size:10.5pt; line-height:1.4; color:#000; background:#fff; padding:0;}
    /* Print the whole document: expand every panel, drop the interactive nav. */
    input.tab-radio, .topbar, .term-max-btn{display:none !important;}
    .panel{display:block !important; margin:1.4rem auto; break-inside:avoid;}
    .terminal{width:100%; margin:1rem 0; transform:none; box-shadow:none;
              -webkit-print-color-adjust:exact; print-color-adjust:exact;}
    .term-body, .term-body.compact{height:auto; max-height:none;}
    h2,h3{break-after:avoid;} figure,table{break-inside:avoid;} @page{margin:2cm;} }
</style></head>
<body>

{% for s in sections %}<input class="tab-radio" type="radio" name="rtab" id="tab-{{ loop.index0 }}"{% if loop.first %} checked{% endif %}>
{% endfor %}
<header class="topbar">
  <div class="topbar-inner">
    <span class="brand">SeeDNAP<span class="brand-sub">{{ marker }} run report</span></span>
    <nav class="tabs" role="tablist" aria-label="Report sections">
{% for s in sections %}      <label for="tab-{{ loop.index0 }}">{{ s.title }}</label>
{% endfor %}    </nav>
  </div>
</header>

{% for s in sections %}<section class="panel" id="panel-{{ loop.index0 }}">
<h2>{{ s.title }}</h2>
{{ s.html }}
</section>
{% endfor %}
</body></html>
"""
)


def _esc(v) -> str:
    return _html.escape(str(v))


class HTMLReportBuilder:
    """Render a self-contained, paper-styled HTML run report."""

    def __init__(
        self,
        marker: str,
        tracking_df: pd.DataFrame,
        warnings: Optional[List[str]] = None,
        summary: Optional[Dict[str, object]] = None,
        steps: Optional[List[str]] = None,
        state: Optional[Dict[str, object]] = None,
        taxonomy_csv: Optional[Union[str, Path]] = None,
        otu_table_full: Optional[Union[str, Path]] = None,
        field_metadata_csv: Optional[Union[str, Path]] = None,
        project_metadata_csv: Optional[Union[str, Path]] = None,
        log_file: Optional[Union[str, Path]] = None,
        max_log_lines: int = 1500,
    ) -> None:
        self.marker = marker
        self.df = tracking_df if tracking_df is not None else pd.DataFrame()
        self.warnings = warnings or []
        self.summary = summary or {}
        self.state = state or {}
        self.taxonomy_csv = Path(taxonomy_csv) if taxonomy_csv else None
        self.otu_table_full = Path(otu_table_full) if otu_table_full else None
        self.field_metadata_csv = Path(field_metadata_csv) if field_metadata_csv else None
        self.project_metadata_csv = Path(project_metadata_csv) if project_metadata_csv else None
        self.log_file = Path(log_file) if log_file else None
        self.max_log_lines = max_log_lines
        if steps:
            self.steps = steps
        else:
            self.steps = [c for c in _DEFAULT_STEPS if c in self.df.columns] or _DEFAULT_STEPS
        self.is_dada2 = "nonchim" in self.steps
        self.final_step = self.steps[-1] if self.steps else "raw"
        self.warn_pct = float(self.summary.get("warn_below_retention_pct", 30.0))
        self._tax_cache: Optional[pd.DataFrame] = None
        self._fig_n = 0
        self._tbl_n = 0

    # ------------------------------------------------------------------ #
    # Optional data sources (lazy, [WARN]-guarded)
    # ------------------------------------------------------------------ #
    def _tax(self) -> Optional[pd.DataFrame]:
        if self._tax_cache is not None:
            return self._tax_cache
        if self.taxonomy_csv is None:
            return None
        if not self.taxonomy_csv.exists():
            logger.warning(f"[WARN] html_report: expected=taxonomy_csv, got=missing "
                           f"({self.taxonomy_csv}), fallback=omit taxonomy section")
            return None
        try:
            self._tax_cache = pd.read_csv(self.taxonomy_csv)
        except (pd.errors.EmptyDataError, OSError) as exc:
            logger.warning(f"[WARN] html_report: expected=readable taxonomy_csv, got=unreadable "
                           f"({self.taxonomy_csv}: {exc}), fallback=omit taxonomy section")
            return None
        return self._tax_cache

    def _otu_full(self) -> Optional[pd.DataFrame]:
        if self.otu_table_full is None:
            return None
        if not self.otu_table_full.exists():
            logger.warning(f"[WARN] html_report: expected=otu_table_full, got=missing "
                           f"({self.otu_table_full}), fallback=omit feature-QC section")
            return None
        try:
            return pd.read_csv(self.otu_table_full)
        except (pd.errors.EmptyDataError, OSError) as exc:
            logger.warning(f"[WARN] html_report: expected=readable otu_table_full, got=unreadable "
                           f"({self.otu_table_full}: {exc}), fallback=omit feature-QC section")
            return None

    def _tax_sample_cols(self, tax: pd.DataFrame) -> List[str]:
        return [c for c in tax.columns if c not in _TAX_META and c not in _RANKS]

    def _richness(self) -> Optional[pd.Series]:
        """Per-sample feature richness (count of features with >0 reads), indexed
        by sample name. Drawn from the taxonomy table, else the full OTU table;
        ``None`` if neither carries per-sample columns."""
        tax = self._tax()
        if tax is not None:
            cols = self._tax_sample_cols(tax)
            src = tax
        else:
            otu = self._otu_full()
            if otu is None:
                return None
            meta = {"OTU", "OTU_ID", "ASV_ID", "total", "length", "chimera", "spread",
                    "sequence", "Sequence"}
            cols = [c for c in otu.columns if c not in meta]
            src = otu
        if not cols:
            return None
        counts = (src[cols].apply(pd.to_numeric, errors="coerce").fillna(0) > 0).sum(axis=0)
        return counts.astype(int)

    def _read_meta(self, path: Optional[Path], kind: str) -> Optional[pd.DataFrame]:
        if path is None:
            return None
        if not path.exists():
            logger.warning(f"[WARN] html_report: expected={kind} metadata CSV, got=missing "
                           f"({path}), fallback=omit it from the dataset section")
            return None
        try:
            # utf-8-sig tolerates a BOM (some lab metadata files carry one).
            return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
        except (pd.errors.EmptyDataError, OSError) as exc:
            logger.warning(f"[WARN] html_report: expected=readable {kind} metadata, got=unreadable "
                           f"({path}: {exc}), fallback=omit it from the dataset section")
            return None

    # ------------------------------------------------------------------ #
    # Numbers / metadata
    # ------------------------------------------------------------------ #
    def _n_controls(self) -> int:
        if "sample" not in self.df.columns:
            return 0
        return int(self.df["sample"].astype(str).str.lower().str.startswith("blank").sum())

    def _run_date(self) -> str:
        for key in ("completed_at", "started_at"):
            v = self.state.get(key) if isinstance(self.state, dict) else None
            if v:
                # Date only, whether the timestamp is "T"- or space-separated.
                return str(v).replace("T", " ").split(" ")[0]
        # No state timestamp -- fall back to build date, said plainly.
        return f"{datetime.now().date().isoformat()} (report build date)"

    def _descriptor(self) -> str:
        n = len(self.df)
        nc = self._n_controls()
        method = "DADA2 ASV path" if self.is_dada2 else "SWARM OTU path"
        bits = [f"{n} sample{'s' if n != 1 else ''}" + (f" ({nc} control{'s' if nc != 1 else ''})" if nc else ""),
                method]
        tax = self._tax()
        if tax is not None:
            bits.append(f"{len(tax):,} {'ASVs' if self.is_dada2 else 'OTUs'}")
        return ", ".join(bits)

    def _abstract(self) -> str:
        n = len(self.df)
        if n == 0:
            return "No samples were found for this run; the read-tracking table is empty."
        nc = self._n_controls()
        method = "DADA2 ASV" if self.is_dada2 else "SWARM OTU"
        feature = "ASV" if self.is_dada2 else "OTU"
        clauses = [
            f"Marker <i>{_esc(self.marker)}</i>, {n} sample{'s' if n != 1 else ''}"
            + (f" ({n - nc} biological, {nc} control{'s' if nc != 1 else ''})" if nc else "")
            + f", {method} path."
        ]
        raw = pd.to_numeric(self.df.get("raw", pd.Series(dtype=float)), errors="coerce")
        fin = pd.to_numeric(self.df.get(self.final_step, pd.Series(dtype=float)), errors="coerce")
        if raw.notna().any() and fin.notna().any() and raw.sum() > 0:
            pct = fin.sum() / raw.sum() * 100
            clauses.append(f"{int(raw.sum()):,} raw read pairs, {pct:.1f}% retained to the "
                           f"{self.final_step} stage.")
        tax = self._tax()
        if tax is not None and "species" in tax.columns and len(tax):
            sp = int((~tax["species"].astype(str).isin(_UNASSIGNED)).sum())
            extra = ""
            if "pident" in tax.columns:
                pid = pd.to_numeric(tax["pident"], errors="coerce").dropna()
                if not pid.empty:
                    extra = f", median best-hit identity {pid.median():.1f}%"
            clauses.append(f"{sp:,} of {len(tax):,} {feature}s ({sp / len(tax) * 100:.1f}%) "
                           f"assigned to species{extra}.")
        nw = len(self.warnings)
        clauses.append("All steps completed. "
                       + ("No retention warnings were raised." if nw == 0
                          else f"{nw} read-tracking warning{'s' if nw != 1 else ''} "
                               f"(see the Read tracking tab)."))
        return " ".join(clauses)

    # ------------------------------------------------------------------ #
    # HTML fragment helpers
    # ------------------------------------------------------------------ #
    def _fig(self, b64: Optional[str], caption: str) -> str:
        if not b64:
            return ""
        self._fig_n += 1
        return (f'<figure><img alt="figure {self._fig_n}" src="data:image/png;base64,{b64}">'
                f'<figcaption><b>Figure {self._fig_n}.</b> {caption}</figcaption></figure>')

    def _table(self, caption: str, headers: List[str], rows: List[List[str]],
               scroll: bool = False) -> str:
        self._tbl_n += 1
        head = "".join(f"<th>{_esc(h)}</th>" for h in headers)
        body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
        tbl = (f'<table><caption><b>Table {self._tbl_n}.</b> {caption}</caption>'
               f'<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>')
        return f'<div class="scroll">{tbl}</div>' if scroll else tbl

    @staticmethod
    def _fmt(v) -> str:
        return "n/a" if v is None or (isinstance(v, float) and pd.isna(v)) else f"{int(v):,}"

    # ------------------------------------------------------------------ #
    # Figures (within scoped rc_context; called once)
    # ------------------------------------------------------------------ #
    def _make_figures(self) -> Dict[str, str]:
        try:
            import matplotlib as mpl
            mpl.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError as exc:
            logger.warning(f"[WARN] html_report: expected=matplotlib for figures, got=missing "
                           f"({exc}), fallback=text+tables only")
            return {}
        figs: Dict[str, str] = {}

        def emit(fig) -> str:
            buf = io.BytesIO()
            fig.savefig(buf, format="png")
            plt.close(fig)
            return base64.b64encode(buf.getvalue()).decode()

        with mpl.rc_context(PAPER_RC):
            # F: read funnel (totals per step), final bar accented + % labels.
            if not self.df.empty:
                totals = [pd.to_numeric(self.df[s], errors="coerce").sum(skipna=True) for s in self.steps]
                if any(t > 0 for t in totals):
                    fig, ax = plt.subplots(figsize=(5.5, 2.7))
                    colors = [GREY] * len(self.steps); colors[-1] = ACCENT
                    bars = ax.bar(self.steps, totals, color=colors, width=.62)
                    base = totals[0] if totals[0] else 1
                    ax.bar_label(bars, labels=[f"{int(t):,}\n{t / base * 100:.0f}%" for t in totals],
                                 fontsize=7.5, color=INK, padding=2)
                    ax.margins(y=.22); ax.set_ylabel("read pairs"); ax.set_title("Reads surviving each step")
                    figs["funnel"] = emit(fig)

                pr = pd.to_numeric(self.df["pct_retained"], errors="coerce").dropna()
                if not pr.empty:
                    if len(pr) <= 50:
                        # Per-sample horizontal bars (readable for modest sample counts).
                        order = pr.sort_values()
                        names = self.df.loc[order.index, "sample"].astype(str).tolist()
                        colors = [ACCENT if v < self.warn_pct else GREY for v in order]
                        fig, ax = plt.subplots(figsize=(5.5, max(2.4, .23 * len(order))))
                        ax.barh(range(len(order)), order.values, color=colors, height=.7)
                        ax.set_yticks(range(len(order))); ax.set_yticklabels(names, fontsize=7)
                        ax.set_title("Per-sample retention")
                    else:
                        # Many samples: a per-sample bar chart would be metres tall, so
                        # show the retention *distribution* instead (the table lists every sample).
                        fig, ax = plt.subplots(figsize=(5.5, 2.7))
                        ax.hist(pr, bins=30, color=GREY, edgecolor="white")
                        ax.set_ylabel("samples"); ax.set_title(f"Retention distribution (n={len(pr)})")
                    ax.axvline(self.warn_pct, color=INK, lw=.8, ls="--")
                    ax.set_xlim(0, 100); ax.set_xlabel("reads retained (%)")
                    figs["retention"] = emit(fig)

            # Per-sample feature richness (OTUs/ASVs detected per sample).
            rich = self._richness()
            if rich is not None and not rich.empty and rich.sum() > 0:
                feat = "ASVs" if self.is_dada2 else "OTUs"
                if len(rich) <= 50:
                    order = rich.sort_values()
                    fig, ax = plt.subplots(figsize=(5.5, max(2.4, .23 * len(order))))
                    ax.barh(range(len(order)), order.values, color=GREY, height=.7)
                    ax.set_yticks(range(len(order)))
                    ax.set_yticklabels([str(s) for s in order.index], fontsize=7)
                    ax.set_xlabel(f"{feat} detected"); ax.set_title(f"Per-sample richness ({feat} detected)")
                    for i, v in enumerate(order.values):
                        ax.text(v, i, f" {int(v):,}", va="center", fontsize=6.5, color=INK)
                else:
                    fig, ax = plt.subplots(figsize=(5.5, 2.7))
                    ax.hist(rich.values, bins=30, color=GREY, edgecolor="white")
                    ax.set_xlabel(f"{feat} detected"); ax.set_ylabel("samples")
                    ax.set_title(f"Richness distribution (n={len(rich)})")
                figs["richness"] = emit(fig)

            tax = self._tax()
            if tax is not None:
                ranks = [r for r in _RANKS if r in tax.columns]
                counts = [int((~tax[r].astype(str).isin(_UNASSIGNED)).sum()) for r in ranks]
                if any(counts):
                    fig, ax = plt.subplots(figsize=(5.5, 2.6))
                    ax.barh(range(len(ranks)), counts, color=GREY, height=.7)
                    ax.set_yticks(range(len(ranks))); ax.set_yticklabels(ranks)
                    ax.invert_yaxis()
                    for i, v in enumerate(counts):
                        ax.text(v, i, f" {v:,}", va="center", fontsize=7.5, color=INK)
                    ax.set_xlabel("features assigned"); ax.set_title("Assignment by rank")
                    ax.grid(axis="x", color=HAIR, lw=.5, ls="--", alpha=.7); ax.set_axisbelow(True)
                    figs["rank"] = emit(fig)
                if "pident" in tax.columns:
                    pid = pd.to_numeric(tax["pident"], errors="coerce").dropna()
                    if not pid.empty:
                        fig, ax = plt.subplots(figsize=(5.5, 2.5))
                        ax.hist(pid, bins=20, color=GREY, edgecolor="white")
                        ax.axvline(pid.median(), color=ACCENT, lw=1.1, ls="--")
                        ax.set_xlabel("BLAST identity (%)"); ax.set_ylabel("features")
                        ax.set_title("Assignment identity")
                        figs["pident"] = emit(fig)

            otu = self._otu_full()
            if otu is not None:
                if "chimera" in otu.columns:
                    vc = otu["chimera"].astype(str).value_counts()
                    seg = [("clean", int(vc.get("N", 0)), GREY), ("chimeric", int(vc.get("Y", 0)), INK),
                           ("borderline", int(vc.get("?", 0)), ACCENT)]
                    seg = [s for s in seg if s[1] > 0]
                    if seg:
                        fig, ax = plt.subplots(figsize=(5.5, 1.5))
                        left = 0
                        for label, val, col in seg:
                            ax.barh([0], [val], left=left, color=col, label=f"{label} ({val:,})", height=.5)
                            left += val
                        ax.set_yticks([]); ax.set_xlabel("OTUs"); ax.set_title("Chimera classification")
                        ax.legend(loc="upper center", bbox_to_anchor=(.5, -.55), ncol=3, fontsize=8)
                        figs["chimera"] = emit(fig)
                if "length" in otu.columns:
                    ln = pd.to_numeric(otu["length"], errors="coerce").dropna()
                    if not ln.empty:
                        fig, ax = plt.subplots(figsize=(5.5, 2.4))
                        ax.hist(ln, bins=30, color=GREY, edgecolor="white")
                        ax.axvline(ln.median(), color=ACCENT, lw=1.1, ls="--")
                        ax.set_xlabel("OTU length (bp)"); ax.set_ylabel("OTUs"); ax.set_title("Sequence-length distribution")
                        figs["length"] = emit(fig)
        return figs

    # ------------------------------------------------------------------ #
    # Section builders
    # ------------------------------------------------------------------ #
    def _section_dataset(self) -> str:
        """Dataset identity + provenance, from config facts + field/project metadata."""
        prov = self.summary.get("provenance") or {}
        field = self._read_meta(self.field_metadata_csv, "field")
        proj = self._read_meta(self.project_metadata_csv, "project")
        rows: List[Tuple[str, str]] = []

        def add(label, value):
            if value not in (None, "", "NA", "nan"):
                rows.append((label, str(value)))

        add("Dataset", prov.get("dataset_name") or self.marker)
        # Prefer the project metadata's marker (authoritative, e.g. "teleo") over the
        # run/dataset name passed on the CLI (e.g. "teleo_rhone").
        add("Marker", self._proj_val(proj, "marker") or prov.get("marker"))
        if prov.get("primer_fwd") or prov.get("primer_rev"):
            add("Primers", f"{prov.get('primer_fwd', '?')} / {prov.get('primer_rev', '?')} (fwd / rev)")

        # --- field metadata: location, dates, sites, institution ---
        if field is not None:
            # biological samples only (drop Blank* controls) for geography
            idc = field.columns[0]
            bio = field[~field[idc].astype(str).str.lower().str.startswith("blank")]
            lat = pd.to_numeric(bio.get("decimalLatitude", pd.Series(dtype=str)), errors="coerce").dropna()
            lon = pd.to_numeric(bio.get("decimalLongitude", pd.Series(dtype=str)), errors="coerce").dropna()
            if not lat.empty and not lon.empty:
                add("Location (lat, lon)",
                    f"{lat.mean():.4f}, {lon.mean():.4f} centroid "
                    f"(lat {lat.min():.3f} to {lat.max():.3f}, lon {lon.min():.3f} to {lon.max():.3f})")
                n_sites = len({(round(a, 4), round(o, 4)) for a, o in zip(lat, lon)})
                add("Distinct sites", f"{n_sites}")
            for col, label in (("site_names", "Site names"), ("area_basin", "Area / basin"),
                               ("ecosystem", "Ecosystem"), ("body", "Water body"),
                               ("env_medium", "Environment")):
                if col in field.columns:
                    vals = sorted({v for v in bio[col].astype(str) if v not in ("", "NA", "nan")})
                    if vals:
                        shown = ", ".join(vals[:6]) + (" …" if len(vals) > 6 else "")
                        add(label, shown)
            dates = sorted({v for v in bio.get("eventDate", pd.Series(dtype=str)).astype(str)
                            if v not in ("", "NA", "nan")})
            if dates:
                add("Sampling dates", f"{dates[0]} to {dates[-1]}" if len(dates) > 1 else dates[0])
            depth = pd.to_numeric(bio.get("depth", pd.Series(dtype=str)), errors="coerce").dropna()
            if not depth.empty:
                add("Depth (m)", f"{depth.min():g} to {depth.max():g}" if depth.min() != depth.max() else f"{depth.min():g}")
            for col, label in (("institution", "Institution"), ("laboratory", "Laboratory")):
                if col in field.columns:
                    vals = sorted({v for v in field[col].astype(str) if v not in ("", "NA", "nan")})
                    if vals:
                        add(label, ", ".join(vals[:3]))

        # --- project metadata: recorder, sequencing, reference DB ---
        add("Recorded by", self._proj_val(proj, "recordedby"))
        add("Sequencing", self._proj_val(proj, "seqmet"))
        add("Reference DB", self._proj_val(proj, "otu_db") or prov.get("reference_db"))
        add("Assignment", self._proj_val(proj, "identificationRemarks"))
        add("Raw data", prov.get("raw_data"))

        if not rows:
            return ("<p>No dataset metadata was provided. Pass <code>--field-metadata</code> and "
                    "<code>--project-metadata</code> (the lab <code>metadata_field_*.csv</code> / "
                    "<code>metadata_proj_*.csv</code>) to embed sampling location, dates, institution "
                    "and sequencing provenance here.</p>")
        body = [[_esc(k), _esc(v)] for k, v in rows]
        note = ""
        if field is None and proj is None:
            note = ("<p>Sampling-location and sequencing metadata were not provided "
                    "(<code>--field-metadata</code> / <code>--project-metadata</code>); the table below "
                    "is limited to pipeline-configuration facts.</p>")
        return note + self._table("Dataset identity and provenance.", ["field", "value"], body)

    @staticmethod
    def _proj_val(proj: Optional[pd.DataFrame], col: str) -> Optional[str]:
        if proj is None or col not in proj.columns or proj.empty:
            return None
        v = str(proj[col].iloc[0]).strip()
        return v or None

    def _section_read_tracking(self, figs) -> str:
        n = len(self.df)
        parts = []
        parts.append(f"<p>Reads were tracked per sample across the {len(self.steps)} stages of the "
                     f"{'DADA2 ASV' if self.is_dada2 else 'SWARM OTU'} path "
                     f"({' &rarr; '.join(self.steps)}).</p>")
        parts.append(self._fig(figs.get("funnel"),
                     f"Total read pairs retained after each step, summed across all {n} samples; "
                     f"labels give absolute counts and the percentage of raw input."))
        pr = pd.to_numeric(self.df.get("pct_retained", pd.Series(dtype=float)), errors="coerce").dropna()
        n_below = int((pr < self.warn_pct).sum()) if not pr.empty else 0
        if len(pr) > 50:
            ret_cap = (f"Distribution of per-sample read retention (raw &rarr; {self.final_step}) across "
                       f"{len(pr)} samples; the dashed line marks the {self.warn_pct:.0f}% threshold "
                       f"({n_below} below).")
        else:
            ret_cap = (f"Per-sample read retention (raw &rarr; {self.final_step}); the dashed line marks "
                       f"the {self.warn_pct:.0f}% threshold ({n_below} sample(s) below).")
        parts.append(self._fig(figs.get("retention"), ret_cap))
        # per-sample table
        rows = []
        for _, r in self.df.iterrows():
            pr = r["pct_retained"]
            low = pd.notna(pr) and pr < self.warn_pct
            name = f'<span class="flag-low">{_esc(r["sample"])}</span>' if low else _esc(r["sample"])
            cells = [name]
            for s in self.steps:
                v = r[s]
                cells.append('<span class="na">NA</span>' if pd.isna(v) else f"{int(v):,}")
            cells.append("NA" if pd.isna(pr) else f"{pr:.1f}%")
            rows.append(cells)
        parts.append(self._table("Per-sample read counts at each pipeline step. Rows flagged with "
                                 "an asterisk fall below the retention threshold.",
                                 ["sample", *self.steps, "% retained"], rows, scroll=True))
        if self.warnings:
            step_loss = float(self.summary.get("warn_step_loss_pct", 70.0))
            nw = len(self.warnings)
            parts.append(
                f"<p>The run raised {nw} read-tracking "
                f"warning{'s' if nw != 1 else ''}. Each line below names a sample (or a "
                f"single step) whose read retention crossed a configured threshold (a "
                f"sample retaining less than {self.warn_pct:.0f}% of its raw reads overall, or "
                f"one step dropping more than {step_loss:.0f}% of a sample's reads). These mark "
                f"where reads were lost; they are not necessarily errors (negative controls are "
                f"<i>expected</i> to retain almost nothing).</p>")
            warn_pre = self._render_warn_terminal(self.warnings)
            if warn_pre is None:  # rich unavailable -- plain escaped fallback
                warn_pre = f'<pre class="runlog">{_esc(chr(10).join(self.warnings))}</pre>'
            parts.append('<p class="warn-head">Read-tracking warnings</p>'
                         + self._terminal_window("read-tracking warnings", warn_pre, compact=True))
        else:
            parts.append('<p class="warn-none">No read-tracking warnings were raised: every '
                         'sample retained reads above the configured thresholds.</p>')
        return "\n".join(p for p in parts if p)

    def _section_per_sample(self, figs) -> str:
        """Per-sample yield: reads retained, feature richness, and % retained."""
        feat = "ASVs" if self.is_dada2 else "OTUs"
        rich = self._richness()
        parts = [f"<p>Sequencing yield per sample: reads retained to the {self.final_step} stage, "
                 f"the number of {feat} detected (features with at least one read), and overall "
                 f"retention. Controls (<code>Blank*</code>) are expected to yield little.</p>"]
        rcap = (f"{feat} detected per sample." if rich is not None
                else f"{feat} richness was unavailable (no per-sample feature table).")
        if rich is not None:
            parts.append(self._fig(figs.get("richness"), rcap))
        else:
            logger.warning("[WARN] html_report: expected=per-sample richness, got=no feature "
                           "table, fallback=reads/retention only")
            parts.append(f"<p>{rcap}</p>")
        rows = []
        for _, r in self.df.iterrows():
            sample = str(r["sample"])
            is_ctrl = sample.lower().startswith("blank")
            name = f'<span class="flag-low">{_esc(sample)}</span>' if is_ctrl else _esc(sample)
            reads = r.get(self.final_step)
            reads_c = '<span class="na">NA</span>' if pd.isna(reads) else f"{int(reads):,}"
            if rich is not None:
                rv = rich.get(sample)
                rich_c = '<span class="na">NA</span>' if rv is None else f"{int(rv):,}"
            else:
                rich_c = '<span class="na">NA</span>'
            pr = r.get("pct_retained")
            pr_c = "NA" if pd.isna(pr) else f"{pr:.1f}%"
            rows.append([name, reads_c, rich_c, pr_c])
        parts.append(self._table(
            f"Per-sample reads, {feat} detected, and retention. Controls are flagged with an asterisk.",
            ["sample", f"{self.final_step} reads", f"{feat} detected", "% retained"],
            rows, scroll=True))
        return "\n".join(p for p in parts if p)

    def _section_taxonomy(self, figs) -> Optional[str]:
        tax = self._tax()
        if tax is None:
            return None
        total = len(tax)
        sp = int((~tax["species"].astype(str).isin(_UNASSIGNED)).sum()) if "species" in tax.columns else 0
        sample_cols = self._tax_sample_cols(tax)
        ctrl = [c for c in sample_cols if str(c).lower().startswith("blank")]
        bio = [c for c in sample_cols if c not in ctrl] or sample_cols
        parts = [f"<p>Of {total:,} features, {sp:,} ({sp / total * 100:.1f}%) were assigned to species. "
                 f"Counts at each rank are shown in Figure {self._fig_n + 1}; the best-hit identity "
                 f"distribution (assigned features only) in the following figure. Features with no "
                 f"assignment are reported as <code>Unassigned</code>, never as a taxon.</p>"]
        parts.append(self._fig(figs.get("rank"),
                     f"Number of features with a non-empty assignment at each taxonomic rank, out of "
                     f"{total:,} total."))
        parts.append(self._fig(figs.get("pident"),
                     "Distribution of best-hit BLAST percent identity over assigned features; "
                     "dashed line marks the median."))
        # full rank-assignment table (exact counts behind the figure)
        rank_rows = []
        for rk in _RANKS:
            if rk in tax.columns:
                a = int((~tax[rk].astype(str).isin(_UNASSIGNED)).sum())
                rank_rows.append([rk, f"{a:,}", f"{a / total * 100:.1f}%"])
        if rank_rows:
            parts.append(self._table("Features assigned at each taxonomic rank.",
                                     ["rank", "features assigned", "% of total"], rank_rows))
        # detected species: reads + occupancy (number of biological samples)
        if "species" in tax.columns and bio:
            bio_num = tax[bio].apply(pd.to_numeric, errors="coerce").fillna(0)
            grp = bio_num.groupby(tax["species"].astype(str)).sum()
            sp_reads = grp.sum(axis=1).sort_values(ascending=False)
            occ = (grp > 0).sum(axis=1)
            rows = []
            for name, val in sp_reads.items():
                if name in _UNASSIGNED:
                    continue
                rows.append([_esc(name.replace("_", " ")), f"{int(val):,}",
                             f"{int(occ.get(name, 0)):,}"])
                if len(rows) >= 20:
                    break
            if rows:
                n_sp = sum(1 for k in sp_reads.index if k not in _UNASSIGNED)
                cap = (f"Top species by total read count across biological samples, with occupancy "
                       f"(number of the {len(bio)} biological samples each was detected in). "
                       f"{n_sp:,} species detected in total.")
                parts.append(self._table(cap, ["species", "reads", "samples"], rows, scroll=True))
        if "genus" in tax.columns:
            gc = tax.loc[~tax["genus"].astype(str).isin(_UNASSIGNED), "genus"].value_counts().head(15)
            if not gc.empty:
                parts.append(self._table("Genera with the most assigned features.",
                                         ["genus", "features"],
                                         [[_esc(str(k).replace("_", " ")), f"{int(v):,}"] for k, v in gc.items()]))
        return "\n".join(p for p in parts if p)

    def _section_feature_qc(self, figs) -> Optional[str]:
        otu = self._otu_full()
        if otu is None:
            return None
        parts = []
        if "chimera" in otu.columns:
            vc = otu["chimera"].astype(str).value_counts()
            parts.append(f"<p>Of {len(otu):,} candidate OTUs, {int(vc.get('Y', 0)):,} were flagged "
                         f"chimeric and removed; {int(vc.get('?', 0)):,} borderline cases were retained.</p>")
            parts.append(self._fig(figs.get("chimera"), "De-novo chimera classification of all candidate OTUs."))
        parts.append(self._fig(figs.get("length"),
                     "OTU sequence-length distribution; dashed line marks the median length."))
        return "\n".join(p for p in parts if p) or None

    def _section_contamination(self) -> str:
        tax = self._tax()
        if tax is None:
            return ("<p>No taxonomy table was available, so contamination screening against "
                    "controls was not performed.</p>")
        sample_cols = self._tax_sample_cols(tax)
        ctrl = [c for c in sample_cols if str(c).lower().startswith("blank")]
        if not ctrl:
            return ("<p>No negative-control columns (<code>Blank*</code>) were found in the taxonomy "
                    "table; contamination screening was not performed.</p>")
        bio = [c for c in sample_cols if c not in ctrl]
        blanks = tax[ctrl].apply(pd.to_numeric, errors="coerce").fillna(0)
        in_blank = blanks.sum(axis=1) > 0
        sub = tax[in_blank]
        if sub.empty:
            return (f"<p>None of the {len(tax):,} features carried reads in the {len(ctrl)} "
                    f"control{'s' if len(ctrl) != 1 else ''} (<code>{', '.join(ctrl)}</code>).</p>")
        sub_sorted = sub.assign(_b=blanks[in_blank].sum(axis=1)).sort_values("_b", ascending=False)
        rows = []
        for _, r in sub_sorted.head(15).iterrows():
            taxon = next((str(r[c]) for c in reversed(_RANKS)
                          if c in tax.columns and str(r[c]) not in _UNASSIGNED), "Unassigned")
            cells = [_esc(taxon.replace("_", " "))]
            cells += [f"{int(pd.to_numeric(r[c], errors='coerce') or 0):,}" for c in ctrl]
            cells.append(f"{int(pd.to_numeric(r[bio], errors='coerce').fillna(0).sum()):,}" if bio else "n/a")
            rows.append(cells)
        intro = (f"<p>{int(in_blank.sum())} feature(s) carried reads in a control. Controls are "
                 f"identified by the <code>Blank*</code> name prefix and contamination is computed "
                 f"from blank read counts (not a precomputed flag).</p>")
        return intro + self._table("Features detected in negative controls, ranked by total control reads.",
                                    ["feature taxon", *ctrl, "total in samples"], rows)

    def _section_timeline(self) -> Optional[str]:
        steps = self.state.get("steps") if isinstance(self.state, dict) else None
        if not isinstance(steps, dict) or not steps:
            return None
        rows = []
        for name, s in steps.items():
            if not isinstance(s, dict):
                continue
            dur = s.get("duration_seconds")
            rows.append([_esc(name), _esc(str(s.get("status", "")).lower()),
                         "n/a" if dur is None else f"{float(dur):.1f} s"])
        if not rows:
            return None
        note = ("<p>Per-step wall-clock timing from the run state. Report and export steps run after "
                "the recorded completion and are not timed here.</p>")
        return note + self._table("Pipeline step status and duration.", ["step", "status", "duration"], rows)

    # ------------------------------------------------------------------ #
    # Run log (colorized console transcript)
    # ------------------------------------------------------------------ #
    @staticmethod
    def _line_level(line: str) -> Optional[str]:
        """Classify a log line by level from the ``| LEVEL |`` field, with a
        text-heuristic fallback for continuation lines (tracebacks, tool output)."""
        parts = line.split(" | ", 3)
        if len(parts) >= 2:
            lvl = parts[1].strip().upper()
            if lvl in _LOG_LEVELS:
                return lvl
        upper = line.upper()
        if "[ERROR]" in upper or "[FAIL]" in upper or "TRACEBACK" in upper or "EXCEPTION" in upper:
            return "ERROR"
        if "[WARN]" in upper or "WARNING" in upper:
            return "WARNING"
        return None

    def _select_log_lines(self, lines: List[str]) -> Tuple[List[Tuple[int, str]], bool, int]:
        """Pick which log lines to embed.

        Short logs are shown whole. Long logs keep run start/end context and
        *every* warning/error line; intervening routine lines are collapsed into
        explicit ``… N omitted …`` markers (never silently dropped, CLAUDE.md
        section 4). Returns ``(items, truncated, total_lines)`` where ``items`` is
        a list of ``(index, text)`` and a marker has index ``-1``.
        """
        n = len(lines)
        if n <= self.max_log_lines:
            return [(i, ln) for i, ln in enumerate(lines)], False, n

        head_n = tail_n = 45
        head = set(range(min(head_n, n)))
        tail = set(range(max(0, n - tail_n), n))
        events = [i for i, ln in enumerate(lines) if self._line_level(ln) in _EVENT_LEVELS]
        budget = max(0, self.max_log_lines - len(head) - len(tail))
        kept_events = set(events[:budget])
        if len(events) > budget:
            logger.warning(
                f"[WARN] html_report: expected=embed all run-log events, "
                f"got={len(events)} events > budget {budget}, "
                f"fallback=keep earliest {budget}; full log on disk",
            )
        logger.warning(
            f"[WARN] html_report: expected=embed full run log, got={n} lines "
            f"> max {self.max_log_lines}, fallback=head/tail + {len(kept_events)} "
            f"event line(s); full log on disk",
        )

        keep = sorted(head | tail | kept_events)
        items: List[Tuple[int, str]] = []
        prev: Optional[int] = None
        for i in keep:
            if prev is not None and i > prev + 1:
                gap = i - prev - 1
                items.append((-1, f"      … {gap:,} routine line(s) omitted …"))
            items.append((i, lines[i]))
            prev = i
        return items, True, n

    @staticmethod
    def _themed_console():
        """A rich recording console with the bright dark-terminal log palette
        (info blue, warning amber, error red). Raises ImportError if rich is
        absent so callers can fall back to plain text."""
        import io as _io

        from rich.console import Console
        from rich.theme import Theme
        theme = Theme({
            "log.time": "#7d8590",
            "log.logger": "#6b7280",
            "logging.level.debug": "#56b6c2",
            "logging.level.info": "#6cb6ff",
            "logging.level.warning": "#e3b341",
            "logging.level.error": "bold #ff6b6b",
            "logging.level.critical": "bold reverse #ff6b6b",
        }, inherit=True)
        return Console(record=True, width=400, file=_io.StringIO(),
                       force_terminal=True, highlight=False, theme=theme)

    @staticmethod
    def _export(con) -> str:
        return con.export_html(inline_styles=True, code_format='<pre class="runlog">{code}</pre>')

    def _render_log_html(self, items: List[Tuple[int, str]]) -> Optional[str]:
        """Render selected log lines to a self-contained, level-colored ``<pre>``
        using rich's own HTML export, so the palette matches the live console
        exactly. Returns ``None`` if rich is unavailable (caller falls back)."""
        try:
            from rich.text import Text
            con = self._themed_console()
        except ImportError as exc:  # pragma: no cover -- rich is a hard dependency
            logger.warning(f"[WARN] html_report: expected=rich for colorized log, got=missing "
                           f"({exc}), fallback=plain monospace log")
            return None

        for idx, text in items:
            if idx == -1:  # omission marker
                con.print(Text(text, style="#7d8590 italic"), soft_wrap=True, highlight=False)
                continue
            parts = text.split(" | ", 3)
            line = Text(no_wrap=True)
            if len(parts) == 4 and parts[1].strip().upper() in _LOG_LEVELS:
                line.append(parts[0] + " | ", style="log.time")
                line.append(parts[1] + " | ", style=f"logging.level.{parts[1].strip().lower()}")
                line.append(parts[2] + " | ", style="log.logger")
                line.append(parts[3])
            else:
                lvl = self._line_level(text)
                line.append(text, style=f"logging.level.{lvl.lower()}" if lvl else "")
            con.print(line, soft_wrap=True, highlight=False)
        return self._export(con)

    def _render_warn_terminal(self, warnings: List[str]) -> Optional[str]:
        """Render read-tracking warnings as a colorized CLI transcript (the same
        terminal palette as the run log). Returns ``None`` if rich is absent."""
        try:
            import re

            from rich.text import Text
            con = self._themed_console()
        except ImportError:
            return None
        pat = re.compile(r"^\[WARN\]\s+(read_tracking)\s+(.+?):\s*(.*)$")
        for w in warnings:
            line = Text(no_wrap=True)
            m = pat.match(w)
            if m:
                line.append("[WARN] ", style="logging.level.warning")
                line.append(m.group(1) + " ", style="log.logger")
                line.append(m.group(2) + ": ", style="logging.level.info")
                line.append(m.group(3))
            else:
                line.append(w, style="logging.level.warning")
            con.print(line, soft_wrap=True, highlight=False)
        return self._export(con)

    def _terminal_window(self, title: str, body_html: str, *,
                         fullscreen: bool = False, compact: bool = False) -> str:
        """Wrap pre-rendered terminal HTML in window chrome (dots + title bar).

        ``fullscreen`` adds a pure-CSS maximize toggle; ``compact`` caps the
        body height (for short transcripts such as the warnings list)."""
        body_cls = "term-body compact" if compact else "term-body"
        toggle, button = "", ""
        if fullscreen:
            toggle = '<input type="checkbox" id="termmax" class="term-max-toggle">'
            button = ('<label for="termmax" class="term-max-btn">'
                      '<span class="lbl-open">Fullscreen</span>'
                      '<span class="lbl-close">Exit fullscreen</span></label>')
        return (
            f'{toggle}<div class="terminal">'
            '<div class="term-bar">'
            '<span class="term-dot r"></span><span class="term-dot y"></span>'
            '<span class="term-dot g"></span>'
            f'<span class="term-title">{_esc(title)}</span>{button}</div>'
            f'<div class="{body_cls}">{body_html}</div></div>'
        )

    def _section_run_log(self) -> Optional[str]:
        """Embed the full pipeline run log, colorized by level like the console.

        A missing/disabled/unreadable log yields an explanatory sentence rather
        than vanishing silently (CLAUDE.md section 4)."""
        if self.log_file is None:
            return ("<p>No run-log file was passed to the report (logging to file may have been "
                    "disabled), so the console transcript is not embedded here.</p>")
        path = self.log_file
        if not path.exists():
            logger.warning(f"[WARN] html_report: expected=run log, got=missing "
                           f"({path}), fallback=note absence in report")
            return (f"<p>The run-log file <code>{_esc(path.name)}</code> was not found, so the "
                    f"console transcript could not be embedded.</p>")
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError as exc:
            logger.warning(f"[WARN] html_report: expected=readable run log, got=unreadable "
                           f"({path}: {exc}), fallback=omit log section")
            return (f"<p>The run-log file <code>{_esc(path.name)}</code> could not be read "
                    f"({_esc(str(exc))}).</p>")
        if not lines:
            return f"<p>The run-log file <code>{_esc(path.name)}</code> is empty.</p>"

        items, truncated, total = self._select_log_lines(lines)
        pre = self._render_log_html(items)
        if pre is None:  # rich unavailable -- plain escaped fallback, still self-contained
            body = _esc("\n".join(t for _, t in items))
            pre = f'<pre class="runlog">{body}</pre>'

        counts: Dict[str, int] = {}
        for ln in lines:
            lvl = self._line_level(ln)
            if lvl:
                counts[lvl] = counts.get(lvl, 0) + 1
        summary_bits = [f"{total:,} lines"]
        for lvl in ("WARNING", "ERROR", "CRITICAL"):
            if counts.get(lvl):
                summary_bits.append(f"{counts[lvl]:,} {lvl.lower()}")

        intro = ("<p>Complete console transcript of the run "
                 f"(<code>{_esc(path.name)}</code>), colorized by log level exactly as the live "
                 "SeeDNAP console renders it: "
                 "<span class=\"term-legend\"><span class=\"lvl-chip lvl-info\">info</span> "
                 "<span class=\"lvl-chip lvl-warning\">warning</span> "
                 "<span class=\"lvl-chip lvl-error\">error</span></span>. ")
        if truncated:
            intro += (f"The log is long ({total:,} lines), so every warning and error is kept "
                      "alongside the run's start and end, and intervening routine lines are "
                      "collapsed (markers state how many). The complete log is on disk at "
                      f"<code>{_esc(str(path))}</code>.</p>")
        else:
            intro += "</p>"

        # Render the transcript inside a real terminal window (chrome + dark
        # body) with a pure-CSS Fullscreen toggle for easier reading.
        meta = f'<p class="runlog-meta">{", ".join(summary_bits)}</p>'
        terminal = self._terminal_window(path.name, pre, fullscreen=True)
        return intro + meta + terminal

    def _summary_table_html(self) -> str:
        df = self.df
        n = len(df)
        rows = [["samples", f"{n:,}"]]
        nc = self._n_controls()
        if nc:
            rows.append(["controls", f"{nc:,}"])
        raw = pd.to_numeric(df.get("raw", pd.Series(dtype=float)), errors="coerce")
        fin = pd.to_numeric(df.get(self.final_step, pd.Series(dtype=float)), errors="coerce")
        if raw.notna().any():
            rows.append(["raw read pairs", self._fmt(raw.sum())])
        if fin.notna().any():
            rows.append([f"{self.final_step} reads", self._fmt(fin.sum())])
        pr = pd.to_numeric(df.get("pct_retained", pd.Series(dtype=float)), errors="coerce")
        if pr.notna().any():
            rows.append(["mean retention", f"{pr.mean():.1f}%"])
        tax = self._tax()
        if tax is not None:
            rows.append([f"{'ASVs' if self.is_dada2 else 'OTUs'} (features)", f"{len(tax):,}"])
            for rank in ("species", "genus", "phylum"):
                if rank in tax.columns and len(tax):
                    nas = int((~tax[rank].astype(str).isin(_UNASSIGNED)).sum())
                    rows.append([f"{rank} assigned", f"{nas:,} ({nas / len(tax) * 100:.1f}%)"])
            if "pident" in tax.columns:
                pid = pd.to_numeric(tax["pident"], errors="coerce").dropna()
                if not pid.empty:
                    rows.append(["median identity", f"{pid.median():.1f}%"])
        otu = self._otu_full()
        if otu is not None and "chimera" in otu.columns:
            rows.append(["chimeras removed", f"{int(otu['chimera'].astype(str).eq('Y').sum()):,}"])
        return self._table("Run summary.", ["quantity", "value"], rows)

    def _methods(self) -> str:
        s = self.summary
        parts = [
            "<p>Controls are identified by the <code>Blank*</code> sample-name prefix; "
            "contamination is computed from blank read counts, not a precomputed flag. "
            "Features with no assignment are reported as <code>Unassigned</code> and never counted "
            "as a taxon; median percent identity is computed over assigned features only.</p>",
            f"<p>Retention thresholds: a sample is flagged below "
            f"{self.warn_pct:.0f}% overall retention; a step is flagged when it drops more than "
            f"{float(s.get('warn_step_loss_pct', 70.0)):.0f}% of a sample's reads.</p>",
        ]
        if s.get("footer"):
            parts.append(f"<p>{_esc(s['footer'])}</p>")
        parts.append("<p>Generated by SeeDNAP. Self-contained report (no external assets).</p>")
        return "".join(parts)

    # ------------------------------------------------------------------ #
    # Render
    # ------------------------------------------------------------------ #
    def render(self) -> str:
        self._fig_n = 0
        self._tbl_n = 0
        figs = self._make_figures()

        # The summary table is Table 1 and leads the Summary tab.
        summary_table = self._summary_table_html()

        sections: List[Dict[str, str]] = []
        sections.append({"title": "Summary", "html": self._section_summary(summary_table)})
        sections.append({"title": "Dataset", "html": self._section_dataset()})
        sections.append({"title": "Read tracking", "html": self._section_read_tracking(figs)})
        sections.append({"title": "Per-sample detail", "html": self._section_per_sample(figs)})
        tax_html = self._section_taxonomy(figs)
        if tax_html:
            sections.append({"title": "Taxonomic assignment", "html": tax_html})
        fqc = self._section_feature_qc(figs)
        if fqc:
            sections.append({"title": "OTU / feature QC", "html": fqc})
        sections.append({"title": "Controls & contamination", "html": self._section_contamination()})
        tl = self._section_timeline()
        if tl:
            sections.append({"title": "Run provenance", "html": tl})
        runlog = self._section_run_log()
        if runlog:
            sections.append({"title": "Run log", "html": runlog})
        sections.append({"title": "Notes & methods", "html": self._methods()})

        return _TEMPLATE.render(
            marker=_esc(self.marker),
            sections=sections,
        )

    def _section_summary(self, summary_table: str) -> str:
        """Lead tab: a one-line descriptor, the run abstract, and Table 1."""
        return (f'<p class="descriptor">{self._descriptor()}. {_esc(self._run_date())}.</p>'
                f'<p class="summary-lead">{self._abstract()}</p>'
                f'{summary_table}')

    def write(self, output_path: Union[str, Path]) -> Path:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.render(), encoding="utf-8")
        logger.info(f"Wrote HTML run report: {path}")
        return path
