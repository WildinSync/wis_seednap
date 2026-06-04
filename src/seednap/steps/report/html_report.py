"""Self-contained HTML run report (opt-in, ``report.html_report``).

Renders a single portable ``.html`` file styled like a typeset scientific
paper: a warm-paper page, serif (Computer Modern) typography, justified text,
CSS-numbered sections, ``Figure N`` / ``Table N`` captions, and restrained
monochrome publication figures with a single SeeDNAP-green accent.

Charts are matplotlib PNGs embedded as base64, so there are no external
assets, no CDN, and no JavaScript. It is dataset-agnostic: every number and
label is derived from the data passed in (read-tracking df, optional taxonomy
CSV, optional ``otu_table_full``, optional state JSON). Optional sources are
``[WARN]``-guarded -- a missing one yields an explanatory sentence rather than
vanishing silently (CLAUDE.md section 4) -- and matplotlib is imported lazily
so the report still renders (text + tables) if it is absent.
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
<title>SeeDNAP run report &mdash; {{ marker }}</title>
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
       background:var(--paper); max-width:var(--measure); margin:0 auto; padding:3rem 1.25rem 5rem;
       text-rendering:optimizeLegibility; font-kerning:normal;
       font-feature-settings:"kern" 1,"liga" 1; counter-reset:section;}
  .title-block{text-align:center; margin:0 0 2.2rem;}
  .title-block .org{font-variant-caps:small-caps; letter-spacing:.12em; font-size:.82rem; color:var(--muted); margin:0 0 .5rem;}
  h1{font-size:1.7rem; font-weight:700; line-height:1.25; margin:0 0 .35rem;}
  .title-block .marker{font-style:italic;}
  .title-block .meta{color:var(--muted); font-size:.95rem; margin:.12rem 0;}
  .title-block .descriptor{font-style:italic; margin-top:.35rem;}
  .title-block::after{content:""; display:block; width:4rem; height:2px; background:var(--accent); margin:1.1rem auto 0;}
  .abstract{margin:0 auto 2.2rem; max-width:94%; font-size:.95rem;}
  .abstract .label{text-align:center; font-weight:700; font-variant-caps:small-caps; letter-spacing:.08em; margin:0 0 .45rem;}
  .abstract p{text-align:justify; hyphens:auto;}
  p{margin:0 0 .85em; text-align:justify; hyphens:auto; -webkit-hyphens:auto; hyphenate-limit-chars:6 3 3;}
  h2,h3{font-family:var(--serif); font-weight:700; line-height:1.2; text-align:left; hyphens:manual;}
  h2{font-size:1.28rem; margin:2.0em 0 .55em; counter-increment:section; border-bottom:1px solid var(--hair); padding-bottom:.2em;}
  h2::before{content:counter(section) "\\00a0\\00a0"; color:var(--accent);}
  h2.nonum{counter-increment:none;} h2.nonum::before{content:"";}
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
  .warn-list{font-family:var(--mono); font-size:.8rem; line-height:1.5; border-left:2px solid var(--accent);
             padding:.3rem 0 .3rem 1rem; margin:.8rem 0; max-height:22rem; overflow:auto;}
  .warn-list div{margin:.12rem 0;}
  .warn-none{font-style:italic; color:var(--muted);}
  code{font-family:var(--mono); font-size:.85em;}
  .methods{font-size:.84rem; margin-top:2.4rem; border-top:1px solid var(--hair); padding-top:1rem; color:var(--ink);}
  .methods h2{font-size:1.02rem;}
  @media print{ body{max-width:100%; font-size:10.5pt; line-height:1.4; color:#000; background:#fff; padding:0;}
    h2,h3{break-after:avoid;} figure,table{break-inside:avoid;} @page{margin:2cm;} }
</style></head>
<body>

<div class="title-block">
  <p class="org">SeeDNAP &middot; eDNA metabarcoding pipeline</p>
  <h1>Run report &mdash; marker <span class="marker">{{ marker }}</span></h1>
  <p class="meta">{{ run_date }}</p>
  <p class="descriptor">{{ descriptor }}</p>
</div>

<div class="abstract">
  <p class="label">Summary</p>
  <p>{{ abstract }}</p>
</div>
{{ summary_table }}

{% for s in sections %}<h2>{{ s.title }}</h2>
{{ s.html }}
{% endfor %}

<div class="methods"><h2 class="nonum">Notes &amp; methods</h2>{{ methods }}</div>
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
    ) -> None:
        self.marker = marker
        self.df = tracking_df if tracking_df is not None else pd.DataFrame()
        self.warnings = warnings or []
        self.summary = summary or {}
        self.state = state or {}
        self.taxonomy_csv = Path(taxonomy_csv) if taxonomy_csv else None
        self.otu_table_full = Path(otu_table_full) if otu_table_full else None
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
                return str(v).split("T")[0].split(".")[0]
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
        return " · ".join(bits)

    def _abstract(self) -> str:
        n = len(self.df)
        if n == 0:
            return "No samples were found for this run; the read-tracking table is empty."
        nc = self._n_controls()
        method = "DADA2 ASV" if self.is_dada2 else "SWARM OTU"
        clauses = [
            f"This report summarizes a SeeDNAP run of marker <i>{_esc(self.marker)}</i> over "
            f"{n} sample{'s' if n != 1 else ''}"
            + (f" ({n - nc} biological, {nc} control{'s' if nc != 1 else ''})" if nc else "")
            + f", processed through the {method} path."
        ]
        raw = pd.to_numeric(self.df.get("raw", pd.Series(dtype=float)), errors="coerce")
        fin = pd.to_numeric(self.df.get(self.final_step, pd.Series(dtype=float)), errors="coerce")
        if raw.notna().any() and fin.notna().any() and raw.sum() > 0:
            pct = fin.sum() / raw.sum() * 100
            clauses.append(f"Of {int(raw.sum()):,} raw read pairs, {pct:.1f}% survived to the "
                           f"final {self.final_step} stage.")
        tax = self._tax()
        if tax is not None and "species" in tax.columns and len(tax):
            sp = int((~tax["species"].astype(str).isin(_UNASSIGNED)).sum())
            extra = ""
            if "pident" in tax.columns:
                pid = pd.to_numeric(tax["pident"], errors="coerce").dropna()
                if not pid.empty:
                    extra = f" at a median best-hit identity of {pid.median():.1f}%"
            clauses.append(f"Taxonomic assignment resolved {sp:,} of {len(tax):,} features "
                           f"({sp / len(tax) * 100:.1f}%) to species{extra}.")
        nw = len(self.warnings)
        clauses.append("All steps completed; "
                       + ("no retention warnings were raised." if nw == 0
                          else f"{nw} read-tracking warning{'s' if nw != 1 else ''} "
                               f"are listed in Section 1."))
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
        return "—" if v is None or (isinstance(v, float) and pd.isna(v)) else f"{int(v):,}"

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
            items = "".join(f"<div>{_esc(w)}</div>" for w in self.warnings)
            parts.append(f'<div class="warn-list">{items}</div>')
        else:
            parts.append('<p class="warn-none">No read-tracking warnings were raised.</p>')
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
                 f"distribution (assigned features only) in the following figure. "
                 f"&lsquo;Unassigned&rsquo; is reported separately, never as a taxon.</p>"]
        parts.append(self._fig(figs.get("rank"),
                     f"Number of features with a non-empty assignment at each taxonomic rank, out of "
                     f"{total:,} total."))
        parts.append(self._fig(figs.get("pident"),
                     "Distribution of best-hit BLAST percent identity over assigned features; "
                     "dashed line marks the median."))
        # top species by reads
        if "species" in tax.columns and bio:
            reads = tax[bio].apply(pd.to_numeric, errors="coerce").sum(axis=1)
            by_sp = reads.groupby(tax["species"].astype(str)).sum().sort_values(ascending=False)
            rows = []
            for name, val in by_sp.items():
                if name in _UNASSIGNED:
                    continue
                rows.append([_esc(name.replace("_", " ")), f"{int(val):,}"])
                if len(rows) >= 10:
                    break
            if rows:
                parts.append(self._table("Top species by total read count (biological samples).",
                                         ["species", "reads"], rows))
        if "genus" in tax.columns:
            gc = tax.loc[~tax["genus"].astype(str).isin(_UNASSIGNED), "genus"].value_counts().head(10)
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
            cells.append(f"{int(pd.to_numeric(r[bio], errors='coerce').fillna(0).sum()):,}" if bio else "—")
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
                         "—" if dur is None else f"{float(dur):.1f} s"])
        if not rows:
            return None
        note = ("<p>Per-step wall-clock timing from the run state. Report and export steps run after "
                "the recorded completion and are not timed here.</p>")
        return note + self._table("Pipeline step status and duration.", ["step", "status", "duration"], rows)

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
            "&lsquo;Unassigned&rsquo; features are reported separately and never counted as a taxon; "
            "median percent identity is computed over assigned features only.</p>",
            f"<p>Retention thresholds: a sample is flagged below "
            f"{self.warn_pct:.0f}% overall retention; a step is flagged when it drops more than "
            f"{float(s.get('warn_step_loss_pct', 70.0)):.0f}% of a sample&rsquo;s reads.</p>",
        ]
        if s.get("footer"):
            parts.append(f"<p>{_esc(s['footer'])}</p>")
        parts.append("<p>Generated by SeeDNAP &middot; self-contained report (no external assets).</p>")
        return "".join(parts)

    # ------------------------------------------------------------------ #
    # Render
    # ------------------------------------------------------------------ #
    def render(self) -> str:
        self._fig_n = 0
        self._tbl_n = 0
        figs = self._make_figures()
        summary_table = self._summary_table_html()   # Table 1

        sections: List[Dict[str, str]] = []
        sections.append({"title": "Read tracking", "html": self._section_read_tracking(figs)})
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

        return _TEMPLATE.render(
            marker=_esc(self.marker),
            run_date=_esc(self._run_date()),
            descriptor=self._descriptor(),
            abstract=self._abstract(),
            summary_table=summary_table,
            sections=sections,
            methods=self._methods(),
        )

    def write(self, output_path: Union[str, Path]) -> Path:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.render(), encoding="utf-8")
        logger.info(f"Wrote HTML run report: {path}")
        return path
