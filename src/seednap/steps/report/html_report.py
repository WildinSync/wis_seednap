"""Self-contained HTML run report (opt-in, ``report.html_report``).

Renders a single portable ``.html`` file in the SeeDNAP brand theme (dark
green-black + neon-green/teal accents). Charts are matplotlib figures embedded
as base64 PNGs, so there are no external assets, no CDN, and no JavaScript.

Beyond the read-tracking funnel it surfaces, when the sources are available,
a run timeline (state JSON), a taxonomy headline (taxonomy CSV) and a chimera
summary (SWARM ``otu_table_full.csv``). Every extra source is optional: a
missing/unreadable one is skipped with a ``[WARN]`` (CLAUDE.md section 4), and
matplotlib is imported lazily so the report still renders (tables only) if it
is absent.
"""

import base64
import io
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
from jinja2 import Template

from seednap.steps.report.read_tracking import DADA2_STEPS
from seednap.utils.logging import get_logger

logger = get_logger(__name__)

_DEFAULT_STEPS = DADA2_STEPS

# Single source of truth for colors -- the CSS :root variables below and the
# matplotlib chart colors are both derived from this (sampled from the SeeDNAP
# teaser: dark green-black bg, neon-green + teal accents).
THEME = {
    "bg": "#041912", "card": "#0c241b", "card2": "#0f2d22",
    "line": "#1b3b2e", "ink": "#eafff4", "ink_soft": "#bfe9d6", "muted": "#6f9d88",
    "accent": "#44f187", "accent2": "#35ecc4", "warn": "#fbbf24", "bad": "#f87171",
}

_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
_TAX_META = {"ASV_ID", "OTU_ID", "pident", "is_contaminant_candidate", "Sequence", "sequence"}
_UNASSIGNED = {"Unassigned", "unassigned", "", "NA", "nan"}

_TEMPLATE = Template(
    """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SeeDNAP run report &mdash; {{ marker }}</title>
<style>
  :root {
    --bg:{{T.bg}}; --card:{{T.card}}; --card2:{{T.card2}}; --line:{{T.line}};
    --ink:{{T.ink}}; --ink-soft:{{T.ink_soft}}; --muted:{{T.muted}};
    --accent:{{T.accent}}; --accent2:{{T.accent2}}; --warn:{{T.warn}}; --bad:{{T.bad}};
    --glow:0 0 0 1px rgba(68,241,135,.16), 0 10px 34px rgba(0,0,0,.5); --radius:14px;
  }
  * { box-sizing:border-box; }
  body { margin:0; color:var(--ink); line-height:1.55;
         font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         background:
           radial-gradient(1100px 380px at 12% -8%, rgba(68,241,135,.10), transparent 60%),
           radial-gradient(900px 360px at 96% 0%, rgba(53,236,196,.09), transparent 55%),
           var(--bg);
         background-attachment:fixed; font-variant-numeric:tabular-nums; }
  .wrap { max-width:1140px; margin:0 auto; padding:30px 24px 70px; }
  .hero { display:flex; align-items:center; gap:18px; border-bottom:1px solid var(--line);
          padding-bottom:22px; margin-bottom:26px; }
  .logo { width:52px; height:52px; border-radius:14px; flex:none; display:grid; place-items:center;
          font-weight:800; font-size:26px; color:#04231a;
          background:linear-gradient(135deg,var(--accent),var(--accent2)); box-shadow:var(--glow); }
  .hero h1 { margin:0; font-size:25px; letter-spacing:.3px; }
  .hero h1 .mk { background:linear-gradient(90deg,var(--accent),var(--accent2));
                 -webkit-background-clip:text; background-clip:text; color:transparent; font-weight:800; }
  .brand { color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:1.4px; }
  .sub { color:var(--ink-soft); font-size:13.5px; margin-top:2px; }
  .pill { margin-left:auto; padding:7px 14px; border-radius:999px; font-size:13px; font-weight:700;
          align-self:center; }
  .pill.good { background:rgba(68,241,135,.14); color:var(--accent); border:1px solid rgba(68,241,135,.4); }
  .pill.warn { background:rgba(251,191,36,.14); color:var(--warn); border:1px solid rgba(251,191,36,.4); }
  .pill.bad  { background:rgba(248,113,113,.14); color:var(--bad); border:1px solid rgba(248,113,113,.4); }
  h2 { font-size:13px; text-transform:uppercase; letter-spacing:1.5px; color:var(--accent2);
       margin:38px 0 14px; display:flex; align-items:center; gap:10px; }
  h2::before { content:attr(data-n); color:var(--muted); font-weight:700; }
  h2::after { content:""; flex:1; height:1px; background:linear-gradient(90deg,var(--line),transparent); }
  .cards { display:grid; grid-template-columns:repeat(auto-fit,minmax(165px,1fr)); gap:13px; }
  .card { position:relative; background:linear-gradient(180deg,var(--card2),var(--card));
          border:1px solid var(--line); border-radius:var(--radius); padding:15px 16px 14px 18px;
          overflow:hidden; }
  .card::before { content:""; position:absolute; left:0; top:0; bottom:0; width:3px;
                  background:linear-gradient(180deg,var(--accent),var(--accent2)); }
  .card .v { font-size:25px; font-weight:750; }
  .card .l { color:var(--muted); font-size:11px; text-transform:uppercase; letter-spacing:.7px; margin-top:3px; }
  .card .s { color:var(--ink-soft); font-size:11.5px; margin-top:2px; }
  .meter { height:6px; border-radius:5px; background:rgba(255,255,255,.07); margin-top:9px; overflow:hidden; }
  .meter > i { display:block; height:100%; border-radius:5px; }
  .grid2 { display:grid; grid-template-columns:repeat(auto-fit,minmax(330px,1fr)); gap:14px; }
  .chart { background:var(--card); border:1px solid var(--line); border-radius:var(--radius); padding:14px; }
  .chart img { width:100%; height:auto; display:block; border-radius:8px; }
  .panel { background:var(--card); border:1px solid var(--line); border-radius:var(--radius); padding:4px 0; }
  .note { color:var(--muted); font-size:12px; margin:8px 2px 14px; font-style:italic; }
  table { width:100%; border-collapse:collapse; font-size:13px; }
  th,td { padding:7px 12px; text-align:right; border-bottom:1px solid var(--line); white-space:nowrap; }
  th:first-child, td:first-child { text-align:left; }
  thead th { color:var(--muted); text-transform:uppercase; font-size:10.5px; letter-spacing:.5px;
             position:sticky; top:0; background:var(--card); }
  .scroll { max-height:560px; overflow:auto; border:1px solid var(--line); border-radius:var(--radius); }
  .na { color:var(--muted); }
  .ret { display:inline-block; min-width:58px; padding:1px 8px; border-radius:999px; font-weight:600; font-size:12px; }
  .ret.g { background:rgba(68,241,135,.13); color:var(--accent); }
  .ret.w { background:rgba(251,191,36,.13); color:var(--warn); }
  .ret.b { background:rgba(248,113,113,.13); color:var(--bad); }
  .timeline { display:flex; gap:10px; flex-wrap:wrap; }
  .tstep { background:var(--card); border:1px solid var(--line); border-radius:12px; padding:11px 15px; min-width:130px; }
  .tstep .n { font-weight:700; text-transform:capitalize; }
  .tstep .d { color:var(--accent2); font-size:13px; }
  .tstep .st { font-size:11px; text-transform:uppercase; letter-spacing:.6px; }
  .tstep .st.completed { color:var(--accent); } .tstep .st.failed { color:var(--bad); }
  .warns { background:var(--card); border:1px solid var(--line); border-left:3px solid var(--warn);
           border-radius:10px; padding:4px 0; max-height:280px; overflow:auto; }
  .warns div { padding:6px 16px; font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
               font-size:12px; color:var(--warn); border-bottom:1px solid var(--line); }
  .warns div:last-child { border-bottom:none; }
  .ok { color:var(--accent); padding:11px 16px; }
  footer { color:var(--muted); font-size:12px; margin-top:44px; border-top:1px solid var(--line); padding-top:14px; }
</style></head>
<body><div class="wrap">

<div class="hero">
  <div class="logo">S</div>
  <div>
    <div class="brand">SeeDNAP &middot; eDNA metabarcoding</div>
    <h1>Run report &mdash; <span class="mk">{{ marker }}</span></h1>
    <div class="sub">{{ subtitle }}</div>
  </div>
  <div class="pill {{ health.cls }}">{{ health.label }}</div>
</div>

<div class="cards">
  {% for c in cards %}<div class="card"><div class="v">{{ c.value }}</div><div class="l">{{ c.label }}</div>
    {% if c.sub %}<div class="s">{{ c.sub }}</div>{% endif %}
    {% if c.meter is not none %}<div class="meter"><i style="width:{{ c.meter }}%;background:{{ c.meter_color }}"></i></div>{% endif %}
  </div>{% endfor %}
</div>

{% if timeline %}<h2 data-n="01">Run timeline</h2>
<div class="timeline">
  {% for t in timeline %}<div class="tstep"><div class="n">{{ t.name }}</div>
    <div class="d">{{ t.duration }}</div><div class="st {{ t.status }}">{{ t.status }}</div></div>{% endfor %}
</div>
{% if timeline_note %}<div class="note">{{ timeline_note }}</div>{% endif %}
{% endif %}

{% if tax_section %}<h2 data-n="02">Taxonomy</h2>
{% if rank_chart %}<div class="chart" style="margin-bottom:14px"><img alt="rank completeness" src="data:image/png;base64,{{ rank_chart }}"></div>{% endif %}
<div class="note">{{ tax_note }}</div>
<div class="grid2">
  {% for tbl in tax_tables %}<div class="panel">
    <table><thead><tr><th>{{ tbl.title }}</th><th>{{ tbl.col }}</th></tr></thead><tbody>
    {% for r in tbl.rows %}<tr><td>{{ r.name }}</td><td>{{ r.value }}</td></tr>{% endfor %}
    </tbody></table></div>{% endfor %}
</div>
{% endif %}

{% if charts %}<h2 data-n="03">Reads &amp; features</h2>
<div class="grid2">
  {% for ch in charts %}<div class="chart"><img alt="{{ ch.alt }}" src="data:image/png;base64,{{ ch.b64 }}"></div>{% endfor %}
</div>
{% endif %}

{% if contamination %}<h2 data-n="04">Controls &amp; contamination</h2>
<div class="note">{{ contamination.note }}</div>
<div class="panel"><div class="scroll"><table><thead><tr><th>feature taxon</th>{% for c in contamination.controls %}<th>{{ c }}</th>{% endfor %}<th>total in samples</th></tr></thead><tbody>
{% for r in contamination.rows %}<tr><td>{{ r.taxon }}</td>{% for v in r.blank_vals %}<td>{{ v }}</td>{% endfor %}<td>{{ r.sample_total }}</td></tr>{% endfor %}
</tbody></table></div></div>
{% endif %}

<h2 data-n="05">Notable events {% if warnings %}({{ warnings|length }}){% endif %}</h2>
{% if warnings %}<div class="warns">{% for w in warnings %}<div>{{ w }}</div>{% endfor %}</div>
{% else %}<div class="ok">No data-loss warnings &mdash; every sample passed the retention thresholds.</div>{% endif %}

<h2 data-n="06">Read tracking</h2>
<div class="scroll"><table><thead><tr><th>sample</th>{% for s in steps %}<th>{{ s }}</th>{% endfor %}<th>% retained</th></tr></thead>
<tbody>
{% for r in rows %}<tr><td>{{ r.sample }}</td>{% for c in r.cells %}<td class="{{ c.cls }}">{{ c.txt }}</td>{% endfor %}<td><span class="ret {{ r.retcls }}">{{ r.pct }}</span></td></tr>{% endfor %}
</tbody></table></div>

<footer>Generated by SeeDNAP &middot; self-contained report (no external assets){% if footer %} &middot; {{ footer }}{% endif %}</footer>
</div></body></html>
"""
)


class HTMLReportBuilder:
    """Render a self-contained, SeeDNAP-themed HTML run report."""

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
        self._tax_cache: Optional[pd.DataFrame] = None

    # ------------------------------------------------------------------
    # Lazy data sources (degrade gracefully + [WARN])
    # ------------------------------------------------------------------

    def _tax(self) -> Optional[pd.DataFrame]:
        if self._tax_cache is not None:
            return self._tax_cache
        if self.taxonomy_csv is None:
            return None
        if not self.taxonomy_csv.exists():
            logger.warning(
                f"[WARN] html_report: expected=taxonomy_csv for taxonomy panels, "
                f"got=missing ({self.taxonomy_csv}), fallback=omit taxonomy section",
            )
            return None
        try:
            self._tax_cache = pd.read_csv(self.taxonomy_csv)
        except (pd.errors.EmptyDataError, OSError) as exc:
            logger.warning(
                f"[WARN] html_report: expected=readable taxonomy_csv, got=unreadable "
                f"({self.taxonomy_csv}: {exc}), fallback=omit taxonomy section",
            )
            return None
        return self._tax_cache

    def _tax_sample_cols(self, tax: pd.DataFrame) -> List[str]:
        return [c for c in tax.columns if c not in _TAX_META and c not in _RANKS]

    def _otu_full(self) -> Optional[pd.DataFrame]:
        if self.otu_table_full is None:
            return None
        if not self.otu_table_full.exists():
            logger.warning(
                f"[WARN] html_report: expected=otu_table_full for chimera summary, "
                f"got=missing ({self.otu_table_full}), fallback=omit chimera panel",
            )
            return None
        try:
            return pd.read_csv(self.otu_table_full)
        except (pd.errors.EmptyDataError, OSError) as exc:
            logger.warning(
                f"[WARN] html_report: expected=readable otu_table_full, got=unreadable "
                f"({self.otu_table_full}: {exc}), fallback=omit chimera panel",
            )
            return None

    # ------------------------------------------------------------------
    # Charts (lazy matplotlib)
    # ------------------------------------------------------------------

    def _mpl(self):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            return plt
        except ImportError as exc:
            logger.warning(
                f"[WARN] html_report: expected=matplotlib for charts, "
                f"got=missing ({exc}), fallback=table-only report",
            )
            return None

    @staticmethod
    def _style(ax):
        ax.set_facecolor(THEME["card"])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(THEME["line"])
        ax.tick_params(colors=THEME["muted"], labelsize=8)
        ax.title.set_color(THEME["ink"])
        ax.xaxis.label.set_color(THEME["ink_soft"]); ax.yaxis.label.set_color(THEME["ink_soft"])
        ax.set_axisbelow(True)
        ax.grid(axis="y", color=THEME["line"], linewidth=.6, alpha=.7, linestyle="--")

    def _emit(self, plt, fig) -> str:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor=THEME["card"])
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode()

    def _rank_chart(self, plt) -> Optional[str]:
        tax = self._tax()
        if tax is None:
            return None
        ranks = [r for r in _RANKS if r in tax.columns]
        counts = [int((~tax[r].astype(str).isin(_UNASSIGNED)).sum()) for r in ranks]
        if not any(counts):
            return None
        fig, ax = plt.subplots(figsize=(7.6, 2.9), facecolor=THEME["card"])
        ax.barh(range(len(ranks)), counts, color=THEME["accent"])
        ax.set_yticks(range(len(ranks))); ax.set_yticklabels(ranks)
        ax.invert_yaxis()
        for i, v in enumerate(counts):
            ax.text(v, i, f" {v:,}", va="center", color=THEME["ink_soft"], fontsize=8)
        ax.set_title("Features assigned at each rank")
        self._style(ax); ax.grid(axis="x", color=THEME["line"], linewidth=.6, alpha=.7, linestyle="--")
        ax.grid(axis="y", visible=False)
        return self._emit(plt, fig)

    def _charts(self) -> List[Dict[str, str]]:
        plt = self._mpl()
        if plt is None:
            return []
        charts: List[Dict[str, str]] = []
        warn_pct = float(self.summary.get("warn_below_retention_pct", 30.0))

        # C1: read funnel (totals per step).
        if not self.df.empty:
            totals = [pd.to_numeric(self.df[s], errors="coerce").sum(skipna=True) for s in self.steps]
            if any(t > 0 for t in totals):
                fig, ax = plt.subplots(figsize=(7.6, 3.0), facecolor=THEME["card"])
                bars = ax.bar(self.steps, totals, color=THEME["accent2"])
                ax.bar_label(bars, labels=[f"{int(t):,}" for t in totals],
                             color=THEME["ink_soft"], fontsize=8, padding=2)
                ax.margins(y=.18); ax.set_title("Total reads surviving each step")
                self._style(ax)
                charts.append({"alt": "Reads surviving each step", "b64": self._emit(plt, fig)})

            # C2: per-sample retention.
            pr = pd.to_numeric(self.df["pct_retained"], errors="coerce").dropna()
            if not pr.empty:
                order = pr.sort_values()
                names = self.df.loc[order.index, "sample"].astype(str).tolist()
                colors = [THEME["bad"] if v < warn_pct else THEME["accent"] for v in order]
                fig, ax = plt.subplots(figsize=(7.6, max(2.6, .26 * len(order))), facecolor=THEME["card"])
                ax.barh(range(len(order)), order.values, color=colors)
                ax.set_yticks(range(len(order))); ax.set_yticklabels(names, fontsize=7)
                ax.axvline(warn_pct, color=THEME["warn"], lw=1, ls="--")
                ax.set_xlim(0, 100); ax.set_title("Per-sample retention (raw → final, %)")
                self._style(ax); ax.grid(axis="y", visible=False)
                charts.append({"alt": "Per-sample retention", "b64": self._emit(plt, fig)})

        # C4: pident histogram (assigned rows only).
        tax = self._tax()
        if tax is not None and "pident" in tax.columns:
            pid = pd.to_numeric(tax["pident"], errors="coerce").dropna()
            if not pid.empty:
                fig, ax = plt.subplots(figsize=(7.6, 2.9), facecolor=THEME["card"])
                ax.hist(pid, bins=20, color=THEME["accent"], edgecolor=THEME["card"])
                ax.set_title(f"BLAST %identity of assigned features (n={len(pid):,})")
                ax.set_xlabel("% identity")
                self._style(ax)
                charts.append({"alt": "pident distribution", "b64": self._emit(plt, fig)})

        # C5: chimera mini-bar.
        otu = self._otu_full()
        if otu is not None and "chimera" in otu.columns:
            vc = otu["chimera"].astype(str).value_counts()
            seg = [("clean", int(vc.get("N", 0)), THEME["accent"]),
                   ("chimeric", int(vc.get("Y", 0)), THEME["bad"]),
                   ("borderline", int(vc.get("?", 0)), THEME["warn"])]
            seg = [s for s in seg if s[1] > 0]
            if seg:
                fig, ax = plt.subplots(figsize=(7.6, 1.5), facecolor=THEME["card"])
                left = 0
                for label, val, col in seg:
                    ax.barh([0], [val], left=left, color=col, label=f"{label} ({val:,})")
                    left += val
                ax.set_yticks([]); ax.set_xlabel("OTUs"); ax.set_title("Chimera detection")
                ax.legend(loc="upper center", bbox_to_anchor=(.5, -.5), ncol=3, frameon=False,
                          labelcolor=THEME["ink_soft"], fontsize=8)
                self._style(ax); ax.grid(visible=False)
                charts.append({"alt": "Chimera detection", "b64": self._emit(plt, fig)})

        # C7: OTU length histogram.
        if otu is not None and "length" in otu.columns:
            ln = pd.to_numeric(otu["length"], errors="coerce").dropna()
            if not ln.empty:
                fig, ax = plt.subplots(figsize=(7.6, 2.7), facecolor=THEME["card"])
                ax.hist(ln, bins=30, color=THEME["accent2"], edgecolor=THEME["card"])
                ax.set_title(f"OTU length distribution (median {int(ln.median())} bp)")
                ax.set_xlabel("length (bp)")
                self._style(ax)
                charts.append({"alt": "OTU length distribution", "b64": self._emit(plt, fig)})

        return charts

    # ------------------------------------------------------------------
    # Cards / sections
    # ------------------------------------------------------------------

    @staticmethod
    def _fmt(v) -> str:
        return "—" if v is None or (isinstance(v, float) and pd.isna(v)) else f"{int(v):,}"

    def _health(self) -> Dict[str, str]:
        pr = pd.to_numeric(self.df.get("pct_retained", pd.Series(dtype=float)), errors="coerce").dropna()
        warn_pct = float(self.summary.get("warn_below_retention_pct", 30.0))
        if pr.empty:
            return {"cls": "warn", "label": "no retention data"}
        if (pr < warn_pct).any():
            n = int((pr < warn_pct).sum())
            return {"cls": "bad", "label": f"{n} low-retention sample(s)"}
        if pr.mean() < 60:
            return {"cls": "warn", "label": f"mean retention {pr.mean():.0f}%"}
        return {"cls": "good", "label": f"healthy · mean {pr.mean():.0f}%"}

    def _cards(self) -> List[Dict[str, object]]:
        df = self.df
        n = len(df)
        warn_pct = float(self.summary.get("warn_below_retention_pct", 30.0))
        final_step = self.steps[-1] if self.steps else "raw"
        raw_total = pd.to_numeric(df["raw"], errors="coerce").sum() if n and "raw" in df else None
        final_total = pd.to_numeric(df[final_step], errors="coerce").sum() if n and final_step in df else None
        pr = pd.to_numeric(df["pct_retained"], errors="coerce") if n else pd.Series(dtype=float)
        mean_ret = pr.mean() if not pr.empty else None

        def card(label, value, sub=None, meter=None, meter_color=None):
            return {"label": label, "value": value, "sub": sub, "meter": meter, "meter_color": meter_color}

        cards = [card("samples", f"{n:,}"),
                 card("raw reads", self._fmt(raw_total)),
                 card(f"{final_step} reads", self._fmt(final_total))]
        if mean_ret is not None and pd.notna(mean_ret):
            mc = THEME["bad"] if mean_ret < warn_pct else (THEME["warn"] if mean_ret < 60 else THEME["accent"])
            cards.append(card("mean retention", f"{mean_ret:.1f}%", meter=min(100, max(0, mean_ret)), meter_color=mc))
        else:
            cards.append(card("mean retention", "—"))

        # Taxonomy headline cards.
        tax = self._tax()
        if tax is not None:
            total = len(tax)
            cards.append(card("OTUs / ASVs", f"{total:,}"))
            for rank in ("species", "genus", "phylum"):
                if rank in tax.columns and total:
                    nas = int((~tax[rank].astype(str).isin(_UNASSIGNED)).sum())
                    cards.append(card(f"{rank} assigned", f"{nas:,}", sub=f"{nas / total * 100:.1f}% of OTUs"))
            if "pident" in tax.columns:
                pid = pd.to_numeric(tax["pident"], errors="coerce").dropna()
                if not pid.empty:
                    cards.append(card("median %id", f"{pid.median():.1f}", sub="assigned features only"))
        # SWARM feature/chimera card.
        otu = self._otu_full()
        if otu is not None and "chimera" in otu.columns:
            vc = otu["chimera"].astype(str).value_counts()
            cards.append(card("chimeras removed", f"{int(vc.get('Y', 0)):,}",
                              sub=f"{int(vc.get('?', 0))} borderline kept"))
        return cards

    def _timeline(self):
        steps = self.state.get("steps") if isinstance(self.state, dict) else None
        if not isinstance(steps, dict) or not steps:
            return [], None
        out = []
        for name, s in steps.items():
            if not isinstance(s, dict):
                continue
            dur = s.get("duration_seconds")
            out.append({
                "name": name,
                "duration": "—" if dur is None else f"{float(dur):.1f} s",
                "status": str(s.get("status", "")).lower(),
            })
        note = ("Steps recorded in the run state; report/export steps run after the "
                "pipeline's recorded completion and are not timed here.")
        return out, note

    def _taxonomy_section(self):
        tax = self._tax()
        if tax is None:
            return None
        sample_cols = self._tax_sample_cols(tax)
        ctrl_cols = [c for c in sample_cols if str(c).lower().startswith("blank")]
        bio_cols = [c for c in sample_cols if c not in ctrl_cols] or sample_cols
        tables = []
        # Top taxa by reads (species; Unassigned listed separately).
        if "species" in tax.columns and bio_cols:
            reads = tax[bio_cols].apply(pd.to_numeric, errors="coerce").sum(axis=1)
            by_sp = reads.groupby(tax["species"].astype(str)).sum().sort_values(ascending=False)
            rows = []
            unassigned = int(by_sp.get("Unassigned", 0))
            for name, val in by_sp.items():
                if name in _UNASSIGNED:
                    continue
                rows.append({"name": name.replace("_", " "), "value": f"{int(val):,}"})
                if len(rows) >= 10:
                    break
            tables.append({"title": "Top species (reads)", "col": "reads", "rows": rows})
            # Top genera by OTU count.
            if "genus" in tax.columns:
                gc = tax.loc[~tax["genus"].astype(str).isin(_UNASSIGNED), "genus"].value_counts().head(10)
                tables.append({"title": "Top genera (OTUs)", "col": "OTUs",
                               "rows": [{"name": str(k).replace("_", " "), "value": f"{int(v):,}"} for k, v in gc.items()]})
        total = len(tax)
        sp = int((~tax["species"].astype(str).isin(_UNASSIGNED)).sum()) if "species" in tax.columns else 0
        unassigned_reads = 0
        if "species" in tax.columns and bio_cols:
            unassigned_reads = int(tax.loc[tax["species"].astype(str).isin(_UNASSIGNED), bio_cols]
                                   .apply(pd.to_numeric, errors="coerce").sum().sum())
        note = (f"{sp:,} of {total:,} OTUs ({sp / total * 100:.1f}%) resolved to species. "
                f"'Unassigned' ({unassigned_reads:,} reads) is shown separately, not as a taxon. "
                f"Rank counts come from the taxonomy CSV; %identity is over assigned features only.")
        return {"tables": tables, "note": note}

    def _contamination(self):
        tax = self._tax()
        if tax is None:
            return None
        sample_cols = self._tax_sample_cols(tax)
        ctrl_cols = [c for c in sample_cols if str(c).lower().startswith("blank")]
        if not ctrl_cols:
            return None
        bio_cols = [c for c in sample_cols if c not in ctrl_cols]
        blanks = tax[ctrl_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        in_blank = blanks.sum(axis=1) > 0
        sub = tax[in_blank]
        if sub.empty:
            return None
        rows = []
        sub_sorted = sub.assign(_b=blanks[in_blank].sum(axis=1)).sort_values("_b", ascending=False)
        for _, r in sub_sorted.head(15).iterrows():
            taxon = next((str(r[c]) for c in reversed(_RANKS)
                          if c in tax.columns and str(r[c]) not in _UNASSIGNED), "Unassigned")
            rows.append({
                "taxon": taxon.replace("_", " "),
                "blank_vals": [f"{int(pd.to_numeric(r[c], errors='coerce') or 0):,}" for c in ctrl_cols],
                "sample_total": f"{int(pd.to_numeric(r[bio_cols], errors='coerce').fillna(0).sum()):,}" if bio_cols else "—",
            })
        note = (f"{int(in_blank.sum())} features carry reads in a control. Controls identified by the "
                f"'Blank*' name prefix; contamination computed from blank read counts, NOT a precomputed "
                f"flag. Review these before trusting low-abundance detections.")
        return {"controls": ctrl_cols, "rows": rows, "note": note}

    def _rows(self) -> List[Dict[str, object]]:
        rows = []
        warn_pct = float(self.summary.get("warn_below_retention_pct", 30.0))
        for _, r in self.df.iterrows():
            cells = []
            for s in self.steps:
                v = r[s]
                cells.append({"txt": "NA", "cls": "na"} if pd.isna(v) else {"txt": f"{int(v):,}", "cls": ""})
            pr = r["pct_retained"]
            if pd.isna(pr):
                retcls, pct = "w", "NA"
            else:
                retcls = "b" if pr < warn_pct else ("w" if pr < 60 else "g")
                pct = f"{pr:.1f}%"
            rows.append({"sample": r["sample"], "cells": cells, "pct": pct, "retcls": retcls})
        return rows

    # ------------------------------------------------------------------
    # Render
    # ------------------------------------------------------------------

    def render(self) -> str:
        timeline, timeline_note = self._timeline()
        tax_section = self._taxonomy_section()
        contamination = self._contamination()
        plt = self._mpl()
        rank_chart = self._rank_chart(plt) if plt is not None else None
        return _TEMPLATE.render(
            T=THEME,
            marker=self.marker,
            subtitle=self.summary.get("subtitle") or f"{len(self.df)} samples",
            health=self._health(),
            cards=self._cards(),
            timeline=timeline, timeline_note=timeline_note,
            tax_section=tax_section is not None,
            tax_tables=tax_section["tables"] if tax_section else [],
            tax_note=tax_section["note"] if tax_section else "",
            rank_chart=rank_chart,
            charts=self._charts(),
            contamination=contamination,
            warnings=self.warnings,
            steps=self.steps,
            rows=self._rows(),
            footer=self.summary.get("footer"),
        )

    def write(self, output_path: Union[str, Path]) -> Path:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.render(), encoding="utf-8")
        logger.info(f"Wrote HTML run report: {path}")
        return path
