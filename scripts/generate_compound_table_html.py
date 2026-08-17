#!/usr/bin/env python3
"""
Turn a compound_summary.tsv (or the all_significant_features_summary.tsv
rollup) into a single self-contained, offline-viewable HTML file: a
sortable/filterable table for scanning a comparison's significant
features and spotting which ones already have a compound identity vs.
still need work. No CDN/network dependency (opened straight from the
HPC filesystem, often via scp/mounted drive) and no build step -- data is
embedded as JSON, sorting/filtering is plain JS.

Design (informed by the dataviz skill + an outside UX consult on
scientist-facing data tables, not a generic sortable-table template):
  - Numeric columns sort as numbers, never as re-parsed text (q-values in
    scientific notation and negative log2FC would sort wrong
    lexicographically) -- the raw numeric value from the TSV is embedded
    directly in the row's JS object and used for both sort and the
    threshold filters; only the rendered <td> text is formatted.
  - Identity is shown as a shape+color glyph (not color alone) in the
    leftmost data column, so "already identified" vs "still needs work"
    reads at a glance while scrolling: filled circle = exact library
    match, filled square = analog library match, hollow diamond = SIRIUS
    structure, hollow circle = SIRIUS formula only, dash = unidentified.
    Same Okabe-Ito palette used throughout this project's plots.
  - Only ~8 columns show in the main table; the other ~15 (library
    cosine/SMILES/InChI, SIRIUS confidence/CANOPUS class, ...) are in a
    per-row expand panel (click the row) rather than a column-visibility
    toggle (no persistence in a static file anyway) or horizontal scroll
    (defeats the point of scanning many rows at once).
  - Three filter affordances, in priority order: click-header sort,
    numeric threshold filters (q-value <=, |log2FC| >=), then a
    best_identity_source category filter (button toggles) and free-text
    search across identity fields last.
  - Plain <table> + DocumentFragment rendering, debounced search input --
    fine up to several thousand rows (checked: this project's tables top
    out around 6,000 for the cross-comparison rollup, well under where
    naive DOM rebuilds start to jank).

Usage:
    python3 scripts/generate_compound_table_html.py [TSV ...]
        (default: every compound_summary.tsv under analysis/differential_features/
        plus the all_significant_features_summary.tsv rollup)
"""
import argparse
import html
import json
import math
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DIFF_ROOT = REPO / "analysis" / "differential_features"
MASTER_TSV = DIFF_ROOT / "all_significant_features_summary.tsv"

# (json_key, header_label, kind) -- kind in {"int", "float", "sci", "pct", "bool", "text"}
PRIMARY_COLS = [
    ("best_identity_source", "ID", "glyph"),
    ("row ID", "row ID", "int"),
    ("row m/z", "m/z", "float"),
    ("row retention time", "RT (min)", "float"),
    ("adduct", "adduct", "text"),
    ("log2FC_a_over_b", "log2FC", "float"),
    ("q_used", "q-value", "sci"),
    ("blocking_robust", "robust", "bool"),
    ("comparison", "comparison", "text"),  # only present in the rollup; dropped if absent
    ("best_identity", "best identity", "text"),
]
SECONDARY_LABELS = {
    "q_value": "Mann-Whitney q-value",
    "q_value_perm_plate": "plate-block-permutation q-value",
    "q_used_source": "q-value source",
    "has_ms2": "has MS2 spectrum",
    "library_NAME": "GNPS exact-match name",
    "library_cosine": "GNPS exact-match cosine",
    "library_matched_peaks": "GNPS exact-match matched peaks",
    "library_SMILES": "GNPS exact-match SMILES",
    "library_INCHI": "GNPS exact-match InChI",
    "library_ORGANISM": "GNPS exact-match organism",
    "analog_NAME": "GNPS analog-match name",
    "analog_cosine": "GNPS analog-match cosine",
    "analog_matched_peaks": "GNPS analog-match matched peaks",
    "analog_SMILES": "GNPS analog-match SMILES",
    "analog_INCHI": "GNPS analog-match InChI",
    "analog_ORGANISM": "GNPS analog-match organism",
    "sirius_formula": "SIRIUS formula",
    "sirius_adduct": "SIRIUS adduct",
    "sirius_structure_name": "SIRIUS structure name",
    "sirius_structure_smiles": "SIRIUS structure SMILES",
    "sirius_structure_confidence": "SIRIUS structure confidence",
    "sirius_npc_pathway": "SIRIUS/CANOPUS NPC pathway",
    "sirius_npc_class": "SIRIUS/CANOPUS NPC class",
    "sirius_classyfire_class": "SIRIUS/CANOPUS ClassyFire class",
    "source_run": "SIRIUS source run(s)",
}

GLYPH_BY_SOURCE = {
    "exact_library_match": {"glyph": "●", "color": "#0072B2", "label": "exact library match"},
    "analog_library_match": {"glyph": "■", "color": "#009E73", "label": "analog library match"},
    "sirius_structure": {"glyph": "◇", "color": "#E69F00", "label": "SIRIUS structure"},
    "sirius_formula_only": {"glyph": "○", "color": "#999999", "label": "SIRIUS formula only"},
    "unidentified": {"glyph": "—", "color": "#BBBBBB", "label": "unidentified"},
}


def clean_value(v):
    if v is None:
        return None
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        # SIRIUS reports ConfidenceScoreExact as -Infinity for some rows
        # (structure found but confidence uncomputable) -- null it out
        # rather than embed a raw Infinity token that would sort as an
        # extreme value and print as literal "-Infinity" in the detail panel.
        return None
    if isinstance(v, bool):
        return v
    return v


def build_html(df: pd.DataFrame, title: str) -> str:
    cols_present = [c for c in PRIMARY_COLS if c[0] in df.columns]
    all_json_cols = list(df.columns)
    secondary_cols = [c for c in all_json_cols if c not in {c0 for c0, *_ in PRIMARY_COLS}]

    records = []
    for _, row in df.iterrows():
        rec = {k: clean_value(row[k]) for k in all_json_cols}
        records.append(rec)

    n_total = len(records)
    n_identified = sum(1 for r in records if r.get("best_identity_source") not in (None, "unidentified"))
    n_robust = sum(1 for r in records if r.get("blocking_robust") is True)

    primary_cols_json = [{"key": k, "label": lbl, "kind": kind} for k, lbl, kind in cols_present]
    secondary_labels_json = {k: SECONDARY_LABELS.get(k, k) for k in secondary_cols}

    return TEMPLATE.format(
        title=html.escape(title),
        n_total=n_total,
        n_identified=n_identified,
        n_robust=n_robust,
        data_json=json.dumps(records),
        primary_cols_json=json.dumps(primary_cols_json),
        secondary_cols_json=json.dumps(secondary_cols),
        secondary_labels_json=json.dumps(secondary_labels_json),
        glyphs_json=json.dumps(GLYPH_BY_SOURCE),
    )


TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  :root {{
    --bg: #ffffff; --surface: #f7f7f8; --border: #e2e2e6;
    --ink: #1a1a1e; --ink-muted: #6b6b74; --accent: #0072B2;
    --row-hover: #f0f4f8; --robust-bg: #eef8f2;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; padding: 24px; background: var(--bg); color: var(--ink);
    font: 14px/1.4 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
  }}
  h1 {{ font-size: 18px; margin: 0 0 4px; }}
  .subtitle {{ color: var(--ink-muted); margin: 0 0 16px; font-size: 13px; }}
  .toolbar {{
    display: flex; flex-wrap: wrap; gap: 10px; align-items: center;
    padding: 12px; background: var(--surface); border: 1px solid var(--border);
    border-radius: 8px; margin-bottom: 12px; position: sticky; top: 0; z-index: 2;
  }}
  .toolbar label {{ font-size: 12px; color: var(--ink-muted); display: flex; align-items: center; gap: 6px; }}
  .toolbar input[type="text"], .toolbar input[type="number"] {{
    font: inherit; padding: 5px 8px; border: 1px solid var(--border); border-radius: 5px; width: 100px;
  }}
  .toolbar input[type="text"] {{ width: 220px; }}
  .chip {{
    font: inherit; font-size: 12px; padding: 5px 10px; border-radius: 999px;
    border: 1px solid var(--border); background: #fff; cursor: pointer; color: var(--ink);
  }}
  .chip[aria-pressed="true"] {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
  .chip .swatch {{ display: inline-block; width: 10px; text-align: center; margin-right: 4px; }}
  #count {{ margin-left: auto; font-size: 12px; color: var(--ink-muted); }}
  table {{ border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }}
  thead th {{
    position: sticky; top: 58px; background: var(--surface); text-align: left;
    padding: 8px 10px; border-bottom: 2px solid var(--border); cursor: pointer;
    white-space: nowrap; user-select: none; font-size: 12px; color: var(--ink-muted);
  }}
  thead th:hover {{ color: var(--ink); }}
  thead th.sorted {{ color: var(--accent); }}
  thead th .arrow {{ font-size: 10px; margin-left: 3px; }}
  tbody td {{ padding: 7px 10px; border-bottom: 1px solid var(--border); white-space: nowrap; }}
  tbody tr.data-row {{ cursor: pointer; }}
  tbody tr.data-row:hover {{ background: var(--row-hover); }}
  tbody tr.data-row.robust {{ background: var(--robust-bg); }}
  tbody tr.data-row.robust:hover {{ background: #e2f2ea; }}
  td.identity {{ text-align: center; font-size: 15px; }}
  td.name-cell {{ max-width: 380px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
  tr.detail-row {{ display: none; }}
  tr.detail-row.open {{ display: table-row; }}
  tr.detail-row td {{ background: var(--surface); white-space: normal; padding: 12px 16px; }}
  .detail-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 6px 20px; }}
  .detail-item dt {{ font-size: 11px; color: var(--ink-muted); margin: 0; }}
  .detail-item dd {{ margin: 0 0 6px; font-size: 13px; word-break: break-word; }}
  .legend {{ display: flex; gap: 14px; flex-wrap: wrap; font-size: 12px; color: var(--ink-muted); margin: 0 0 12px; }}
  .legend span.g {{ margin-right: 3px; }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg: #16161a; --surface: #1e1e23; --border: #303038;
      --ink: #eaeaef; --ink-muted: #9a9aa5; --row-hover: #26262c; --robust-bg: #17251d;
    }}
    .chip {{ background: #222228; }}
    tr.detail-row.robust:hover {{ background: #1c2f24; }}
  }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="subtitle">{n_total} significant features &middot; {n_identified} already identified (library/SIRIUS) &middot; {n_robust} survive plate-block permutation testing</p>
<p class="legend" id="legend"></p>
<div class="toolbar">
  <label>Search <input type="text" id="search" placeholder="identity, formula, adduct..."></label>
  <label>q-value &le; <input type="number" id="qmax" step="0.001" placeholder="0.05"></label>
  <label>|log2FC| &ge; <input type="number" id="fcmin" step="0.1" placeholder="0"></label>
  <label><input type="checkbox" id="robustOnly"> robust only</label>
  <span id="sourceChips"></span>
  <span id="count"></span>
</div>
<table>
  <thead><tr id="headRow"></tr></thead>
  <tbody id="body"></tbody>
</table>

<script>
const DATA = {data_json};
const PRIMARY_COLS = {primary_cols_json};
const SECONDARY_COLS = {secondary_cols_json};
const SECONDARY_LABELS = {secondary_labels_json};
const GLYPHS = {glyphs_json};

let sortKey = "q_used", sortDir = 1;
let openRow = null;

function fmt(kind, v) {{
  if (v === null || v === undefined) return "–";
  if (kind === "sci") return Number(v).toExponential(2);
  if (kind === "bool") return v ? "✓" : "";
  return String(v);
}}
function fmtFloat(v, digits) {{
  if (v === null || v === undefined) return "–";
  return Number(v).toFixed(digits);
}}

function glyphCell(source) {{
  const g = GLYPHS[source] || GLYPHS["unidentified"];
  return `<span title="${{g.label}}" style="color:${{g.color}}">${{g.glyph}}</span>`;
}}

function buildLegend() {{
  const el = document.getElementById("legend");
  el.innerHTML = Object.values(GLYPHS).map(g =>
    `<span><span class="g" style="color:${{g.color}}">${{g.glyph}}</span>${{g.label}}</span>`
  ).join("");
}}

function buildHead() {{
  const tr = document.getElementById("headRow");
  tr.innerHTML = "";
  PRIMARY_COLS.forEach(col => {{
    const th = document.createElement("th");
    th.textContent = col.label;
    th.dataset.key = col.key;
    th.dataset.kind = col.kind;
    const arrow = document.createElement("span");
    arrow.className = "arrow";
    th.appendChild(arrow);
    th.addEventListener("click", () => {{
      if (sortKey === col.key) sortDir *= -1; else {{ sortKey = col.key; sortDir = col.kind === "text" || col.kind === "glyph" ? 1 : -1; }}
      render();
    }});
    tr.appendChild(th);
  }});
  const th = document.createElement("th");
  tr.appendChild(th); // spacer for detail toggle
}}

let sourceFilter = new Set(Object.keys(GLYPHS));
function buildSourceChips() {{
  const el = document.getElementById("sourceChips");
  el.innerHTML = "";
  Object.entries(GLYPHS).forEach(([key, g]) => {{
    const btn = document.createElement("button");
    btn.className = "chip";
    btn.setAttribute("aria-pressed", "true");
    btn.innerHTML = `<span class="swatch" style="color:${{g.color}}">${{g.glyph}}</span>${{g.label}}`;
    btn.addEventListener("click", () => {{
      if (sourceFilter.has(key)) {{ sourceFilter.delete(key); btn.setAttribute("aria-pressed", "false"); }}
      else {{ sourceFilter.add(key); btn.setAttribute("aria-pressed", "true"); }}
      render();
    }});
    el.appendChild(btn);
  }});
}}

function currentFilters() {{
  return {{
    search: document.getElementById("search").value.trim().toLowerCase(),
    qmax: parseFloat(document.getElementById("qmax").value),
    fcmin: parseFloat(document.getElementById("fcmin").value),
    robustOnly: document.getElementById("robustOnly").checked,
  }};
}}

function matches(row, f) {{
  if (!sourceFilter.has(row.best_identity_source || "unidentified")) return false;
  if (f.robustOnly && row.blocking_robust !== true) return false;
  if (!isNaN(f.qmax) && !(row.q_used <= f.qmax)) return false;
  if (!isNaN(f.fcmin) && !(Math.abs(row.log2FC_a_over_b) >= f.fcmin)) return false;
  if (f.search) {{
    const hay = [row["row ID"], row.best_identity, row.library_NAME, row.analog_NAME, row.sirius_structure_name,
                 row.sirius_formula, row.adduct].filter(Boolean).join(" ").toLowerCase();
    if (!hay.includes(f.search)) return false;
  }}
  return true;
}}

function detailPanel(row) {{
  const items = SECONDARY_COLS.map(k => {{
    const v = row[k];
    if (v === null || v === undefined || v === "") return "";
    return `<div class="detail-item"><dt>${{SECONDARY_LABELS[k] || k}}</dt><dd>${{String(v)}}</dd></div>`;
  }}).join("");
  return `<dl class="detail-grid">${{items}}</dl>`;
}}

let debounceTimer;
function scheduleRender() {{
  clearTimeout(debounceTimer);
  debounceTimer = setTimeout(render, 150);
}}

function render() {{
  const f = currentFilters();
  let rows = DATA.filter(r => matches(r, f));

  document.querySelectorAll("#headRow th[data-key]").forEach(th => {{
    th.classList.toggle("sorted", th.dataset.key === sortKey);
    const arrow = th.querySelector(".arrow");
    arrow.textContent = th.dataset.key === sortKey ? (sortDir === 1 ? "▲" : "▼") : "";
  }});

  const kind = (PRIMARY_COLS.find(c => c.key === sortKey) || {{}}).kind;
  rows.sort((a, b) => {{
    let av = a[sortKey], bv = b[sortKey];
    if (kind === "text" || kind === "glyph") {{
      av = (av || "").toString(); bv = (bv || "").toString();
      return sortDir * av.localeCompare(bv);
    }}
    av = av === null || av === undefined ? -Infinity : Number(av);
    bv = bv === null || bv === undefined ? -Infinity : Number(bv);
    return sortDir * (av - bv);
  }});

  document.getElementById("count").textContent = `${{rows.length}} / ${{DATA.length}} rows`;

  const tbody = document.getElementById("body");
  tbody.innerHTML = "";
  const frag = document.createDocumentFragment();

  rows.forEach((row, i) => {{
    const tr = document.createElement("tr");
    tr.className = "data-row" + (row.blocking_robust ? " robust" : "");
    PRIMARY_COLS.forEach(col => {{
      const td = document.createElement("td");
      if (col.key === "best_identity_source") {{
        td.className = "identity";
        td.innerHTML = glyphCell(row.best_identity_source);
      }} else if (col.key === "best_identity") {{
        td.className = "name-cell";
        td.textContent = row.best_identity || "–";
      }} else if (col.kind === "float") {{
        td.textContent = fmtFloat(row[col.key], col.key.includes("retention") ? 2 : 4);
      }} else {{
        td.textContent = fmt(col.kind, row[col.key]);
      }}
      tr.appendChild(td);
    }});
    const toggleTd = document.createElement("td");
    toggleTd.textContent = "›";
    toggleTd.style.color = "var(--ink-muted)";
    tr.appendChild(toggleTd);

    const detailTr = document.createElement("tr");
    detailTr.className = "detail-row";
    const detailTd = document.createElement("td");
    detailTd.colSpan = PRIMARY_COLS.length + 1;
    detailTd.innerHTML = detailPanel(row);
    detailTr.appendChild(detailTd);

    tr.addEventListener("click", () => {{
      const wasOpen = detailTr.classList.contains("open");
      document.querySelectorAll("tr.detail-row.open").forEach(el => el.classList.remove("open"));
      if (!wasOpen) detailTr.classList.add("open");
    }});

    frag.appendChild(tr);
    frag.appendChild(detailTr);
  }});
  tbody.appendChild(frag);
}}

buildLegend();
buildHead();
buildSourceChips();
["search"].forEach(id => document.getElementById(id).addEventListener("input", scheduleRender));
["qmax", "fcmin"].forEach(id => document.getElementById(id).addEventListener("input", scheduleRender));
document.getElementById("robustOnly").addEventListener("change", render);
render();
</script>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsvs", nargs="*", type=Path)
    args = ap.parse_args()

    if args.tsvs:
        tsv_paths = args.tsvs
    else:
        tsv_paths = sorted(DIFF_ROOT.glob("*/compound_summary.tsv"))
        if MASTER_TSV.exists():
            tsv_paths.append(MASTER_TSV)

    for tsv_path in tsv_paths:
        if not tsv_path.exists():
            print(f"skip: {tsv_path} not found", file=sys.stderr)
            continue
        df = pd.read_csv(tsv_path, sep="\t")
        if df.empty:
            print(f"skip: {tsv_path} has no rows", file=sys.stderr)
            continue
        title = tsv_path.parent.name if tsv_path.name == "compound_summary.tsv" else tsv_path.stem
        out_path = tsv_path.with_suffix(".html")
        out_path.write_text(build_html(df, title))
        print(f"{tsv_path} ({len(df)} rows) -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
