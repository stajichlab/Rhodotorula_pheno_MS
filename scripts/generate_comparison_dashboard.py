#!/usr/bin/env python3
"""
Generate a per-comparison dashboard.html next to each
compound_summary.tsv under analysis/differential_features/.

Each dashboard embeds the volcano and top-features PNG plots,
links to their PDFs, and links to the compound_summary.html
interactive table -- bringing all the pairwise comparison
outputs onto a single page.

A README.md index is also written to analysis/differential_features/
linking to every dashboard and to the all_significant_features_summary
rollup, for navigation from a folder.

Usage:
    python3 scripts/generate_comparison_dashboard.py
"""
from __future__ import annotations

import base64
import html
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DIFF_ROOT = REPO / "analysis" / "differential_features"
MASTER_TSV = DIFF_ROOT / "all_significant_features_summary.tsv"


def _label(name: str) -> str:
    """cell_diobovata_vs_mucilaginosa -> 'R. diobovata vs R. mucilaginosa'"""
    if name.startswith("cell_"):
        frac = "Cell pellet"
        pair = name[len("cell_"):]
    elif name.startswith("supernatant_"):
        frac = "Supernatant"
        pair = name[len("supernatant_"):]
    else:
        return name.replace("_", " ")

    def sp(s: str) -> str:
        if s == "sp_clade_I":
            return "sp. clade I"
        return s

    parts = pair.split("_vs_")
    if len(parts) != 2:
        return name.replace("_", " ")
    a = parts[0].split("_")
    b = parts[1].split("_")
    a_name = " ".join(sp(x) for x in a)
    b_name = " ".join(sp(x) for x in b)
    return f"{frac}: R. {a_name} vs R. {b_name}"


def _img_data(path: Path) -> str | None:
    if not path.exists():
        return None
    return base64.b64encode(path.read_bytes()).decode("ascii")


DASHBOARD_CSS = """
  :root {
    --bg: #ffffff; --surface: #f7f7f8; --border: #e2e2e6;
    --ink: #1a1a1e; --ink-muted: #6b6b74; --accent: #0072B2;
  }
  @media (prefers-color-scheme: dark) {
    :root {
      --bg: #16161a; --surface: #1e1e23; --border: #303038;
      --ink: #eaeaef; --ink-muted: #9a9aa5; --accent: #56B4E9;
    }
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; padding: 24px; background: var(--bg); color: var(--ink);
    font: 14px/1.5 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
  }
  .header { display: flex; align-items: center; gap: 12px; margin-bottom: 16px; }
  .header h1 { font-size: 18px; margin: 0; }
  .back { font-size: 13px; color: var(--accent); text-decoration: none; }
  .stats {
    display: flex; gap: 20px; flex-wrap: wrap; margin-bottom: 20px;
    padding: 12px 16px; background: var(--surface); border: 1px solid var(--border);
    border-radius: 8px;
  }
  .stat .num { font-size: 20px; font-weight: 600; color: var(--accent); }
  .stat .label { font-size: 11px; color: var(--ink-muted); text-transform: uppercase; letter-spacing: 0.5px; }
  .plots { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px; }
  @media (max-width: 900px) { .plots { grid-template-columns: 1fr; } }
  .plot-card { background: var(--surface); border: 1px solid var(--border); border-radius: 8px; padding: 16px; }
  .plot-card h2 { font-size: 13px; color: var(--ink-muted); text-transform: uppercase; letter-spacing: 0.5px; margin: 0 0 8px; }
  .plot-card img { width: 100%; border-radius: 4px; }
  .plot-links { margin-top: 8px; font-size: 12px; }
  .plot-links a { color: var(--accent); text-decoration: none; margin-right: 12px; }
  .table-link {
    display: inline-block; padding: 10px 20px; background: var(--accent); color: #fff;
    border-radius: 6px; text-decoration: none; font-size: 14px; font-weight: 500;
  }
  .table-link:hover { opacity: 0.9; }
  .section { margin-bottom: 20px; }
  .section h2 { font-size: 13px; color: var(--ink-muted); text-transform: uppercase; letter-spacing: 0.5px; }
"""


def build_dashboard(diff_dir: Path) -> str | None:
    """Generate dashboard.html for a comparison directory. Returns None if no data."""
    tsv_path = diff_dir / "compound_summary.tsv"
    if not tsv_path.exists():
        return None

    df = pd.read_csv(tsv_path, sep="\t")
    if df.empty:
        return None

    name = diff_dir.name
    label = _label(name)
    n_total = len(df)
    n_identified = int(df["best_identity"].notna().sum()) if "best_identity" in df.columns else 0
    n_robust = int(df["blocking_robust"].sum()) if "blocking_robust" in df.columns else 0

    volcano_png = _img_data(diff_dir / "volcano.png")
    top_png = _img_data(diff_dir / "top_features.png")

    volcano_html = ""
    if volcano_png:
        volcano_html = f'''
    <div class="plot-card">
      <h2>Volcano plot</h2>
      <img src="data:image/png;base64,{volcano_png}" alt="Volcano plot">
      <div class="plot-links">
        <a href="volcano.pdf">PDF</a>
        <a href="volcano.png">PNG</a>
      </div>
    </div>'''
    else:
        volcano_html = '<div class="plot-card"><h2>Volcano plot</h2><p>Not available</p></div>'

    top_html = ""
    if top_png:
        top_html = f'''
    <div class="plot-card">
      <h2>Top features</h2>
      <img src="data:image/png;base64,{top_png}" alt="Top features plot">
      <div class="plot-links">
        <a href="top_features.pdf">PDF</a>
        <a href="top_features.png">PNG</a>
      </div>
    </div>'''
    else:
        top_html = '<div class="plot-card"><h2>Top features</h2><p>Not available</p></div>'

    has_table_html = (diff_dir / "compound_summary.html").exists()

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(label)}</title>
<style>{DASHBOARD_CSS}</style>
</head>
<body>
<div class="header">
  <h1>{html.escape(label)}</h1>
</div>

<div class="stats">
  <div class="stat"><div class="num">{n_total:,}</div><div class="label">Significant features</div></div>
  <div class="stat"><div class="num">{n_identified:,}</div><div class="label">Identified</div></div>
  <div class="stat"><div class="num">{n_robust:,}</div><div class="label">Plate-block robust</div></div>
</div>

<div class="plots">{volcano_html}
{top_html}
</div>

<div class="section">
  <h2>Compound summary table</h2>
  {'<a class="table-link" href="compound_summary.html">Open interactive table &rsaquo;</a>' if has_table_html else '<p>Table not generated</p>'}
</div>

</body>
</html>"""


INDEX_CSS = """
  :root {
    --bg: #ffffff; --surface: #f7f7f8; --border: #e2e2e6;
    --ink: #1a1a1e; --ink-muted: #6b6b74; --accent: #0072B2;
  }
  @media (prefers-color-scheme: dark) {
    :root {
      --bg: #16161a; --surface: #1e1e23; --border: #303038;
      --ink: #eaeaef; --ink-muted: #9a9aa5; --accent: #56B4E9;
    }
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; padding: 32px; background: var(--bg); color: var(--ink);
    font: 14px/1.5 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
  }
  h1 { font-size: 22px; margin: 0 0 4px; }
  h2 { font-size: 14px; color: var(--ink-muted); text-transform: uppercase; letter-spacing: 0.5px; margin: 24px 0 10px; }
  .subtitle { color: var(--ink-muted); margin: 0 0 20px; font-size: 13px; }
  .stats { display: flex; gap: 24px; flex-wrap: wrap; margin-bottom: 8px; }
  .stat .num { font-size: 22px; font-weight: 600; color: var(--accent); }
  .stat .label { font-size: 11px; color: var(--ink-muted); text-transform: uppercase; letter-spacing: 0.5px; }
  table { border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }
  thead th { background: var(--surface); text-align: left; padding: 8px 12px; border-bottom: 2px solid var(--border); font-size: 12px; color: var(--ink-muted); text-transform: uppercase; letter-spacing: 0.5px; }
  tbody td { padding: 8px 12px; border-bottom: 1px solid var(--border); }
  tbody tr:hover { background: var(--surface); }
  a { color: var(--accent); text-decoration: none; }
  a:hover { text-decoration: underline; }
  .frac { font-size: 11px; color: var(--ink-muted); text-transform: uppercase; }
  code { background: var(--surface); padding: 1px 4px; border-radius: 3px; font-size: 12px; }
"""


def _row(name: str, n_total: int, n_identified: int) -> str:
    return (
        f'<tr><td>{html.escape(_label(name))}</td>'
        f'<td>{n_total:,}</td><td>{n_identified:,}</td>'
        f'<td><a href="./{name}/dashboard.html">dashboard</a> &middot; '
        f'<a href="./{name}/compound_summary.html">table</a></td></tr>'
    )


def build_readme(dashboards: list[tuple[str, int, int, int]]) -> str:
    """Generate README.md index."""
    lines = [
        "# Differential Features Analysis",
        "",
        "Pairwise species comparisons of MS2 features, annotated with",
        "GNPS library matches and SIRIUS/CANOPUS predictions.",
        "",
        f"- **{sum(d[1] for d in dashboards):,}** significant features (FDR < 5%)",
        f"- **{len(dashboards)}** comparisons with significant features (out of 110 total)",
        f"- **790** SIRIUS annotations (663 with structure predictions)",
        "",
        "## Cross-comparison rollup",
        "",
        "| View | Features | Link |",
        "|------|----------|------|",
        f"| All significant features (every comparison concatenated) | {sum(d[1] for d in dashboards):,} | [open](all_significant_features_summary.html) |",
        "",
    ]

    cell_dirs = [d for d in dashboards if d[0].startswith("cell_")]
    sup_dirs = [d for d in dashboards if d[0].startswith("supernatant_")]

    if cell_dirs:
        lines += ["## Cell pellet comparisons", "", "| Comparison | Features | Identified | Dashboard |", "|------------|----------|------------|-----------|"]
        for name, n_total, n_identified, n_robust in sorted(cell_dirs, key=lambda x: -x[1]):
            lines.append(f"| {_label(name)} | {n_total:,} | {n_identified:,} | [dashboard](./{name}/dashboard.html) |")
        lines.append("")

    if sup_dirs:
        lines += ["## Supernatant comparisons", "", "| Comparison | Features | Identified | Dashboard |", "|------------|----------|------------|-----------|"]
        for name, n_total, n_identified, n_robust in sorted(sup_dirs, key=lambda x: -x[1]):
            lines.append(f"| {_label(name)} | {n_total:,} | {n_identified:,} | [dashboard](./{name}/dashboard.html) |")
        lines.append("")

    lines += [
        "## Individual comparison tables",
        "",
        "Each comparison directory also contains:",
        "- `compound_summary.html` -- sortable/filterable table of significant features",
        "- `volcano.pdf` / `volcano.png` -- volcano plot",
        "- `top_features.pdf` / `top_features.png` -- top features plot",
        "- `differential_features.csv.gz` -- full differential features table",
        "",
    ]

    return "\n".join(lines)


def build_readme_html(dashboards: list[tuple[str, int, int, int]]) -> str:
    """Generate README.html — same content as README.md but rendered as
    a self-contained HTML page, since .nojekyll disables Jekyll markdown
    processing on GitHub Pages."""
    total = sum(d[1] for d in dashboards)
    cell_dirs = sorted([d for d in dashboards if d[0].startswith("cell_")], key=lambda x: -x[1])
    sup_dirs = sorted([d for d in dashboards if d[0].startswith("supernatant_")], key=lambda x: -x[1])

    rows_html = []
    if cell_dirs:
        rows_html.append("<h2>Cell pellet comparisons</h2>")
        rows_html.append("<table><thead><tr><th>Comparison</th><th>Features</th><th>Identified</th><th>Links</th></tr></thead><tbody>")
        for name, n_total, n_identified, _ in cell_dirs:
            rows_html.append(_row(name, n_total, n_identified))
        rows_html.append("</tbody></table>")
    if sup_dirs:
        rows_html.append("<h2>Supernatant comparisons</h2>")
        rows_html.append("<table><thead><tr><th>Comparison</th><th>Features</th><th>Identified</th><th>Links</th></tr></thead><tbody>")
        for name, n_total, n_identified, _ in sup_dirs:
            rows_html.append(_row(name, n_total, n_identified))
        rows_html.append("</tbody></table>")

    body_rows = "\n".join(rows_html)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Differential Features Analysis</title>
<style>{INDEX_CSS}</style>
</head>
<body>
<h1>Differential Features Analysis</h1>
<p class="subtitle">Pairwise species comparisons of MS2 features, annotated with
GNPS library matches and SIRIUS/CANOPUS predictions.</p>

<div class="stats">
  <div class="stat"><div class="num">{total:,}</div><div class="label">Significant features</div></div>
  <div class="stat"><div class="num">{len(dashboards)}</div><div class="label">Comparisons with hits</div></div>
  <div class="stat"><div class="num">110</div><div class="label">Total comparisons</div></div>
  <div class="stat"><div class="num">790</div><div class="label">SIRIUS annotations</div></div>
</div>

<h2>Cross-comparison rollup</h2>
<table>
<thead><tr><th>View</th><th>Features</th><th>Link</th></tr></thead>
<tbody>
<tr><td>All significant features (every comparison concatenated)</td><td>{total:,}</td>
<td><a href="all_significant_features_summary.html">open</a></td></tr>
</tbody>
</table>

{body_rows}

<h2>Individual comparison files</h2>
<p>Each comparison directory also contains:</p>
<ul>
<li><code>compound_summary.html</code> — sortable/filterable table of significant features</li>
<li><code>volcano.pdf</code> / <code>volcano.png</code> — volcano plot</li>
<li><code>top_features.pdf</code> / <code>top_features.png</code> — top features plot</li>
<li><code>differential_features.csv.gz</code> — full differential features table</li>
</ul>

</body>
</html>"""


def main():
    diff_dirs = sorted(p.parent for p in DIFF_ROOT.glob("*/compound_summary.tsv"))

    dashboards = []
    for diff_dir in diff_dirs:
        tsv_path = diff_dir / "compound_summary.tsv"
        df = pd.read_csv(tsv_path, sep="\t")
        if df.empty:
            continue

        name = diff_dir.name
        n_total = len(df)
        n_identified = int(df["best_identity"].notna().sum()) if "best_identity" in df.columns else 0
        n_robust = int(df["blocking_robust"].sum()) if "blocking_robust" in df.columns else 0

        html_content = build_dashboard(diff_dir)
        if html_content:
            out_path = diff_dir / "dashboard.html"
            out_path.write_text(html_content)
            print(f"{name}: dashboard -> {out_path}", file=sys.stderr)

        dashboards.append((name, n_total, n_identified, n_robust))

    readme_md = build_readme(dashboards)
    readme_html = build_readme_html(dashboards)
    md_path = DIFF_ROOT / "README.md"
    html_path = DIFF_ROOT / "README.html"
    md_path.write_text(readme_md)
    html_path.write_text(readme_html)
    print(f"README index -> {md_path} + {html_path} ({len(dashboards)} comparisons)", file=sys.stderr)


if __name__ == "__main__":
    main()
