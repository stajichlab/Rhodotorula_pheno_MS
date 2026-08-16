#!/usr/bin/env python3
"""
Extends compare_phenotype_sources.py (YPD2 vs. the single old
"control_late" file) now that the upstream Rhodotorula_phenotypes project
split that file into three explicit imaging-timepoint windows (2026-08-15):
control_70_80 (passes 75/78h), control_80_90 (87/90h), control_90_110
(105/108h, latest/most-plateaued -- this is what "control_late" used to
mean). See .../control_late_timepoint_phenotype/CONTROL_LATE_PHENOTYPE.md
in the sibling project for the full method.

Two questions, two figures:
  1. timepoint_progression_scatter.png -- do the three control-media windows
     agree with EACH OTHER? (70_80 vs 80_90, 80_90 vs 90_110, 70_80 vs
     90_110). If growth has plateaued by 70-80h, these should all sit near
     y=x; systematic drift would mean the phenotype is still developing and
     which window you pick matters for downstream analysis.
  2. ypd2_vs_control_windows_scatter.png -- does EACH window agree with the
     original YPD2 table? (YPD2 vs 70_80, vs 80_90, vs 90_110).
Both cover color (L*, a*, b*, C*) and colony size (area, YPD2's
Median_Shape_Area vs. the control tables' area_median -- same px units).

Matching: strain codes normalized (strip underscore/dash/space, uppercase)
and merged, same convention as compare_phenotype_sources.py and
scripts/sync_bfd_symlinks.py -- comparisons only use strains whose codes
match across all sources involved in that panel (a real per-strain repeat
measurement, not an approximate join).

Usage:
    python3 analysis/examine_phenotype_calling/scripts/compare_timepoint_sources.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "scripts"))
from pcoa_ms_features import canonical_species_order, full_species_color_map, savefig_multi

REPO = Path(__file__).resolve().parent.parent.parent.parent
CONTROL_DIR = Path(
    "/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/analysis/"
    "control_late_timepoint_phenotype/results"
)
OUT_DIR = REPO / "analysis" / "examine_phenotype_calling"

SOURCES = {
    "ypd2": dict(
        path=REPO / "data" / "metadata" / "EXFAB_UCR-005" / "YPD2_phenotypic.20260702.fixed.csv.gz",
        strain_col="Strain",
        species_col="Species",
        cols={"L*": "Median_ColorLab_L*Mean", "a*": "Median_ColorLab_a*Mean", "b*": "Median_ColorLab_b*Mean", "area": "Median_Shape_Area"},
    ),
    "control_70_80": dict(
        path=CONTROL_DIR / "phenotype_control_timepoint_70_80.csv",
        strain_col="strain_code", species_col="species",
        cols={"L*": "l_median", "a*": "a_median", "b*": "b_median", "area": "area_median"},
    ),
    "control_80_90": dict(
        path=CONTROL_DIR / "phenotype_control_timepoint_80_90.csv",
        strain_col="strain_code", species_col="species",
        cols={"L*": "l_median", "a*": "a_median", "b*": "b_median", "area": "area_median"},
    ),
    "control_90_110": dict(
        path=CONTROL_DIR / "phenotype_control_timepoint_90_110.csv",
        strain_col="strain_code", species_col="species",
        cols={"L*": "l_median", "a*": "a_median", "b*": "b_median", "area": "area_median"},
    ),
}
TRAITS = ["L*", "a*", "b*", "C*", "area"]
LABELS = {"ypd2": "YPD2", "control_70_80": "control 70-80h", "control_80_90": "control 80-90h", "control_90_110": "control 90-110h"}


def norm(s: str) -> str:
    return re.sub(r"[_\-\s]", "", str(s)).upper()


def load_source(name: str) -> pd.DataFrame:
    src = SOURCES[name]
    df = pd.read_csv(src["path"])
    value_cols = list(src["cols"].values())
    df = df.dropna(subset=value_cols)
    agg = {v: "mean" for v in value_cols}
    agg[src["species_col"]] = "first"
    g = df.groupby(src["strain_col"], as_index=False).agg(agg).rename(
        columns={src["strain_col"]: "strain_id", src["species_col"]: "species", **{v: k for k, v in src["cols"].items()}}
    )
    g["C*"] = np.sqrt(g["a*"] ** 2 + g["b*"] ** 2)
    g["norm_id"] = g["strain_id"].map(norm)
    return g


def scatter_grid(pairs: list[tuple[str, str]], out_path: Path, title: str, all_data: dict) -> list[str]:
    order, markers = canonical_species_order()
    colors = full_species_color_map(order)

    fig, axes = plt.subplots(len(TRAITS), len(pairs), figsize=(5.2 * len(pairs), 4.3 * len(TRAITS)))
    if len(pairs) == 1:
        axes = axes.reshape(-1, 1)
    summary_lines = []

    for col, (src_x, src_y) in enumerate(pairs):
        paired = all_data[src_x].merge(all_data[src_y], on="norm_id", suffixes=("_x", "_y"))
        paired["species"] = paired.get("species_y", paired.get("species")).fillna(paired.get("species_x"))
        paired["species_label"] = paired["species"].fillna("Unknown")
        summary_lines.append(f"\n{LABELS[src_x]} vs {LABELS[src_y]} (n={len(paired)} matched strains):")

        for row, trait in enumerate(TRAITS):
            ax = axes[row, col]
            x = paired[f"{trait}_x"].to_numpy(dtype=float)
            y = paired[f"{trait}_y"].to_numpy(dtype=float)
            r, _ = pearsonr(x, y)
            rho, _ = spearmanr(x, y)
            summary_lines.append(f"  {trait}: Pearson r={r:.4f}  Spearman rho={rho:.4f}  n={len(x)}")

            for sp in order:
                sub = paired[paired["species_label"] == sp]
                if sub.empty:
                    continue
                ax.scatter(
                    sub[f"{trait}_x"], sub[f"{trait}_y"],
                    marker=markers[sp], s=22, color=colors[sp], edgecolor="white", linewidth=0.3, zorder=2,
                )
            lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
            pad = 0.03 * (hi - lo) if hi > lo else 1.0
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#999999", linewidth=1, linestyle="--", zorder=1)
            ax.set_xlim(lo - pad, hi + pad)
            ax.set_ylim(lo - pad, hi + pad)
            ax.set_aspect("equal", adjustable="box")
            if row == 0:
                ax.set_title(f"{LABELS[src_x]}  vs.  {LABELS[src_y]}", fontsize=10)
            ax.set_xlabel(f"{trait} ({LABELS[src_x]})", fontsize=8)
            ax.set_ylabel(f"{trait} ({LABELS[src_y]})", fontsize=8)
            ax.text(
                0.03, 0.97, f"r={r:.3f}\nρ={rho:.3f}\nn={len(x)}",
                transform=ax.transAxes, fontsize=7.5, va="top", ha="left",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="none"),
            )
            ax.spines[["top", "right"]].set_visible(False)

    handles = [
        plt.Line2D([0], [0], marker=markers[sp], linestyle="none", markerfacecolor=colors[sp],
                   markeredgecolor="white", markersize=7, label=sp)
        for sp in order
    ]
    fig.legend(handles=handles, title="Species", frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 0.90, 0.97])
    savefig_multi(fig, out_path)
    plt.close(fig)
    return summary_lines


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_data = {name: load_source(name) for name in SOURCES}
    for name, df in all_data.items():
        print(f"{name}: {len(df)} strains with usable color+area")

    progression_pairs = [
        ("control_70_80", "control_80_90"),
        ("control_80_90", "control_90_110"),
        ("control_70_80", "control_90_110"),
    ]
    ypd2_pairs = [("ypd2", "control_70_80"), ("ypd2", "control_80_90"), ("ypd2", "control_90_110")]

    summary = []
    summary.append("=== Timepoint progression: do the 3 control-media windows agree with each other? ===")
    summary += scatter_grid(
        progression_pairs, OUT_DIR / "timepoint_progression_scatter.png",
        "Timepoint progression: control-media color & colony size across imaging windows (dashed = y=x)",
        all_data,
    )
    summary.append("\n\n=== YPD2 agreement: does each control-media window agree with the original YPD2 table? ===")
    summary += scatter_grid(
        ypd2_pairs, OUT_DIR / "ypd2_vs_control_windows_scatter.png",
        "YPD2 vs. control-media timepoint windows: color & colony size (dashed = y=x)",
        all_data,
    )

    summary_path = OUT_DIR / "timepoint_comparison_summary.txt"
    with summary_path.open("w") as fh:
        fh.write("Phenotype-source correlation summary (color: L*/a*/b*/C*, colony size: area)\n")
        fh.write("=" * 78 + "\n")
        fh.write("\n".join(summary) + "\n")

    print(f"\nWrote {OUT_DIR / 'timepoint_progression_scatter.png'} (+ .pdf)")
    print(f"Wrote {OUT_DIR / 'ypd2_vs_control_windows_scatter.png'} (+ .pdf)")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
