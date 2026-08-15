#!/usr/bin/env python3
"""
Do the two color-phenotype tables agree on the strains they share?

Compares per-strain L*, a*, b*, C* between:
  - YPD2 (data/metadata/EXFAB_UCR-005/YPD2_phenotypic.20260702.fixed.csv.gz,
    318 rows, 15 strains missing a Species call as of 2026-08)
  - control_late (/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/
    analysis/control_late_timepoint_phenotype/results/
    phenotype_control_late_timepoint.csv, 303 rows, 1 missing Species,
    corrected strain codes/species for several strains -- see
    analysis/examine_phenotype_calling/strain_membership_diff.csv from the
    membership-diff step of this same investigation)

Both are read straight from source (not from build_strain_phenotype_table.py's
output, which only holds whichever --source was run last) and independently
collapsed to one row per strain (mean across replicate rows), so this script
is a standalone check, not a dependency of the Phase 1 pipeline.

Matching: strain codes are matched after normalizing (strip underscore/dash/
space, uppercase), same convention as analysis/copper/scripts/common.py and
scripts/sync_bfd_symlinks.py. A handful of strains only match after this
normalization because control_late corrected typos in the strain code (e.g.
YPD2's TFCN_223D-8 == control_late's TFCN_223A-8 by *Strain ID*, not by
normalized code -- those do NOT get merged here, since matching by numeric
Strain ID vs. normalized code text are different, non-interchangeable join
keys; see the membership-diff writeup for that nuance). This script only
compares strains whose *codes* match, which is a same-strain guarantee.

If two independent phenotyping runs of the same strain agree well
(high r, points hugging y=x), that strengthens confidence in using either
table. Points that are strong outliers from y=x are worth a manual look --
either a strain-code collision (two different physical isolates matched to
the same normalized code) or a real batch/plate effect between runs.

Usage:
    python3 analysis/examine_phenotype_calling/scripts/compare_phenotype_sources.py
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
YPD2_PATH = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "YPD2_phenotypic.20260702.fixed.csv.gz"
CONTROL_LATE_PATH = Path(
    "/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/analysis/"
    "control_late_timepoint_phenotype/results/phenotype_control_late_timepoint.csv"
)
OUT_DIR = REPO / "analysis" / "examine_phenotype_calling"
PLOT_PATH = OUT_DIR / "ypd2_vs_control_late_scatter.png"
PAIRED_CSV = OUT_DIR / "ypd2_vs_control_late_paired_values.csv"
SUMMARY_TXT = OUT_DIR / "ypd2_vs_control_late_summary.txt"

TRAITS = ["L*", "a*", "b*", "C*"]


def norm(s: str) -> str:
    return re.sub(r"[_\-\s]", "", str(s)).upper()


def load_ypd2() -> pd.DataFrame:
    df = pd.read_csv(YPD2_PATH)
    lab_cols = {"L*": "Median_ColorLab_L*Mean", "a*": "Median_ColorLab_a*Mean", "b*": "Median_ColorLab_b*Mean"}
    df = df.dropna(subset=list(lab_cols.values()))
    agg = {v: "mean" for v in lab_cols.values()}
    agg["Species"] = "first"
    g = df.groupby("Strain", as_index=False).agg(agg).rename(
        columns={"Strain": "strain_id", "Species": "species", **{v: k for k, v in lab_cols.items()}}
    )
    g["C*"] = np.sqrt(g["a*"] ** 2 + g["b*"] ** 2)
    g["norm_id"] = g["strain_id"].map(norm)
    return g


def load_control_late() -> pd.DataFrame:
    df = pd.read_csv(CONTROL_LATE_PATH)
    lab_cols = {"L*": "l_median", "a*": "a_median", "b*": "b_median"}
    df = df.dropna(subset=list(lab_cols.values()))
    agg = {v: "mean" for v in lab_cols.values()}
    agg["species"] = "first"
    g = df.groupby("strain_code", as_index=False).agg(agg).rename(
        columns={"strain_code": "strain_id", **{v: k for k, v in lab_cols.items()}}
    )
    g["C*"] = np.sqrt(g["a*"] ** 2 + g["b*"] ** 2)
    g["norm_id"] = g["strain_id"].map(norm)
    return g


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ypd2 = load_ypd2()
    ctrl = load_control_late()

    paired = ypd2.merge(ctrl, on="norm_id", suffixes=("_ypd2", "_control_late"))
    # prefer control_late's (more complete/corrected) species label; fall back to YPD2's
    paired["species"] = paired["species_control_late"].fillna(paired["species_ypd2"])
    n_paired = len(paired)
    print(f"{n_paired} strains matched by normalized code between YPD2 (n={len(ypd2)}) and control_late (n={len(ctrl)})")

    order, markers = canonical_species_order()
    colors = full_species_color_map(order)
    paired["species_label"] = paired["species"].fillna("Unknown")

    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    stats = {}
    for ax, trait in zip(axes.flat, TRAITS):
        x = paired[f"{trait}_ypd2"].to_numpy(dtype=float)
        y = paired[f"{trait}_control_late"].to_numpy(dtype=float)
        r, r_p = pearsonr(x, y)
        rho, rho_p = spearmanr(x, y)
        stats[trait] = dict(pearson_r=r, pearson_p=r_p, spearman_rho=rho, spearman_p=rho_p, n=len(x))

        for sp in order:
            sub = paired[paired["species_label"] == sp]
            if sub.empty:
                continue
            ax.scatter(
                sub[f"{trait}_ypd2"], sub[f"{trait}_control_late"],
                marker=markers[sp], s=28, color=colors[sp], edgecolor="white", linewidth=0.3, zorder=2,
            )
        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        pad = 0.03 * (hi - lo)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#999999", linewidth=1, linestyle="--", zorder=1)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"{trait} (YPD2)")
        ax.set_ylabel(f"{trait} (control_late)")
        ax.set_title(f"{trait}  (Pearson r={r:.3f}, Spearman ρ={rho:.3f}, n={len(x)})", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)

    handles = [
        plt.Line2D(
            [0], [0], marker=markers[sp], linestyle="none", markerfacecolor=colors[sp],
            markeredgecolor="white", markersize=7, label=sp,
        )
        for sp in order
        if (paired["species_label"] == sp).any()
    ]
    fig.legend(
        handles=handles, title="Species", frameon=False, fontsize=7,
        loc="center left", bbox_to_anchor=(1.0, 0.5), labelspacing=0.5,
    )
    fig.suptitle("YPD2 vs. control_late per-strain color phenotype (dashed line = y=x)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 0.86, 0.96])
    savefig_multi(fig, PLOT_PATH)
    plt.close(fig)

    out_cols = ["strain_id_ypd2", "strain_id_control_late", "species"] + [
        f"{t}_{src}" for t in TRAITS for src in ("ypd2", "control_late")
    ]
    paired.rename(columns={"strain_id_ypd2": "strain_id_ypd2", "strain_id_control_late": "strain_id_control_late"})[
        out_cols
    ].to_csv(PAIRED_CSV, index=False)

    # flag the largest per-trait disagreements as a quick outlier list
    outlier_lines = []
    for trait in TRAITS:
        d = (paired[f"{trait}_control_late"] - paired[f"{trait}_ypd2"]).abs()
        top = d.sort_values(ascending=False).head(5)
        outlier_lines.append(f"\nLargest |control_late - YPD2| for {trait}:")
        for idx in top.index:
            row = paired.loc[idx]
            outlier_lines.append(
                f"  {row['strain_id_ypd2']} ({row['species']}): "
                f"YPD2={row[f'{trait}_ypd2']:.2f}  control_late={row[f'{trait}_control_late']:.2f}  "
                f"diff={d[idx]:.2f}"
            )

    with SUMMARY_TXT.open("w") as fh:
        fh.write("YPD2 vs. control_late color-phenotype agreement\n")
        fh.write("=" * 50 + "\n")
        fh.write(f"YPD2 strains (usable L*/a*/b*): {len(ypd2)}\n")
        fh.write(f"control_late strains (usable L*/a*/b*): {len(ctrl)}\n")
        fh.write(f"Matched by normalized strain code: {n_paired}\n\n")
        for trait, s in stats.items():
            fh.write(
                f"{trait}: Pearson r={s['pearson_r']:.4f} (p={s['pearson_p']:.2e}), "
                f"Spearman rho={s['spearman_rho']:.4f} (p={s['spearman_p']:.2e}), n={s['n']}\n"
            )
        fh.write("\n" + "\n".join(outlier_lines) + "\n")

    print(f"Wrote {PLOT_PATH} (+ .pdf)")
    print(f"Wrote {PAIRED_CSV}")
    print(f"Wrote {SUMMARY_TXT}")
    for trait, s in stats.items():
        print(f"  {trait}: r={s['pearson_r']:.3f}  rho={s['spearman_rho']:.3f}")


if __name__ == "__main__":
    main()
