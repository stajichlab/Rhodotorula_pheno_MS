#!/usr/bin/env python3
"""
AHL autoinducer / colony morphology hypothesis: figures showing how much
the colony-texture morphology proxy (strain_texture_table.csv, from
salinity_texture_baseline -- see AHL_AUTOINDUCER_SEARCH.md for full
caveats: this is a texture PROXY, not a validated smooth/rough call)
varies across strains and species.

Figures (PNG for embedding, PDF for vector close-up), written to
analysis/ahl_autoinducer_search/figures/:

  fig1_smoothness_by_species   boxplot + strip of smoothness_z per species,
                                sorted by median, n>=3 species only (species
                                with n<3 plotted separately as points -- a
                                box/quartile is not meaningful below n=3)
  fig2_smoothness_distribution histogram of smoothness_z across all 298
                                strains, with median/IQR annotated
  fig3_raw_texture_components  2x2 small multiples of the 4 raw Haralick
                                metrics feeding smoothness_z (Contrast,
                                Entropy, AngularSecondMoment,
                                InverseDifferenceMoment) by species, so the
                                composite's variance isn't taken on faith

Design choices (dataviz skill): single consistent hue throughout (species
identity is carried by axis position/labels, not color, so no categorical
palette is needed here); recessive gridlines; direct annotation instead of
a legend where there is only one series per panel.

Usage:
    python3 analysis/scripts/ahl_morphology_variance_figures.py
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
TEXTURE_TABLE = REPO / "analysis" / "ahl_autoinducer_search" / "strain_texture_table.csv"
FIG_DIR = REPO / "analysis" / "ahl_autoinducer_search" / "figures"

HUE = "#4393c3"        # single consistent hue, matches network_components figure palette
HUE_DARK = "#2166ac"
POINT_GRAY = "#555555"
GRID_GRAY = "#e0e0e0"
MIN_N_FOR_BOX = 3

TEXTURE_METRICS = [
    ("AngularSecondMoment_median", "Angular Second Moment (energy/uniformity; higher = smoother-reading)"),
    ("Contrast_median", "Contrast (local intensity variation; higher = rougher-reading)"),
    ("Correlation_median", "Correlation (linear dependency between neighbor pixels)"),
    ("HaralickVariance_median", "Haralick Variance (gray-level pair dispersion)"),
    ("InverseDifferenceMoment_median", "Inverse Difference Moment (local homogeneity; higher = smoother-reading)"),
    ("SumAverage_median", "Sum Average"),
    ("SumVariance_median", "Sum Variance"),
    ("SumEntropy_median", "Sum Entropy"),
    ("Entropy_median", "Entropy (textural randomness; higher = rougher-reading)"),
    ("DiffVariance_median", "Difference Variance"),
    ("DiffEntropy_median", "Difference Entropy"),
    ("InfoCorrelation1_median", "Information Measure of Correlation 1 (Haralick f12)"),
    ("InfoCorrelation2_median", "Information Measure of Correlation 2 (Haralick f13)"),
]


def zscore(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=0)


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color=GRID_GRAY, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(TEXTURE_TABLE)

    df["smooth_asm_z"] = zscore(df["AngularSecondMoment_median"])
    df["smooth_idm_z"] = zscore(df["InverseDifferenceMoment_median"])
    df["smooth_contrast_z"] = zscore(df["Contrast_median"])
    df["smooth_entropy_z"] = zscore(df["Entropy_median"])
    df["smoothness_z"] = (
        df["smooth_asm_z"] + df["smooth_idm_z"] - df["smooth_contrast_z"] - df["smooth_entropy_z"]
    )

    # ---- fig1: smoothness_z by species, sorted by median ----
    counts = df.groupby("Species")["smoothness_z"].count()
    big_species = counts[counts >= MIN_N_FOR_BOX].index.tolist()
    small_species = counts[counts < MIN_N_FOR_BOX].index.tolist()
    order = df[df["Species"].isin(big_species)].groupby("Species")["smoothness_z"].median().sort_values().index.tolist()

    # ONE pathological outlier (a single strain with an extreme raw Contrast
    # value, see fig4) sits at smoothness_z=-43 against a population range of
    # roughly -6..+9 for everyone else -- it alone would squash the axis for
    # all 298 strains. Excluded from the plotted axis range only (a fixed
    # floor, not a symmetric percentile, so the legitimate -6..+9 spread,
    # including the single-strain diamonds at the extremes, is NOT clipped
    # along with it) and named explicitly rather than silently dropped.
    OUTLIER_FLOOR = -10.0
    clipped = df[df["smoothness_z"] < OUTLIER_FLOOR][["strain_code", "smoothness_z"]]
    visible = df.loc[df["smoothness_z"] >= OUTLIER_FLOOR, "smoothness_z"]
    lo, hi = visible.min(), visible.max()
    pad = 0.08 * (hi - lo)
    xlim = (lo - pad, hi + pad)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.32 * (len(order) + len(small_species)))))
    box_data = [df.loc[df["Species"] == sp, "smoothness_z"].dropna().values for sp in order]
    bp = ax.boxplot(
        box_data, vert=False, positions=range(len(order)), widths=0.6,
        patch_artist=True, showfliers=False,
        boxprops=dict(facecolor=HUE, alpha=0.35, edgecolor=HUE_DARK, linewidth=1.2),
        medianprops=dict(color=HUE_DARK, linewidth=1.8),
        whiskerprops=dict(color=HUE_DARK, linewidth=1.0),
        capprops=dict(color=HUE_DARK, linewidth=1.0),
    )
    rng = np.random.default_rng(0)
    for i, sp in enumerate(order):
        vals = df.loc[df["Species"] == sp, "smoothness_z"].dropna().values
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(vals, np.full(len(vals), i) + jitter, s=14, color=POINT_GRAY, alpha=0.5, zorder=3, linewidths=0)

    y_offset = len(order)
    small_labels = []
    for j, sp in enumerate(sorted(small_species, key=lambda s: df.loc[df["Species"] == s, "smoothness_z"].median())):
        vals = df.loc[df["Species"] == sp, "smoothness_z"].dropna().values
        y = y_offset + j
        ax.scatter(vals, np.full(len(vals), y), s=22, color=HUE_DARK, alpha=0.8, zorder=3, marker="D")
        small_labels.append(f"{sp} (n={len(vals)})")

    all_labels = [f"{sp} (n={int(counts[sp])})" for sp in order] + small_labels
    ax.set_yticks(range(len(all_labels)))
    ax.set_yticklabels(all_labels, fontsize=8)
    ax.axvline(0, color=GRID_GRAY, linewidth=1.0, zorder=0)
    ax.set_xlim(*xlim)
    xlabel = "\n".join(textwrap.wrap(
        "smoothness_z  (composite: +AngularSecondMoment +InverseDifferenceMoment "
        "-Contrast -Entropy; higher = smoother-reading texture)", width=70))
    ax.set_xlabel(xlabel, fontsize=8)
    title_lines = [
        "Colony-texture morphology proxy: variance by species",
        "(diamonds = species with n<3, boxes not meaningful)",
    ]
    if len(clipped):
        note = "; ".join(f"{r.strain_code}={r.smoothness_z:.1f}" for r in clipped.itertuples())
        title_lines.append(f"axis excludes outlier(s) below {OUTLIER_FLOOR:.0f}: {note}")
    ax.set_title("\n".join(title_lines), fontsize=10)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_smoothness_by_species.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig1_smoothness_by_species.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---- fig2: overall distribution across strains ----
    # same single-outlier exclusion as fig1 (TFCN_43A-4, smoothness_z=-43),
    # named explicitly rather than silently dropped from the underlying data
    # (which is untouched -- this only affects the plotted axis/bins).
    all_vals = df["smoothness_z"].dropna()
    vals = all_vals[all_vals >= OUTLIER_FLOOR].values
    n_excluded = len(all_vals) - len(vals)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(vals, bins=30, color=HUE, alpha=0.6, edgecolor=HUE_DARK, linewidth=0.8, zorder=2)
    med, q1, q3 = np.median(vals), np.percentile(vals, 25), np.percentile(vals, 75)
    ax.axvline(med, color=HUE_DARK, linewidth=1.6, linestyle="-", zorder=3)
    ax.axvline(q1, color=HUE_DARK, linewidth=1.0, linestyle="--", zorder=3)
    ax.axvline(q3, color=HUE_DARK, linewidth=1.0, linestyle="--", zorder=3)
    ax.text(med, ax.get_ylim()[1] * 0.95, f" median={med:.2f}", fontsize=8, color=HUE_DARK, va="top")
    ax.set_xlabel("smoothness_z")
    ax.set_ylabel("n strains")
    title_lines = [f"Colony-texture morphology proxy (smoothness_z):", f"strain-level distribution (n={len(vals)})"]
    if n_excluded:
        title_lines.append(f"axis excludes {n_excluded} outlier(s) below {OUTLIER_FLOOR:.0f} (out of {len(all_vals)} total)")
    ax.set_title("\n".join(title_lines), fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=GRID_GRAY, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_smoothness_distribution.png", dpi=200)
    fig.savefig(FIG_DIR / "fig2_smoothness_distribution.pdf")
    plt.close(fig)

    # ---- fig3: raw component small multiples by species (n>=3 species only) ----
    components = [
        ("Contrast_median", "Contrast (higher = rougher-reading)"),
        ("Entropy_median", "Entropy (higher = rougher-reading)"),
        ("AngularSecondMoment_median", "Angular Second Moment (higher = smoother-reading)"),
        ("InverseDifferenceMoment_median", "Inverse Difference Moment (higher = smoother-reading)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, max(6, 0.22 * len(order) * 2)))
    for ax, (col, label) in zip(axes.flat, components):
        col_order = df[df["Species"].isin(big_species)].groupby("Species")[col].median().sort_values().index.tolist()
        data = [df.loc[df["Species"] == sp, col].dropna().values for sp in col_order]
        ax.boxplot(
            data, vert=False, positions=range(len(col_order)), widths=0.6,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor=HUE, alpha=0.35, edgecolor=HUE_DARK, linewidth=1.0),
            medianprops=dict(color=HUE_DARK, linewidth=1.4),
            whiskerprops=dict(color=HUE_DARK, linewidth=0.8),
            capprops=dict(color=HUE_DARK, linewidth=0.8),
        )
        ax.set_yticks(range(len(col_order)))
        ax.set_yticklabels(col_order, fontsize=6)
        ax.set_xlabel(label, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="x", color=GRID_GRAY, linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)
    fig.suptitle(f"Raw Haralick texture components by species (n>={MIN_N_FOR_BOX} strains only, {len(order)} species)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / "fig3_raw_texture_components_by_species.png", dpi=200)
    fig.savefig(FIG_DIR / "fig3_raw_texture_components_by_species.pdf")
    plt.close(fig)

    # ---- fig4: distribution of EACH of the 13 raw Haralick metrics across all strains ----
    # (per PI request: "several measures of smooth/rough morphology ... how they
    # distribute across the samples" -- one histogram per metric, all strains
    # pooled, independent of species, to see each measure's own spread/shape
    # before any of them get combined into the smoothness_z composite.)
    fig, axes = plt.subplots(4, 4, figsize=(14, 11))
    for ax, (col, label) in zip(axes.flat, TEXTURE_METRICS):
        vals = df[col].dropna().values
        ax.hist(vals, bins=25, color=HUE, alpha=0.6, edgecolor=HUE_DARK, linewidth=0.6, zorder=2)
        med = np.median(vals)
        ax.axvline(med, color=HUE_DARK, linewidth=1.2, zorder=3)
        ax.set_title("\n".join(textwrap.wrap(label, width=32)), fontsize=7.5)
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color=GRID_GRAY, linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)
    for ax in axes.flat[len(TEXTURE_METRICS):]:
        ax.axis("off")
    fig.suptitle(f"Distribution of each raw GLCM/Haralick texture measure across all strains (n={len(df)})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / "fig4_all_texture_metrics_distribution.png", dpi=200)
    fig.savefig(FIG_DIR / "fig4_all_texture_metrics_distribution.pdf")
    plt.close(fig)

    print(f"Wrote figures to {FIG_DIR}")
    for f in sorted(FIG_DIR.glob("*.png")):
        print(f" - {f.name}")


if __name__ == "__main__":
    main()
