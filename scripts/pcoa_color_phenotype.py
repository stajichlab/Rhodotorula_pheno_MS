#!/usr/bin/env python3
"""
PCoA of the CIELAB color phenotype (L*, a*, b*) across all 318 phenotyped
strains, to see how they cluster in color space.

Distance metric: Euclidean, not Bray-Curtis. CIELAB was specifically
designed so that Euclidean distance between two points approximates
perceived color difference (DeltaE); a*/b* are signed (green-red,
blue-yellow) so a compositional/non-negative metric like Bray-Curtis
doesn't apply here the way it does to MS peak areas. With only 3 input
variables, classical PCoA on a Euclidean distance matrix is equivalent to
PCA on the raw (unstandardized) L*/a*/b* values -- unstandardized is the
right choice here because the whole point of CIELAB is that its units are
already perceptually comparable across axes.

Species shapes/colors come from pcoa_ms_features.canonical_species_order(),
computed once from the full 318-strain phenotype table and shared with
the MS feature ordinations, so a given species draws the same marker
everywhere in this project rather than each plot re-ranking its own
sample subset.

Three plots are written:
  - pcoa_color_by_species.png/.pdf : PCoA points by species, one shape
    per species (no "Other" bucket) with a cycled color as a secondary
    channel.
  - pcoa_color_swatches.png/.pdf   : PCoA points colored by the strain's
    own measured color (Lab -> sRGB, D65, clipped to gamut) -- lets you
    see directly whether nearby points in the ordination actually look
    alike.
  - ab_plane_swatches.png/.pdf     : NOT a PCoA -- the raw a*/b*
    chromaticity plane (a*=green-red, b*=blue-yellow; this pair *is* a
    natural, uncomputed 2D space, unlike the PCoA axes above), point
    fill = the strain's own measured color, marker *shape* = species (so
    color is free to encode the actual color, per request).

Usage:
    python3 scripts/pcoa_color_phenotype.py
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pcoa_ms_features import (
    canonical_species_order,
    classical_pcoa,
    full_species_color_map,
    savefig_multi,
)

REPO = Path(__file__).resolve().parent.parent
YPD2_FIXED = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "YPD2_phenotypic.20260702.fixed.csv.gz"
OUT_DIR = REPO / "analysis" / "color_phenotype_ordination"
LAB_COLS = ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
# This dataset's a*/b* chroma is genuinely low (a* mostly 0-14, b* -1 to 9,
# L*~74), so the true-color swatches read as near-monochrome pale pink --
# real biology, not a bug, but it hides structure visually. The companion
# "chroma-amplified" plots scale a*/b* (not L*, and not the point's actual
# plotted position) by this factor before Lab->sRGB conversion, purely to
# make hue/chroma differences between strains easier to see; they are
# clearly labeled as not true color.
CHROMA_GAIN = 3.0


def lab_to_srgb(lab: np.ndarray) -> np.ndarray:
    """CIELAB (D65) -> sRGB, clipped to [0, 1] for any out-of-gamut colors."""
    L, a, b = lab[:, 0], lab[:, 1], lab[:, 2]
    xn, yn, zn = 95.047, 100.0, 108.883

    fy = (L + 16) / 116
    fx = fy + a / 500
    fz = fy - b / 200

    def finv(t):
        t3 = t**3
        return np.where(t3 > 0.008856, t3, (t - 16 / 116) / 7.787)

    X = xn * finv(fx) / 100
    Y = yn * finv(fy) / 100
    Z = zn * finv(fz) / 100

    R = 3.2406 * X - 1.5372 * Y - 0.4986 * Z
    G = -0.9689 * X + 1.8758 * Y + 0.0415 * Z
    B = 0.0557 * X - 0.2040 * Y + 1.0570 * Z
    rgb = np.stack([R, G, B], axis=1)

    def gamma(c):
        return np.where(c <= 0.0031308, 12.92 * c, 1.055 * np.clip(c, 0, None) ** (1 / 2.4) - 0.055)

    return np.clip(gamma(rgb), 0, 1)


def plot_by_species(df, prop, out_path: Path):
    order, markers = canonical_species_order()
    colors = full_species_color_map(order)
    df = df.assign(species_label=df["Species"].fillna("Unknown"))

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    for sp in order:
        sub = df[df["species_label"] == sp]
        if sub.empty:
            continue
        ax.scatter(
            sub["PCoA1"], sub["PCoA2"], marker=markers[sp], s=32, color=colors[sp],
            edgecolor="white", linewidth=0.4,
        )
    handles = [
        Line2D(
            [0], [0], marker=markers[sp], linestyle="none", markerfacecolor=colors[sp],
            markeredgecolor="white", markersize=8,
            label=f"{sp} (n={(df['species_label'] == sp).sum()})",
        )
        for sp in order
        if (df["species_label"] == sp).any()
    ]
    ax.set_xlabel(f"PCoA1 ({prop[0]:.1%})")
    ax.set_ylabel(f"PCoA2 ({prop[1]:.1%})")
    ax.set_title(
        "PCoA of CIELAB color phenotype (L*, a*, b*) -- by species\n"
        "(star-family shapes = non-Rhodotorula genera)"
    )
    ax.legend(
        handles=handles, title="Species", frameon=False, fontsize=7,
        loc="center left", bbox_to_anchor=(1.0, 0.5), labelspacing=0.6,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def plot_swatches(df, prop, out_path: Path, swatch_col="swatch", chroma_gain=None):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        df["PCoA1"], df["PCoA2"], s=60, c=df[swatch_col].tolist(),
        edgecolor="#333333", linewidth=0.4,
    )
    ax.set_xlabel(f"PCoA1 ({prop[0]:.1%})")
    ax.set_ylabel(f"PCoA2 ({prop[1]:.1%})")
    title = "PCoA of CIELAB color phenotype (L*, a*, b*) -- point = strain's own color"
    if chroma_gain:
        title += f"\n(a*/b* amplified {chroma_gain:.0f}x for visibility -- NOT true color)"
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def plot_ab_plane(df, out_path: Path, swatch_col="swatch", chroma_gain=None):
    order, markers = canonical_species_order()
    df = df.assign(species_label=df["Species"].fillna("Unknown"))

    fig, ax = plt.subplots(figsize=(8.5, 7))
    ax.axhline(0, color="#DDDDDD", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#DDDDDD", linewidth=0.8, zorder=0)
    for sp in order:
        sub = df[df["species_label"] == sp]
        if sub.empty:
            continue
        ax.scatter(
            sub["a"], sub["b"], marker=markers[sp], s=55,
            c=sub[swatch_col].tolist(), edgecolor="#333333", linewidth=0.5, zorder=2,
        )

    handles = [
        Line2D(
            [0], [0], marker=markers[sp], linestyle="none", markerfacecolor="#AAAAAA",
            markeredgecolor="#333333", markersize=8,
            label=f"{sp} (n={(df['species_label'] == sp).sum()})",
        )
        for sp in order
        if (df["species_label"] == sp).any()
    ]
    ax.legend(
        handles=handles, title="Species (shape)", frameon=False, fontsize=7,
        loc="center left", bbox_to_anchor=(1.0, 0.5), labelspacing=0.6,
    )
    ax.set_xlabel("a*  (green- / red+)")
    ax.set_ylabel("b*  (blue- / yellow+)")
    title = (
        "CIELAB chromaticity plane -- fill = strain's own color, shape = species\n"
        "(star-family shapes = non-Rhodotorula genera)"
    )
    if chroma_gain:
        title += f"\n(fill a*/b* amplified {chroma_gain:.0f}x for visibility -- point position is still true a*/b*)"
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ypd2 = pd.read_csv(YPD2_FIXED)
    ypd2 = ypd2.dropna(subset=LAB_COLS).reset_index(drop=True)

    lab = ypd2[LAB_COLS].to_numpy(dtype=float)
    dist = squareform(pdist(lab, metric="euclidean"))
    coords, prop, eigvals = classical_pcoa(dist)
    n_negative = int((eigvals < -1e-8).sum())
    print(
        f"{len(ypd2)} strains, axis1={prop[0]:.1%} axis2={prop[1]:.1%} of "
        f"positive-eigenvalue variance ({n_negative} negative eigenvalues -- "
        "expected to be ~0 since Euclidean distance in <=3 dims is exactly "
        "embeddable)",
        file=sys.stderr,
    )

    df = pd.DataFrame(
        {
            "Strain ID": ypd2["Strain ID"],
            "Strain": ypd2["Strain"],
            "Species": ypd2["Species"],
            "PCoA1": coords[:, 0],
            "PCoA2": coords[:, 1],
            "L": lab[:, 0],
            "a": lab[:, 1],
            "b": lab[:, 2],
        }
    )
    df["swatch"] = list(lab_to_srgb(lab))
    lab_amplified = lab.copy()
    lab_amplified[:, 1:] *= CHROMA_GAIN  # scale a*/b* only, leave L* alone
    df["swatch_chroma_amp"] = list(lab_to_srgb(lab_amplified))

    df.drop(columns=["swatch", "swatch_chroma_amp"]).to_csv(
        OUT_DIR / "pcoa_axes_color_phenotype.csv", index=False
    )
    plot_by_species(df, prop, OUT_DIR / "pcoa_color_by_species.png")
    plot_swatches(df, prop, OUT_DIR / "pcoa_color_swatches.png")
    plot_swatches(
        df, prop, OUT_DIR / "pcoa_color_swatches_chroma_amplified.png",
        swatch_col="swatch_chroma_amp", chroma_gain=CHROMA_GAIN,
    )
    plot_ab_plane(df, OUT_DIR / "ab_plane_swatches.png")
    plot_ab_plane(
        df, OUT_DIR / "ab_plane_swatches_chroma_amplified.png",
        swatch_col="swatch_chroma_amp", chroma_gain=CHROMA_GAIN,
    )
    print(f"wrote plots + axis coordinates to {OUT_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
