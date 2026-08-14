#!/usr/bin/env python3
"""
PCoA (Bray-Curtis) ordination of the MS2 feature-abundance data, using the
linked tables from build_analysis_table.py.

Three ordinations, per the project's request:
  1. combined  -- cell + supernatant samples together, colored by fraction,
                  with a thin line joining each strain's cell/supernatant
                  pair (AGENTS.md: the cell-vs-supernatant contrast for a
                  given strain is a specific interest).
  2. cell      -- cell-pellet samples only, colored by species.
  3. supernatant -- supernatant samples only, colored by species.

Preprocessing per ordination (features are independently filtered/
normalized within each sample subset, since prevalence depends on which
samples are included):
  1. Prevalence filter: keep features detected (peak area > 0) in >= 10%
     of the samples in that subset.
  2. Total-sum-scaling (TSS): divide each sample's feature vector by its
     own total, so every sample's kept-feature abundances sum to 1.
  3. Bray-Curtis distance between all sample pairs.
  4. Classical (Torgerson/Gower) PCoA: double-center the squared distance
     matrix and eigendecompose; report axes 1-2 and the fraction of
     variance they explain (of the *positive* eigenvalues only -- Bray-
     Curtis is non-Euclidean, so a classical PCoA typically has some
     small negative eigenvalues; the standard convention is to exclude
     them from the variance-explained denominator, which is noted on
     each plot).

Usage:
    python3 scripts/pcoa_ms_features.py
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

REPO = Path(__file__).resolve().parent.parent
LINKED = REPO / "analysis" / "linked_data"
OUT_DIR = REPO / "analysis" / "ms_feature_ordination"
PREVALENCE_MIN = 0.10

# Okabe-Ito colorblind-safe categorical palette (fixed assignment, not cycled).
FRACTION_COLORS = {"cell": "#0072B2", "supernatant": "#D55E00"}
SPECIES_PALETTE = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # bluish green
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
]
OTHER_COLOR = "#F0E442"  # yellow, reserved for the "Other species" bucket
UNKNOWN_COLOR = "#999999"  # neutral gray, reserved for missing species


def prep_matrix(feature_df: pd.DataFrame, sample_cols: list[str]):
    """Prevalence-filter then TSS-normalize. Returns (samples x features
    ndarray, n_features_kept, n_features_total)."""
    mat = feature_df[sample_cols].to_numpy(dtype=float)
    prevalence = (mat > 0).mean(axis=1)
    keep = prevalence >= PREVALENCE_MIN
    mat = mat[keep]
    col_sums = mat.sum(axis=0)
    if (col_sums == 0).any():
        empty = [s for s, ok in zip(sample_cols, col_sums == 0) if ok]
        sys.exit(f"sample(s) with zero total abundance after filtering: {empty}")
    mat_norm = mat / col_sums
    return mat_norm.T, int(keep.sum()), len(keep)


def classical_pcoa(dist_matrix: np.ndarray):
    """Torgerson/Gower classical PCoA. Returns (coords [n x k_pos],
    proportion_explained [k_pos], all_eigenvalues)."""
    n = dist_matrix.shape[0]
    d2 = dist_matrix**2
    centering = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * centering @ d2 @ centering
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-8
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    prop = eigvals[pos] / eigvals[pos].sum()
    return coords, prop, eigvals


def run_ordination(name: str, sample_ids: list[str], feature_df: pd.DataFrame):
    mat, n_kept, n_total = prep_matrix(feature_df, sample_ids)
    dist = squareform(pdist(mat, metric="braycurtis"))
    coords, prop, eigvals = classical_pcoa(dist)
    n_negative = int((eigvals < -1e-8).sum())
    print(
        f"[{name}] {len(sample_ids)} samples, {n_kept}/{n_total} features kept "
        f"(prevalence >= {PREVALENCE_MIN:.0%}), axis1={prop[0]:.1%} axis2={prop[1]:.1%} "
        f"of positive-eigenvalue variance ({n_negative} negative eigenvalues excluded)",
        file=sys.stderr,
    )
    axes = pd.DataFrame(
        {"sample_id": sample_ids, "PCoA1": coords[:, 0], "PCoA2": coords[:, 1]}
    )
    return axes, prop


def bucket_species(species: pd.Series) -> tuple[pd.Series, list[str]]:
    top = species.value_counts().index[: len(SPECIES_PALETTE)].tolist()
    order = top + ["Other", "Unknown"]

    def label(sp):
        if pd.isna(sp):
            return "Unknown"
        return sp if sp in top else "Other"

    return species.map(label), order


def species_color_map(order: list[str]) -> dict:
    colors = {}
    for sp, c in zip(order, SPECIES_PALETTE):
        colors[sp] = c
    colors["Other"] = OTHER_COLOR
    colors["Unknown"] = UNKNOWN_COLOR
    return colors


def savefig_multi(fig, png_path: Path):
    """Save both a PNG (raster, dpi=150) and a PDF (vector) of the figure,
    same basename."""
    fig.savefig(png_path, dpi=150)
    fig.savefig(png_path.with_suffix(".pdf"))


def plot_combined(axes: pd.DataFrame, meta: pd.DataFrame, prop, out_path: Path):
    df = axes.merge(meta[["sample_id", "fraction", "canonical_strain"]], on="sample_id")
    fig, ax = plt.subplots(figsize=(7, 6))

    # thin gray line joining each strain's cell/supernatant pair
    for _, g in df.groupby("canonical_strain"):
        if len(g) == 2:
            ax.plot(g["PCoA1"], g["PCoA2"], color="#BBBBBB", linewidth=0.5, zorder=1)

    for fraction in ["cell", "supernatant"]:
        sub = df[df["fraction"] == fraction]
        ax.scatter(
            sub["PCoA1"],
            sub["PCoA2"],
            s=22,
            color=FRACTION_COLORS[fraction],
            edgecolor="white",
            linewidth=0.4,
            label=fraction,
            zorder=2,
        )

    ax.set_xlabel(f"PCoA1 ({prop[0]:.1%})")
    ax.set_ylabel(f"PCoA2 ({prop[1]:.1%})")
    ax.set_title("Bray-Curtis PCoA of MS2 features -- cell + supernatant")
    ax.legend(title="Fraction", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def plot_by_species(axes: pd.DataFrame, meta: pd.DataFrame, prop, title: str, out_path: Path):
    df = axes.merge(meta[["sample_id", "Species"]], on="sample_id")
    labels, order = bucket_species(df["Species"])
    df = df.assign(species_label=labels)
    colors = species_color_map(order)

    fig, ax = plt.subplots(figsize=(7.5, 6))
    for sp in order:
        sub = df[df["species_label"] == sp]
        if sub.empty:
            continue
        ax.scatter(
            sub["PCoA1"],
            sub["PCoA2"],
            s=22,
            color=colors[sp],
            edgecolor="white",
            linewidth=0.4,
            label=f"{sp} (n={len(sub)})",
        )

    ax.set_xlabel(f"PCoA1 ({prop[0]:.1%})")
    ax.set_ylabel(f"PCoA2 ({prop[1]:.1%})")
    ax.set_title(title)
    ax.legend(title="Species", frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(LINKED / "sample_metadata.csv.gz")
    feature_df = pd.read_csv(LINKED / "feature_abundance_matrix.csv.gz")

    cell_ids = meta.loc[meta["fraction"] == "cell", "sample_id"].tolist()
    sup_ids = meta.loc[meta["fraction"] == "supernatant", "sample_id"].tolist()
    all_ids = meta["sample_id"].tolist()

    axes_all, prop_all = run_ordination("combined", all_ids, feature_df)
    axes_all.to_csv(OUT_DIR / "MS_pcoa_axes_combined.csv", index=False)
    plot_combined(axes_all, meta, prop_all, OUT_DIR / "MS_pcoa_combined.png")

    axes_cell, prop_cell = run_ordination("cell", cell_ids, feature_df)
    axes_cell.to_csv(OUT_DIR / "MS_pcoa_axes_cell.csv", index=False)
    plot_by_species(
        axes_cell, meta, prop_cell,
        "Bray-Curtis PCoA of MS2 features -- cell pellet only",
        OUT_DIR / "MS_pcoa_cell.png",
    )

    axes_sup, prop_sup = run_ordination("supernatant", sup_ids, feature_df)
    axes_sup.to_csv(OUT_DIR / "MS_pcoa_axes_supernatant.csv", index=False)
    plot_by_species(
        axes_sup, meta, prop_sup,
        "Bray-Curtis PCoA of MS2 features -- supernatant only",
        OUT_DIR / "MS_pcoa_supernatant.png",
    )

    print(f"wrote plots + axis coordinates to {OUT_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
