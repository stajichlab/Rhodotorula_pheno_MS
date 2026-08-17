#!/usr/bin/env python3
"""
PCoA (Bray-Curtis) ordination of the MS2 feature-abundance data, using the
linked tables from build_analysis_table.py.

Three ordinations, per the project's request:
  1. combined  -- cell + supernatant samples together, colored by fraction,
                  with a thin line joining each strain's cell/supernatant
                  pair (AGENTS.md: the cell-vs-supernatant contrast for a
                  given strain is a specific interest).
  2. cell      -- cell-pellet samples only, one marker shape per species
                  (no "Other" bucket) with a cycled color as a secondary
                  channel.
  3. supernatant -- same, supernatant samples only.

Preprocessing per ordination (features are independently filtered/
normalized within each sample subset, since prevalence depends on which
samples are included):
  1. Prevalence filter: keep features detected (peak area > 0) in >= 10%
     of the samples in that subset.
  2. Total-sum-scaling (TSS): divide each sample's feature vector by its
     own total, so every sample's kept-feature abundances sum to 1.
  3. Fourth-root power transform (variance-stabilizing): without it, the
     top ~20 of ~13,500 kept features carry ~28% of the TSS-normalized
     mass, so Bray-Curtis is effectively ordinating a few hundred
     high-intensity ions rather than the whole feature set. A log
     transform would work too but would produce negative values (TSS
     proportions are all < 1) that Bray-Curtis can't take; a power
     transform stays non-negative while still compressing the dynamic
     range.
  4. Bray-Curtis distance between all sample pairs.
  5. Classical (Torgerson/Gower) PCoA: double-center the squared distance
     matrix and eigendecompose; report axes 1-2 and the fraction of
     variance they explain (of the *positive* eigenvalues only -- Bray-
     Curtis is non-Euclidean, so a classical PCoA typically has some
     small negative eigenvalues; the standard convention is to exclude
     them from the variance-explained denominator, which is noted on
     each plot). This differs from scikit-bio's `pcoa()`, which divides
     by the sum of *all* eigenvalues including negative ones -- ours
     reads slightly higher for the same distance matrix; both are
     defensible, just don't mix the two conventions when comparing.

Species shapes/colors are assigned once from the *full* 318-strain YPD2
phenotype table (canonical_species_order(), memoized) rather than
per-subset, so a given species draws the same marker and color in every
plot in this project (MS cell/supernatant here, color-phenotype PCoA and
ab_plane in pcoa_color_phenotype.py) -- ranking a subset's own
value_counts() would let two plots disagree on which shape means which
species whenever their sample composition differs.

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
from matplotlib.lines import Line2D
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
UNKNOWN_COLOR = "#999999"  # neutral gray, reserved for missing species

# One shape per species, no "Other" bucket -- color alone caps out at the
# 8-hue colorblind-safe palette, but shape has no such limit, so every
# species gets its own marker. Rhodotorula species draw from the plain-
# polygon set (most-common species first); non-Rhodotorula genera get
# the star-family shapes so a different genus reads as visually distinct
# at a glance; "Unknown" (missing Species) is always a bold X. Ordered so
# that stepping through by color (mod len(SPECIES_PALETTE) = 6) never
# lands two of the roundish, easy-to-confuse shapes (h/8/H/heptagon/
# nonagon) on the same color -- they're deliberately spread >=6 apart.
RHODOTORULA_MARKERS = [
    "o", "^", "s", "D", "v", "P", "<", ">", "d", "p",
    "h", (7, 0, 0), "8", "H", (9, 0, 0),
]
OTHER_GENUS_MARKERS = ["*", (6, 1, 0), (5, 1, 0), (4, 1, 0)]
UNKNOWN_MARKER = "X"

_canonical_order_cache = None


def full_species_order(species: pd.Series) -> tuple[list[str], dict]:
    """All species individually (no 'Other' bucket), Rhodotorula sorted by
    descending count (ties broken alphabetically, so the order is
    deterministic) first (drawing from RHODOTORULA_MARKERS), then other
    genera the same way (drawing from OTHER_GENUS_MARKERS), then 'Unknown'
    for missing Species. Returns (order, marker_map).

    Raises if a subset has more species in either group than there are
    markers for it, rather than silently truncating (the old `zip`-based
    version would drop the excess species from `markers` and KeyError
    later when the plot code looked one up)."""
    counts = species.value_counts()

    def ranked(names):
        return sorted(names, key=lambda s: (-counts[s], s))

    rhodo = ranked([s for s in counts.index if s.startswith("Rhodotorula")])
    other_genus = ranked([s for s in counts.index if not s.startswith("Rhodotorula")])
    if len(rhodo) > len(RHODOTORULA_MARKERS):
        raise ValueError(f"{len(rhodo)} Rhodotorula species but only {len(RHODOTORULA_MARKERS)} markers")
    if len(other_genus) > len(OTHER_GENUS_MARKERS):
        raise ValueError(f"{len(other_genus)} non-Rhodotorula genera but only {len(OTHER_GENUS_MARKERS)} markers")

    order = rhodo + other_genus + (["Unknown"] if species.isna().any() else [])
    markers = dict(zip(rhodo, RHODOTORULA_MARKERS))
    markers.update(dict(zip(other_genus, OTHER_GENUS_MARKERS)))
    markers["Unknown"] = UNKNOWN_MARKER
    return order, markers


def canonical_species_order() -> tuple[list[str], dict]:
    """The full_species_order() computed once from every phenotyped strain
    (YPD2, 318 rows, one per strain) and memoized, so every plot in this
    project draws a given species with the *same* shape regardless of
    which sample subset (MS cell-only, supernatant-only, all 318
    color-phenotyped strains, ...) it happens to be plotting -- computing
    the order fresh from each subset's own value_counts() would let two
    plots disagree whenever their species composition differs."""
    global _canonical_order_cache
    if _canonical_order_cache is None:
        ypd2_species_path = (
            REPO / "data" / "metadata" / "EXFAB_UCR-005" / "YPD2_phenotypic.20260702.fixed.csv.gz"
        )
        species = pd.read_csv(ypd2_species_path, usecols=["Species"])["Species"]
        _canonical_order_cache = full_species_order(species)
    return _canonical_order_cache


def full_species_color_map(order: list[str]) -> dict:
    """Color is a secondary channel here (shape already makes every
    species unique) -- cycle the 6-color categorical palette across
    Rhodotorula species (spaced so shape-confusable species never also
    share a color, see RHODOTORULA_MARKERS), give the non-Rhodotorula
    genera a fixed black to match their star shapes' 'stands apart'
    reading, and Unknown its usual neutral gray."""
    colors = {}
    rhodo = [s for s in order if s != "Unknown" and s.startswith("Rhodotorula")]
    other_genus = [s for s in order if s != "Unknown" and not s.startswith("Rhodotorula")]
    for i, sp in enumerate(rhodo):
        colors[sp] = SPECIES_PALETTE[i % len(SPECIES_PALETTE)]
    for sp in other_genus:
        colors[sp] = "#000000"
    colors["Unknown"] = UNKNOWN_COLOR
    return colors


def prep_matrix(feature_df: pd.DataFrame, sample_cols: list[str]):
    """Prevalence-filter, TSS-normalize, then fourth-root transform.
    Returns (samples x features ndarray, n_features_kept, n_features_total)."""
    mat = feature_df[sample_cols].to_numpy(dtype=float)
    prevalence = (mat > 0).mean(axis=1)
    keep = prevalence >= PREVALENCE_MIN
    mat = mat[keep]
    col_sums = mat.sum(axis=0)
    if (col_sums == 0).any():
        empty = [s for s, ok in zip(sample_cols, col_sums == 0) if ok]
        sys.exit(f"sample(s) with zero total abundance after filtering: {empty}")
    mat_norm = mat / col_sums
    mat_norm = mat_norm**0.25
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
    pos = eigvals > 1e-8 * eigvals[0]  # relative, not absolute -- scale-independent
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


def savefig_multi(fig, png_path: Path):
    """Save both a PNG (raster, dpi=150) and a PDF (vector) of the figure,
    same basename."""
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(png_path.with_suffix(".pdf"), bbox_inches="tight")


def plot_combined(axes: pd.DataFrame, meta: pd.DataFrame, prop, out_path: Path):
    df = axes.merge(meta[["sample_id", "fraction", "Strain ID"]], on="sample_id")
    fig, ax = plt.subplots(figsize=(7, 6))

    # Thin gray line joining each strain's cell/supernatant pair. Grouped by
    # the numeric Strain ID, not canonical_strain -- a strain *code* can be
    # legitimately reused across two different Strain IDs (e.g. DBVPG_3853
    # at ids 192 and 195; see build_analysis_table.py), which would silently
    # 4-way-group and skip those pairs (or worse, cross-wire them) under a
    # text-based groupby.
    for _, g in df.groupby("Strain ID"):
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
    order, markers = canonical_species_order()
    colors = full_species_color_map(order)
    df = df.assign(species_label=df["Species"].fillna("Unknown"))

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    for sp in order:
        sub = df[df["species_label"] == sp]
        if sub.empty:
            continue
        ax.scatter(
            sub["PCoA1"], sub["PCoA2"], marker=markers[sp], s=28, color=colors[sp],
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
    ax.set_title(f"{title}\n(star-family shapes = non-Rhodotorula genera)")
    ax.legend(
        handles=handles, title="Species", frameon=False, fontsize=7,
        loc="center left", bbox_to_anchor=(1.0, 0.5), labelspacing=0.6,
    )
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
