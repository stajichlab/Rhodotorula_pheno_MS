#!/usr/bin/env python3
"""
Find MS2 features that differ in abundance between two species (or any two
groups defined by the 'Species' column), within one sample fraction --
e.g. what's picking R. diobovata out from the R. mucilaginosa background
in the cell-pellet PCoA (MS_pcoa_cell.png).

Method: for each feature independently --
  1. Same prevalence filter + TSS normalization as the PCoA scripts
     (prep_matrix in pcoa_ms_features.py), but *without* the 4th-root
     power transform used there -- that transform is for stabilizing a
     Bray-Curtis distance matrix, and would only rescale (not reorder)
     ranks here, so it's dropped in favor of reporting fold-change on the
     directly-interpretable TSS proportions.
  2. Mann-Whitney U test (rank-based, robust to the ~9-vs-209 sample-size
     imbalance and to the heavy right-skew typical of peak-area data;
     does not assume normality).
  3. Effect size = log2 fold-change of group medians (pseudocount = the
     smallest nonzero TSS proportion observed for that feature, so the
     pseudocount scales with the feature rather than being one constant
     that swamps low-abundance features).
  4. Benjamini-Hochberg FDR correction across all tested features.

Output (analysis/differential_features/<fraction>_<groupA>_vs_<groupB>/):
  - differential_features.csv.gz : full ranked table, one row per tested
    feature (row ID, m/z, RT, adduct, group medians, log2FC, U-stat,
    p-value, FDR q-value), sorted by q-value then |log2FC|.
  - volcano.png/.pdf              : log2FC vs -log10(p), FDR-significant
    features colored by direction.
  - top_features.png/.pdf         : per-sample TSS abundance (small
    multiples) for the top 12 features by q-value, so you can eyeball
    that the statistical hits are actually clean separations and not
    driven by one outlier sample.

Usage:
    python3 scripts/differential_features_by_species.py \
        [--species-a "Rhodotorula diobovata"] \
        [--species-b "Rhodotorula mucilaginosa"] \
        [--fraction cell|supernatant] [--prevalence-min 0.10] [--fdr 0.05]
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pcoa_ms_features import savefig_multi

REPO = Path(__file__).resolve().parent.parent
LINKED = REPO / "analysis" / "linked_data"
FEATURE_ANNOTATION_COLS = [
    "row ID", "row m/z", "row retention time", "adduct",
    "is_default_adduct", "has_ms2", "detection_count", "detection_rate",
]
UP_COLOR = "#D55E00"  # vermillion -- higher in species A
DOWN_COLOR = "#0072B2"  # blue -- higher in species B
NS_COLOR = "#BBBBBB"


def slugify(name: str) -> str:
    return name.replace("Rhodotorula ", "").replace(" ", "_").replace(".", "")


def tss_normalize(feature_df: pd.DataFrame, sample_cols: list[str], prevalence_min: float):
    """Prevalence-filter then TSS-normalize (no power transform -- see
    module docstring). Returns (kept feature_df subset with a fresh index,
    samples x features ndarray)."""
    mat = feature_df[sample_cols].to_numpy(dtype=float)
    prevalence = (mat > 0).mean(axis=1)
    keep = prevalence >= prevalence_min
    kept_annot = feature_df.loc[keep, FEATURE_ANNOTATION_COLS].reset_index(drop=True)
    mat = mat[keep]
    col_sums = mat.sum(axis=0)
    if (col_sums == 0).any():
        empty = [s for s, ok in zip(sample_cols, col_sums == 0) if ok]
        sys.exit(f"sample(s) with zero total abundance after filtering: {empty}")
    mat_norm = (mat / col_sums).T  # samples x features
    return kept_annot, mat_norm


def test_features(mat_a: np.ndarray, mat_b: np.ndarray) -> pd.DataFrame:
    n_features = mat_a.shape[1]
    pvals = np.empty(n_features)
    log2fc = np.empty(n_features)
    median_a = np.empty(n_features)
    median_b = np.empty(n_features)
    ustat = np.empty(n_features)

    for i in range(n_features):
        a, b = mat_a[:, i], mat_b[:, i]
        median_a[i], median_b[i] = np.median(a), np.median(b)
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        ustat[i], pvals[i] = u, p
        pseudocount = min(x[x > 0].min() if (x > 0).any() else 1e-12 for x in (a, b)) / 2
        log2fc[i] = np.log2((median_a[i] + pseudocount) / (median_b[i] + pseudocount))

    _, qvals, _, _ = multipletests(pvals, method="fdr_bh")
    return pd.DataFrame(
        {
            "median_a": median_a, "median_b": median_b, "log2FC_a_over_b": log2fc,
            "U_stat": ustat, "p_value": pvals, "q_value": qvals,
        }
    )


def plot_volcano(df: pd.DataFrame, species_a: str, species_b: str, fdr: float, out_path: Path):
    sig_up = (df["q_value"] < fdr) & (df["log2FC_a_over_b"] > 0)
    sig_down = (df["q_value"] < fdr) & (df["log2FC_a_over_b"] < 0)
    ns = ~(sig_up | sig_down)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        df.loc[ns, "log2FC_a_over_b"], -np.log10(df.loc[ns, "p_value"]),
        s=8, color=NS_COLOR, edgecolor="none", alpha=0.5, label=f"not significant (n={ns.sum()})",
    )
    ax.scatter(
        df.loc[sig_up, "log2FC_a_over_b"], -np.log10(df.loc[sig_up, "p_value"]),
        s=14, color=UP_COLOR, edgecolor="white", linewidth=0.3,
        label=f"higher in {species_a} (n={sig_up.sum()})",
    )
    ax.scatter(
        df.loc[sig_down, "log2FC_a_over_b"], -np.log10(df.loc[sig_down, "p_value"]),
        s=14, color=DOWN_COLOR, edgecolor="white", linewidth=0.3,
        label=f"higher in {species_b} (n={sig_down.sum()})",
    )
    ax.set_xlabel(f"log2 fold-change (median {species_a} / median {species_b})")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"Differential MS2 features: {species_a} vs {species_b}\n(FDR < {fdr:.0%} colored)")
    ax.axvline(0, color="#DDDDDD", linewidth=0.8, zorder=0)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def plot_top_features(
    df: pd.DataFrame, mat_a: np.ndarray, mat_b: np.ndarray,
    species_a: str, species_b: str, out_path: Path, top_n: int = 12,
):
    """df must already be sorted by significance (most significant first)
    and have a fresh 0..n-1 RangeIndex matching mat_a/mat_b's columns."""
    top = df.head(top_n)

    ncols = 4
    nrows = -(-top_n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, (_, row) in zip(axes, top.iterrows()):
        i = int(row.name)
        rng = np.random.default_rng(0)
        xa = 0 + rng.uniform(-0.12, 0.12, size=mat_a.shape[0])
        xb = 1 + rng.uniform(-0.12, 0.12, size=mat_b.shape[0])
        ax.scatter(xa, mat_a[:, i], s=14, color=UP_COLOR, edgecolor="none", alpha=0.8)
        ax.scatter(xb, mat_b[:, i], s=14, color=DOWN_COLOR, edgecolor="none", alpha=0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([species_a.split()[-1], species_b.split()[-1]], fontsize=7)
        ax.set_yscale("symlog", linthresh=1e-6)
        ax.set_title(
            f"row {row['row ID']:.0f}  m/z {row['row m/z']:.2f}\nq={row['q_value']:.1e}  log2FC={row['log2FC_a_over_b']:.2f}",
            fontsize=7,
        )
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[len(top):]:
        ax.set_visible(False)

    fig.suptitle(f"Top {top_n} differential features: {species_a} vs {species_b} (TSS proportion, symlog y)", fontsize=10)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-a", default="Rhodotorula diobovata")
    ap.add_argument("--species-b", default="Rhodotorula mucilaginosa")
    ap.add_argument("--fraction", choices=["cell", "supernatant"], default="cell")
    ap.add_argument("--prevalence-min", type=float, default=0.10)
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()

    meta = pd.read_csv(LINKED / "sample_metadata.csv.gz")
    feature_df = pd.read_csv(LINKED / "feature_abundance_matrix.csv.gz")

    sub_meta = meta[meta["fraction"] == args.fraction]
    ids_a = sub_meta.loc[sub_meta["Species"] == args.species_a, "sample_id"].tolist()
    ids_b = sub_meta.loc[sub_meta["Species"] == args.species_b, "sample_id"].tolist()
    if not ids_a or not ids_b:
        sys.exit(f"no {args.fraction} samples for '{args.species_a}' and/or '{args.species_b}'")

    kept_annot, mat_all = tss_normalize(feature_df, ids_a + ids_b, args.prevalence_min)
    mat_a, mat_b = mat_all[: len(ids_a)], mat_all[len(ids_a) :]

    stats = test_features(mat_a, mat_b)
    result = pd.concat([kept_annot, stats], axis=1)
    result = result.assign(_abs_log2fc=result["log2FC_a_over_b"].abs()).sort_values(
        ["q_value", "_abs_log2fc"], ascending=[True, False]
    ).drop(columns="_abs_log2fc")

    out_dir = REPO / "analysis" / "differential_features" / f"{args.fraction}_{slugify(args.species_a)}_vs_{slugify(args.species_b)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_dir / "differential_features.csv.gz", index=False)
    plot_volcano(stats, args.species_a, args.species_b, args.fdr, out_dir / "volcano.png")
    # result's index still holds each row's original position in mat_a/mat_b's
    # columns (concat + sort_values preserves index labels, only reorders rows)
    # -- plot_top_features relies on that to pull the right column per feature.
    plot_top_features(result, mat_a, mat_b, args.species_a, args.species_b, out_dir / "top_features.png")

    n_sig = (stats["q_value"] < args.fdr).sum()
    print(
        f"{args.fraction}: {len(ids_a)} {args.species_a} vs {len(ids_b)} {args.species_b} samples, "
        f"{len(stats)}/{len(feature_df)} features tested (prevalence >= {args.prevalence_min:.0%}), "
        f"{n_sig} significant at FDR < {args.fdr:.0%}",
        file=sys.stderr,
    )
    print(f"wrote outputs to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
