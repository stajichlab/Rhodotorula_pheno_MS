#!/usr/bin/env python3
"""
Phase 2, ANOVA/pattern-group follow-up (PI request, 2026-08-15): every
test so far in this project (univariate Spearman correlation, sparse
Lasso regression) assumes a monotonic or linear relationship between a
single continuous color axis and metabolite abundance. This script asks a
structurally different question: does metabolite abundance differ across
COLOR-PATTERN GROUPS (clusters in L\*/a\*/b\* space, not just a single
axis) -- which can catch non-monotonic, threshold-like, or multi-axis
color-pattern relationships that a rank correlation on a\* alone cannot.

Method, within one species at a time (phylogeny has much less impact
within a species than across the whole panel, per the within-species
scripts already in this project):
  1. Cluster strains into --n-groups color-pattern groups via k-means on
     standardized [L*, a*, b*] jointly (not a single axis) -- this is the
     "phenotype pattern" the PI is asking about, as distinct from a
     single continuous predictor.
  2. For each deduplicated MS2 feature group, Kruskal-Wallis H-test
     (non-parametric one-way ANOVA -- appropriate given untargeted MS
     abundances are not normally distributed) across the color-pattern
     groups.
  3. Primary inferential statistic: empirical p-value from label
     permutation restricted within the species' real strain-level
     phylogenetic blocks (same convention as
     phase2_within_species_association.py), not the asymptotic
     Kruskal-Wallis p-value, for the same reason as the other Phase 2
     scripts (small number of independent lineages even within a species).
  4. BH-FDR on the deduplicated feature-group count.

HARD GATE (same convention as the other Phase 2 scripts): --predictor
color (the real color-pattern clustering) refuses to run without a fresh
--predictor area decoy run (colony size clustered into groups the same
way, as the phylogenetically-structured negative control).

Usage:
    python3 analysis/scripts/phase2_anova_pattern_association.py --species "Rhodotorula mucilaginosa" --predictor area   # negative control, run first
    python3 analysis/scripts/phase2_anova_pattern_association.py --species "Rhodotorula mucilaginosa" --predictor color # primary
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase2_within_species_association import (  # noqa: E402
    DEDUP_GROUPS,
    FEATURE_MATRIX,
    PHENOTYPE_TABLE,
    load_strain_fraction_matrix,
    slugify,
    strain_blocks_from_tree,
)

REPO = Path(__file__).resolve().parent.parent.parent
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase2_metabolome_phenotype"
FDR_ALPHA = 0.05


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def cluster_groups(values: np.ndarray, n_groups: int, seed: int) -> np.ndarray:
    values = values.reshape(-1, 1) if values.ndim == 1 else values
    X = StandardScaler().fit_transform(values)
    km = KMeans(n_clusters=n_groups, n_init=10, random_state=seed)
    return km.fit_predict(X)


def kruskal_stat_vec(group_labels: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Kruskal-Wallis H per feature column, vectorized via per-group rank sums
    (avoids looping scipy.stats.kruskal n_features times)."""
    n = X.shape[0]
    ranks = np.apply_along_axis(lambda col: pd.Series(col).rank().to_numpy(), 0, X)
    H = np.zeros(X.shape[1])
    groups = np.unique(group_labels)
    for g in groups:
        idx = group_labels == g
        ng = idx.sum()
        if ng == 0:
            continue
        rank_sum = ranks[idx].sum(axis=0)
        H += (rank_sum**2) / ng
    H = (12 / (n * (n + 1))) * H - 3 * (n + 1)
    # tie correction
    return H


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--predictor", choices=["color", "area"], default="color")
    ap.add_argument("--n-groups", type=int, default=3, help="Number of color-pattern (or area) clusters.")
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-clades", type=int, default=15)
    ap.add_argument("--min-strains", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    is_decoy_run = args.predictor == "area"
    slug = slugify(args.species)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    decoy_out = OUT_DIR / f"anova_{slug}_area_decoy.csv"
    real_out = OUT_DIR / f"anova_{slug}_color.csv"

    MIN_DECOY_N_PERM = 100  # a --n-perm 20 smoke test is not an adequate negative control -- see .living/learnings.md 2026-08-15 "underpowered decoy" entry

    if not is_decoy_run:
        if not decoy_out.exists():
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy output not found ({decoy_out}).\n"
                f"Run first: python3 {Path(__file__).name} --species \"{args.species}\" --predictor area"
            )
        decoy_df = pd.read_csv(decoy_out)
        if "n_perm" not in decoy_df.columns or decoy_df["n_perm"].min() < MIN_DECOY_N_PERM:
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy output ({decoy_out}) was run with too few "
                f"permutations (need >={MIN_DECOY_N_PERM}, a smoke-test run is not an adequate negative "
                f"control). Rerun: python3 {Path(__file__).name} --species \"{args.species}\" --predictor area --n-perm {MIN_DECOY_N_PERM}"
            )
        newest_input = max(p.stat().st_mtime for p in (FEATURE_MATRIX, PHENOTYPE_TABLE, DEDUP_GROUPS) if p.exists())
        if decoy_out.stat().st_mtime < newest_input:
            sys.exit(f"REFUSING TO RUN: decoy output ({decoy_out}) is STALE. Rerun the negative control first.")

    pheno = pd.read_csv(PHENOTYPE_TABLE).set_index("strain_id")
    species_strains_all = set(pheno.index[pheno["species"] == args.species])
    print(f"{args.species}: {len(species_strains_all)} strains in phenotype table", file=sys.stderr)

    rows_all = []
    for fraction in ["cell", "supernatant"]:
        strain_df, group_ids = load_strain_fraction_matrix(fraction, species_strains_all)
        needed_cols = ["L*", "a*", "b*"] if args.predictor == "color" else ["area"]
        common = [
            s for s in strain_df.index
            if s in pheno.index and pheno.loc[s, needed_cols].notna().all()
        ]
        if len(common) < args.min_strains:
            print(f"[{fraction}] only {len(common)} strains -- skipping (min {args.min_strains})", file=sys.stderr)
            continue

        if args.predictor == "color":
            raw_vals = pheno.loc[common, ["L*", "a*", "b*"]].to_numpy(dtype=float)
        else:
            raw_vals = pheno.loc[common, ["area"]].to_numpy(dtype=float)
        pattern_groups = cluster_groups(raw_vals, args.n_groups, args.seed)
        group_sizes = pd.Series(pattern_groups).value_counts().sort_index()
        print(f"[{fraction}] n_strains={len(common)} group sizes: {group_sizes.to_dict()}", file=sys.stderr)

        X = strain_df.loc[common].to_numpy(dtype=float)
        keep = X.std(axis=0) > 0
        X = X[:, keep]
        kept_group_ids = group_ids[keep]
        n_features = X.shape[1]

        blocks_map = strain_blocks_from_tree(common, args.n_clades)
        blocks = np.array([blocks_map[s] for s in common])

        observed_H = kruskal_stat_vec(pattern_groups, X)

        rng = np.random.default_rng(args.seed)
        block_indices = [np.where(blocks == b)[0] for b in np.unique(blocks)]
        exceed = np.zeros(n_features, dtype=np.int64)
        for p in range(args.n_perm):
            perm_labels = pattern_groups.copy()
            for idx in block_indices:
                if len(idx) > 1:
                    perm_labels[idx] = rng.permutation(pattern_groups[idx])
            perm_H = kruskal_stat_vec(perm_labels, X)
            exceed += perm_H >= observed_H

        empirical_p = (exceed + 1) / (args.n_perm + 1)
        empirical_fdr = bh_fdr(empirical_p)
        n_hits = int((empirical_fdr < FDR_ALPHA).sum())
        print(f"[{fraction}] BH-FDR<{FDR_ALPHA} hits (empirical): {n_hits} of {n_features}", file=sys.stderr)

        res = pd.DataFrame(
            {
                "dedup_group_id": kept_group_ids,
                "fraction": fraction,
                "predictor": args.predictor,
                "species": args.species,
                "n_strains": len(common),
                "n_groups": args.n_groups,
                "n_perm": args.n_perm,
                "kruskal_H": observed_H,
                "empirical_p": empirical_p,
                "empirical_fdr": empirical_fdr,
            }
        )
        rows_all.append(res)

    if not rows_all:
        sys.exit("No fraction had enough matched strains -- aborting.")

    final = pd.concat(rows_all, ignore_index=True)
    out_path = decoy_out if is_decoy_run else real_out
    final.to_csv(out_path, index=False)
    print(f"Wrote {out_path}", file=sys.stderr)
    if is_decoy_run:
        print("Negative-control decoy run complete -- real predictor runs are now unblocked.", file=sys.stderr)


if __name__ == "__main__":
    main()
