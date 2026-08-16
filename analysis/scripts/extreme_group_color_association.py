#!/usr/bin/env python3
"""
PI follow-up (2026-08-16): instead of a continuous Spearman correlation
across the whole panel (phase2_color_metabolome_association.py, null),
test whether metabolite features differ in abundance between STRAINS AT
THE COLOR EXTREMES -- top-quantile ("high orange/red") vs bottom-quantile
("low orange/red") on a*, using a Mann-Whitney-U-equivalent rank-sum
statistic (rank-biserial effect size). Extreme-group contrasts can have
more power to detect threshold/nonlinear effects than a continuous
correlation at the same n, at the cost of discarding the middle of the
distribution.

PI explicitly acknowledged this design lumps species together and
ignores phylogeny (dominated by *R. mucilaginosa*, 216/303 strains, whose
own a* range spans much of the whole panel -- see the species-composition
diagnostic printed per run). This is NOT a substitute for the
phylogenetically-corrected Phase 2 tests; it is a different, higher-power
QUESTION (is there a compound whose abundance tracks color extremes at
all, regardless of species membership) run alongside them -- and it still
uses species-tree-block-restricted label permutation (not a fully naive
shuffle), so it isn't completely blind to phylogenetic structure either.

Reuses phase2_color_metabolome_association.py's TSS-normalization,
dedup-group, and species-block machinery unchanged. The rank-sum
statistic is computed with a single rankdata() call per feature (fixed
data ranks) and varies only which samples are labeled high/low per
permutation -- vectorized the same way as that script's Spearman rho, so
it stays fast even at --n-perm 500 across ~10,000 features.

Predictor / negative-control conventions carried over unchanged:
  --predictor a     PRIMARY. Extreme groups on a* (CIELAB red-green axis).
  --predictor C     SECONDARY. Extreme groups on C* (chroma).
  --predictor area  DECOY. Extreme groups on colony area -- same
                     hard-gated negative-control design as
                     phase2_color_metabolome_association.py.

Cell and supernatant fractions are ALWAYS tested and reported separately
(never pooled) -- same as every other script in this project.

Usage:
    python3 analysis/scripts/extreme_group_color_association.py --predictor area   # negative control, run first
    python3 analysis/scripts/extreme_group_color_association.py --predictor a      # primary
    python3 analysis/scripts/extreme_group_color_association.py --predictor C      # secondary
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase2_color_metabolome_association import (  # noqa: E402
    PHENOTYPE_TABLE, SPECIES_TREE, PREDICTOR_COL, FDR_ALPHA,
    bh_fdr, species_blocks_from_tree, load_strain_fraction_matrix,
)

REPO = Path(__file__).resolve().parent.parent.parent
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase2_metabolome_phenotype"


def rank_sum_effect(is_high: np.ndarray, X_rank: np.ndarray) -> np.ndarray:
    """Vectorized rank-biserial effect size (Mann-Whitney U equivalent) per
    feature column. X_rank: [n_samples x n_features] ranks computed ONCE on
    the fixed data (ties-averaged); is_high: boolean group-membership mask
    that varies per permutation. Positive = high group ranks higher."""
    n1 = int(is_high.sum())
    n2 = len(is_high) - n1
    R1 = X_rank[is_high].sum(axis=0)
    U1 = R1 - n1 * (n1 + 1) / 2
    return 1 - (2 * U1) / (n1 * n2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictor", choices=["a", "C", "area"], default="a")
    ap.add_argument("--quantile", type=float, default=0.25,
                     help="Top/bottom quantile defining the 'high'/'low' extreme groups (default 0.25 = quartiles).")
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--n-clades", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    predictor_col = PREDICTOR_COL[args.predictor]
    is_decoy_run = args.predictor == "area"
    rng = np.random.default_rng(args.seed)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    decoy_out = OUT_DIR / "extreme_group_association_area_decoy.csv"
    real_out = OUT_DIR / f"extreme_group_association_{args.predictor}.csv"

    if not is_decoy_run:
        if not decoy_out.exists():
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy output not found ({decoy_out}).\n"
                f"Run first: python3 {Path(__file__).name} --predictor area"
            )
        if decoy_out.stat().st_mtime < PHENOTYPE_TABLE.stat().st_mtime:
            sys.exit(f"REFUSING TO RUN: decoy output ({decoy_out}) is STALE. Rerun: python3 {Path(__file__).name} --predictor area")

    pheno = pd.read_csv(PHENOTYPE_TABLE).set_index("strain_id")
    blocks = species_blocks_from_tree(SPECIES_TREE, args.n_clades)
    pheno["block"] = pheno["species"].map(blocks).fillna("clade_unmatched")

    lo_q, hi_q = pheno[predictor_col].quantile([args.quantile, 1 - args.quantile])
    pheno["group"] = np.select(
        [pheno[predictor_col] <= lo_q, pheno[predictor_col] >= hi_q],
        ["low", "high"], default="mid",
    )
    print(f"Predictor {predictor_col}: low group <= {lo_q:.3f}, high group >= {hi_q:.3f}", file=sys.stderr)
    print("Species composition of extreme groups (n strains):", file=sys.stderr)
    print(pheno[pheno["group"] != "mid"].groupby(["group", "species"]).size().to_string(), file=sys.stderr)

    all_results = []
    for fraction in ["cell", "supernatant"]:
        strain_df, dedup_ids, _ = load_strain_fraction_matrix(fraction)
        common = [s for s in strain_df.index if s in pheno.index and pheno.loc[s, "group"] != "mid"]
        sub_pheno = pheno.loc[common]
        n_high = int((sub_pheno["group"] == "high").sum())
        n_low = int((sub_pheno["group"] == "low").sum())
        if n_high < 3 or n_low < 3:
            print(f"[{fraction}] only {n_high} high / {n_low} low strains -- skipping", file=sys.stderr)
            continue

        X = strain_df.loc[common].to_numpy(dtype=float)
        is_high = (sub_pheno["group"] == "high").to_numpy()
        block = sub_pheno["block"].to_numpy()
        n_features = X.shape[1]

        X_rank = np.apply_along_axis(rankdata, 0, X)  # fixed, computed once
        constant = np.array([len(np.unique(X[:, fi])) < 2 for fi in range(n_features)])
        valid = ~constant
        n_valid = int(valid.sum())
        print(f"[{fraction}] n_high={n_high} n_low={n_low} n_features={n_features} predictor={predictor_col} "
              f"({n_features - n_valid} constant-abundance feature(s) excluded)", file=sys.stderr)

        observed_effect = rank_sum_effect(is_high, X_rank)

        block_indices = [np.where(block == b)[0] for b in np.unique(block)]
        exceed = np.zeros(n_valid, dtype=np.int64)
        null_hit_counts = np.zeros(args.n_perm, dtype=np.int64)
        for p in range(args.n_perm):
            perm_is_high = is_high.copy()
            for idx in block_indices:
                if len(idx) > 1:
                    perm_is_high[idx] = rng.permutation(is_high[idx])
            if perm_is_high.sum() < 2 or (~perm_is_high).sum() < 2:
                continue
            perm_effect = rank_sum_effect(perm_is_high, X_rank)[valid]
            exceed += np.abs(perm_effect) >= np.abs(observed_effect[valid])
            perm_rank_of_effect = rankdata(-np.abs(perm_effect)) / n_valid  # rough per-perm pseudo-p, null-hit-rate calibration only
            null_hit_counts[p] = int((bh_fdr(perm_rank_of_effect) < FDR_ALPHA).sum())

        empirical_p = np.full(n_features, np.nan)
        empirical_p[valid] = (exceed + 1) / (args.n_perm + 1)
        empirical_fdr = np.full(n_features, np.nan)
        empirical_fdr[valid] = bh_fdr(empirical_p[valid])

        res = pd.DataFrame({
            "dedup_group_id": dedup_ids,
            "fraction": fraction,
            "predictor": predictor_col,
            "n_high": n_high, "n_low": n_low,
            "quantile": args.quantile,
            "rank_biserial": observed_effect,
            "empirical_p": empirical_p,
            "empirical_fdr": empirical_fdr,
        })
        all_results.append(res)
        print(f"[{fraction}] BH-FDR<{FDR_ALPHA} hits (empirical): {int((empirical_fdr < FDR_ALPHA).sum())} "
              f"| null mean hits/perm: {null_hit_counts.mean():.1f} (sd {null_hit_counts.std():.1f})", file=sys.stderr)

    if not all_results:
        sys.exit("No fraction had enough matched strains -- aborting.")

    final = pd.concat(all_results, ignore_index=True)
    out_path = decoy_out if is_decoy_run else real_out
    final.to_csv(out_path, index=False)
    print(f"Wrote {out_path}", file=sys.stderr)
    if is_decoy_run:
        print("Negative-control decoy run complete -- real predictor runs are now unblocked.", file=sys.stderr)


if __name__ == "__main__":
    main()
