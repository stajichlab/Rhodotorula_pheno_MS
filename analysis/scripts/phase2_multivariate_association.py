#!/usr/bin/env python3
"""
Phase 2, multivariate follow-up (PHASE2_SUMMARY.md / WITHIN_SPECIES_MUCILAGINOSA.md
recommended next step #2): per-feature univariate tests (both the
whole-panel scan and the R. mucilaginosa within-species follow-up) found
zero significant color-metabolome associations. A per-feature test with
FDR correction is specifically bad at detecting a signal that is real but
DIFFUSE -- spread thinly and weakly across many features rather than
concentrated in one or two strong ones. This script asks a different
question: can color be predicted from the metabolome at all, using a
sparse multivariate model that can combine many weak features?

Method: Lasso regression (L1-regularized linear regression -- performs
its own feature selection, appropriate for p >> n) predicting the color
predictor from all deduplicated feature-group abundances jointly.
  1. Fit LassoCV with GroupKFold(5) (groups = phylogenetic blocks from the
     real strain-level genome tree, same blocking as
     phase2_within_species_association.py -- prevents closely related
     strains from leaking between train/test folds) to pick alpha and
     get real out-of-fold predictions -> observed cross-validated R^2
     (Q^2).
  2. NEGATIVE CONTROL / permutation test (same hard-gate convention as the
     univariate scripts): shuffle the color predictor within the same
     phylogenetic blocks, refit with the SAME alpha and SAME fold
     assignment (not researched per permutation -- keeps this
     computationally tractable), recompute Q^2, repeat --n-perm times.
     Empirical p-value = fraction of permutations with Q^2 >= observed.
  3. Report which features survive in the final model (nonzero Lasso
     coefficients, fit on the full data at the chosen alpha) as the
     candidate multivariate signature, IF the permutation test says the
     overall model is doing better than chance -- do not report/interpret
     individual coefficients if the omnibus test is not significant, that
     would be double-dipping on a null result.

HARD GATE (same convention as the univariate scripts): --predictor a or C
refuse to run final output without a fresh --predictor area decoy run.

Usage:
    python3 analysis/scripts/phase2_multivariate_association.py --species "Rhodotorula mucilaginosa" --predictor area   # negative control, run first
    python3 analysis/scripts/phase2_multivariate_association.py --species "Rhodotorula mucilaginosa" --predictor a
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LassoCV
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase2_within_species_association import (  # noqa: E402
    DEDUP_GROUPS,
    FEATURE_MATRIX,
    PHENOTYPE_TABLE,
    PREDICTOR_COL,
    load_strain_fraction_matrix,
    slugify,
    strain_blocks_from_tree,
)

REPO = Path(__file__).resolve().parent.parent.parent
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase2_metabolome_phenotype"


def cv_r2(y: np.ndarray, X: np.ndarray, groups: np.ndarray, alpha: float, n_splits: int) -> float:
    n_splits_eff = min(n_splits, len(np.unique(groups)))
    if n_splits_eff < 2:
        return np.nan
    gkf = GroupKFold(n_splits=n_splits_eff)
    preds = cross_val_predict(Lasso(alpha=alpha, max_iter=20000), X, y, cv=gkf, groups=groups)
    ss_res = np.sum((y - preds) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--predictor", choices=["a", "C", "area"], default="a")
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-clades", type=int, default=15)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    predictor_col = PREDICTOR_COL[args.predictor]
    is_decoy_run = args.predictor == "area"
    slug = slugify(args.species)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    decoy_out = OUT_DIR / f"multivariate_{slug}_area_decoy.csv"
    real_out = OUT_DIR / f"multivariate_{slug}_{args.predictor}.csv"

    if not is_decoy_run:
        if not decoy_out.exists():
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy output not found ({decoy_out}).\n"
                f"Run first: python3 {Path(__file__).name} --species \"{args.species}\" --predictor area"
            )
        newest_input = max(p.stat().st_mtime for p in (FEATURE_MATRIX, PHENOTYPE_TABLE, DEDUP_GROUPS) if p.exists())
        if decoy_out.stat().st_mtime < newest_input:
            sys.exit(f"REFUSING TO RUN: decoy output ({decoy_out}) is STALE. Rerun the negative control first.")

    pheno = pd.read_csv(PHENOTYPE_TABLE).set_index("strain_id")
    species_strains_all = set(pheno.index[pheno["species"] == args.species])
    print(f"{args.species}: {len(species_strains_all)} strains in phenotype table", file=sys.stderr)

    rows = []
    for fraction in ["cell", "supernatant"]:
        strain_df, group_ids = load_strain_fraction_matrix(fraction, species_strains_all)
        common = [s for s in strain_df.index if s in pheno.index and pd.notna(pheno.loc[s, predictor_col])]
        if len(common) < 30:
            print(f"[{fraction}] only {len(common)} strains -- skipping (need >=30 for a stable multivariate CV)", file=sys.stderr)
            continue

        y = pheno.loc[common, predictor_col].to_numpy(dtype=float)
        X = strain_df.loc[common].to_numpy(dtype=float)
        # drop constant features (zero variance -> no information, and breaks standardization)
        keep = X.std(axis=0) > 0
        X = X[:, keep]
        kept_group_ids = group_ids[keep]
        X = StandardScaler().fit_transform(X)

        blocks_map = strain_blocks_from_tree(common, args.n_clades)
        groups = np.array([blocks_map[s] for s in common])
        n_groups = len(np.unique(groups))
        print(f"[{fraction}] n_strains={len(common)} n_features={X.shape[1]} n_blocks={n_groups}", file=sys.stderr)

        # Pick alpha via a manual group-aware CV grid search (LassoCV's built-in
        # `groups=` routing to GroupKFold is not supported in this sklearn
        # version without enabling metadata routing globally -- doing the grid
        # search by hand avoids that and keeps the CV split explicitly
        # phylogeny-aware, which is the whole point).
        alpha_grid = np.logspace(-3, 1, 15) * np.std(y)
        alpha_scores = [cv_r2(y, X, groups, a, args.n_splits) for a in alpha_grid]
        alpha = alpha_grid[int(np.nanargmax(alpha_scores))]
        observed_r2 = cv_r2(y, X, groups, alpha, args.n_splits)

        rng = np.random.default_rng(args.seed)
        block_indices = [np.where(groups == g)[0] for g in np.unique(groups)]
        null_r2 = np.zeros(args.n_perm)
        for p in range(args.n_perm):
            y_perm = y.copy()
            for idx in block_indices:
                if len(idx) > 1:
                    y_perm[idx] = rng.permutation(y[idx])
            null_r2[p] = cv_r2(y_perm, X, groups, alpha, args.n_splits)

        valid_null = null_r2[~np.isnan(null_r2)]
        empirical_p = (np.sum(valid_null >= observed_r2) + 1) / (len(valid_null) + 1)

        # final full-data fit at the chosen alpha, for reporting selected features
        final_model = Lasso(alpha=alpha, max_iter=20000).fit(X, y)
        n_nonzero = int(np.sum(final_model.coef_ != 0))
        top_idx = np.argsort(-np.abs(final_model.coef_))[:20]
        top_features = kept_group_ids[top_idx][final_model.coef_[top_idx] != 0]
        top_coefs = final_model.coef_[top_idx][final_model.coef_[top_idx] != 0]

        print(
            f"[{fraction}] alpha={alpha:.4g} observed Q2={observed_r2:.4f} "
            f"null Q2 mean={valid_null.mean():.4f} (sd {valid_null.std():.4f}) "
            f"empirical_p={empirical_p:.4f} | {n_nonzero} nonzero coefficients",
            file=sys.stderr,
        )

        rows.append(
            dict(
                species=args.species, fraction=fraction, predictor=predictor_col,
                n_strains=len(common), n_features=X.shape[1], alpha=alpha,
                observed_q2=observed_r2, null_q2_mean=valid_null.mean(), null_q2_sd=valid_null.std(),
                empirical_p=empirical_p, n_nonzero_coefficients=n_nonzero,
                top_features=";".join(f"{g}:{c:.4g}" for g, c in zip(top_features, top_coefs)),
            )
        )

    if not rows:
        sys.exit("No fraction had enough matched strains -- aborting.")

    res = pd.DataFrame(rows)
    out_path = decoy_out if is_decoy_run else real_out
    res.to_csv(out_path, index=False)
    print(f"Wrote {out_path}", file=sys.stderr)
    if is_decoy_run:
        print("Negative-control decoy run complete -- real predictor runs are now unblocked.", file=sys.stderr)


if __name__ == "__main__":
    main()
