#!/usr/bin/env python3
"""
Phase 2 (analysis/INTEGRATED_ANALYSIS_STRATEGY.md): whole-panel color <->
metabolome association. For each deduplicated MS2 feature group, test
whether abundance across strains correlates with color -- using the
entire strain panel (~275 strains with both phenotype and MS data) as the
sample, not a single-clade contrast.

Predictor (PI-confirmed, 2026-08-15 grilling session):
  --predictor a     PRIMARY, default. CIELAB a* (green-red axis) from the
                     canonical control_90_110 phenotype table -- the most
                     direct, mechanistically-motivated axis for
                     carotenoid-driven color.
  --predictor C     SECONDARY. Chroma.
  --predictor area  The phylogenetically-structured DECOY TRAIT for the
                     negative-control design below -- not a color axis at
                     all (colony size), used only to calibrate whether a
                     "hit" reflects color specifically or just general
                     phylogenetic structure in the metabolome.
(orange_score is demoted to exploratory-only per the same PI decision and
is not wired into this script -- rerun manually against
strain_phenotype_table.csv's orange_score column if wanted, clearly
labeled exploratory.)

Statistical design:
  - Spearman correlation (rank-based) between predictor and TSS-normalized
    feature abundance, computed on the DEDUPLICATED representative feature
    set (ms_feature_dedup_groups.csv, is_group_representative==True --
    16,332 raw features collapse to 10,949 groups, see
    analysis/scripts/dedupe_ms_features.py). BH-FDR is applied on this
    10,949-group count, not the raw 16,332.
  - Primary inferential statistic: an EMPIRICAL p-value from
    phylogenetic-block-restricted permutation (shuffle each strain's
    predictor value only among strains in the same species-tree-derived
    block, recompute correlation, repeat --n-perm times). This is
    preferred over the asymptotic Spearman p-value because with only
    ~17-18 species-level blocks, the correlation structure is
    phylogenetically non-random in a way the asymptotic null doesn't
    model -- see "Phylogenetic correction" in the strategy doc. The
    asymptotic p-value is still reported as a secondary/diagnostic column.
  - The SAME permutation replicates used for per-feature empirical
    p-values also give the negative-control (a) calibration "for free":
    the average number of BH-FDR<0.05 "hits" a random permuted dataset
    would produce is reported alongside the real hit count, without a
    separate costly outer-loop run.
  - Leave-one/few-species-out sensitivity check on the top hits (Fable
    review): recompute the naive correlation after dropping each of the
    top-3 orange species (R. taiwanensis, R. sphaerocarpa, R. glutinis,
    per Phase 1's ranking) individually and together.

HARD GATE (PI-confirmed, 2026-08-15 grilling session -- negative controls
are enforced, not just documented): running with --predictor a or
--predictor C (a "real" scientific run) refuses to write final output
unless a --predictor area run (the decoy-trait negative control) already
exists and is newer than the input data. Run the decoy first:
    python3 analysis/scripts/phase2_color_metabolome_association.py --predictor area

Usage:
    python3 analysis/scripts/phase2_color_metabolome_association.py --predictor area   # negative control, run first
    python3 analysis/scripts/phase2_color_metabolome_association.py --predictor a      # primary
    python3 analysis/scripts/phase2_color_metabolome_association.py --predictor C      # secondary
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

REPO = Path(__file__).resolve().parent.parent.parent
FEATURE_MATRIX = REPO / "analysis" / "linked_data" / "feature_abundance_matrix.csv.gz"
DEDUP_GROUPS = REPO / "analysis" / "linked_data" / "ms_feature_dedup_groups.csv"
SAMPLE_METADATA = REPO / "analysis" / "linked_data" / "sample_metadata.csv.gz"
PHENOTYPE_TABLE = REPO / "analysis" / "integrated_analysis" / "phase1_phenotype" / "strain_phenotype_table.csv"
SPECIES_TREE = REPO / "analysis" / "integrated_analysis" / "phase1_phenotype" / "species_tree.nwk"

OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase2_metabolome_phenotype"

PREDICTOR_COL = {"a": "a*", "C": "C*", "area": "area"}
TOP3_ORANGE_SPECIES = ["Rhodotorula taiwanensis", "Rhodotorula sphaerocarpa", "Rhodotorula glutinis"]
FDR_ALPHA = 0.05


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def species_blocks_from_tree(tree_path: Path, n_clades: int) -> dict[str, str]:
    """species -> 'clade_<k>', cut from the species-level tree's patristic
    distances via average-linkage hierarchical clustering. Operates
    directly on species tip names (this repo's species_tree.nwk), not
    strain-level tips -- avoids depending on the stale
    BFD/strain_coverage_summary.tsv strain<->tip mapping that
    block_permutation.py's tree path uses."""
    from Bio import Phylo

    tree = Phylo.read(str(tree_path), "newick")
    terminals = tree.get_terminals()
    names = [t.name for t in terminals]
    n = len(names)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = tree.distance(terminals[i], terminals[j])
            dist[i, j] = dist[j, i] = d
    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    cluster_ids = fcluster(z, t=n_clades, criterion="maxclust")
    # tree tip names have spaces->underscores (ape/Bio.Phylo Newick round-trip); normalize back
    return {name.replace("_", " "): f"clade_{c}" for name, c in zip(names, cluster_ids)}


def load_strain_fraction_matrix(fraction: str) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """Returns (strain-indexed metadata-ish df with predictor cols, feature
    abundance matrix [n_strains x n_features] TSS-normalized, dedup_group_id list)."""
    feat = pd.read_csv(FEATURE_MATRIX)
    dedup = pd.read_csv(DEDUP_GROUPS)
    dedup_rep = dedup[dedup["is_group_representative"]]
    feat = feat.merge(dedup_rep[["row ID", "dedup_group_id"]], on="row ID", how="inner")

    meta = pd.read_csv(SAMPLE_METADATA)
    meta = meta[meta["fraction"] == fraction]

    sample_cols = [c for c in feat.columns if c in set(meta["sample_id"])]
    mat = feat[sample_cols].to_numpy(dtype=float)
    # TSS normalize per sample (column)
    col_sums = mat.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    mat = mat / col_sums

    strain_by_sample = meta.set_index("sample_id")["canonical_strain"].to_dict()
    strains = [strain_by_sample[c] for c in sample_cols]

    # average replicate samples of the same strain within this fraction (rare, 2 strains x2)
    strain_df = pd.DataFrame(mat.T, index=strains)
    strain_df = strain_df.groupby(level=0).mean()

    return strain_df, feat["dedup_group_id"].to_numpy(), sample_cols


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictor", choices=["a", "C", "area"], default="a")
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--n-clades", type=int, default=6, help="Number of species-tree blocks for restricted permutation.")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    predictor_col = PREDICTOR_COL[args.predictor]
    is_decoy_run = args.predictor == "area"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    decoy_out = OUT_DIR / "color_metabolome_association_area_decoy.csv"
    real_out = OUT_DIR / f"color_metabolome_association_{args.predictor}.csv"

    # --- hard gate: real predictor runs require a fresh decoy run first ---
    if not is_decoy_run:
        if not decoy_out.exists():
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy-trait output not found ({decoy_out}).\n"
                f"Run the negative control first:\n"
                f"  python3 {Path(__file__).name} --predictor area\n"
                f"(PI-confirmed hard gate, 2026-08-15 grilling session -- see .living/decisions.md)"
            )
        newest_input = max(p.stat().st_mtime for p in (FEATURE_MATRIX, PHENOTYPE_TABLE, DEDUP_GROUPS) if p.exists())
        if decoy_out.stat().st_mtime < newest_input:
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy output ({decoy_out}) is STALE "
                f"(older than input data). Rerun the negative control first:\n"
                f"  python3 {Path(__file__).name} --predictor area"
            )

    pheno = pd.read_csv(PHENOTYPE_TABLE).set_index("strain_id")
    species_by_strain = pheno["species"]
    blocks_by_species = species_blocks_from_tree(SPECIES_TREE, args.n_clades)

    all_results = []
    for fraction in ["cell", "supernatant"]:
        strain_df, group_ids, _ = load_strain_fraction_matrix(fraction)

        common = strain_df.index.intersection(pheno.index)
        common = [s for s in common if pd.notna(pheno.loc[s, predictor_col])]
        if len(common) < 10:
            print(f"[{fraction}] only {len(common)} strains with both MS and {predictor_col} data -- skipping", file=sys.stderr)
            continue

        y = pheno.loc[common, predictor_col].to_numpy(dtype=float)
        X = strain_df.loc[common].to_numpy(dtype=float)  # n_strains x n_features
        n_strains, n_features = X.shape
        species = species_by_strain.loc[common].to_numpy()
        blocks = np.array([blocks_by_species.get(sp, "clade_unknown") for sp in species])

        print(f"[{fraction}] n_strains={n_strains} n_features={n_features} predictor={predictor_col}", file=sys.stderr)

        y_rank = rankdata(y)
        X_rank = np.apply_along_axis(rankdata, 0, X)  # rank each feature column across strains

        def spearman_vec(yr, Xr):
            yc = yr - yr.mean()
            Xc = Xr - Xr.mean(axis=0)
            num = yc @ Xc
            den = np.sqrt((yc**2).sum()) * np.sqrt((Xc**2).sum(axis=0))
            den[den == 0] = np.nan
            return num / den

        observed_rho = spearman_vec(y_rank, X_rank)
        # constant-abundance features (zero variance across strains) give a
        # NaN rho (0/0). NaN comparisons (e.g. `nan >= x`) silently evaluate
        # False in numpy, which would make these features falsely "never
        # exceeded" by any permutation -> a spuriously tiny empirical p.
        # Exclude them entirely from the permutation test and BH-FDR instead.
        valid = ~np.isnan(observed_rho)
        n_valid = int(valid.sum())
        if n_valid < n_features:
            print(f"[{fraction}] {n_features - n_valid} constant-abundance feature(s) excluded from testing", file=sys.stderr)

        rng = np.random.default_rng(args.seed)
        block_indices = [np.where(blocks == b)[0] for b in np.unique(blocks)]
        exceed = np.zeros(n_valid, dtype=np.int64)
        null_hit_counts = np.zeros(args.n_perm, dtype=np.int64)
        for p in range(args.n_perm):
            perm_rank = y_rank.copy()
            for idx in block_indices:
                if len(idx) > 1:
                    perm_rank[idx] = rng.permutation(y_rank[idx])
            perm_rho = spearman_vec(perm_rank, X_rank)[valid]
            exceed += np.abs(perm_rho) >= np.abs(observed_rho[valid])
            perm_rank_of_rho = rankdata(-np.abs(perm_rho)) / n_valid  # rough per-perm rank-based p, for null-hit-count calibration only
            null_hit_counts[p] = int((bh_fdr(perm_rank_of_rho) < FDR_ALPHA).sum())

        empirical_p = np.full(n_features, np.nan)
        empirical_p[valid] = (exceed + 1) / (args.n_perm + 1)
        empirical_fdr = np.full(n_features, np.nan)
        empirical_fdr[valid] = bh_fdr(empirical_p[valid])

        from scipy.stats import spearmanr

        asymptotic_p = np.full(n_features, np.nan)
        for i in np.where(valid)[0]:
            asymptotic_p[i] = spearmanr(y, X[:, i]).pvalue
        asymptotic_fdr = np.full(n_features, np.nan)
        asymptotic_fdr[valid] = bh_fdr(asymptotic_p[valid])

        res = pd.DataFrame(
            {
                "dedup_group_id": group_ids,
                "fraction": fraction,
                "predictor": predictor_col,
                "n_strains": n_strains,
                "spearman_rho": observed_rho,
                "empirical_p": empirical_p,
                "empirical_fdr": empirical_fdr,
                "asymptotic_p": asymptotic_p,
                "asymptotic_fdr": asymptotic_fdr,
            }
        )

        # leave-one/few-species-out sensitivity on top hits by empirical_fdr
        top_hits = res.sort_values("empirical_p").head(20)["dedup_group_id"].tolist()
        sens_cols = {}
        drop_sets = {sp: [sp] for sp in TOP3_ORANGE_SPECIES}
        drop_sets["ALL_TOP3"] = TOP3_ORANGE_SPECIES
        for label, drop_species in drop_sets.items():
            keep_mask = ~np.isin(species, drop_species)
            if keep_mask.sum() < 10:
                continue
            yk = y_rank[keep_mask] if False else rankdata(y[keep_mask])
            col_name = f"rho_drop_{label.replace(' ', '_')}"
            vals = {}
            for gid in top_hits:
                fi = np.where(group_ids == gid)[0]
                if len(fi) == 0:
                    continue
                fi = fi[0]
                xk = rankdata(X[keep_mask, fi])
                if len(np.unique(xk)) < 2 or len(np.unique(yk)) < 2:
                    vals[gid] = np.nan
                    continue
                vals[gid] = np.corrcoef(yk, xk)[0, 1]
            sens_cols[col_name] = vals
        for col_name, vals in sens_cols.items():
            res[col_name] = res["dedup_group_id"].map(vals)

        print(
            f"[{fraction}] observed BH-FDR<{FDR_ALPHA} hits (empirical): {int((empirical_fdr < FDR_ALPHA).sum())} "
            f"| null mean hits per permutation: {null_hit_counts.mean():.1f} (sd {null_hit_counts.std():.1f})",
            file=sys.stderr,
        )

        all_results.append(res)

    if not all_results:
        sys.exit("No fraction had enough matched strains -- aborting.")

    final = pd.concat(all_results, ignore_index=True)
    out_path = decoy_out if is_decoy_run else real_out
    final.to_csv(out_path, index=False)
    print(f"Wrote {out_path}", file=sys.stderr)
    if is_decoy_run:
        print("Negative-control decoy run complete -- real predictor runs (--predictor a / C) are now unblocked.", file=sys.stderr)


if __name__ == "__main__":
    main()
