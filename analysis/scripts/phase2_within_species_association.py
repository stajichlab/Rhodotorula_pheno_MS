#!/usr/bin/env python3
"""
Phase 2, within-species variant (analysis/INTEGRATED_ANALYSIS_STRATEGY.md
"Species-Level Collapse" Step 6, and PHASE2_SUMMARY.md's recommended next
step after the whole-panel test came back null): does within-species
variation in color track within-species variation in the metabolome?

Why this is a different, better-powered question than the whole-panel
test in phase2_color_metabolome_association.py: that test's effective
sample size for phylogenetic correction is bounded by the number of
independent SPECIES-level lineages (~17-18), not the ~275 strains --
Phase 2's null result is consistent with that power ceiling
(.living/findings/phenotype-metabolome-association-statistical-power.md).
A single species with many strains sidesteps this almost entirely: all
strains share a recent common ancestor, so there's no small-number-of-
independent-lineages tax, and *R. mucilaginosa* (216 phenotyped strains,
206 with both color + MS data, 201 with genome data) is by far the best
candidate in this panel for this design.

This is explicitly a genuinely different test, not a rerun of Phase 2 on
a subset -- it asks "does color variation *within* one species correlate
with metabolome variation *within* that species," which is a
phylogeny-free (or nearly so) question, answerable with much higher power
here than the whole-panel test.

Predictor / negative-control / hard-gate conventions carried over
unchanged from phase2_color_metabolome_association.py (a* primary, C*
secondary, colony-area decoy trait, hard-gated negative control -- see
that script's docstring and .living/decisions.md's 2026-08-15 grilling
session entries for why).

Phylogenetic blocking here uses the ACTUAL strain-level genome tree
(BFD/results/phyling_pep/...treefile) pruned to the target species' tips
and hierarchically clustered -- not the species-level tree used by the
whole-panel script (which only has one tip per species and can't resolve
within-species structure at all). Strain<->tree-tip matching reuses the
suffix-matching convention from analysis/copper/scripts/common.py.

Usage:
    python3 analysis/scripts/phase2_within_species_association.py --species "Rhodotorula mucilaginosa" --predictor area   # negative control, run first
    python3 analysis/scripts/phase2_within_species_association.py --species "Rhodotorula mucilaginosa" --predictor a      # primary
    python3 analysis/scripts/phase2_within_species_association.py --species "Rhodotorula mucilaginosa" --predictor C      # secondary
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr

REPO = Path(__file__).resolve().parent.parent.parent
FEATURE_MATRIX = REPO / "analysis" / "linked_data" / "feature_abundance_matrix.csv.gz"
DEDUP_GROUPS = REPO / "analysis" / "linked_data" / "ms_feature_dedup_groups.csv"
SAMPLE_METADATA = REPO / "analysis" / "linked_data" / "sample_metadata.csv.gz"
PHENOTYPE_TABLE = REPO / "analysis" / "integrated_analysis" / "phase1_phenotype" / "strain_phenotype_table.csv"
STRAIN_TREE = (
    REPO
    / "BFD"
    / "results"
    / "phyling_pep"
    / "protein"
    / "buildtree"
    / "fungi_odb10"
    / "fasttree"
    / "protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile"
)
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase2_metabolome_phenotype"

PREDICTOR_COL = {"a": "a*", "C": "C*", "area": "area"}
FDR_ALPHA = 0.05


def norm(s: str) -> str:
    return re.sub(r"[_\-\s]", "", str(s)).upper()


def slugify(species: str) -> str:
    return species.replace(" ", "_").replace(".", "")


def match_to_tree(key: str, tip_norm_map: dict) -> str | None:
    """Suffix-matching convention from analysis/copper/scripts/common.py."""
    nk = norm(key)
    hits = []
    for tnorm, tip in tip_norm_map.items():
        stripped = re.sub(r"\.?PROTEINS$", "", tnorm)
        if stripped.endswith(nk):
            hits.append(tip)
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        return sorted(hits, key=len)[0]
    return None


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def strain_blocks_from_tree(strains: list[str], n_clades: int) -> dict[str, str]:
    """canonical_strain -> 'clade_<k>' using the real strain-level genome
    tree, pruned (implicitly, via distance matrix subset) to the given
    strains' tips and hierarchically clustered. Strains with no tree tip
    match are returned mapped to 'clade_unmatched' (their own pseudo-block
    -- excluded from block-restricted shuffling since a block of size 1
    contributes nothing to the permutation, same convention as
    phase2_color_metabolome_association.py's block_indices logic)."""
    from Bio import Phylo

    tree = Phylo.read(str(STRAIN_TREE), "newick")
    terminals = tree.get_terminals()
    tip_norm_map = {norm(t.name): t.name for t in terminals}

    strain_to_tip = {}
    for s in strains:
        tip = match_to_tree(s, tip_norm_map)
        if tip is not None:
            strain_to_tip[s] = tip
    n_missing = len(strains) - len(strain_to_tip)
    if n_missing:
        print(f"strain_blocks_from_tree: {n_missing}/{len(strains)} strains have no tree tip -> own singleton block", file=sys.stderr)

    matched_strains = list(strain_to_tip.keys())
    tips = [strain_to_tip[s] for s in matched_strains]
    terminal_by_name = {t.name: t for t in terminals}
    n = len(tips)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = tree.distance(terminal_by_name[tips[i]], terminal_by_name[tips[j]])
            dist[i, j] = dist[j, i] = d
    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    cluster_ids = fcluster(z, t=min(n_clades, max(n - 1, 1)), criterion="maxclust")

    blocks = {s: f"clade_{c}" for s, c in zip(matched_strains, cluster_ids)}
    for s in strains:
        if s not in blocks:
            blocks[s] = f"singleton_{s}"
    return blocks


def load_strain_fraction_matrix(fraction: str, species_strains: set[str]) -> tuple[pd.DataFrame, np.ndarray]:
    feat = pd.read_csv(FEATURE_MATRIX)
    dedup = pd.read_csv(DEDUP_GROUPS)
    dedup_rep = dedup[dedup["is_group_representative"]]
    feat = feat.merge(dedup_rep[["row ID", "dedup_group_id"]], on="row ID", how="inner")

    meta = pd.read_csv(SAMPLE_METADATA)
    meta = meta[(meta["fraction"] == fraction) & (meta["canonical_strain"].isin(species_strains))]

    sample_cols = [c for c in feat.columns if c in set(meta["sample_id"])]
    mat = feat[sample_cols].to_numpy(dtype=float)
    col_sums = mat.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    mat = mat / col_sums

    strain_by_sample = meta.set_index("sample_id")["canonical_strain"].to_dict()
    strains = [strain_by_sample[c] for c in sample_cols]
    strain_df = pd.DataFrame(mat.T, index=strains)
    strain_df = strain_df.groupby(level=0).mean()

    return strain_df, feat["dedup_group_id"].to_numpy()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True, help='e.g. "Rhodotorula mucilaginosa"')
    ap.add_argument("--predictor", choices=["a", "C", "area"], default="a")
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--n-clades", type=int, default=15, help="Number of strain-tree blocks for restricted permutation.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--min-strains", type=int, default=20,
        help="Minimum matched strains per fraction to run at all. Species below this are explicitly "
             "underpowered for formal inference -- lower this deliberately for a labeled quick/exploratory "
             "robustness check (e.g. --min-strains 8), never silently.",
    )
    args = ap.parse_args()
    if args.min_strains < 20:
        print(
            f"NOTE: --min-strains={args.min_strains} < 20 -- this run is an EXPLORATORY robustness check, "
            f"not a formally-powered test. Treat results accordingly.",
            file=sys.stderr,
        )

    predictor_col = PREDICTOR_COL[args.predictor]
    is_decoy_run = args.predictor == "area"
    slug = slugify(args.species)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    decoy_out = OUT_DIR / f"within_species_{slug}_association_area_decoy.csv"
    real_out = OUT_DIR / f"within_species_{slug}_association_{args.predictor}.csv"

    MIN_DECOY_N_PERM = 100  # a --n-perm 20 smoke test is not an adequate negative control -- see .living/learnings.md 2026-08-15 "underpowered decoy" entry

    if not is_decoy_run:
        if not decoy_out.exists():
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy output not found ({decoy_out}).\n"
                f"Run first: python3 {Path(__file__).name} --species \"{args.species}\" --predictor area"
            )
        newest_input = max(p.stat().st_mtime for p in (FEATURE_MATRIX, PHENOTYPE_TABLE, DEDUP_GROUPS) if p.exists())
        if decoy_out.stat().st_mtime < newest_input:
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy output ({decoy_out}) is STALE.\n"
                f"Rerun: python3 {Path(__file__).name} --species \"{args.species}\" --predictor area"
            )
        decoy_df = pd.read_csv(decoy_out)
        if "n_perm" not in decoy_df.columns or decoy_df["n_perm"].min() < MIN_DECOY_N_PERM:
            sys.exit(
                f"REFUSING TO RUN: negative-control decoy output ({decoy_out}) was run with too few "
                f"permutations (need >={MIN_DECOY_N_PERM}). Rerun: python3 {Path(__file__).name} "
                f"--species \"{args.species}\" --predictor area --n-perm {MIN_DECOY_N_PERM}"
            )

    pheno = pd.read_csv(PHENOTYPE_TABLE).set_index("strain_id")
    species_strains_all = set(pheno.index[pheno["species"] == args.species])
    print(f"{args.species}: {len(species_strains_all)} strains in phenotype table", file=sys.stderr)

    all_results = []
    for fraction in ["cell", "supernatant"]:
        strain_df, group_ids = load_strain_fraction_matrix(fraction, species_strains_all)

        common = [s for s in strain_df.index if s in pheno.index and pd.notna(pheno.loc[s, predictor_col])]
        if len(common) < args.min_strains:
            print(f"[{fraction}] only {len(common)} {args.species} strains with MS + {predictor_col} -- skipping", file=sys.stderr)
            continue

        y = pheno.loc[common, predictor_col].to_numpy(dtype=float)
        X = strain_df.loc[common].to_numpy(dtype=float)
        n_strains, n_features = X.shape
        blocks_map = strain_blocks_from_tree(common, args.n_clades)
        blocks = np.array([blocks_map[s] for s in common])

        print(f"[{fraction}] n_strains={n_strains} n_features={n_features} predictor={predictor_col} n_blocks={len(set(blocks))}", file=sys.stderr)

        y_rank = rankdata(y)
        X_rank = np.apply_along_axis(rankdata, 0, X)

        def spearman_vec(yr, Xr):
            yc = yr - yr.mean()
            Xc = Xr - Xr.mean(axis=0)
            num = yc @ Xc
            den = np.sqrt((yc**2).sum()) * np.sqrt((Xc**2).sum(axis=0))
            den[den == 0] = np.nan
            return num / den

        observed_rho = spearman_vec(y_rank, X_rank)
        valid = ~np.isnan(observed_rho)
        n_valid = int(valid.sum())
        if n_valid < n_features:
            print(f"[{fraction}] {n_features - n_valid} constant-abundance feature(s) excluded", file=sys.stderr)

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
            perm_rank_of_rho = rankdata(-np.abs(perm_rho)) / n_valid
            null_hit_counts[p] = int((bh_fdr(perm_rank_of_rho) < FDR_ALPHA).sum())

        empirical_p = np.full(n_features, np.nan)
        empirical_p[valid] = (exceed + 1) / (args.n_perm + 1)
        empirical_fdr = np.full(n_features, np.nan)
        empirical_fdr[valid] = bh_fdr(empirical_p[valid])

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
                "species": args.species,
                "n_strains": n_strains,
                "n_perm": args.n_perm,
                "spearman_rho": observed_rho,
                "empirical_p": empirical_p,
                "empirical_fdr": empirical_fdr,
                "asymptotic_p": asymptotic_p,
                "asymptotic_fdr": asymptotic_fdr,
            }
        )
        print(
            f"[{fraction}] BH-FDR<{FDR_ALPHA} hits (empirical): {int((empirical_fdr < FDR_ALPHA).sum())} "
            f"| null mean hits/perm: {null_hit_counts.mean():.1f} (sd {null_hit_counts.std():.1f})",
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
        print("Negative-control decoy run complete -- real predictor runs are now unblocked.", file=sys.stderr)


if __name__ == "__main__":
    main()
