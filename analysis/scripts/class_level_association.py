#!/usr/bin/env python3
"""
Class-level metabolome <-> phenotype association (strategy 1 + 2 + 3 of the
"combine features by class" plan; see analysis/class_level_aggregation/).

Motivation: the project's feature-level analyses (Phase 2 and its follow-ups)
test ~10,949 deduplicated features individually. When features are correlated
within a chemical class and each feature's effect is small, no individual
feature survives BH-FDR even though the class as a whole is informative. This
script collapses features to classes and re-tests the phenotype at class
resolution, cutting the multiple-testing load 18- to 1,500-fold.

Pipeline (all three strategies implemented):
  Strategy 1 -- collapse redundant features. Reuses the existing
      dedup grouping (analysis/linked_data/ms_feature_dedup_groups.csv:
      16,332 raw features -> 10,949 adduct/isotopologue/ISF-collapsed
      groups). Only the representative feature of each group enters the
      per-class aggregation below, so a compound that fails to appear in a
      class as several redundant peaks is counted once.
  Strategy 2 -- class-summary aggregation. For each class (NPC pathway /
      NPC class / ClassyFire class), build a per-strain score = mean of the
      per-feature z-scores (features standardized across strains within a
      fraction, so each member contributes equal weight -- PI-confirmed).
      Test Spearman(class score, phenotype) with phylogenetically block-
      restricted permutation + BH-FDR across classes.
  Strategy 3 -- set-based enrichment. Rank the annotated features by their
      individual association with the phenotype (signed Spearman rho) and
      ask whether a class's members sit systematically at the top of that
      ranking. Primary statistic: GSEA-style Kolmogorov-Smirnov enrichment
      score (ES); also reported: class mean signed rho. Both tested with the
      same block-restricted permutation null, BH-FDR across classes.

Phenotypes (PI-confirmed, same vocabulary as Phase 2):
  --predictor a     PRIMARY, default. CIELAB a* (green-red)
  --predictor C     SECONDARY. Chroma
  --predictor area  The phylogenetically-structured DECOY TRAIT (colony
                    area), the negative-control calibrator.

HARD GATE (same as Phase 2): --predictor a / C refuse to run unless a fresh
area (decoy) run already exists and is newer than the input data.

Class ontologies: NPC pathway (7 classes), NPC class (181), ClassyFire class
(66), all reported in one output. IMPORTANT annotation-coverage caveat: only
~29% of the 10,949 dedup groups carry a SIRIUS CANOPUS class call, so both
strategies operate on the annotated subset (3,217 groups for NPC, 3,131 for
ClassyFire). Unannotated features are dropped from both the aggregation and
the enrichment background -- enrichment is conditional on annotation.

Usage:
    python3 analysis/scripts/class_level_association.py --predictor area   # negative control first
    python3 analysis/scripts/class_level_association.py --predictor a      # primary
    python3 analysis/scripts/class_level_association.py --predictor C      # secondary
"""
from __future__ import annotations

import argparse
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
SPECIES_TREE = REPO / "analysis" / "integrated_analysis" / "phase1_phenotype" / "species_tree.nwk"
SIRIUS_ANNOT = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"

OUT_DIR = REPO / "analysis" / "class_level_aggregation" / "outputs"

PREDICTOR_COL = {"a": "a*", "C": "C*", "area": "area"}
ONTOLOGY_COLS = {
    "npc_pathway": "sirius_npc_pathway",
    "npc_class": "sirius_npc_class",
    "classyfire_class": "sirius_classyfire_class",
}
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
    """species -> 'clade_<k>', cut from the species-level tree (same logic as
    Phase 2's phase2_color_metabolome_association.py)."""
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
    return {name.replace("_", " "): f"clade_{c}" for name, c in zip(names, cluster_ids)}


def _centered_unit(a: np.ndarray) -> np.ndarray:
    """Columnwise center and L2-normalize (for Pearson-on-ranks = Spearman)."""
    a = a - a.mean(axis=0)
    norm = np.sqrt((a**2).sum(axis=0))
    norm[norm == 0] = np.nan
    return a / norm


def build_matrix(fraction: str, ontology_col: str) -> tuple[pd.DataFrame, np.ndarray, dict]:
    """Return (strain-indexed z-scored abundance df of annotated reps,
    1/0 class-membership array [n_features x n_classes], class_id list)."""
    feat = pd.read_csv(FEATURE_MATRIX)
    dedup = pd.read_csv(DEDUP_GROUPS)
    annot = pd.read_csv(SIRIUS_ANNOT, sep="\t")
    meta = pd.read_csv(SAMPLE_METADATA)

    rep = dedup[dedup["is_group_representative"]]
    feat = feat.merge(rep[["row ID", "dedup_group_id"]], on="row ID", how="inner")
    feat = feat.merge(annot[["row ID", ontology_col]], on="row ID", how="left")
    feat = feat[feat[ontology_col].notna()]  # annotated subset (the exposed universe)

    meta = meta[meta["fraction"] == fraction]
    sample_cols = [c for c in feat.columns if c in set(meta["sample_id"])]

    mat = feat[sample_cols].to_numpy(dtype=float)
    col_sums = mat.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    mat = mat / col_sums  # TSS per sample

    strain_by_sample = meta.set_index("sample_id")["canonical_strain"].to_dict()
    strains = [strain_by_sample[c] for c in sample_cols]
    strain_df = pd.DataFrame(mat.T, index=strains)
    strain_df = strain_df.groupby(level=0).mean()  # collapse replicate samples

    z = (strain_df - strain_df.mean()) / strain_df.std(ddof=0)
    names = feat[ontology_col].tolist()
    class_ids = np.unique(names)
    members = np.stack([np.array(names) == c for c in class_ids]).astype(float)  # n_class x n_feat

    return z, members, {i: c for i, c in enumerate(class_ids)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictor", choices=["a", "C", "area"], default="a")
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--n-clades", type=int, default=6, help="Species-tree blocks for restricted permutation.")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    predictor_col = PREDICTOR_COL[args.predictor]
    is_decoy_run = args.predictor == "area"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    decoy_out = OUT_DIR / f"class_association_area.csv"
    agg_out = OUT_DIR / f"class_association_{args.predictor}.csv"
    enr_out = OUT_DIR / f"class_enrichment_{args.predictor}.csv"

    if not is_decoy_run:
        for gate in (decoy_out,):
            if not gate.exists():
                sys.exit(f"REFUSING TO RUN: negative-control decoy output not found ({gate}).\nRun first:\n  python3 {Path(__file__).name} --predictor area")
        newest_input = max(p.stat().st_mtime for p in (FEATURE_MATRIX, PHENOTYPE_TABLE, DEDUP_GROUPS, SIRIUS_ANNOT) if p.exists())
        if decoy_out.stat().st_mtime < newest_input:
            sys.exit(f"REFUSING TO RUN: decoy output ({decoy_out}) is STALE. Rerun:\n  python3 {Path(__file__).name} --predictor area")

    pheno = pd.read_csv(PHENOTYPE_TABLE).set_index("strain_id")
    species_by_strain = pheno["species"]
    blocks_by_species = species_blocks_from_tree(SPECIES_TREE, args.n_clades)

    agg_rows, enr_rows = [], []
    for fraction in ["cell", "supernatant"]:
        for ont_name, ont_col in ONTOLOGY_COLS.items():
            z, members, id2class = build_matrix(fraction, ont_col)
            n_features, n_classes = z.shape[1], len(id2class)

            common = [s for s in z.index if s in pheno.index and pd.notna(pheno.loc[s, predictor_col])]
            if len(common) < 10:
                print(f"[{fraction}/{ont_name}] only {len(common)} matched strains -- skipping", file=sys.stderr)
                continue

            zt = z.loc[common]
            y_rank = rankdata(pheno.loc[common, predictor_col].to_numpy(dtype=float))
            species = species_by_strain.loc[common].to_numpy()
            blocks = np.array([blocks_by_species.get(sp, "clade_unknown") for sp in species])
            n_strains = len(common)

            za = zt.to_numpy(dtype=float)
            have_value = ~np.isnan(za).any(axis=0)
            za = za[:, have_value]
            members_hv = members[:, have_value]  # n_class x n_feat (non-NA features)
            col_span = np.nanmax(za, axis=0) - np.nanmin(za, axis=0)
            scale = np.maximum.reduce([np.abs(np.nanmin(za, axis=0)), np.abs(np.nanmax(za, axis=0)), np.full(za.shape[1], 1.0)])
            valid_feat = col_span > np.finfo(float).eps * scale * 100
            za = za[:, valid_feat]
            members_v = members_hv[:, valid_feat]  # n_class x n_feat, only usable features counted
            n_members = members_v.sum(axis=1).astype(int)
            keep_class = n_members >= 1
            n_classes_use = int(keep_class.sum())

            # --- strategy 2: class score = mean member z; Spearman vs phenotype ---
            score_rank = rankdata(za @ members_v[keep_class].T / n_members[keep_class], axis=0)  # n_strain x n_class
            yc = _centered_unit(y_rank[:, None])
            Xc = _centered_unit(score_rank.astype(float))
            obs_agg = (Xc.T @ yc).ravel()

            # feature-level signed rho (ranking statistic for strategy 3)
            Xf = _centered_unit(rankdata(za, axis=0).astype(float))
            obs_rho = (Xf.T @ yc).ravel()

            rng = np.random.default_rng(args.seed)
            block_indices = [np.where(blocks == b)[0] for b in np.unique(blocks)]
            perm_yc = np.zeros((args.n_perm, n_strains))
            for p in range(args.n_perm):
                perm_rank = y_rank.copy()
                for idx in block_indices:
                    if len(idx) > 1:
                        perm_rank[idx] = rng.permutation(perm_rank[idx])
                perm_yc[p] = _centered_unit(perm_rank[:, None]).ravel()
            n_perm = args.n_perm

            perm_agg = (Xc.T @ perm_yc.T)  # n_class x n_perm
            agg_exceed = (np.abs(perm_agg) >= np.abs(obs_agg)[:, None]).sum(axis=1)
            agg_p = (agg_exceed + 1) / (n_perm + 1)

            perm_rho = (Xf.T @ perm_yc.T)  # n_feat x n_perm
            nobs = np.where(keep_class)[0]
            mem_use = members_v[keep_class]  # n_class x n_feat

            def es_over_ordering(order_desc: np.ndarray, mem: np.ndarray) -> np.ndarray:
                nm = mem.sum(axis=1)
                n_all = order_desc.size
                pref = np.cumsum(mem[:, order_desc], axis=1)  # n_class x n_feat
                pos = np.arange(n_all)
                run = pref / nm[:, None] - (pos[None, :] - pref) / (n_all - nm)[:, None]
                return np.nanmax(np.abs(run), axis=1)

            enr_es_obs = es_over_ordering(np.argsort(obs_rho)[::-1], mem_use)
            enr_es_perm = np.zeros((len(nobs), n_perm))
            for p in range(n_perm):
                enr_es_perm[:, p] = es_over_ordering(np.argsort(perm_rho[:, p])[::-1], mem_use)
            enr_mean_obs = (mem_use @ obs_rho) / n_members[keep_class]
            enr_mean_perm = (mem_use @ perm_rho) / n_members[keep_class][:, None]

            enr_mean_perm_abs = np.abs(enr_mean_perm)
            enr_mean_exceed = (enr_mean_perm_abs >= np.abs(enr_mean_obs)[:, None]).sum(axis=1)
            enr_mean_p = (enr_mean_exceed + 1) / (n_perm + 1)

            enr_es_exceed = (enr_es_perm >= enr_es_obs[:, None]).sum(axis=1)
            enr_es_p = (enr_es_exceed + 1) / (n_perm + 1)

            # BH-FDR per level/fraction
            agg_fdr = np.full(n_classes_use, np.nan)
            enr_mean_fdr = np.full(n_classes_use, np.nan)
            enr_es_fdr = np.full(n_classes_use, np.nan)
            if n_classes_use:
                agg_fdr = bh_fdr(agg_p)
                enr_mean_fdr = bh_fdr(enr_mean_p)
                enr_es_fdr = bh_fdr(enr_es_p)

            class_ids_use = [id2class[i] for i in np.where(keep_class)[0]]

            # asymptotic p for strategy 2 (diagnostic)
            for k in range(n_classes_use):
                agg_rows.append(
                    {
                        "ont_level": ont_name,
                        "class_id": class_ids_use[k],
                        "fraction": fraction,
                        "predictor": predictor_col,
                        "n_members": int(n_members[np.where(keep_class)[0][k]]),
                        "spearman_rho": float(obs_agg[k]),
                        "empirical_p": float(agg_p[k]),
                        "empirical_fdr": float(agg_fdr[k]),
                        "n_strains": n_strains,
                    }
                )
            for k in range(n_classes_use):
                enr_rows.append(
                    {
                        "ont_level": ont_name,
                        "class_id": class_ids_use[k],
                        "fraction": fraction,
                        "predictor": predictor_col,
                        "n_members": int(n_members[np.where(keep_class)[0][k]]),
                        "mean_signed_rho": float(enr_mean_obs[k]),
                        "mean_signed_rho_empirical_p": float(enr_mean_p[k]),
                        "mean_signed_rho_empirical_fdr": float(enr_mean_fdr[k]),
                        "gsea_es": float(enr_es_obs[k]),
                        "gsea_es_empirical_p": float(enr_es_p[k]),
                        "gsea_es_empirical_fdr": float(enr_es_fdr[k]),
                    }
                )

            n_hits_agg = int((agg_fdr < FDR_ALPHA).sum())
            n_hits_mean = int((enr_mean_fdr < FDR_ALPHA).sum())
            n_hits_es = int((enr_es_fdr < FDR_ALPHA).sum())
            null_agg = (agg_fdr[np.newaxis, :] < FDR_ALPHA).sum() / 1 if False else None
            print(
                f"[{fraction}/{ont_name}] {n_classes_use} classes, {n_features} annotated features, n={n_strains} "
                f"| Str2 hits(FDR<0.05): {n_hits_agg}/{n_classes_use} | Str3 mean-rho hits: {n_hits_mean} | ES hits: {n_hits_es}",
                file=sys.stderr,
            )

    if not agg_rows:
        sys.exit("No class-level results produced -- aborting.")

    pd.DataFrame(agg_rows).to_csv(agg_out, index=False)
    pd.DataFrame(enr_rows).to_csv(enr_out, index=False)
    print(f"Wrote {agg_out} (strategy 2, class aggregation)", file=sys.stderr)
    print(f"Wrote {enr_out} (strategy 3, enrichment)", file=sys.stderr)


if __name__ == "__main__":
    main()
