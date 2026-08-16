"""Idea 1 follow-up: naive (uncorrected) Spearman quick-check of targeted
mass-search candidates against copper-resistance growth rate (mean_auc_rate).

Not a validated pipeline -- no phylogenetic block permutation, no negative
control, no BH-FDR. Purpose is only to screen whether a candidate is worth
promoting to the full rigor (phase2_within_species_association.py-style) test.
Uses the same TSS-normalization / dedup-group-representative convention as
phase2_color_metabolome_association.py's load_strain_fraction_matrix().

Usage:
    python3 analysis/scripts/idea1_auc_quickcheck.py --row-id 846 9852 6682 35014
"""
import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

LINKED = Path("analysis/linked_data")
FEATURE_MATRIX = LINKED / "feature_abundance_matrix.csv.gz"
DEDUP_GROUPS = LINKED / "ms_feature_dedup_groups.csv"
SAMPLE_METADATA = LINKED / "sample_metadata.csv.gz"


def load_strain_fraction_matrix(fraction: str, feat: pd.DataFrame, dedup: pd.DataFrame) -> pd.DataFrame:
    dedup_rep = dedup[dedup["is_group_representative"]]
    f = feat.merge(dedup_rep[["row ID", "dedup_group_id"]], on="row ID", how="inner")

    meta = pd.read_csv(SAMPLE_METADATA)
    meta = meta[meta["fraction"] == fraction]

    sample_cols = [c for c in f.columns if c in set(meta["sample_id"])]
    mat = f[sample_cols].to_numpy(dtype=float)
    col_sums = mat.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    mat = mat / col_sums

    strain_by_sample = meta.set_index("sample_id")["canonical_strain"].to_dict()
    strains = [strain_by_sample[c] for c in sample_cols]

    strain_df = pd.DataFrame(mat.T, index=strains, columns=f["dedup_group_id"].to_numpy())
    strain_df = strain_df.groupby(level=0).mean()
    strain_df = strain_df.groupby(level=0, axis=1).first()  # collapse duplicate dedup_group_id cols if any
    return strain_df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--row-id", type=int, nargs="+", required=True,
                     help="Raw feature row IDs to test (need not be group representatives; "
                          "resolved to their dedup group's representative automatically).")
    args = ap.parse_args()

    feat = pd.read_csv(FEATURE_MATRIX)
    dedup = pd.read_csv(DEDUP_GROUPS)
    meta = pd.read_csv(SAMPLE_METADATA)
    auc_by_strain = (
        meta.dropna(subset=["mean_auc_rate"])
        .groupby("canonical_strain")["mean_auc_rate"].mean()
    )

    # resolve each requested row ID to its dedup_group_id
    group_of = dedup.set_index("row ID")["dedup_group_id"].to_dict()
    rep_of_group = dedup[dedup["is_group_representative"]].set_index("dedup_group_id")["row ID"].to_dict()

    rows = []
    for rid in args.row_id:
        if rid not in group_of:
            print(f"row {rid}: not found in dedup groups, skipping")
            continue
        gid = group_of[rid]
        rep_rid = rep_of_group.get(gid, rid)
        for fraction in ["cell", "supernatant"]:
            strain_df = load_strain_fraction_matrix(fraction, feat, dedup)
            if gid not in strain_df.columns:
                print(f"row {rid} (group {gid}, rep {rep_rid}): group not present in {fraction} fraction, skipping")
                continue
            joined = pd.DataFrame({"abund": strain_df[gid]}).join(auc_by_strain, how="inner").dropna()
            rho, p = spearmanr(joined["abund"], joined["mean_auc_rate"])
            rows.append({
                "row_id": rid, "dedup_group_id": gid, "rep_row_id": rep_rid,
                "fraction": fraction, "n": len(joined), "auc_spearman_rho": rho, "auc_p_naive": p,
            })
            print(f"row {rid} (group {gid}, rep {rep_rid}) [{fraction}]: "
                  f"n={len(joined)} rho={rho:.3f} p={p:.4g}")

    out = pd.DataFrame(rows)
    out_path = Path("analysis/integrated_analysis/phase3_metabolome_phenotype_idea1/sterol_cluster_auc_quickcheck.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
