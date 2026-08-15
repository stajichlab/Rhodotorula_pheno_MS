"""Join Cu_AUC copper-resistance phenotype to per-strain proteome-wide
amino-acid composition, restricted to strains present in the phyling
protein tree.

Output: analysis/copper/outputs/copper_aa_master_table.csv
        one row per strain with mean_auc_rate, 20 AA frequencies, and the
        matched tree tip label (join key for downstream PGLS).
Also writes analysis/copper/outputs/join_diagnostics.txt summarizing how
many Cu_AUC strains were dropped for lacking genome data, and why.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from common import AA20, REPO, load_aa_freq_table, load_cu_auc, load_tree_tips, match_to_tree, norm

OUT_DIR = os.path.join(REPO, "analysis/copper/outputs")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    tips = load_tree_tips()
    tip_norm_map = {norm(t): t for t in tips}

    aa_table = load_aa_freq_table()
    # map aa_freq prefix -> tree tip
    aa_prefix_to_tip = {}
    for prefix in aa_table:
        tip = match_to_tree(prefix, tip_norm_map)
        if tip is not None:
            aa_prefix_to_tip[prefix] = tip
    tip_to_aa_prefix = {v: k for k, v in aa_prefix_to_tip.items()}

    cu_rows = load_cu_auc()

    matched = []
    unmatched_no_tip = []
    unmatched_no_aa = []
    unmatched_bad_auc = []
    for row in cu_rows:
        sample = row["SAMPLE_NAME"]
        tip = match_to_tree(sample, tip_norm_map)
        if tip is None:
            unmatched_no_tip.append(sample)
            continue
        aa_prefix = tip_to_aa_prefix.get(tip)
        if aa_prefix is None:
            unmatched_no_aa.append(sample)
            continue
        freqs = aa_table[aa_prefix]
        try:
            auc = float(row["mean_auc_rate"])
        except (TypeError, ValueError):
            unmatched_bad_auc.append(sample)
            continue
        out = {
            "sample_name": sample,
            "species": row["SPECIES"],
            "tree_tip": tip,
            "aa_freq_prefix": aa_prefix,
            "mean_auc_rate": auc,
        }
        for aa in AA20:
            out[aa] = freqs.get(aa, "")
        matched.append(out)

    out_csv = os.path.join(OUT_DIR, "copper_aa_master_table.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample_name", "species", "tree_tip", "aa_freq_prefix", "mean_auc_rate"] + AA20)
        w.writeheader()
        w.writerows(matched)

    diag = os.path.join(OUT_DIR, "join_diagnostics.txt")
    with open(diag, "w") as fh:
        fh.write(f"Cu_AUC total rows: {len(cu_rows)}\n")
        fh.write(f"Matched to tree tip + aa_freq: {len(matched)}\n")
        fh.write(f"No matching tree tip (no genome/not in phyling tree): {len(unmatched_no_tip)}\n")
        fh.write(f"Tree tip matched but no aa_freq file: {len(unmatched_no_aa)}\n")
        fh.write(f"Matched but mean_auc_rate missing/non-numeric: {len(unmatched_bad_auc)}\n")
        fh.write("\nUnique species among matched strains:\n")
        from collections import Counter

        c = Counter(r["species"] for r in matched)
        for sp, n in sorted(c.items(), key=lambda x: -x[1]):
            fh.write(f"  {sp}: {n}\n")

    print(f"Wrote {len(matched)} rows to {out_csv}")
    print(f"Diagnostics: {diag}")


if __name__ == "__main__":
    main()
