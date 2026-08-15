"""Join YPD2 non-stress phenotypes (colony color L*a*b*, colony Shape_Area
as a growth-size proxy) to per-strain proteome-wide amino-acid composition,
restricted to strains present in the phyling protein tree. Mirrors
analysis/copper/scripts/00_build_master_table.py; run as an independent
baseline (non-stress condition, so not subject to copper-specific
confounds) against which the copper-condition PGLS hits can be compared.

A strain with multiple YPD2 rows (different plates/replicates) is
collapsed to its mean across rows before joining.
"""
import csv
import os
import sys
from collections import defaultdict

REPO_ROOT = "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS"
sys.path.insert(0, os.path.join(REPO_ROOT, "analysis/copper/scripts"))
from common import AA20, REPO, load_aa_freq_table, load_tree_tips, load_ypd2, match_to_tree, norm

OUT_DIR = os.path.join(REPO, "analysis/YPD/color_shape_growth/outputs")

PHENO_COLS = [
    "Mean_Shape_Area",
    "Mean_ColorLab_L*Mean",
    "Mean_ColorLab_a*Mean",
    "Mean_ColorLab_b*Mean",
]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    tips = load_tree_tips()
    tip_norm_map = {norm(t): t for t in tips}

    aa_table = load_aa_freq_table()
    aa_prefix_to_tip = {}
    for prefix in aa_table:
        tip = match_to_tree(prefix, tip_norm_map)
        if tip is not None:
            aa_prefix_to_tip[prefix] = tip
    tip_to_aa_prefix = {v: k for k, v in aa_prefix_to_tip.items()}

    ypd_rows = load_ypd2()
    by_strain = defaultdict(list)
    for row in ypd_rows:
        by_strain[row["Strain"]].append(row)

    matched = []
    unmatched_no_tip = 0
    unmatched_no_aa = 0
    unmatched_missing_pheno = 0
    for strain, rows in by_strain.items():
        tip = match_to_tree(strain, tip_norm_map)
        if tip is None:
            unmatched_no_tip += 1
            continue
        aa_prefix = tip_to_aa_prefix.get(tip)
        if aa_prefix is None:
            unmatched_no_aa += 1
            continue
        freqs = aa_table[aa_prefix]

        vals = {c: [] for c in PHENO_COLS}
        ok = True
        for r in rows:
            for c in PHENO_COLS:
                try:
                    vals[c].append(float(r[c]))
                except (TypeError, ValueError):
                    ok = False
        if not ok or any(len(vals[c]) == 0 for c in PHENO_COLS):
            unmatched_missing_pheno += 1
            continue

        out = {
            "strain": strain,
            "species": rows[0]["Species"],
            "tree_tip": tip,
            "aa_freq_prefix": aa_prefix,
            "n_replicates": len(rows),
        }
        for c in PHENO_COLS:
            out[c] = sum(vals[c]) / len(vals[c])
        for aa in AA20:
            out[aa] = freqs.get(aa, "")
        matched.append(out)

    out_csv = os.path.join(OUT_DIR, "ypd_aa_master_table.csv")
    fieldnames = ["strain", "species", "tree_tip", "aa_freq_prefix", "n_replicates"] + PHENO_COLS + AA20
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(matched)

    diag = os.path.join(OUT_DIR, "join_diagnostics.txt")
    with open(diag, "w") as fh:
        fh.write(f"YPD2 unique strains: {len(by_strain)}\n")
        fh.write(f"Matched to tree tip + aa_freq + complete phenotypes: {len(matched)}\n")
        fh.write(f"No matching tree tip: {unmatched_no_tip}\n")
        fh.write(f"Tree tip matched but no aa_freq file: {unmatched_no_aa}\n")
        fh.write(f"Matched but missing/non-numeric phenotype values: {unmatched_missing_pheno}\n")
        from collections import Counter

        fh.write("\nUnique species among matched strains:\n")
        c = Counter(r["species"] for r in matched)
        for sp, n in sorted(c.items(), key=lambda x: -x[1]):
            fh.write(f"  {sp}: {n}\n")

    print(f"Wrote {len(matched)} rows to {out_csv}")
    print(f"Diagnostics: {diag}")


if __name__ == "__main__":
    main()
