"""Naive correlation of proteome-wide AA frequency against each of 4
non-stress YPD phenotypes (colony size proxy for growth, L*/a*/b* color).
Mirrors analysis/copper/scripts/01_naive_correlation.py. BH-FDR is applied
within each phenotype across the 20 AAs (4 separate families of 20 tests,
not pooled across phenotypes).
"""
import csv
import os
import sys

from scipy import stats

REPO = "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS"
sys.path.insert(0, os.path.join(REPO, "analysis/copper/scripts"))
from common import AA20

IN_CSV = os.path.join(REPO, "analysis/YPD/color_shape_growth/outputs/ypd_aa_master_table.csv")
OUT_CSV = os.path.join(REPO, "analysis/YPD/color_shape_growth/outputs/naive_correlation_results.csv")

PHENOTYPES = ["Mean_Shape_Area", "Mean_ColorLab_L*Mean", "Mean_ColorLab_a*Mean", "Mean_ColorLab_b*Mean"]


def bh_fdr(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    ranked = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = order[n - rank]
        q = pvals[i] * n / (n - rank + 1)
        prev = min(prev, q)
        ranked[i] = prev
    return ranked


def main():
    with open(IN_CSV) as fh:
        rows = list(csv.DictReader(fh))

    all_results = []
    for pheno in PHENOTYPES:
        y = [float(r[pheno]) for r in rows]
        results = []
        for aa in AA20:
            vals = [float(r[aa]) for r in rows]
            pear_r, pear_p = stats.pearsonr(vals, y)
            results.append({"phenotype": pheno, "amino_acid": aa, "n": len(vals),
                             "pearson_r": pear_r, "pearson_p": pear_p})
        qvals = bh_fdr([r["pearson_p"] for r in results])
        for r, q in zip(results, qvals):
            r["pearson_p_fdr_bh"] = q
        all_results.extend(results)

    all_results.sort(key=lambda r: r["pearson_p"])
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_results[0].keys()))
        w.writeheader()
        w.writerows(all_results)

    print(f"n strains = {len(rows)}")
    print("Naive hits at BH q<0.05, by phenotype:")
    for pheno in PHENOTYPES:
        hits = [r for r in all_results if r["phenotype"] == pheno and r["pearson_p_fdr_bh"] < 0.05]
        hits.sort(key=lambda r: r["pearson_p"])
        print(f"  {pheno}: {[h['amino_acid'] for h in hits]}")
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
