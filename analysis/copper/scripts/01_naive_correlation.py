"""Naive (non-phylogenetic) correlation of each proteome-wide amino-acid
frequency against copper-resistance AUC, treating each strain as an
independent observation (i.e. ignoring shared ancestry). This is the
'first principles' baseline that the PGLS step (02) is designed to
stress-test: any AA whose naive correlation collapses once phylogeny is
accounted for is a candidate lineage artifact rather than a trait-level
association with copper resistance.
"""
import csv
import os
import sys

from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from common import AA20, REPO

IN_CSV = os.path.join(REPO, "analysis/copper/outputs/copper_aa_master_table.csv")
OUT_CSV = os.path.join(REPO, "analysis/copper/outputs/naive_correlation_results.csv")

# Pre-registered candidates: residues with known metal-binding / redox roles
# relevant to copper handling (thiol Cys, imidazole His, carboxylate Asp/Glu,
# thioether Met). Tested and reported separately from the exploratory
# all-20-AA scan to avoid conflating a directed hypothesis with a fishing
# expedition.
PRIMARY_AA = ["C", "H", "D", "E", "M"]


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
    auc = [float(r["mean_auc_rate"]) for r in rows]

    results = []
    for aa in AA20:
        vals = [float(r[aa]) for r in rows]
        pear_r, pear_p = stats.pearsonr(vals, auc)
        spear_r, spear_p = stats.spearmanr(vals, auc)
        results.append({
            "amino_acid": aa,
            "n": len(vals),
            "pearson_r": pear_r,
            "pearson_p": pear_p,
            "spearman_r": spear_r,
            "spearman_p": spear_p,
            "primary_hypothesis": aa in PRIMARY_AA,
        })

    pvals = [r["pearson_p"] for r in results]
    qvals = bh_fdr(pvals)
    for r, q in zip(results, qvals):
        r["pearson_p_fdr_bh"] = q

    results.sort(key=lambda r: r["pearson_p"])

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)

    print(f"n strains = {len(rows)}")
    print("Top naive associations (by Pearson p, uncorrected):")
    for r in results[:6]:
        flag = "*" if r["primary_hypothesis"] else " "
        print(f"  {flag} {r['amino_acid']}: pearson_r={r['pearson_r']:.3f} p={r['pearson_p']:.4g} "
              f"q(BH)={r['pearson_p_fdr_bh']:.4g}")
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
