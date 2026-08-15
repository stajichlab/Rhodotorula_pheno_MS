"""Merge naive and PGLS results into one comparison table: which AA
associations survive phylogenetic correction vs which were lineage
artifacts of the naive (non-phylogenetic) test."""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from common import REPO

NAIVE_CSV = os.path.join(REPO, "analysis/copper/outputs/naive_correlation_results.csv")
PGLS_CSV = os.path.join(REPO, "analysis/copper/outputs/pgls_correlation_results.csv")
OUT_CSV = os.path.join(REPO, "analysis/copper/outputs/naive_vs_pgls_comparison.csv")

Q_THRESH = 0.05


def main():
    with open(NAIVE_CSV) as fh:
        naive = {r["amino_acid"]: r for r in csv.DictReader(fh)}
    with open(PGLS_CSV) as fh:
        pgls = {r["amino_acid"]: r for r in csv.DictReader(fh)}

    rows = []
    for aa in naive:
        n, p = naive[aa], pgls[aa]
        naive_sig = float(n["pearson_p_fdr_bh"]) < Q_THRESH
        pgls_sig = p["p_value_fdr_bh"] not in ("", "NA") and float(p["p_value_fdr_bh"]) < Q_THRESH
        rows.append({
            "amino_acid": aa,
            "primary_hypothesis": n["primary_hypothesis"],
            "naive_pearson_r": n["pearson_r"],
            "naive_q_bh": n["pearson_p_fdr_bh"],
            "naive_significant": naive_sig,
            "pgls_slope": p["slope"],
            "pgls_lambda": p["lambda"],
            "pgls_q_bh": p["p_value_fdr_bh"],
            "pgls_significant": pgls_sig,
            "verdict": (
                "survives_phylo_correction" if naive_sig and pgls_sig else
                "naive_only_lineage_artifact" if naive_sig and not pgls_sig else
                "pgls_only" if pgls_sig and not naive_sig else
                "not_significant"
            ),
        })

    rows.sort(key=lambda r: (r["verdict"] != "survives_phylo_correction", r["pgls_q_bh"]))

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {OUT_CSV}\n")
    print(f"{'AA':<3}{'primary':<8}{'naive_q':<10}{'pgls_lambda':<12}{'pgls_q':<10}verdict")
    for r in rows:
        print(f"{r['amino_acid']:<3}{r['primary_hypothesis']:<8}{float(r['naive_q_bh']):<10.3g}"
              f"{float(r['pgls_lambda']):<12.3f}{float(r['pgls_q_bh']):<10.3g}{r['verdict']}")


if __name__ == "__main__":
    main()
