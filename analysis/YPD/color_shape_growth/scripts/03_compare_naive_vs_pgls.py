"""Merge naive and PGLS results per YPD phenotype into one comparison
table, classifying each AA x phenotype pair as surviving phylogenetic
correction, a naive-only lineage artifact, PGLS-only (masked in naive),
or not significant. Mirrors analysis/copper/scripts/03_compare_naive_vs_pgls.py.
"""
import csv
import os

REPO = "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS"
NAIVE_CSV = os.path.join(REPO, "analysis/YPD/color_shape_growth/outputs/naive_correlation_results.csv")
PGLS_CSV = os.path.join(REPO, "analysis/YPD/color_shape_growth/outputs/pgls_correlation_results.csv")
OUT_CSV = os.path.join(REPO, "analysis/YPD/color_shape_growth/outputs/naive_vs_pgls_comparison.csv")

# PGLS phenotype column names get mangled by R's read.csv (L*Mean -> L.Mean etc).
PHENO_MAP = {
    "Mean_Shape_Area": "Mean_Shape_Area",
    "Mean_ColorLab_L*Mean": "Mean_ColorLab_L.Mean",
    "Mean_ColorLab_a*Mean": "Mean_ColorLab_a.Mean",
    "Mean_ColorLab_b*Mean": "Mean_ColorLab_b.Mean",
}

Q_THRESH = 0.05


def main():
    with open(NAIVE_CSV) as fh:
        naive = list(csv.DictReader(fh))
    with open(PGLS_CSV) as fh:
        pgls = {(r["phenotype"], r["amino_acid"]): r for r in csv.DictReader(fh)}

    rows = []
    for n in naive:
        pheno_naive, aa = n["phenotype"], n["amino_acid"]
        pgls_key = (PHENO_MAP[pheno_naive], aa)
        p = pgls.get(pgls_key)
        naive_sig = float(n["pearson_p_fdr_bh"]) < Q_THRESH
        pgls_sig = p is not None and p["p_value_fdr_bh"] not in ("", "NA") and float(p["p_value_fdr_bh"]) < Q_THRESH
        rows.append({
            "phenotype": pheno_naive,
            "amino_acid": aa,
            "naive_pearson_r": n["pearson_r"],
            "naive_q_bh": n["pearson_p_fdr_bh"],
            "naive_significant": naive_sig,
            "pgls_lambda": p["lambda"] if p else "",
            "pgls_q_bh": p["p_value_fdr_bh"] if p else "",
            "pgls_significant": pgls_sig,
            "verdict": (
                "survives_phylo_correction" if naive_sig and pgls_sig else
                "naive_only_lineage_artifact" if naive_sig and not pgls_sig else
                "pgls_only" if pgls_sig and not naive_sig else
                "not_significant"
            ),
        })

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {OUT_CSV}\n")
    for pheno in PHENO_MAP:
        sub = [r for r in rows if r["phenotype"] == pheno]
        survive = [r["amino_acid"] for r in sub if r["verdict"] == "survives_phylo_correction"]
        n_naive = sum(1 for r in sub if r["naive_significant"])
        print(f"{pheno}: {n_naive} naive hits -> {len(survive)} survive PGLS: {survive}")


if __name__ == "__main__":
    main()
