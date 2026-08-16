#!/usr/bin/env python3
"""
Real ortholog search for the rhodotorulic-acid NRPS, superseding the
coarse Pfam-domain stand-in in siderophore_nrps_pfam_screen.py (which
could not resolve a 2-module-NRPS architecture call from fragmented draft
gene models). PI supplied the actual reference sequence
(tmpin/RA_NRPS.fa -> analysis/integrated_analysis/phase_siderophore/reference/RA_NRPS.fa):
protein F2DD6D01_006956-T1 from *Rhodotorula kratochvilovae* Y14 (external
reference genome, NOT one of the 278 BFD-panel strains -- though 3
DBVPG-strain *R. kratochvilovae* ARE in the panel, giving a useful
self-hit sanity check). Confirmed by antiSMASH
(.../Rhodotorula_kratochvilovae_Y14/antismash_local/JAFEUJ010000019.1.region001.gbk,
product="NRPS") to sit in a real biosynthetic gene cluster, with a
"biosynthetic-additional" smCOG gene (F2DD6D01_006955, plausibly the
ornithine N5-hydroxylase partner) immediately adjacent -- consistent
cluster architecture for a siderophore BGC, not an isolated domain hit.

Runs `diamond blastp` (query = single reference protein, target = all
2,188,032 BFD-panel proteins) and calls per-strain best-hit orthology.
The identity distribution is cleanly bimodal (see RESULTS.md) -- 303
hits below 30% identity (generic NRPS/AMP-binding-domain cross-hits, same
noise the coarse Pfam screen couldn't filter out) vs. 305 hits at >=60%
identity with high query coverage (real orthologs, including the 3
in-panel *R. kratochvilovae* strains at 74-99.8% identity as the
positive control). Default ortholog threshold: pident>=50 AND
qcovhsp>=70 (falls in the empty gap between the two modes).

Usage:
    python3 analysis/scripts/siderophore_nrps_diamond_search.py
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import duckdb
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
REFERENCE = REPO / "analysis" / "integrated_analysis" / "phase_siderophore" / "reference" / "RA_NRPS.fa"
BFD_PROTEOMES_DIR = REPO / "BFD" / "input" / "pep"
BFD_DUCKDB = REPO / "BFD" / "db" / "BFD.duckdb"
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase_siderophore" / "outputs"

DIAMOND_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                "qlen", "slen", "qcovhsp", "scovhsp"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pident-min", type=float, default=50.0)
    ap.add_argument("--qcovhsp-min", type=float, default=70.0)
    ap.add_argument("--rebuild-db", action="store_true", help="Rebuild the diamond database from BFD/input/pep even if it already exists.")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    db_path = OUT_DIR / "bfd_proteomes.dmnd"
    all_proteomes_fa = OUT_DIR / "_all_proteomes_tmp.fa"

    if args.rebuild_db or not db_path.exists():
        with open(all_proteomes_fa, "wb") as out:
            for fa in sorted(BFD_PROTEOMES_DIR.glob("*.proteins.fa")):
                out.write(fa.read_bytes())
        subprocess.run(["diamond", "makedb", "--in", str(all_proteomes_fa), "-d", str(db_path.with_suffix(""))], check=True)
        all_proteomes_fa.unlink()

    hits_path = OUT_DIR / "RA_NRPS_diamond_hits.tsv"
    subprocess.run([
        "diamond", "blastp", "-q", str(REFERENCE), "-d", str(db_path.with_suffix("")),
        "-o", str(hits_path), "--outfmt", "6", *DIAMOND_COLS,
        "--max-target-seqs", "1000", "--evalue", "1e-5", "--sensitive", "--threads", "4",
    ], check=True)

    df = pd.read_csv(hits_path, sep="\t", names=DIAMOND_COLS)
    best = df.sort_values("bitscore", ascending=False).drop_duplicates("sseqid")
    best["species_prefix"] = best["sseqid"].str.split("_").str[0]
    best["is_ortholog"] = (best["pident"] >= args.pident_min) & (best["qcovhsp"] >= args.qcovhsp_min)

    con = duckdb.connect(str(BFD_DUCKDB), read_only=True)
    species = con.execute("SELECT LOCUSTAG, GENUS, SPECIES, STRAIN FROM species").df()
    con.close()

    best = best.merge(species, left_on="species_prefix", right_on="LOCUSTAG", how="left")
    best_out = best[["sseqid", "species_prefix", "GENUS", "SPECIES", "STRAIN", "pident", "qcovhsp",
                      "scovhsp", "evalue", "bitscore", "is_ortholog"]].sort_values("bitscore", ascending=False)
    best_out.to_csv(OUT_DIR / "RA_NRPS_diamond_best_hit_per_protein.csv", index=False)

    # per-strain: best (highest-bitscore) hit among all proteins in that strain
    per_strain_best = best.sort_values("bitscore", ascending=False).drop_duplicates("species_prefix")
    per_strain_best = per_strain_best[["species_prefix", "GENUS", "SPECIES", "STRAIN", "sseqid",
                                        "pident", "qcovhsp", "scovhsp", "evalue", "bitscore", "is_ortholog"]]
    summary = species.copy()
    summary = summary.merge(per_strain_best.drop(columns=["GENUS", "SPECIES", "STRAIN"]),
                             left_on="LOCUSTAG", right_on="species_prefix", how="left")
    summary["has_hit"] = summary["sseqid"].notna()
    summary["has_ortholog"] = summary["is_ortholog"].fillna(False)
    summary.to_csv(OUT_DIR / "RA_NRPS_strain_summary.csv", index=False)

    print(f"Total unique target proteins hit: {len(best)}")
    print(f"Identity distribution (all hits):\n{best['pident'].describe().to_string()}")
    print(f"\nStrains with >=1 hit at all (any e-value<1e-5): {int(summary['has_hit'].sum())}/278")
    print(f"Strains with a confirmed ortholog (pident>={args.pident_min}, qcovhsp>={args.qcovhsp_min}): "
          f"{int(summary['has_ortholog'].sum())}/278")
    print("\nPer-species ortholog presence:")
    print(summary.groupby(["GENUS", "SPECIES"])["has_ortholog"].agg(["sum", "count"]).to_string())
    missing = summary.loc[~summary["has_ortholog"], ["GENUS", "SPECIES", "STRAIN"]]
    print(f"\n{len(missing)} strains WITHOUT a confirmed ortholog:")
    print(missing.to_string(index=False))

    print(f"\nWrote {OUT_DIR / 'RA_NRPS_diamond_best_hit_per_protein.csv'}")
    print(f"Wrote {OUT_DIR / 'RA_NRPS_strain_summary.csv'}")


if __name__ == "__main__":
    main()
