#!/usr/bin/env python3
"""
Build a multifasta of the best-candidate rhodotorulic-acid NRPS ortholog
per strain (from siderophore_nrps_diamond_search.py's per-strain best-hit
table), for downstream alignment (mafft) and phylogeny (e.g. IQ-TREE/
FastTree). Includes the PI-supplied reference protein
(F2DD6D01_006956-T1, *R. kratochvilovae* Y14) as an anchor/outgroup.

Filtering (per PI, 2026-08-16): confirmed ortholog (pident/qcovhsp hit
the threshold used by siderophore_nrps_diamond_search.py) AND BUSCO
genome completeness >=90% -- excludes the 2 low-completeness
*R. mucilaginosa* strains (DBVPG_3236: 90.4%... note: 90.4 passes a
>=90 cutoff; DBVPG_3855: 81.0%, excluded) whose apparent gene absence is
attributed to incomplete assembly (per PI decision) rather than pursued
further as biology. Strains already lacking a confirmed ortholog (the
*Cystobasidium* outgroup, and low-completeness *R. mucilaginosa*) are
naturally excluded either way.

Usage:
    python3 analysis/scripts/siderophore_nrps_build_multifasta.py --min-busco 90
"""
from __future__ import annotations

import argparse
from pathlib import Path

import duckdb
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
REFERENCE = REPO / "analysis" / "integrated_analysis" / "phase_siderophore" / "reference" / "RA_NRPS.fa"
BFD_DUCKDB = REPO / "BFD" / "db" / "BFD.duckdb"
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase_siderophore" / "outputs"
STRAIN_SUMMARY = OUT_DIR / "RA_NRPS_strain_summary.csv"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-busco", type=float, default=90.0)
    args = ap.parse_args()

    summary = pd.read_csv(STRAIN_SUMMARY)
    con = duckdb.connect(str(BFD_DUCKDB), read_only=True)
    busco = con.execute("SELECT LOCUSTAG, complete_pct FROM busco_genome").df()
    summary = summary.merge(busco, on="LOCUSTAG", how="left")

    n_total = len(summary)
    n_ortholog = int(summary["has_ortholog"].sum())
    keep = summary[summary["has_ortholog"].fillna(False) & (summary["complete_pct"] >= args.min_busco)]
    dropped_low_busco = summary[summary["has_ortholog"].fillna(False) & (summary["complete_pct"] < args.min_busco)]
    print(f"{n_total} strains total, {n_ortholog} confirmed ortholog, "
          f"{len(dropped_low_busco)} dropped for BUSCO<{args.min_busco} "
          f"({', '.join(dropped_low_busco['STRAIN'].tolist())}), {len(keep)} retained")

    protein_ids = keep["sseqid"].tolist()
    placeholders = ",".join(f"'{p}'" for p in protein_ids)
    seqs = con.execute(f"SELECT protein_id, peptide FROM gene_proteins WHERE protein_id IN ({placeholders})").df()
    con.close()

    seqs = seqs.merge(keep[["sseqid", "GENUS", "SPECIES", "STRAIN", "pident", "qcovhsp", "bitscore"]],
                       left_on="protein_id", right_on="sseqid", how="left")
    missing = set(protein_ids) - set(seqs["protein_id"])
    if missing:
        print(f"WARNING: {len(missing)} protein_id(s) not found in gene_proteins table: {missing}")

    out_path = OUT_DIR / "RA_NRPS_candidates.faa"
    with open(out_path, "w") as fh:
        ref_lines = REFERENCE.read_text().splitlines()
        ref_header = ref_lines[0].split()[0][1:]
        ref_seq = "".join(ref_lines[1:])
        fh.write(f">{ref_header} REFERENCE Rhodotorula_kratochvilovae_Y14 (PI-supplied, antiSMASH-confirmed NRPS cluster gene)\n")
        for i in range(0, len(ref_seq), 80):
            fh.write(ref_seq[i:i + 80] + "\n")

        for _, row in seqs.sort_values(["GENUS", "SPECIES", "STRAIN"]).iterrows():
            species = str(row["SPECIES"]).replace(" ", "_")
            strain = str(row["STRAIN"]).replace(" ", "_")
            header = f">{row['protein_id']} {species}_{strain} pident={row['pident']:.1f} qcovhsp={row['qcovhsp']:.1f}"
            fh.write(header + "\n")
            pep = row["peptide"]
            for i in range(0, len(pep), 80):
                fh.write(pep[i:i + 80] + "\n")

    print(f"Wrote {out_path} ({len(seqs) + 1} sequences: {len(seqs)} candidates + 1 reference)")


if __name__ == "__main__":
    main()
