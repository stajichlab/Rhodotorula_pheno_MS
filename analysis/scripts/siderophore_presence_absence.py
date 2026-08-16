#!/usr/bin/env python3
"""
PI request (2026-08-16), step 2 of the siderophore investigation: for the
candidate features found by siderophore_mass_remining.py, determine
per-strain presence/absence (detected in >=1 replicate sample of a
fraction, using a nonzero raw-abundance threshold -- these are raw
(non-TSS-normalized) counts, since presence/absence should not depend on
what else was in the same run) and summarize which strains/species show
the compound at all. This is a detectability/occurrence survey, not a
statistical association test -- no permutation, no FDR; that's a
follow-up once (a) a specific candidate is chosen as the best rhodotorulic
acid signal and (b) it's cross-referenced against the NRPS genome screen
(siderophore_nrps_pfam_screen.py).

Usage:
    python3 analysis/scripts/siderophore_presence_absence.py --row-id 2190 1621 5985 7524 3827
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
FEATURE_MATRIX = REPO / "analysis" / "linked_data" / "feature_abundance_matrix.csv.gz"
DEDUP_GROUPS = REPO / "analysis" / "linked_data" / "ms_feature_dedup_groups.csv"
SAMPLE_METADATA = REPO / "analysis" / "linked_data" / "sample_metadata.csv.gz"
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase_siderophore"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--row-id", type=int, nargs="+", required=True)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    feat = pd.read_csv(FEATURE_MATRIX)
    dedup = pd.read_csv(DEDUP_GROUPS)
    meta = pd.read_csv(SAMPLE_METADATA)

    group_of = dedup.set_index("row ID")["dedup_group_id"].to_dict()

    strain_rows = []
    overview_rows = []
    for rid in args.row_id:
        if rid not in group_of:
            print(f"row {rid}: not in dedup groups, skipping")
            continue
        gid = group_of[rid]
        feat_row = feat.loc[feat["row ID"] == rid]
        if feat_row.empty:
            print(f"row {rid}: not in feature matrix, skipping")
            continue
        feat_row = feat_row.iloc[0]

        for fraction in ["cell", "supernatant"]:
            fmeta = meta[meta["fraction"] == fraction]
            sample_cols = [c for c in feat.columns if c in set(fmeta["sample_id"])]
            vals = feat_row[sample_cols].astype(float)
            present = vals > 0
            strain_by_sample = fmeta.set_index("sample_id")["canonical_strain"].to_dict()
            species_by_sample = fmeta.set_index("sample_id")["Species"].to_dict() if "Species" in fmeta.columns else {}
            df = pd.DataFrame({
                "sample_id": sample_cols,
                "abundance": vals.to_numpy(),
                "present": present.to_numpy(),
                "strain": [strain_by_sample.get(c) for c in sample_cols],
                "species": [species_by_sample.get(c) for c in sample_cols],
            })
            per_strain = df.groupby("strain").agg(
                species=("species", "first"),
                n_samples=("present", "size"),
                n_present=("present", "sum"),
                max_abundance=("abundance", "max"),
            ).reset_index()
            per_strain["row_id"] = rid
            per_strain["dedup_group_id"] = gid
            per_strain["fraction"] = fraction
            strain_rows.append(per_strain)

            n_strains_present = (per_strain["n_present"] > 0).sum()
            n_strains_total = len(per_strain)
            species_present = sorted(per_strain.loc[per_strain["n_present"] > 0, "species"].dropna().unique().tolist())
            overview_rows.append(dict(
                row_id=rid, dedup_group_id=gid, fraction=fraction,
                n_strains_present=n_strains_present, n_strains_total=n_strains_total,
                pct_strains_present=100 * n_strains_present / n_strains_total if n_strains_total else float("nan"),
                n_species_present=len(species_present),
                species_present=";".join(species_present),
            ))

    strain_detail = pd.concat(strain_rows, ignore_index=True)
    strain_detail.to_csv(OUT_DIR / "siderophore_presence_by_strain.csv", index=False)
    overview = pd.DataFrame(overview_rows)
    overview.to_csv(OUT_DIR / "siderophore_presence_overview.csv", index=False)

    print(overview.to_string(index=False))
    print(f"\nWrote {OUT_DIR / 'siderophore_presence_by_strain.csv'}")
    print(f"Wrote {OUT_DIR / 'siderophore_presence_overview.csv'}")


if __name__ == "__main__":
    main()
