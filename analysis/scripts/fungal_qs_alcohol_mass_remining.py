#!/usr/bin/env python3
"""
Fungal quorum-sensing alcohol search (2026-09-16, PI-directed follow-up
to the biofilm-signaling question). Unlike the AHL/aminolipid/ceramide
searches, these are the mechanistically best-established fungal
density-dependent signals for the exact phenomenon of interest
(coordinating yeast-hyphal/biofilm morphogenesis) -- verified this
session (background research, not assumed): farnesol represses the
yeast-to-hypha transition in Candida albicans (Hornby et al. 2001, Appl.
Environ. Microbiol., PMID 11425711); tyrosol does the reverse, promoting
germ-tube formation, and biofilm cells secrete substantially more of it
than planktonic cells. Every citation found for this mechanism is
Ascomycota (Candida-clade) -- NO Basidiomycota/Rhodotorula-specific
literature was found. This search tests, rather than assumes, whether
the chemistry is present in this genus at all.

Compounds (small molecules, not combinatorial homolog series like the
earlier searches -- these are specific, named natural products, not a
systematic sweep of a formula family):
  farnesol       C15H26O   (sesquiterpenoid alcohol; the primary QS molecule)
  farnesoic acid C15H24O2  (farnesol pathway oxidation product, also QS-active
                             in some reports; included as a closely related
                             bonus target, not independently verified this
                             session to the same standard as farnesol/tyrosol)
  tyrosol        C8H10O2   (aromatic alcohol; opposes farnesol)
  tryptophol     C10H11NO  (aromatic alcohol, indole-derived; same pathway family)
  phenylethanol  C8H10O    (aromatic alcohol, same pathway family)

Adducts: [M+H]+, [M+Na]+, [M+NH4]+ (positive mode only -- this project's
feature table has no negative-mode data, same finding as the earlier
N-acyl amino acid search).

Usage:
    python3 analysis/scripts/fungal_qs_alcohol_mass_remining.py --ppm 20
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
FULL_FEATURES_CSV = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9" / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output" / "feature_finding" / "feature_finding_results" / "aligned_features_ms2.csv"
)
DEDUP_GROUPS = REPO / "analysis" / "linked_data" / "ms_feature_dedup_groups.csv"
SIRIUS_ANNOTATIONS = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"
OUT_DIR = REPO / "analysis" / "ahl_autoinducer_search"

MASS = {"C": 12.0, "H": 1.0078250319, "O": 15.9949146221, "N": 14.0030740052, "e": 0.00054858}
PROTON = MASS["H"] - MASS["e"]
NA = 22.98976928 - MASS["e"]
NH4 = MASS["N"] + 4 * MASS["H"] - MASS["e"]

COMPOUNDS = {
    # name: (C, H, N, O)
    "farnesol": (15, 26, 0, 1),
    "farnesoic_acid": (15, 24, 0, 2),
    "tyrosol": (8, 10, 0, 2),
    "tryptophol": (10, 11, 1, 1),
    "phenylethanol": (8, 10, 0, 1),
}

ADDUCTS = {
    "[M+H]+": lambda m: m + PROTON,
    "[M+Na]+": lambda m: m + NA,
    "[M+NH4]+": lambda m: m + NH4,
}


def formula_mass(c: int, h: int, n: int, o: int) -> float:
    return c * MASS["C"] + h * MASS["H"] + n * MASS["N"] + o * MASS["O"]


def parse_formula(f):
    if not isinstance(f, str):
        return None
    d = {}
    for el, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        if not el:
            continue
        d[el] = d.get(el, 0) + (int(cnt) if cnt else 1)
    return d


def build_target_list() -> pd.DataFrame:
    rows = []
    for name, (c, h, n, o) in COMPOUNDS.items():
        mm = formula_mass(c, h, n, o)
        formula_str = f"C{c}H{h}" + (f"N{n}" if n else "") + f"O{o}"
        for adduct_name, fn in ADDUCTS.items():
            rows.append(dict(compound=name, formula=formula_str, neutral_mass=mm,
                              adduct=adduct_name, target_mz=fn(mm)))
    return pd.DataFrame(rows).sort_values("target_mz").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppm", type=float, default=20.0)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = build_target_list()
    targets.to_csv(OUT_DIR / "fungal_qs_alcohol_target_list.csv", index=False)
    print(f"Built {len(targets)} targets ({len(COMPOUNDS)} compounds x {len(ADDUCTS)} adducts)")
    print(targets[["compound", "formula", "neutral_mass"]].drop_duplicates().to_string(index=False))

    feat = pd.read_csv(FULL_FEATURES_CSV, usecols=[
        "row ID", "row m/z", "row retention time", "adduct", "is_default_adduct",
        "is_isf", "isf_parent_id", "adduct_source_id", "has_ms2", "total_scans",
    ])
    dedup = pd.read_csv(DEDUP_GROUPS) if DEDUP_GROUPS.exists() else None
    sirius = pd.read_csv(SIRIUS_ANNOTATIONS, sep="\t") if SIRIUS_ANNOTATIONS.exists() else None

    matches = []
    for _, t in targets.iterrows():
        tol = t["target_mz"] * args.ppm / 1e6
        hits = feat[(feat["row m/z"] >= t["target_mz"] - tol) & (feat["row m/z"] <= t["target_mz"] + tol)]
        for _, h in hits.iterrows():
            ppm_error = (h["row m/z"] - t["target_mz"]) / t["target_mz"] * 1e6
            row = dict(
                compound=t["compound"], formula=t["formula"], adduct=t["adduct"], target_mz=t["target_mz"],
                row_id=h["row ID"], observed_mz=h["row m/z"], ppm_error=ppm_error,
                rt_min=h["row retention time"], has_ms2=h["has_ms2"], is_isf=h["is_isf"],
                total_scans=h["total_scans"],
            )
            if dedup is not None:
                dgrp = dedup.loc[dedup["row ID"] == h["row ID"]]
                if len(dgrp):
                    row["dedup_group_id"] = dgrp.iloc[0]["dedup_group_id"]
            if sirius is not None:
                sann = sirius.loc[sirius["row ID"] == h["row ID"]]
                if len(sann):
                    row["sirius_formula"] = sann.iloc[0].get("sirius_formula")
                    row["sirius_structure_name"] = sann.iloc[0].get("sirius_structure_name")
                    row["sirius_npc_class"] = sann.iloc[0].get("sirius_npc_class")
            matches.append(row)

    matches_df = pd.DataFrame(matches)
    out_path = OUT_DIR / "fungal_qs_alcohol_mass_matches.csv"
    matches_df.to_csv(out_path, index=False)
    print(f"\n{len(matches_df)} raw feature match(es), {matches_df['row_id'].nunique() if len(matches_df) else 0} distinct rows")
    if len(matches_df):
        print(matches_df[["compound", "adduct", "row_id", "observed_mz", "ppm_error", "rt_min", "has_ms2", "total_scans"]].to_string(index=False))
        matches_df["_t"] = matches_df["formula"].apply(parse_formula)
        matches_df["_s"] = matches_df.get("sirius_formula", pd.Series(dtype=object)).apply(parse_formula)
        exact = matches_df[matches_df.apply(lambda r: r["_t"] == r["_s"] if r["_s"] else False, axis=1)].drop_duplicates("row_id")
        print(f"\nRows with SIRIUS formula EXACTLY matching target (element-count normalized): {len(exact)}")
        if len(exact):
            print(exact[["row_id", "compound", "adduct", "sirius_npc_class", "sirius_structure_name"]].to_string(index=False))
    else:
        print("No matches at all.")
    print(f"\nWrote {OUT_DIR / 'fungal_qs_alcohol_target_list.csv'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
