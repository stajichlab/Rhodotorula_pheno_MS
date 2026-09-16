#!/usr/bin/env python3
"""
Fungal ceramide (N-acyl sphingoid base) search (2026-09-16 expansion #2
of 3, PI-directed). Unlike the aminolipid searches (documented as
bacterial), N-acyl-sphingoid-base lipids (ceramides) are well-documented,
textbook FUNGAL membrane chemistry -- basidiomycetes including
Rhodotorula are known to build sphingolipids on phytosphingosine-type
long-chain bases (LCBs). This checks whether any of the broader
"N-acyl amines" SIRIUS-class signal in this dataset is actually
legitimate fungal ceramide chemistry rather than being read only through
a bacterial-aminolipid lens.

Formula derivation (first principles): ceramide = LCB + fatty acid - H2O
(amide bond at the LCB's own amino group).
  LCBs (both real, common fungal C18 long-chain bases):
    phytosphingosine   C18H39NO3 (4-hydroxysphinganine)
    dihydrosphingosine C18H39NO2 (sphinganine)
  Fatty acid, two variants (fungal ceramides commonly carry an
  alpha-hydroxylated fatty acid, not just a plain one):
    plain Cn:      CnH(2n)O2
    2-hydroxy Cn:  CnH(2n)O3

Acyl chain range C14-C26 (even chain lengths only -- fungal/plant
ceramide fatty acids are overwhelmingly even-chain, very-long-chain
fatty acids [VLCFAs] up to C26 are common in this lipid class, unlike
the C4-C18 range appropriate for the AHL/aminolipid searches). This is a
documented, motivated range choice, not an arbitrary extension of the
earlier searches' range.

Usage:
    python3 analysis/scripts/fungal_ceramide_mass_remining.py --ppm 20
"""
from __future__ import annotations

import argparse
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

LCBS = {
    "phytosphingosine": (18, 39, 1, 3),
    "dihydrosphingosine": (18, 39, 1, 2),
}
CHAIN_RANGE = [14, 16, 18, 20, 22, 24, 26]
FA_VARIANTS = {"plain": 2, "2-hydroxy": 3}  # extra O count for the fatty acid

ADDUCTS = {
    "[M+H]+": lambda m: m + PROTON,
    "[M+Na]+": lambda m: m + NA,
    "[M+NH4]+": lambda m: m + NH4,
}


def formula_mass(c: int, h: int, n: int, o: int) -> float:
    return c * MASS["C"] + h * MASS["H"] + n * MASS["N"] + o * MASS["O"]


def build_target_list() -> pd.DataFrame:
    rows = []
    for lcb_name, (lc, lh, ln, lo) in LCBS.items():
        for fa_name, fa_o in FA_VARIANTS.items():
            for n in CHAIN_RANGE:
                c = lc + n
                h = lh + 2 * n - 2
                o = lo + fa_o - 1
                mm = formula_mass(c=c, h=h, n=ln, o=o)
                formula_str = f"C{c}H{h}N{ln}O{o}"
                compound = f"{lcb_name}/{fa_name}-C{n}:0-ceramide"
                for adduct_name, fn in ADDUCTS.items():
                    rows.append(dict(
                        compound=compound, lcb=lcb_name, fa_variant=fa_name, acyl_n=n,
                        formula=formula_str, neutral_mass=mm, adduct=adduct_name, target_mz=fn(mm),
                    ))
    return pd.DataFrame(rows).sort_values("target_mz").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppm", type=float, default=20.0)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = build_target_list()
    targets.to_csv(OUT_DIR / "fungal_ceramide_target_list.csv", index=False)
    print(f"Built {len(targets)} targets ({len(LCBS)} LCBs x {len(FA_VARIANTS)} FA variants x "
          f"{len(CHAIN_RANGE)} chain lengths x {len(ADDUCTS)} adducts)")

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
                compound=t["compound"], lcb=t["lcb"], fa_variant=t["fa_variant"], acyl_n=t["acyl_n"],
                formula=t["formula"], adduct=t["adduct"], target_mz=t["target_mz"],
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
    out_path = OUT_DIR / "fungal_ceramide_mass_matches.csv"
    matches_df.to_csv(out_path, index=False)
    print(f"\n{len(matches_df)} raw feature match(es), {matches_df['row_id'].nunique() if len(matches_df) else 0} distinct rows")
    if len(matches_df):
        import re

        def parse_formula(f):
            if not isinstance(f, str):
                return None
            d = {}
            for el, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", f):
                if not el:
                    continue
                d[el] = d.get(el, 0) + (int(cnt) if cnt else 1)
            return d
        matches_df["_t"] = matches_df["formula"].apply(parse_formula)
        matches_df["_s"] = matches_df.get("sirius_formula", pd.Series(dtype=object)).apply(parse_formula)
        exact = matches_df[matches_df.apply(lambda r: r["_t"] == r["_s"] if r["_s"] else False, axis=1)].drop_duplicates("row_id")
        print(f"Rows with SIRIUS formula EXACTLY matching target (element-count normalized): {len(exact)}")
        if len(exact):
            print(exact[["row_id", "compound", "adduct", "sirius_npc_class", "sirius_structure_name"]].drop_duplicates("row_id").to_string(index=False))
    print(f"\nWrote {OUT_DIR / 'fungal_ceramide_target_list.csv'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
