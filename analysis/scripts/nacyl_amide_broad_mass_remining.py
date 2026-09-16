#!/usr/bin/env python3
"""
Broader N-acyl amide family search (2026-09-16 expansion #3 of 3,
PI-directed). The N-acyl-amino-acid literature review (Cohen et al. 2017,
Nature -- commensal-bacteria-derived GPCR ligands) documents a wider
"N-acyl amide" natural-product superfamily beyond amino acids: fatty
acyl conjugates of ethanolamine, taurine, GABA, and biogenic amines
(dopamine, serotonin, tryptamine). This is the least specifically
motivated of the three 2026-09-16 expansions (no specific biosynthetic
pathway or prior hit in this dataset points here the way the diacyl-
aminolipid and ceramide searches were motivated) -- a systematic,
unbiased sweep of a documented related chemical space, not a targeted
follow-up on an existing lead.

Formula derivation (first principles, same amide-condensation method as
the other searches): headgroup (C_a H_h N_n O_o [S_s]) + fatty acid
(Cm H(2m) O2) - H2O.

Headgroup neutral formulas:
  ethanolamine  C2H7NO
  taurine       C2H7NO3S
  GABA          C4H9NO2
  dopamine      C8H11NO2
  serotonin     C10H12N2O
  tryptamine    C10H12N2

Acyl chain range C8-C22 (wider than the AHL/aminolipid searches, per the
gut-microbiome N-acyl-amide literature's typical reported chain-length
range).

Usage:
    python3 analysis/scripts/nacyl_amide_broad_mass_remining.py --ppm 20
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

MASS = {"C": 12.0, "H": 1.0078250319, "O": 15.9949146221, "N": 14.0030740052, "S": 31.97207069, "e": 0.00054858}
PROTON = MASS["H"] - MASS["e"]
NA = 22.98976928 - MASS["e"]
NH4 = MASS["N"] + 4 * MASS["H"] - MASS["e"]

HEADGROUPS = {
    # name: (C, H, N, O, S)
    "ethanolamine": (2, 7, 1, 1, 0),
    "taurine": (2, 7, 1, 3, 1),
    "GABA": (4, 9, 1, 2, 0),
    "dopamine": (8, 11, 1, 2, 0),
    "serotonin": (10, 12, 2, 1, 0),
    "tryptamine": (10, 12, 2, 0, 0),
}
CHAIN_RANGE = list(range(8, 23))  # C8-C22

ADDUCTS = {
    "[M+H]+": lambda m: m + PROTON,
    "[M+Na]+": lambda m: m + NA,
    "[M+NH4]+": lambda m: m + NH4,
}


def formula_mass(c: int, h: int, n: int, o: int, s: int = 0) -> float:
    return c * MASS["C"] + h * MASS["H"] + n * MASS["N"] + o * MASS["O"] + s * MASS["S"]


def build_target_list() -> pd.DataFrame:
    rows = []
    for hg_name, (a, h_hg, n_hg, o_hg, s_hg) in HEADGROUPS.items():
        for n in CHAIN_RANGE:
            c = a + n
            h = h_hg + 2 * n - 2
            o = o_hg + 1
            mm = formula_mass(c=c, h=h, n=n_hg, o=o, s=s_hg)
            formula_str = f"C{c}H{h}N{n_hg}O{o}" + (f"S{s_hg}" if s_hg else "")
            compound = f"N-C{n}:0-{hg_name}"
            for adduct_name, fn in ADDUCTS.items():
                rows.append(dict(
                    compound=compound, headgroup=hg_name, acyl_n=n,
                    formula=formula_str, neutral_mass=mm, adduct=adduct_name, target_mz=fn(mm),
                ))
    return pd.DataFrame(rows).sort_values("target_mz").reset_index(drop=True)


def parse_formula(f):
    import re
    if not isinstance(f, str):
        return None
    d = {}
    for el, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        if not el:
            continue
        d[el] = d.get(el, 0) + (int(cnt) if cnt else 1)
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppm", type=float, default=20.0)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = build_target_list()
    targets.to_csv(OUT_DIR / "nacyl_amide_broad_target_list.csv", index=False)
    print(f"Built {len(targets)} targets ({len(HEADGROUPS)} headgroups x {len(CHAIN_RANGE)} chain "
          f"lengths x {len(ADDUCTS)} adducts)")

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
                compound=t["compound"], headgroup=t["headgroup"], acyl_n=t["acyl_n"],
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
    out_path = OUT_DIR / "nacyl_amide_broad_mass_matches.csv"
    matches_df.to_csv(out_path, index=False)
    print(f"\n{len(matches_df)} raw feature match(es), {matches_df['row_id'].nunique() if len(matches_df) else 0} distinct rows")
    if len(matches_df):
        matches_df["target_parsed"] = matches_df["formula"].apply(parse_formula)
        matches_df["sirius_parsed"] = matches_df["sirius_formula"].apply(parse_formula) if "sirius_formula" in matches_df else None
        matches_df["formula_match"] = matches_df.apply(
            lambda r: r["target_parsed"] == r["sirius_parsed"] if r.get("sirius_parsed") else False, axis=1)
        exact = matches_df[matches_df["formula_match"]].drop_duplicates("row_id")
        print(f"Rows with SIRIUS formula EXACTLY matching target (element-count normalized): {len(exact)}")
        if len(exact):
            print(exact[["row_id", "compound", "adduct", "sirius_npc_class", "sirius_structure_name"]].to_string(index=False))
    print(f"\nWrote {OUT_DIR / 'nacyl_amide_broad_target_list.csv'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
