#!/usr/bin/env python3
"""
Diacylated ornithine/lysine lipid search (2026-09-16 expansion #1 of 3,
PI-directed follow-up to the N-acyl amino acid search).

Motivation: the mono-acyl N-acyl-lysine/N-oleoyl-ornithine candidates
already found (NACYL_AMINO_ACID_SEARCH.md) match the mono-acyl
INTERMEDIATE of the documented bacterial ornithine-lipid/lysine-lipid
(OL/LL) biosynthetic pathway (OlsB -> OlsA). The MATURE lipid in that
pathway is diacylated: OlsB makes a 3-hydroxy-fatty-acid amide on the
headgroup's alpha-amino group; OlsA then esterifies a SECOND fatty acid
onto that 3-hydroxyl. This search targets that mature, two-chain form,
not yet searched for.

Formula derivation (first principles, same method as
nacyl_amino_acid_mass_remining.py):
  headgroup (ornithine C5H12N2O2 or lysine C6H14N2O2)
  + 3-hydroxy fatty acid 1 (Cn, saturated: CnH(2n)O3 -- one more O than a
    plain fatty acid, for the 3-OH)
  - H2O (amide bond to headgroup's alpha-amino)
  + fatty acid 2 (Cm, saturated: CmH(2m)O2)
  - H2O (ester bond to fatty-acid-1's 3-OH)
  = C(a+n+m) H(h+2n+2m-4) N(1) O(o+3)
where headgroup = C(a)H(h)N(1)O(o).

Both acyl chain lengths (n, m) swept independently over a systematic,
not hand-picked, range (C10-C18, common bacterial membrane-lipid acyl
range) to avoid presupposing which chain-length combination is present.

Usage:
    python3 analysis/scripts/diacyl_aminolipid_mass_remining.py --ppm 20
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

HEADGROUPS = {
    "ornithine": (5, 12, 2, 2),   # C, H, N, O
    "lysine": (6, 14, 2, 2),
}
CHAIN_RANGE = [10, 12, 14, 15, 16, 17, 18]  # both acyl positions, systematic sweep

ADDUCTS = {
    "[M+H]+": lambda m: m + PROTON,
    "[M+Na]+": lambda m: m + NA,
    "[M+NH4]+": lambda m: m + NH4,
}


def formula_mass(c: int, h: int, n: int, o: int) -> float:
    return c * MASS["C"] + h * MASS["H"] + n * MASS["N"] + o * MASS["O"]


def build_target_list() -> pd.DataFrame:
    rows = []
    for hg_name, (a, h_hg, n_hg, o_hg) in HEADGROUPS.items():
        for n in CHAIN_RANGE:
            for m in CHAIN_RANGE:
                c = a + n + m
                h = h_hg + 2 * n + 2 * m - 4
                o = o_hg + 3
                mm = formula_mass(c=c, h=h, n=n_hg, o=o)
                formula_str = f"C{c}H{h}N{n_hg}O{o}"
                compound = f"diacyl-{hg_name}-3OH-C{n}:0/C{m}:0"
                for adduct_name, fn in ADDUCTS.items():
                    rows.append(dict(
                        compound=compound, headgroup=hg_name, chain1_n=n, chain2_n=m,
                        formula=formula_str, neutral_mass=mm, adduct=adduct_name, target_mz=fn(mm),
                    ))
    return pd.DataFrame(rows).sort_values("target_mz").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppm", type=float, default=20.0)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = build_target_list()
    targets.to_csv(OUT_DIR / "diacyl_aminolipid_target_list.csv", index=False)
    print(f"Built {len(targets)} targets ({len(HEADGROUPS)} headgroups x {len(CHAIN_RANGE)}x{len(CHAIN_RANGE)} "
          f"chain-length combos x {len(ADDUCTS)} adducts)")

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
                compound=t["compound"], headgroup=t["headgroup"], chain1_n=t["chain1_n"], chain2_n=t["chain2_n"],
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
    out_path = OUT_DIR / "diacyl_aminolipid_mass_matches.csv"
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
    print(f"\nWrote {OUT_DIR / 'diacyl_aminolipid_target_list.csv'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
