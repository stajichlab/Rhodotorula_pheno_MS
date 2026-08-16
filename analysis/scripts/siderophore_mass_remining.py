#!/usr/bin/env python3
"""
PI request (2026-08-16): targeted exact-mass search for rhodotorulic acid
and related iron-sequestration (siderophore) chemistry in the existing raw
MS2 feature table -- same design as idea1_targeted_mass_remining.py
(carotenoid pathway), applied to a different compound class. Step 1 of
2: "what can be detected at all" (this script). Step 2 (not done here):
once candidates are found, check presence/absence per strain and
cross-reference against the NRPS/ornithine-hydroxylase genome screen
(see siderophore_nrps_pfam_screen.py).

Compound rationale:
  - rhodotorulic acid (C20H32N4O6): the genus-defining extracellular
    siderophore for *Rhodotorula*/*Rhodosporidium* -- a cyclic
    dihydroxamate diketopiperazine built from two N5-acetyl-N5-
    hydroxy-L-ornithine units. This is the primary target.
  - N5-acetyl-N5-hydroxy-L-ornithine (C7H14N2O4): the monomeric
    biosynthetic precursor/building block of rhodotorulic acid --
    smaller and more polar, plausibly easier to ionize by ESI than the
    cyclic dimer; also useful as a pathway-intermediate signal even if
    the mature product isn't detected.
  - dimerumic acid (C22H37N3O9): a related fungal hydroxamate
    siderophore reported in some basidiomycetous yeasts -- included as a
    broader net, NOT strain-confirmed for *Rhodotorula* specifically
    (same caveat status as idea1's "additional_carotenoid" category).
  - ferrichrome (C27H42N9O12, cyclic hexapeptide trihydroxamate): the
    classic fungal siderophore, common across Ascomycota/Basidiomycota
    generally -- included as a genus-agnostic reference point, not
    because it's specifically expected in *Rhodotorula*.

NOT searched: Fe(III)-bound complexes. Iron-siderophore stoichiometry is
not a simple "replace 3H with Fe" substitution for these ligand
geometries (e.g. rhodotorulic acid forms a 3:2 ligand:Fe complex, not
1:1) -- computing the correct target masses needs the actual coordination
chemistry, not a quick formula edit, so this pass is apo (metal-free)
forms only. A real Fe-bound-complex search is a legitimate follow-up, not
done here.

Adducts searched: [M+H]+, [M+Na]+, [M+NH4]+, [M+K]+ (siderophores are
often observed as Na+/K+ adducts in ESI+ due to their multiple
carbonyl/hydroxamate oxygen donors) -- confirmed the raw feature table is
ESI positive-mode only (no [M-H]- or other negative-mode adducts present),
so [M-H]- was not searched.

Usage:
    python3 analysis/scripts/siderophore_mass_remining.py --ppm 20
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
FULL_FEATURES_CSV = (
    REPO
    / "data"
    / "processed"
    / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output"
    / "feature_finding"
    / "feature_finding_results"
    / "aligned_features_ms2.csv"
)
DEDUP_GROUPS = REPO / "analysis" / "linked_data" / "ms_feature_dedup_groups.csv"
SIRIUS_ANNOTATIONS = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"
SAMPLE_METADATA = REPO / "analysis" / "linked_data" / "sample_metadata.csv.gz"
FEATURE_MATRIX = REPO / "analysis" / "linked_data" / "feature_abundance_matrix.csv.gz"
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase_siderophore"

MASS = {"C": 12.0, "H": 1.0078250319, "O": 15.9949146221, "N": 14.0030740052, "K": 38.9637064864, "e": 0.00054858}
PROTON = MASS["H"] - MASS["e"]
NA = 22.98976928 - MASS["e"]
NH4 = MASS["N"] + 4 * MASS["H"] - MASS["e"]
K = MASS["K"] - MASS["e"]


def formula_mass(c: int, h: int, n: int = 0, o: int = 0) -> float:
    return c * MASS["C"] + h * MASS["H"] + n * MASS["N"] + o * MASS["O"]


COMPOUNDS = {
    "rhodotorulic_acid": ("primary_target", dict(c=20, h=32, n=4, o=6)),
    "N5-acetyl-N5-hydroxyornithine": ("precursor", dict(c=7, h=14, n=2, o=4)),
    "dimerumic_acid": ("related_siderophore", dict(c=22, h=37, n=3, o=9)),
    "ferrichrome": ("related_siderophore", dict(c=27, h=42, n=9, o=12)),
}

ADDUCTS = {
    "[M+H]+": lambda m: m + PROTON,
    "[M+Na]+": lambda m: m + NA,
    "[M+NH4]+": lambda m: m + NH4,
    "[M+K]+": lambda m: m + K,
}


def build_target_list() -> pd.DataFrame:
    rows = []
    for name, (category, f) in COMPOUNDS.items():
        m = formula_mass(**f)
        formula_str = f"C{f['c']}H{f['h']}N{f.get('n', 0)}" + (f"O{f['o']}" if f.get("o") else "")
        for adduct_name, fn in ADDUCTS.items():
            rows.append(dict(compound=name, category=category, formula=formula_str,
                              neutral_mass=m, adduct=adduct_name, target_mz=fn(m)))
    return pd.DataFrame(rows).sort_values("target_mz").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppm", type=float, default=20.0)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = build_target_list()
    targets.to_csv(OUT_DIR / "siderophore_target_list.csv", index=False)
    print(f"Built {len(targets)} targets ({len(COMPOUNDS)} compounds x {len(ADDUCTS)} adducts)")

    if not FULL_FEATURES_CSV.exists():
        raise SystemExit(f"Full feature table not found: {FULL_FEATURES_CSV}")
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
                compound=t["compound"], category=t["category"], adduct=t["adduct"], target_mz=t["target_mz"],
                row_id=h["row ID"], observed_mz=h["row m/z"], ppm_error=ppm_error,
                rt_min=h["row retention time"], has_ms2=h["has_ms2"], is_isf=h["is_isf"],
                is_default_adduct=h["is_default_adduct"], total_scans=h["total_scans"],
            )
            if dedup is not None:
                dgrp = dedup.loc[dedup["row ID"] == h["row ID"]]
                if len(dgrp):
                    row["dedup_group_id"] = dgrp.iloc[0]["dedup_group_id"]
                    row["is_group_representative"] = dgrp.iloc[0]["is_group_representative"]
            if sirius is not None and "row ID" in sirius.columns:
                sann = sirius.loc[sirius["row ID"] == h["row ID"]]
                if len(sann):
                    row["sirius_npc_class"] = sann.iloc[0].get("sirius_npc_class")
                    row["sirius_structure_name"] = sann.iloc[0].get("sirius_structure_name")
            matches.append(row)

    matches_df = pd.DataFrame(matches)
    out_path = OUT_DIR / "siderophore_mass_matches.csv"
    matches_df.to_csv(out_path, index=False)

    print(f"\n{len(matches_df)} raw feature(s) matched within {args.ppm} ppm across {len(COMPOUNDS)} compounds x {len(ADDUCTS)} adducts")
    if len(matches_df):
        print(matches_df.sort_values(["category", "compound"])[
            ["category", "compound", "adduct", "row_id", "observed_mz", "ppm_error", "rt_min", "total_scans", "has_ms2", "is_isf"]
        ].to_string(index=False))
        print("\nMatches per compound:")
        print(matches_df.groupby("compound").size().to_string())
    else:
        print("No matches -- none of these 4 compounds' target masses were detected at this ppm tolerance.")
    print(f"\nWrote {OUT_DIR / 'siderophore_target_list.csv'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
