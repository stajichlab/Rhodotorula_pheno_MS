#!/usr/bin/env python3
"""
N-acyl amino acid targeted mass re-mining -- expansion of the AHL
autoinducer search (analysis/ahl_autoinducer_search/) to a chemically
related signaling-molecule class, per PI request (2026-09-16).

N-acyl amino acids (a fatty acyl chain condensed onto an amino acid's free
amine, releasing H2O) are made by the SAME acyl-CoA/acyl-ACP donor
chemistry as AHLs, and some LuxI-family acyltransferases have been shown
to make N-acyl amino acids instead of (or alongside) AHLs -- e.g. the
Brady lab's soil-metagenome-derived "N-acyl amino acid synthases," and
gut-commensal-bacteria-derived N-acyl amino acids that act as host GPCR
ligands (Cohen et al. 2017, Cell). This is a real, distinct literature
from the AHL literature -- N-acyl amino acids are not AHLs and are not
claimed to be quorum-sensing molecules in the same LuxI/LuxR sense, but
they are a chemically adjacent, plausibly co-produced signaling-molecule
family worth searching for given the AHL search's null result.

Same search design as ahl_targeted_mass_remining.py: exact-mass re-mining
of the raw EB feature table (aligned_features_ms2.csv), independent of
SIRIUS's own class calls, at 20 ppm tolerance.

Target list construction (systematic, not cherry-picked): all 20 standard
proteinogenic amino acids, PLUS homoserine (the AHL lactone's own
non-cyclized backbone -- included because this project already
established R. mucilaginosa has lactonase activity (see
AHL_AUTOINDUCER_SEARCH.md's literature-check section); if that activity
acts on any bacterial AHL contaminant, the open-chain N-acyl-homoserine
hydrolysis product is a directly motivated additional target, not merely
a symmetric extension of the amino-acid list). Every amino acid x every
acyl chain length (saturated, Cn n=4-18, same range as the AHL search) x
3 adducts.

N-acyl amino acid neutral formula, derived from first principles (amide
condensation of a saturated fatty acid CnH(2n)O2 with an amino acid
C(a)H(h)N(1)O(o), losing H2O):
    C(n+a) H(2n+h-2) N(1) O(o+1)
Verified against a known reference compound: N-lauroylglycine (n=12,
glycine C2H5NO2) predicts C14H27NO3, matching the published formula for
this compound (MW 257.37).

Adducts: [M+H]+, [M+Na]+, [M+NH4]+ -- same positive-mode set as the AHL
search. NOT searched: negative mode ([M-H]-, which N-acyl amino acids'
free carboxylic acid group would be expected to favor in ESI) -- checked
and confirmed this project's raw feature table (aligned_features_ms2.csv)
contains ONLY positive-mode adducts (verified via the `adduct` column,
2026-09-16); negative-mode LC-MS data does not exist for this sample set,
so this is a real, unclosable gap with the current data, not an oversight.
NOT searched: unsaturated or hydroxylated acyl chains, or non-standard
amino acids beyond homoserine -- a documented scope choice, not
exhaustive coverage of the N-acyl amino acid structural space.

Usage:
    python3 analysis/scripts/nacyl_amino_acid_mass_remining.py --ppm 20
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
OUT_DIR = REPO / "analysis" / "ahl_autoinducer_search"

# monoisotopic atomic masses
MASS = {"C": 12.0, "H": 1.0078250319, "O": 15.9949146221, "N": 14.0030740052, "S": 31.97207069, "e": 0.00054858}
PROTON = MASS["H"] - MASS["e"]
NA = 22.98976928 - MASS["e"]
NH4 = MASS["N"] + 4 * MASS["H"] - MASS["e"]

N_RANGE = range(4, 19)  # acyl chain length Cn, n=4..18 -- same range as the AHL search
LONG_CHAIN_MIN_N = 12   # kept consistent with the AHL search's long_chain threshold

# Standard proteinogenic amino acid neutral monoisotopic formulas (C, H, N, O, S),
# free amino acid form. Plus homoserine (non-proteinogenic; see module docstring
# for why it's included).
AMINO_ACIDS = {
    # name: (C, H, N, O, S)
    "glycine": (2, 5, 1, 2, 0),
    "alanine": (3, 7, 1, 2, 0),
    "serine": (3, 7, 1, 3, 0),
    "cysteine": (3, 7, 1, 2, 1),
    "valine": (5, 11, 1, 2, 0),
    "leucine": (6, 13, 1, 2, 0),
    "isoleucine": (6, 13, 1, 2, 0),
    "proline": (5, 9, 1, 2, 0),
    "methionine": (5, 11, 1, 2, 1),
    "threonine": (4, 9, 1, 3, 0),
    "phenylalanine": (9, 11, 1, 2, 0),
    "tyrosine": (9, 11, 1, 3, 0),
    "tryptophan": (11, 12, 2, 2, 0),
    "aspartate": (4, 7, 1, 4, 0),
    "glutamate": (5, 9, 1, 4, 0),
    "asparagine": (4, 8, 2, 3, 0),
    "glutamine": (5, 10, 2, 3, 0),
    "lysine": (6, 14, 2, 2, 0),
    "arginine": (6, 14, 4, 2, 0),
    "histidine": (6, 9, 3, 2, 0),
    "homoserine": (4, 9, 1, 3, 0),  # non-proteinogenic; see docstring
}

ADDUCTS = {
    "[M+H]+": lambda m: m + PROTON,
    "[M+Na]+": lambda m: m + NA,
    "[M+NH4]+": lambda m: m + NH4,
}


def formula_mass(c: int, h: int, n: int, o: int, s: int = 0) -> float:
    return c * MASS["C"] + h * MASS["H"] + n * MASS["N"] + o * MASS["O"] + s * MASS["S"]


def build_target_list() -> pd.DataFrame:
    rows = []
    for aa_name, (aa_c, aa_h, aa_n, aa_o, aa_s) in AMINO_ACIDS.items():
        for n in N_RANGE:
            chain_tag = "long_chain" if n >= LONG_CHAIN_MIN_N else "short_medium_chain"
            c = n + aa_c
            h = 2 * n + aa_h - 2
            o = aa_o + 1
            m = formula_mass(c=c, h=h, n=aa_n, o=o, s=aa_s)
            formula_str = f"C{c}H{h}N{aa_n}O{o}" + (f"S{aa_s}" if aa_s else "")
            compound = f"N-C{n}:0-{aa_name}"
            for adduct_name, fn in ADDUCTS.items():
                rows.append(dict(
                    compound=compound, amino_acid=aa_name, acyl_n=n, category=chain_tag,
                    formula=formula_str, neutral_mass=m, adduct=adduct_name, target_mz=fn(m),
                ))
    return pd.DataFrame(rows).sort_values("target_mz").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppm", type=float, default=20.0, help="Mass tolerance in ppm for the m/z match.")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = build_target_list()
    targets.to_csv(OUT_DIR / "nacyl_amino_acid_target_list.csv", index=False)
    n_compounds = len(AMINO_ACIDS) * len(N_RANGE)
    print(f"Built {len(targets)} targets ({n_compounds} amino-acid x chain-length combinations x "
          f"{len(ADDUCTS)} adducts, {len(AMINO_ACIDS)} amino acids, acyl chain lengths "
          f"C{min(N_RANGE)}-C{max(N_RANGE)})")

    if not FULL_FEATURES_CSV.exists():
        raise SystemExit(f"Full feature table not found: {FULL_FEATURES_CSV}")
    feat = pd.read_csv(FULL_FEATURES_CSV, usecols=[
        "row ID", "row m/z", "row retention time", "adduct", "is_default_adduct",
        "is_isf", "isf_parent_id", "adduct_source_id", "has_ms2", "total_scans",
    ])

    dedup = pd.read_csv(DEDUP_GROUPS) if DEDUP_GROUPS.exists() else None

    sirius = None
    if SIRIUS_ANNOTATIONS.exists():
        sirius = pd.read_csv(SIRIUS_ANNOTATIONS, sep="\t")

    matches = []
    for _, t in targets.iterrows():
        tol = t["target_mz"] * args.ppm / 1e6
        hits = feat[(feat["row m/z"] >= t["target_mz"] - tol) & (feat["row m/z"] <= t["target_mz"] + tol)]
        for _, h in hits.iterrows():
            ppm_error = (h["row m/z"] - t["target_mz"]) / t["target_mz"] * 1e6
            row = dict(
                compound=t["compound"], amino_acid=t["amino_acid"], acyl_n=t["acyl_n"],
                category=t["category"], formula=t["formula"], adduct=t["adduct"], target_mz=t["target_mz"],
                row_id=h["row ID"], observed_mz=h["row m/z"], ppm_error=ppm_error,
                rt_min=h["row retention time"], has_ms2=h["has_ms2"], is_isf=h["is_isf"],
                is_default_adduct=h["is_default_adduct"], total_scans=h["total_scans"],
            )
            if dedup is not None:
                dgrp = dedup.loc[dedup["row ID"] == h["row ID"]]
                if len(dgrp):
                    row["dedup_group_id"] = dgrp.iloc[0]["dedup_group_id"]
                    row["is_group_representative"] = dgrp.iloc[0]["is_group_representative"]
            if sirius is not None:
                sann = sirius.loc[sirius["row ID"] == h["row ID"]] if "row ID" in sirius.columns else None
                if sann is not None and len(sann):
                    row["sirius_npc_class"] = sann.iloc[0].get("sirius_npc_class")
                    row["sirius_structure_name"] = sann.iloc[0].get("sirius_structure_name")
                    row["sirius_formula"] = sann.iloc[0].get("sirius_formula")
            matches.append(row)

    matches_df = pd.DataFrame(matches)
    out_path = OUT_DIR / "nacyl_amino_acid_mass_matches.csv"
    matches_df.to_csv(out_path, index=False)

    print(f"\n{len(matches_df)} raw feature(s) matched within {args.ppm} ppm across "
          f"{n_compounds} amino-acid/chain-length combinations x {len(ADDUCTS)} adducts")
    if len(matches_df):
        print("\nMatches per amino acid:")
        print(matches_df.groupby("amino_acid").size().sort_values(ascending=False).to_string())
        print("\nMatches per category:")
        print(matches_df.groupby("category").size().to_string())
        n_distinct_rows = matches_df["row_id"].nunique()
        print(f"\nDistinct feature rows matched (a row can match >1 amino acid/adduct combo): {n_distinct_rows}")
    else:
        print("No matches.")
    print(f"\nWrote {OUT_DIR / 'nacyl_amino_acid_target_list.csv'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
