#!/usr/bin/env python3
"""
AHL autoinducer targeted mass re-mining (PI hypothesis: long-chain
N-acyl-L-homoserine lactones (AHLs) vary across Rhodotorula strains/
species and correlate with colony morphology).

Same design as analysis/scripts/idea1_targeted_mass_remining.py: an
exact-mass search of the EXISTING raw EB feature table (no new data, no
SIRIUS dependency), because SIRIUS's own class-level calls could miss or
mis-annotate a low-abundance signaling molecule the same way it
mis-annotated a carotenoid-mass feature in that earlier pass.

AHLs are NOT an endogenously-reported fungal metabolite class -- they are
canonical bacterial (Gram-negative, LuxI/LuxR) quorum-sensing signals. Any
hit here in axenic fungal-strain extracts would be either (a) a real,
previously-unreported case of fungal AHL production/mimicry, (b) trace
bacterial contamination of the culture/extract, or (c) a coincidental
isobaric match to an unrelated fungal metabolite -- this search cannot by
itself distinguish those; MS2 spectral confirmation against an authentic
AHL fragmentation pattern (loss of the homoserine lactone ring, m/z 102.055
[C4H8NO2]+ for the acylium-independent core fragment, plus the diagnostic
lactone immonium ion) is required before treating any hit as biological.

Target list is generated programmatically from the general AHL molecular
formula (not hardcoded per-compound, unlike the pigment script) because
the homologous series is long and regular:

    unsubstituted N-acyl-HSL   (Cn-HSL):        C(n+4) H(2n+5) N O3
    3-oxo-N-acyl-HSL           (3-oxo-Cn-HSL):  C(n+4) H(2n+3) N O4
    3-hydroxy-N-acyl-HSL       (3-hydroxy-Cn-HSL): C(n+4) H(2n+5) N O4

where n = number of carbons in the acyl chain (the "Cn" in standard AHL
nomenclature, e.g. C12-HSL = N-dodecanoyl-L-homoserine lactone). Formulas
verified against known reference compounds: C4-HSL C8H13NO3 (171.09 Da),
3-oxo-C6-HSL C10H15NO4 (213.10 Da, Vibrio fischeri autoinducer),
3-oxo-C12-HSL C16H27NO4 (297.20 Da, Pseudomonas aeruginosa LasI signal),
3-hydroxy-C14-HSL C18H33NO4 (327.24 Da, Burkholderia BHL).

n swept 4-18 (saturated acyl chains only; does not cover the
mono-unsaturated long-chain AHLs reported in some Rhizobiaceae, e.g.
cis-11-C16:1-HSL -- a documented gap, not an oversight). "long_chain" is
tagged for n>=12 per the PI's stated hypothesis; n<12 is retained as an
in-search negative range for comparison, not because short-chain AHLs are
expected to be absent.

Usage:
    python3 analysis/scripts/ahl_targeted_mass_remining.py --ppm 20
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
MASS = {"C": 12.0, "H": 1.0078250319, "O": 15.9949146221, "N": 14.0030740052, "e": 0.00054858}
PROTON = MASS["H"] - MASS["e"]
NA = 22.98976928 - MASS["e"]
NH4 = MASS["N"] + 4 * MASS["H"] - MASS["e"]

N_RANGE = range(4, 19)  # acyl chain length Cn, n=4..18
LONG_CHAIN_MIN_N = 12  # PI hypothesis threshold; n<12 kept in-search as a comparison range

ADDUCTS = {
    "[M+H]+": lambda m: m + PROTON,
    "[M+Na]+": lambda m: m + NA,
    "[M+NH4]+": lambda m: m + NH4,
}


def formula_mass(c: int, h: int, n: int = 1, o: int = 3) -> float:
    return c * MASS["C"] + h * MASS["H"] + n * MASS["N"] + o * MASS["O"]


def build_compounds() -> dict:
    compounds = {}
    for n in N_RANGE:
        chain_tag = "long_chain" if n >= LONG_CHAIN_MIN_N else "short_medium_chain"
        c = n + 4
        # unsubstituted
        h = 2 * n + 5
        compounds[f"C{n}-HSL"] = (chain_tag, dict(c=c, h=h, o=3))
        # 3-oxo
        h_oxo = 2 * n + 3
        compounds[f"3-oxo-C{n}-HSL"] = (chain_tag, dict(c=c, h=h_oxo, o=4))
        # 3-hydroxy
        compounds[f"3-hydroxy-C{n}-HSL"] = (chain_tag, dict(c=c, h=h, o=4))
    return compounds


COMPOUNDS = build_compounds()


def build_target_list() -> pd.DataFrame:
    rows = []
    for name, (category, f) in COMPOUNDS.items():
        m = formula_mass(**f)
        formula_str = f"C{f['c']}H{f['h']}NO{f['o']}"
        for adduct_name, fn in ADDUCTS.items():
            rows.append(dict(compound=name, category=category, formula=formula_str,
                              neutral_mass=m, adduct=adduct_name, target_mz=fn(m)))
    return pd.DataFrame(rows).sort_values("target_mz").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppm", type=float, default=20.0, help="Mass tolerance in ppm for the m/z match.")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = build_target_list()
    targets.to_csv(OUT_DIR / "ahl_target_list.csv", index=False)
    n_compounds = len(COMPOUNDS)
    print(f"Built {len(targets)} targets ({n_compounds} compounds x {len(ADDUCTS)} adducts, "
          f"acyl chain lengths C{min(N_RANGE)}-C{max(N_RANGE)})")

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
                compound=t["compound"], category=t["category"], formula=t["formula"], adduct=t["adduct"],
                target_mz=t["target_mz"], row_id=h["row ID"], observed_mz=h["row m/z"], ppm_error=ppm_error,
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
    out_path = OUT_DIR / "ahl_mass_matches.csv"
    matches_df.to_csv(out_path, index=False)

    print(f"\n{len(matches_df)} raw feature(s) matched within {args.ppm} ppm across {n_compounds} compounds x {len(ADDUCTS)} adducts")
    if len(matches_df):
        print(matches_df.sort_values(["category", "compound"])[
            ["category", "compound", "formula", "adduct", "row_id", "observed_mz", "ppm_error", "rt_min", "total_scans", "has_ms2", "is_isf"]
        ].to_string(index=False))
        print("\nMatches per category:")
        print(matches_df.groupby("category").size().to_string())
    else:
        print("No matches -- see console output above for target m/z range checked.")
    print(f"\nWrote {OUT_DIR / 'ahl_target_list.csv'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
