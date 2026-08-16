#!/usr/bin/env python3
"""
Idea 1, Phase 1 (analysis/ideas/2026-08-15-color-metabolome-genome-null-brainstorm/
IDEA1_CHEMIST_FRAMEWORK.md): predicted-mass targeted re-mining of the
EXISTING raw MS2 feature table -- no new data, no SIRIUS dependency. This
searches the fuller EB feature table directly by exact mass, not by the
untargeted pipeline's own (SIRIUS-based) compound calling, so it can
catch a real pigment/terpenoid feature that got mis-annotated, orphaned
as an unlabeled in-source fragment, or split across adducts the
untargeted dedup left unmerged.

EXPANDED 2026-08-16 (PI request, after the first pass's 2 standout
candidates turned out inconclusive on MS2 inspection -- see
MS2_FRAGMENTATION_CHECK.md) beyond the core carotenogenesis pathway to
three additional compound categories, each chosen for a specific reason:

1. carotenoid_pathway (original 9, unchanged): the linear desaturation
   series (each step removes 2H) + cyclization products + the
   Rhodotorula-specific torulene/torularhodin tailoring step.
       phytoene -> phytofluene -> zeta-carotene -> neurosporene -> lycopene
       -> [cyclization, same mass] gamma-carotene / beta-carotene
       -> torulene -> torularhodin

2. apocarotenoids: oxidative CLEAVAGE fragments of the full C40
   carotenoids -- smaller, more polar, and therefore more likely to
   ionize well by ESI than the intact backbones (this project's SIRIUS
   annotations already called one existing feature "Apocarotenoids(ε-)",
   so this category has independent supporting evidence already).
   Curated to 3 well-characterized, literature-standard short-chain
   cleavage products (NOT an exhaustive combinatorial cleavage-position
   enumeration, which would need dedicated curation beyond this pass):
   beta-ionone, beta-cyclocitral, dihydroactinidiolide.

3. additional_carotenoids: oxygenated carotenoid variants reported in
   other carotenogenic basidiomycetous yeasts (astaxanthin, canthaxanthin,
   echinenone, zeaxanthin) -- included as a broader net in case
   Rhodotorula makes a minor/alternative oxygenated product not on the
   canonical torulene/torularhodin branch. NOT strain-confirmed for
   Rhodotorula specifically -- treat any hit here as needing a literature
   check on whether it's actually plausible for this genus, not just
   chemically detectable.

4. sterol_pathway: ergosterol, squalene, lanosterol -- NOT pigments, but
   share the same upstream isoprenoid/MVA precursor pool as carotenoids
   (HMG-CoA reductase, GGPPS -- see DEVELOPMENT_PLAN.md Part B's
   candidate gene table). Included to test the "precursor-competition"
   hypothesis: if a strain diverts more flux to sterols, less may be
   available for carotenoids, independent of anything downstream. For
   sterols specifically, [M+H-H2O]+ (loss of the 3-OH) is often the
   DOMINANT ion, not [M+H]+ -- worth weighting accordingly when reading
   results.

Adducts searched: [M+H]+, [M+Na]+, [M+NH4]+, [M+H-H2O]+. NOT searched:
[M]+* (radical cation, the classic APCI-preferred carotenoid ion -- this
pipeline is ESI-based) and negative mode (relevant for the carboxylic
acids: torularhodin, and none of the others here). These remain real gaps
in this search, not oversights to ignore when interpreting a null result.

Usage:
    python3 analysis/scripts/idea1_targeted_mass_remining.py --ppm 20
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
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase3_metabolome_phenotype_idea1"

# monoisotopic atomic masses
MASS = {"C": 12.0, "H": 1.0078250319, "O": 15.9949146221, "N": 14.0030740052, "e": 0.00054858}
PROTON = MASS["H"] - MASS["e"]
NA = 22.98976928 - MASS["e"]
NH4 = MASS["N"] + 4 * MASS["H"] - MASS["e"]


def formula_mass(c: int, h: int, o: int = 0) -> float:
    return c * MASS["C"] + h * MASS["H"] + o * MASS["O"]


COMPOUNDS = {
    # category, formula (c, h, o)
    "phytoene": ("carotenoid_pathway", dict(c=40, h=64, o=0)),
    "phytofluene": ("carotenoid_pathway", dict(c=40, h=62, o=0)),
    "zeta-carotene": ("carotenoid_pathway", dict(c=40, h=60, o=0)),
    "neurosporene": ("carotenoid_pathway", dict(c=40, h=58, o=0)),
    "lycopene": ("carotenoid_pathway", dict(c=40, h=56, o=0)),
    "gamma-carotene": ("carotenoid_pathway", dict(c=40, h=56, o=0)),  # isomer of lycopene, same mass
    "beta-carotene": ("carotenoid_pathway", dict(c=40, h=56, o=0)),  # isomer of lycopene, same mass
    "torulene": ("carotenoid_pathway", dict(c=40, h=54, o=0)),
    "torularhodin": ("carotenoid_pathway", dict(c=40, h=52, o=2)),
    # apocarotenoids: curated short-chain cleavage products (not exhaustive)
    "beta-ionone": ("apocarotenoid", dict(c=13, h=20, o=1)),
    "beta-cyclocitral": ("apocarotenoid", dict(c=10, h=16, o=1)),
    "dihydroactinidiolide": ("apocarotenoid", dict(c=11, h=16, o=2)),
    # additional oxygenated carotenoids reported in other carotenogenic yeasts
    # (not strain-confirmed for Rhodotorula -- treat hits as leads, not expectations)
    "astaxanthin": ("additional_carotenoid", dict(c=40, h=52, o=4)),
    "canthaxanthin": ("additional_carotenoid", dict(c=40, h=52, o=2)),  # isobaric with torularhodin
    "echinenone": ("additional_carotenoid", dict(c=40, h=54, o=1)),
    "zeaxanthin": ("additional_carotenoid", dict(c=40, h=56, o=2)),
    # sterol/MVA-pathway precursor-competition markers (not pigments)
    "ergosterol": ("sterol_pathway", dict(c=28, h=44, o=1)),
    "squalene": ("sterol_pathway", dict(c=30, h=50, o=0)),
    "lanosterol": ("sterol_pathway", dict(c=30, h=50, o=1)),
}

ADDUCTS = {
    "[M+H]+": lambda m: m + PROTON,
    "[M+Na]+": lambda m: m + NA,
    "[M+NH4]+": lambda m: m + NH4,
    "[M+H-H2O]+": lambda m: m + PROTON - (2 * MASS["H"] + MASS["O"]),
}


def build_target_list() -> pd.DataFrame:
    rows = []
    for name, (category, f) in COMPOUNDS.items():
        m = formula_mass(**f)
        formula_str = f"C{f['c']}H{f['h']}" + (f"O{f['o']}" if f["o"] else "")
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
    targets.to_csv(OUT_DIR / "carotenoid_pathway_target_list.csv", index=False)
    n_compounds = len(COMPOUNDS)
    print(f"Built {len(targets)} targets ({n_compounds} compounds x {len(ADDUCTS)} adducts, "
          f"{len(set(c for c, _ in COMPOUNDS.values()))} categories)")

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
            if sirius is not None:
                sann = sirius.loc[sirius["row ID"] == h["row ID"]] if "row ID" in sirius.columns else None
                if sann is not None and len(sann):
                    row["sirius_npc_class"] = sann.iloc[0].get("sirius_npc_class")
                    row["sirius_structure_name"] = sann.iloc[0].get("sirius_structure_name")
            matches.append(row)

    matches_df = pd.DataFrame(matches)
    out_path = OUT_DIR / "carotenoid_pathway_mass_matches.csv"
    matches_df.to_csv(out_path, index=False)

    print(f"\n{len(matches_df)} raw feature(s) matched within {args.ppm} ppm across {n_compounds} compounds x {len(ADDUCTS)} adducts")
    if len(matches_df):
        print(matches_df.sort_values(["category", "compound"])[
            ["category", "compound", "adduct", "row_id", "observed_mz", "ppm_error", "rt_min", "total_scans", "has_ms2", "is_isf"]
        ].to_string(index=False))
        print("\nMatches per category:")
        print(matches_df.groupby("category").size().to_string())
    print(f"\nWrote {OUT_DIR / 'carotenoid_pathway_target_list.csv'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
