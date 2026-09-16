#!/usr/bin/env python3
"""
Build the species-level table needed for the phylogenetic-signal test on
the N-acyl amino acid candidates (2026-09-16, PI request).

Response: mean log10(peak_area + 1) per species, CELL FRACTION ONLY --
the compartment analysis (nacyl_amino_acid_compartment_species.csv)
established these compounds are essentially absent from the supernatant,
so a species-level signal test on the supernatant fraction would mostly
be testing noise around zero.

Averages over strains within a species (not a phylogenetically-aware
collapse -- this matches the convention already used for
species_phenotype_table.csv in analysis/scripts/build_species_level_tables.py,
Phase 1). Species with only 1 strain (Cystobasidium sp., Pseudomicrostroma
phylloplanum, R. araucariae, R. sp. clade XIII) contribute a single
observation each -- flagged in the output, not a computed mean.

Usage:
    python3 analysis/scripts/nacyl_amino_acid_build_species_table.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
COMPARTMENT_CSV = REPO / "analysis" / "ahl_autoinducer_search" / "nacyl_amino_acid_compartment_species.csv"
OUT_CSV = REPO / "analysis" / "ahl_autoinducer_search" / "nacyl_amino_acid_species_table.csv"

COMPOUND_LABELS = {
    4109: "row_4109_palmitoyl_arginine",
    51126: "row_51126_myristoyl_arginine",
    51152: "row_51152_palmitoyl_lysine",
    24998: "row_24998_histidine_MH",
    40740: "row_40740_histidine_MNa",
}


def main():
    d = pd.read_csv(COMPARTMENT_CSV)
    d = d[d.fraction == "cell"].copy()
    d = d[d.species.notna()]
    # normalize "Rhodotorula sp_clade_X" (Cu_AUC crosswalk convention) to
    # "Rhodotorula sp. clade X" (species_tree.nwk tip-label convention, per
    # phylogenetic_signal.R's own species-name handling)
    d["species"] = d["species"].str.replace(
        r"^(Rhodotorula) sp_clade_(\w+)$", r"\1 sp. clade \2", regex=True
    )
    d["log_peak"] = np.log10(d["peak_area"] + 1)
    d["compound"] = d["row_id"].map(COMPOUND_LABELS)

    n_strains = d.groupby("species")["strain_code"].nunique().rename("n_strains")
    sp = d.groupby(["species", "compound"])["log_peak"].mean().unstack()
    sp = sp.join(n_strains)
    sp = sp.reset_index()
    sp.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV} ({len(sp)} species)")
    print(sp.to_string(index=False))


if __name__ == "__main__":
    main()
