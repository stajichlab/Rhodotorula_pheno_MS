#!/usr/bin/env python3
"""
N-acyl amino acid candidates -- compartment (cell vs. supernatant) and
strain/species variation analysis (2026-09-16, PI request).

Follows up on the 3 structurally-corroborated candidates from
NACYL_AMINO_ACID_SEARCH.md (SIRIUS's own structure call names the exact
target compound, not just a matching formula):
  - row 4109  -- N-C16:0-arginine  (SIRIUS: "Palmitoyl arginine")
  - row 51126 -- N-C14:0-arginine  (SIRIUS: "N2-(1-Oxotetradecyl)-L-arginine")
  - row 51152 -- N-C16:0-lysine    (SIRIUS: "N6-Palmitoyl lysine")
Also carries the 2 weaker histidine candidates (class-level "N-acyl
amines" match, no specific structure name) as supplementary, clearly
separated from the 3 strong candidates in all output.

For each candidate row:
1. Pulls raw per-sample peak areas from the raw EB feature table.
2. Maps MS sample columns to strain_code + fraction (cell/supernatant)
   via the Cu_AUC crosswalk (ID lookup only, same convention as
   ahl_strain_detection_vs_morphology.py -- NOT a phenotype source).
3. Cell vs. supernatant: paired Wilcoxon signed-rank test on
   log10(peak_area + 1), paired by strain (each strain contributes one
   cell value and one supernatant value) -- appropriate because cell and
   supernatant are two fractions of the SAME strain/culture, not
   independent samples.
4. Strain/species variation: reports per-species summary statistics
   (median, IQR, n) of log10(peak_area + 1) in each fraction, plus a
   simple Kruskal-Wallis test across species (descriptive, NOT
   phylogenetically corrected -- same caveat as this project's other
   descriptive/exploratory passes, e.g. ahl_strain_detection_vs_morphology.py).

Peak areas are raw (not TIC/abundance-normalized). A "not detected"
value (peak area == 0, i.e. below this feature's own blank-sample floor,
same definition as ahl_strain_detection_vs_morphology.py) is retained as
0 -- NOT dropped or imputed -- since presence/absence is itself
informative for a candidate signaling molecule search.

Usage:
    python3 analysis/scripts/nacyl_amino_acid_compartment_species_analysis.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
FULL_FEATURES_CSV = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9" / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output" / "feature_finding" / "feature_finding_results" / "aligned_features_ms2.csv"
)
CU_AUC_CROSSWALK = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "Cu_AUC.20260811.fixed.csv.gz"
OUT_DIR = REPO / "analysis" / "ahl_autoinducer_search"
OUT_CSV = OUT_DIR / "nacyl_amino_acid_compartment_species.csv"
DIAG_TXT = OUT_DIR / "nacyl_amino_acid_compartment_species_diagnostics.txt"

STRONG_CANDIDATES = {
    4109: "N-C16:0-arginine (Palmitoyl arginine)",
    51126: "N-C14:0-arginine (N-myristoyl-arginine)",
    51152: "N-C16:0-lysine (N6-Palmitoyl lysine)",
}
SUPPLEMENTARY_CANDIDATES = {
    24998: "N-C16:0-histidine [M+H]+ (class match only, no named structure)",
    40740: "N-C16:0-histidine [M+Na]+ (class match only, no named structure)",
}
ALL_ROW_IDS = {**STRONG_CANDIDATES, **SUPPLEMENTARY_CANDIDATES}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log_lines = []

    def log(msg):
        print(msg)
        log_lines.append(msg)

    header = pd.read_csv(FULL_FEATURES_CSV, nrows=0).columns
    sample_cols = [c for c in header if c.endswith(".mzML Peak area")]
    strain_sample_cols = [c for c in sample_cols if c.startswith("C_") or c.startswith("SUP_")]

    feat = pd.read_csv(FULL_FEATURES_CSV, usecols=["row ID"] + strain_sample_cols)
    feat = feat[feat["row ID"].isin(ALL_ROW_IDS)].set_index("row ID")
    missing = set(ALL_ROW_IDS) - set(feat.index)
    assert not missing, f"candidate row IDs not found in feature table: {missing}"

    # crosswalk: sample column -> strain_code + fraction ("No MS2 Data" placeholder excluded)
    cw = pd.read_csv(CU_AUC_CROSSWALK)
    cw_cell = cw[cw["MS2_SAMPLE_Cell"] != "No MS2 Data"]
    cw_sup = cw[cw["MS2_SAMPLE_Supernatant"] != "No MS2 Data"]
    cell_map = cw_cell.set_index(cw_cell["MS2_SAMPLE_Cell"] + ".mzML Peak area")[["SAMPLE_NAME", "SPECIES"]]
    sup_map = cw_sup.set_index(cw_sup["MS2_SAMPLE_Supernatant"] + ".mzML Peak area")[["SAMPLE_NAME", "SPECIES"]]

    long_rows = []
    for row_id in ALL_ROW_IDS:
        for col in strain_sample_cols:
            val = feat.loc[row_id, col]
            if col in cell_map.index:
                strain, species, fraction = cell_map.loc[col, "SAMPLE_NAME"], cell_map.loc[col, "SPECIES"], "cell"
            elif col in sup_map.index:
                strain, species, fraction = sup_map.loc[col, "SAMPLE_NAME"], sup_map.loc[col, "SPECIES"], "supernatant"
            else:
                continue
            long_rows.append(dict(row_id=row_id, strain_code=strain, species=species, fraction=fraction, peak_area=val))

    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(OUT_CSV, index=False)
    log(f"Wrote {OUT_CSV} ({len(long_df)} strain x fraction x candidate rows)")

    for row_id, label in ALL_ROW_IDS.items():
        tier = "STRONG" if row_id in STRONG_CANDIDATES else "supplementary"
        log(f"\n{'='*70}\nrow {row_id} [{tier}] -- {label}\n{'='*70}")
        d = long_df[long_df.row_id == row_id]

        # --- cell vs supernatant, paired by strain ---
        wide = d.pivot_table(index="strain_code", columns="fraction", values="peak_area")
        paired = wide.dropna(subset=["cell", "supernatant"]) if {"cell", "supernatant"}.issubset(wide.columns) else wide.iloc[0:0]
        n_pos_cell = (wide["cell"] > 0).sum() if "cell" in wide else 0
        n_pos_sup = (wide["supernatant"] > 0).sum() if "supernatant" in wide else 0
        log(f"Detected (peak area > 0): cell {n_pos_cell}/{wide['cell'].notna().sum() if 'cell' in wide else 0}, "
            f"supernatant {n_pos_sup}/{wide['supernatant'].notna().sum() if 'supernatant' in wide else 0}")
        if len(paired) >= 5:
            log_cell = np.log10(paired["cell"] + 1)
            log_sup = np.log10(paired["supernatant"] + 1)
            stat, p = stats.wilcoxon(log_cell, log_sup)
            direction = "cell > supernatant" if log_cell.median() > log_sup.median() else "supernatant > cell"
            log(f"Paired Wilcoxon signed-rank (n={len(paired)} strains with both fractions): "
                f"stat={stat:.1f}, p={p:.4f}, median log10(peak+1): cell={log_cell.median():.2f}, "
                f"sup={log_sup.median():.2f} ({direction})")
        else:
            log(f"Too few strains with both fractions detected (n={len(paired)}) for a paired test.")

        # --- species-level variation (descriptive, not phylogenetically corrected) ---
        for fraction in ["cell", "supernatant"]:
            dd = d[(d.fraction == fraction) & d.species.notna()].copy()
            dd["log_peak"] = np.log10(dd["peak_area"] + 1)
            counts = dd.groupby("species").size()
            big = counts[counts >= 3].index
            if len(big) < 2:
                log(f"[{fraction}] too few species with n>=3 for a species comparison.")
                continue
            groups = [dd.loc[dd.species == sp, "log_peak"].values for sp in big]
            h, p = stats.kruskal(*groups)
            log(f"[{fraction}] Kruskal-Wallis across {len(big)} species (n>=3 each): H={h:.2f}, p={p:.4f}")
            summary = dd[dd.species.isin(big)].groupby("species")["log_peak"].agg(["median", "count"]).sort_values("median", ascending=False)
            log(summary.to_string())

    with open(DIAG_TXT, "w") as fh:
        fh.write("\n".join(log_lines) + "\n")
    print(f"\nWrote {DIAG_TXT}")


if __name__ == "__main__":
    main()
