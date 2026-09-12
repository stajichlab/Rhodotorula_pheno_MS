#!/usr/bin/env python3
"""
AHL autoinducer / colony morphology hypothesis -- Phase 2: per-strain
detection of the candidate AHL mass-search hits (Phase 1,
ahl_targeted_mass_remining.py) cross-referenced against the colony-texture
morphology proxy (strain_texture_table.csv, from
salinity_texture_baseline / Salinity=0% table -- NOT the Cu_AUC table,
which is used here only as an ID crosswalk, see below).

Per PI direction (2026-09-11): this pass looks for POSITIVE HITS -- which
strains actually show detectable signal for the candidate AHL features,
and whether those skew toward smoother colony texture -- rather than
running this project's full phylogenetically-aware block-permutation
association framework (used elsewhere, e.g.
phase2_color_metabolome_association.py). That heavier framework is
explicitly deferred, not skipped because it's unnecessary --  any
apparent pattern here is descriptive/exploratory and has NOT been
checked for phylogenetic confounding (a few large clades, e.g.
R. mucilaginosa, dominate this panel and could produce a naive-looking
"hit" that is really one clade's texture + one clade's chemistry).

Steps:
1. Load the 103 raw AHL mass-search hits (ahl_mass_matches.csv) -- both
   long_chain (n>=12, the PI's stated hypothesis) and short_medium_chain
   (kept as an in-search comparison group).
2. Pull per-sample peak areas for those row IDs from the raw EB feature
   table (aligned_features_ms2.csv).
3. Map MS sample columns (C_<id>.mzML / SUP_<id>.mzML) to strain_code via
   data/metadata/EXFAB_UCR-005/Cu_AUC.20260811.fixed.csv.gz's
   Strain ID <-> MS2_SAMPLE_Cell/MS2_SAMPLE_Supernatant <-> SAMPLE_NAME
   columns -- an ID CROSSWALK ONLY, not a phenotype source.
4. Define "detected" per (row_id, sample) as peak area > that row's own
   max Blank_* peak area (a per-feature blank floor), not a fixed
   arbitrary cutoff or raw >0 -- a feature that is also present in the
   solvent blanks at comparable intensity is not a clean detection.
5. Build a smoothness composite from strain_texture_table.csv, the same
   z-score-sum pattern this project already uses for color
   (build_strain_phenotype_table.py's orange_score): higher
   AngularSecondMoment (energy/uniformity) and InverseDifferenceMoment
   (homogeneity) plus lower Contrast and Entropy all indicate a smoother-
   reading surface at this pixel scale, so
     smoothness_z = z(ASM) + z(IDM) - z(Contrast) - z(Entropy)
   Component z-scores kept as columns so the composite's weighting is
   auditable, same convention as orange_score.
6. Report, per strain: whether any long-chain AHL candidate was detected
   (cell and/or supernatant), how many distinct candidate features, and
   the strain's smoothness_z -- sorted to surface strains that are BOTH
   AHL-detected and high-smoothness (the pattern the hypothesis predicts)
   at the top.
7. One simple, transparent, NON-phylogenetic cross-tab (Fisher's exact,
   AHL-long-chain-detected yes/no x above/below-median smoothness) is
   reported as descriptive context only -- explicitly not a substitute
   for a real association test, and not adjusted for species/clade
   structure. Species composition of any apparent pattern is reported
   alongside it so an all-one-clade artifact is visible, not hidden.

Usage:
    python3 analysis/scripts/ahl_strain_detection_vs_morphology.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
AHL_MATCHES = REPO / "analysis" / "ahl_autoinducer_search" / "ahl_mass_matches.csv"
FULL_FEATURES_CSV = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9" / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output" / "feature_finding" / "feature_finding_results" / "aligned_features_ms2.csv"
)
CU_AUC_CROSSWALK = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "Cu_AUC.20260811.fixed.csv.gz"
TEXTURE_TABLE = REPO / "analysis" / "ahl_autoinducer_search" / "strain_texture_table.csv"
OUT_DIR = REPO / "analysis" / "ahl_autoinducer_search"
OUT_CSV = OUT_DIR / "ahl_strain_detection_vs_morphology.csv"
DIAG_TXT = OUT_DIR / "ahl_strain_detection_vs_morphology_diagnostics.txt"


def zscore(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=0)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log_lines = []

    def log(msg):
        print(msg)
        log_lines.append(msg)

    hits = pd.read_csv(AHL_MATCHES)
    row_ids = sorted(hits["row_id"].unique())
    log(f"AHL candidate row IDs: {len(row_ids)} ({(hits['category']=='long_chain').sum()} long_chain rows, "
        f"{hits['row_id'][hits['category']=='long_chain'].nunique()} distinct long_chain row IDs)")

    header = pd.read_csv(FULL_FEATURES_CSV, nrows=0).columns
    sample_cols = [c for c in header if c.endswith(".mzML Peak area")]
    blank_cols = [c for c in sample_cols if c.startswith("Blank_")]
    strain_sample_cols = [c for c in sample_cols if c.startswith("C_") or c.startswith("SUP_")]
    log(f"Sample columns: {len(strain_sample_cols)} strain (C_/SUP_) + {len(blank_cols)} blank")

    feat = pd.read_csv(FULL_FEATURES_CSV, usecols=["row ID"] + sample_cols)
    feat = feat[feat["row ID"].isin(row_ids)].copy()
    assert len(feat) == len(row_ids), (
        f"Expected {len(row_ids)} matched feature rows in the full table, found {len(feat)} -- "
        "row IDs from ahl_mass_matches.csv not all present in aligned_features_ms2.csv?"
    )

    # per-feature blank floor
    blank_floor = feat[blank_cols].max(axis=1)
    feat = feat.set_index("row ID")
    blank_floor.index = feat.index

    long_ids = sorted(set(hits.loc[hits["category"] == "long_chain", "row_id"]))

    # boolean detection matrix (rows=candidate feature row_id, cols=sample), long-chain rows only
    detected_long = feat.loc[long_ids, strain_sample_cols].gt(blank_floor.loc[long_ids], axis=0)
    n_detected_long_per_sample = detected_long.sum(axis=0)  # how many of the 41 long-chain candidates detected, per sample
    # per-sample MAX raw peak area across the 41 long-chain candidates -- a continuous
    # alternative to the binary call, needed because the binary call (below) turns out
    # to saturate at ~100% and carries no discriminating information (see results).
    max_peak_area_per_sample = feat.loc[long_ids, strain_sample_cols].max(axis=0)

    # crosswalk: sample column name -> strain_code + fraction
    # ("No MS2 Data" is a literal placeholder string for strains lacking an MS2
    # sample in either fraction -- 9-10 strains as of 2026-09; must be excluded
    # before building the map or those placeholder values collide as duplicate keys)
    cw = pd.read_csv(CU_AUC_CROSSWALK)
    cw_cell = cw[cw["MS2_SAMPLE_Cell"] != "No MS2 Data"]
    cw_sup = cw[cw["MS2_SAMPLE_Supernatant"] != "No MS2 Data"]
    cell_map = cw_cell.set_index(cw_cell["MS2_SAMPLE_Cell"] + ".mzML Peak area")["SAMPLE_NAME"]
    sup_map = cw_sup.set_index(cw_sup["MS2_SAMPLE_Supernatant"] + ".mzML Peak area")["SAMPLE_NAME"]
    sample_to_strain = pd.concat([cell_map, sup_map])
    sample_to_fraction = pd.Series(
        {c: ("cell" if c in cell_map.index else "supernatant") for c in sample_to_strain.index}
    )

    samp_df = pd.DataFrame({
        "sample_col": strain_sample_cols,
        "n_ahl_long_detected": n_detected_long_per_sample.reindex(strain_sample_cols).values,
        "max_peak_area": max_peak_area_per_sample.reindex(strain_sample_cols).values,
    })
    samp_df["strain_code"] = samp_df["sample_col"].map(sample_to_strain)
    samp_df["fraction"] = samp_df["sample_col"].map(sample_to_fraction)
    n_unmapped = samp_df["strain_code"].isna().sum()
    log(f"Sample columns with no strain_code in Cu_AUC crosswalk: {n_unmapped} (dropped)")
    samp_df = samp_df.dropna(subset=["strain_code"])

    # per-fraction detection flags (any of the 41 long-chain candidate features above its blank floor)
    cell_g = samp_df[samp_df["fraction"] == "cell"].groupby("strain_code")
    sup_g = samp_df[samp_df["fraction"] == "supernatant"].groupby("strain_code")
    strain_detect = pd.DataFrame({
        "ahl_long_detected_cell": cell_g["n_ahl_long_detected"].max().gt(0),
        "ahl_long_detected_sup": sup_g["n_ahl_long_detected"].max().gt(0),
        "n_ahl_long_features_cell": cell_g["n_ahl_long_detected"].max(),
        "n_ahl_long_features_sup": sup_g["n_ahl_long_detected"].max(),
        "max_peak_area_cell": cell_g["max_peak_area"].max(),
        "max_peak_area_sup": sup_g["max_peak_area"].max(),
    })
    # a strain is "tested" if it has an MS sample in at least one fraction; the
    # detected_{cell,sup} flags are NaN (not False) for the fraction it lacks,
    # so ahl_long_detected_any below stays a real tri-state (True/False/NaN),
    # not silently collapsed to False for untested strains.
    strain_detect["has_ms_data"] = strain_detect[["ahl_long_detected_cell", "ahl_long_detected_sup"]].notna().any(axis=1)
    strain_detect["ahl_long_detected_any"] = strain_detect[["ahl_long_detected_cell", "ahl_long_detected_sup"]].max(axis=1, skipna=True)
    strain_detect["max_peak_area_long_chain"] = strain_detect[["max_peak_area_cell", "max_peak_area_sup"]].max(axis=1, skipna=True)
    strain_detect = strain_detect.reset_index()

    n_tested = int(strain_detect["has_ms_data"].sum())
    n_detected = int(strain_detect.loc[strain_detect["has_ms_data"], "ahl_long_detected_any"].sum())
    log(f"Strains with MS data (either fraction): {n_tested}")
    log(f"Of those, strains with >=1 long-chain AHL candidate above its own blank floor: "
        f"{n_detected} / {n_tested} ({100*n_detected/n_tested:.1f}%)")
    if n_detected / n_tested > 0.95:
        log("*** Binary detection saturates near 100% -- this call carries essentially NO "
            "discriminating information. A blank-floor-only threshold is too permissive for a "
            "16k+-feature untargeted table: almost every sample has SOME peak above blank noise "
            "somewhere in these 41 narrow mass windows, real signal or not. The binary detection "
            "columns below are retained for transparency but should NOT be read as evidence of "
            "AHL production -- see the continuous max-intensity check instead. ***")

    # texture / smoothness composite
    tex = pd.read_csv(TEXTURE_TABLE)
    tex["smooth_asm_z"] = zscore(tex["AngularSecondMoment_median"])
    tex["smooth_idm_z"] = zscore(tex["InverseDifferenceMoment_median"])
    tex["smooth_contrast_z"] = zscore(tex["Contrast_median"])
    tex["smooth_entropy_z"] = zscore(tex["Entropy_median"])
    tex["smoothness_z"] = (
        tex["smooth_asm_z"] + tex["smooth_idm_z"] - tex["smooth_contrast_z"] - tex["smooth_entropy_z"]
    )

    out = tex.merge(strain_detect, on="strain_code", how="left")
    n_no_ms = (~out["has_ms_data"].fillna(False)).sum()
    log(f"\nStrains with texture data but no MS sample match at all (excluded from all MS-side stats, "
        f"kept in output table as NaN): {n_no_ms}")

    out.to_csv(OUT_CSV, index=False)
    log(f"\nWrote {OUT_CSV}")

    tested = out[out["has_ms_data"].fillna(False)].copy()

    # continuous check: does raw AHL-candidate intensity correlate with smoothness?
    # NOTE: max_peak_area_long_chain is raw, un-normalized peak area (no TIC/total-
    # abundance normalization applied here) -- a strain-level abundance-scaling
    # confound (documented elsewhere in this project, see
    # .living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md)
    # has NOT been ruled out for this quantity.
    rho, p_rho = stats.spearmanr(tested["max_peak_area_long_chain"], tested["smoothness_z"])
    log(f"\n=== Continuous check (descriptive, NOT phylogenetically controlled) ===")
    log(f"Spearman rho(max raw peak area of strongest long-chain AHL-mass candidate, smoothness_z) "
        f"= {rho:.3f}, p = {p_rho:.4f}, n = {len(tested)}")
    log("(raw peak area, not abundance-normalized -- a known confound in this project's other "
        "abundance-based comparisons; do not read a positive rho here as confirmed biology)")

    log(f"\n=== Top 20 strains by max long-chain AHL-candidate peak area, with smoothness_z ===")
    top20 = tested.sort_values("max_peak_area_long_chain", ascending=False).head(20)
    print(top20[["strain_code", "Species", "max_peak_area_long_chain", "smoothness_z", "n_colonies_used"]].to_string(index=False))
    log(f"median smoothness_z, top-20-by-intensity strains: {top20['smoothness_z'].median():.3f}")
    log(f"median smoothness_z, all tested strains: {tested['smoothness_z'].median():.3f}")
    log("\nSpecies composition of top-20-by-intensity strains "
        "(check for single-clade dominance before reading any pattern as hypothesis-consistent):")
    log(str(top20["Species"].value_counts()))

    # binary cross-tab retained for transparency, explicitly flagged as uninformative given saturation
    med = tested["smoothness_z"].median()
    above = tested["smoothness_z"] > med
    table = pd.crosstab(tested["ahl_long_detected_any"], above)
    log(f"\n=== Binary 2x2 (retained for transparency -- SEE SATURATION WARNING ABOVE, not informative) ===")
    log(str(table))
    if table.shape == (2, 2):
        odds, p = stats.fisher_exact(table)
        log(f"Fisher's exact: odds ratio={odds:.3f}, p={p:.4f}")

    with open(DIAG_TXT, "w") as fh:
        fh.write("\n".join(log_lines) + "\n")
    print(f"\nWrote {DIAG_TXT}")


if __name__ == "__main__":
    main()
