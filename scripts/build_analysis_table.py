#!/usr/bin/env python3
"""
Link phenotype metadata, strain identity, and MS2 feature abundances into
a single analysis-ready pair of tables:

- analysis/linked_data/sample_metadata.csv.gz
    One row per MS2 sample (C_<id> or SUP_<id>), with fraction, strain
    code, species, and the YPD2 color/AUC phenotypes joined in.
- analysis/linked_data/feature_abundance_matrix.csv.gz
    Feature (row) x sample (column) peak-area matrix, subset to only the
    samples in sample_metadata.csv.gz, with feature annotation columns
    (m/z, retention time, adduct, etc.) kept alongside as the first
    columns.

Join logic: EB_..._merged_metadata.fixed.tsv.gz's canonical_strain column
(built by fix_metadata_strain_ids.py from the numeric Strain ID <-> C_/SUP_
filename correspondence) is matched exactly against YPD2's fixed Strain
column -- both were derived from the same corrected source, so this is an
exact string join, not fuzzy matching. Rows with no canonical_strain
(controls, and the one unresolved JBEI-13807 sample) are dropped.

Usage:
    python3 scripts/build_analysis_table.py
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
YPD2_FIXED = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "YPD2_phenotypic.20260702.fixed.csv.gz"
CU_AUC_FIXED = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "Cu_AUC.20260811.fixed.csv.gz"
EB_MERGED_FIXED = (
    REPO
    / "data"
    / "processed"
    / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "Exfab_-_20260130_ExFAB_Rhodo_--b773ffa18c2b41e5a3484526293a54f9-merged_metadata.fixed.tsv.gz"
)
EB_FEATURES = (
    REPO
    / "data"
    / "processed"
    / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "aligned_features_ms2.csv.zst"
)
OUT_DIR = REPO / "analysis" / "linked_data"
FEATURE_ANNOTATION_COLS = [
    "row ID",
    "row m/z",
    "row retention time",
    "adduct",
    "is_default_adduct",
    "has_ms2",
    "detection_count",
    "detection_rate",
]


def build_sample_metadata() -> pd.DataFrame:
    eb = pd.read_csv(EB_MERGED_FIXED, sep="\t", dtype=str)
    eb = eb[(eb["ATTRIBUTE_TYPE"] == "sample") & (eb["canonical_strain"].fillna("") != "")].copy()
    eb["sample_id"] = eb["filename"].str.removesuffix(".mzML")
    eb["fraction"] = eb["sample_id"].str.split("_").str[0].map(
        {"C": "cell", "SUP": "supernatant"}
    )
    # Join key is the numeric id shared between "C_<id>"/"SUP_<id>" and YPD2's
    # "Strain ID", NOT the canonical_strain text: a couple of strain *codes*
    # are legitimately reused across two distinct Strain ID rows (e.g.
    # DBVPG_3853 at ids 192 and 195), so a text join fans out one EB sample
    # into two rows. The numeric id is always unique per sample.
    eb["Strain ID"] = eb["sample_id"].str.split("_").str[1]
    eb = eb[["sample_id", "fraction", "Strain ID", "canonical_strain", "ATTRIBUTE_POSITION_PLATE"]]
    eb = eb.rename(columns={"ATTRIBUTE_POSITION_PLATE": "position_plate"})

    ypd2 = pd.read_csv(YPD2_FIXED, dtype=str)
    pheno_cols = [
        "Strain ID",
        "Strain",
        "Species",
        "Origin",
        "Environment",
        "Library Plate",
        "Media",
        "Median_ColorLab_L*Mean",
        "Median_ColorLab_a*Mean",
        "Median_ColorLab_b*Mean",
        "Mean_ColorLab_L*Mean",
        "Mean_ColorLab_a*Mean",
        "Mean_ColorLab_b*Mean",
        "Median_ColorHSV_HueMean",
        "Median_ColorHSV_SaturationMean",
        "Median_ColorHSV_BrightnessMean",
    ]
    ypd2 = ypd2[pheno_cols]
    for col in pheno_cols:
        if col not in ("Strain ID", "Strain", "Species", "Origin", "Environment", "Media"):
            ypd2[col] = pd.to_numeric(ypd2[col], errors="coerce")

    cu_auc = pd.read_csv(CU_AUC_FIXED, dtype=str)[["Strain ID", "mean_auc_rate"]]
    cu_auc["mean_auc_rate"] = pd.to_numeric(cu_auc["mean_auc_rate"], errors="coerce")

    merged = eb.merge(ypd2, on="Strain ID", how="left", indicator=True)
    n_unmatched = (merged["_merge"] == "left_only").sum()
    merged = merged.drop(columns=["_merge"])
    merged = merged.merge(cu_auc, on="Strain ID", how="left")

    if n_unmatched:
        print(
            f"WARNING: {n_unmatched} EB sample rows had no phenotype match "
            "(unexpected given prior congruence fix -- investigate)",
            file=sys.stderr,
        )
    n_blank_species = merged["Species"].isna().sum()
    if n_blank_species:
        print(
            f"NOTE: {n_blank_species} sample rows joined to a YPD2 strain "
            "with blank Species/Origin/Environment (incompletely-catalogued "
            "record, not a join failure)",
            file=sys.stderr,
        )
    disagree = merged[merged["canonical_strain"] != merged["Strain"]]
    if len(disagree):
        print(
            f"WARNING: {len(disagree)} rows where canonical_strain != YPD2 Strain "
            "for the same Strain ID -- investigate",
            file=sys.stderr,
        )

    merged = merged.drop(columns=["Strain"]).sort_values(["fraction", "sample_id"])
    return merged


def build_feature_matrix(sample_ids: list[str]) -> pd.DataFrame:
    features = pd.read_csv(EB_FEATURES, compression="zstd")
    peak_cols = {f"{s}.mzML Peak area": s for s in sample_ids}
    missing = [s for s in sample_ids if f"{s}.mzML Peak area" not in features.columns]
    if missing:
        print(f"WARNING: {len(missing)} sample(s) missing from feature table: {missing}", file=sys.stderr)
    keep_cols = FEATURE_ANNOTATION_COLS + [c for c in peak_cols if c in features.columns]
    subset = features[keep_cols].rename(columns=peak_cols)
    return subset


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sample_meta = build_sample_metadata()
    sample_meta_path = OUT_DIR / "sample_metadata.csv.gz"
    sample_meta.to_csv(sample_meta_path, index=False)
    print(
        f"sample_metadata: {len(sample_meta)} samples "
        f"({(sample_meta['fraction'] == 'cell').sum()} cell, "
        f"{(sample_meta['fraction'] == 'supernatant').sum()} supernatant), "
        f"{sample_meta['canonical_strain'].nunique()} distinct strains -> {sample_meta_path}",
        file=sys.stderr,
    )

    feature_matrix = build_feature_matrix(sample_meta["sample_id"].tolist())
    feature_matrix_path = OUT_DIR / "feature_abundance_matrix.csv.gz"
    feature_matrix.to_csv(feature_matrix_path, index=False)
    n_sample_cols = feature_matrix.shape[1] - len(FEATURE_ANNOTATION_COLS)
    print(
        f"feature_abundance_matrix: {len(feature_matrix)} features x "
        f"{n_sample_cols} samples -> {feature_matrix_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
