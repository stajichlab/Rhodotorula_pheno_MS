#!/usr/bin/env python3
"""Summarize abundance of rhodotorulic acid / siderophore-related MS2 features
across species and strains, using existing SIRIUS structure annotations and
the aligned MS2 feature quant matrix (row ID = feature ID, shared key).

Inputs (already produced by prior work in this repo -- not regenerated here):
  - analysis/sirius_annotation/sirius_annotations.tsv
      SIRIUS/CANOPUS structure predictions per feature (row ID).
  - data/processed/EB_20260130_ExFAB_Rhodo_Sup_and_Pellet/<hash>/aligned_features_ms2.csv.zst
      MZmine/GNPS-style aligned feature table: row ID x per-sample peak area.
  - data/processed/EB_20260130_ExFAB_Rhodo_Sup_and_Pellet/<hash>/*merged_metadata.fixed.tsv.gz
      Sample (filename) -> species / strain / source (cell_pellett vs supernatant) crosswalk.

Output:
  - outputs/siderophore_feature_annotations.tsv   (1 row per candidate feature)
  - outputs/siderophore_abundance_by_species.tsv  (feature x species, mean/median/detection-rate, split by source)
  - outputs/siderophore_abundance_by_strain.tsv   (feature x strain, mean peak area, split by source)
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
SIRIUS_TSV = REPO / "analysis/sirius_annotation/sirius_annotations.tsv"
PROC_DIR = REPO / "data/processed/EB_20260130_ExFAB_Rhodo_Sup_and_Pellet/b773ffa18c2b41e5a3484526293a54f9"
FEATURES_ZST = PROC_DIR / "aligned_features_ms2.csv.zst"
METADATA_GZ = PROC_DIR / "Exfab_-_20260130_ExFAB_Rhodo_--b773ffa18c2b41e5a3484526293a54f9-merged_metadata.fixed.tsv.gz"
OUT_DIR = Path(__file__).resolve().parents[1] / "outputs"

# Keyword search is intentionally broad -- siderophore chemistry in fungi is
# dominated by hydroxamate siderophores (desferrioxamine/ferrioxamine family,
# named after the actinobacterial molecules SIRIUS structure DB uses as the
# closest reference structures) plus rhodotorulate-class siderophores.
SIDEROPHORE_KEYWORDS = [
    "rhodotorulic",
    "rhodotorulate",
    "siderophore",
    "ferrichrome",
    "ferrioxamine",
    "desferrioxamine",
    "deferoxamine",
    "coprogen",
    "fusarinine",
    "dimerum",
    "basidiochrome",
    "rhizoferrin",
    "hydroxamate",
]


def load_sirius_annotations() -> pd.DataFrame:
    df = pd.read_csv(SIRIUS_TSV, sep="\t")
    if "row ID" not in df.columns:
        raise ValueError(f"unexpected sirius_annotations.tsv schema: {df.columns.tolist()}")
    # ANALYSIS_OK[missingness]: SIRIUS leaves name/pathway/class blank when it made
    # no structure/class call for a feature; fill with "" purely so the keyword
    # regex search below has a string to test, not to impute a real value. No
    # downstream numeric aggregation touches these columns.
    name_col = df["sirius_structure_name"].fillna("")  # ANALYSIS_OK[missingness]: blank = no SIRIUS structure call; "" is a search-safe sentinel, not an imputed value
    class_cols = (
        df["sirius_npc_pathway"].fillna("")  # ANALYSIS_OK[missingness]: blank = no NPC pathway call; same sentinel rationale as above
        + " "
        + df["sirius_npc_class"].fillna("")  # ANALYSIS_OK[missingness]: blank = no NPC class call; same sentinel rationale as above
        + " "
        + df["sirius_classyfire_class"].fillna("")  # ANALYSIS_OK[missingness]: blank = no ClassyFire call; same sentinel rationale as above
    )
    pattern = "|".join(SIDEROPHORE_KEYWORDS)
    hit = name_col.str.contains(pattern, case=False, regex=True) | class_cols.str.contains(
        pattern, case=False, regex=True
    )
    hits = df.loc[hit].copy()
    hits = hits.sort_values("sirius_structure_confidence", ascending=False)
    return hits


def load_feature_quant() -> pd.DataFrame:
    # aligned_features_ms2.csv.zst is a zstd-compressed CSV; decompress via
    # `zstd -dc` per HPCC storage convention (internal pipeline intermediate).
    proc = subprocess.run(
        ["zstd", "-dc", str(FEATURES_ZST)], capture_output=True, check=True
    )
    from io import BytesIO

    df = pd.read_csv(BytesIO(proc.stdout), low_memory=False)
    if "row ID" not in df.columns:
        raise ValueError(f"unexpected aligned_features_ms2.csv schema: {df.columns.tolist()[:10]}")
    return df


def load_metadata() -> pd.DataFrame:
    meta = pd.read_csv(METADATA_GZ, sep="\t")
    meta = meta[meta["ATTRIBUTE_TYPE"] == "sample"].copy()
    meta["peak_col"] = meta["filename"] + " Peak area"
    return meta


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    hits = load_sirius_annotations()
    if hits.empty:
        print("No siderophore-related SIRIUS annotations found.", file=sys.stderr)
        sys.exit(1)
    hits.to_csv(OUT_DIR / "siderophore_feature_annotations.tsv", sep="\t", index=False)
    print(f"[annotations] {len(hits)} candidate feature(s) matched siderophore keyword search")
    print(hits[["row ID", "sirius_structure_name", "sirius_structure_confidence"]].to_string(index=False))

    quant = load_feature_quant()
    meta = load_metadata()

    peak_cols_present = [c for c in meta["peak_col"] if c in quant.columns]
    missing = set(meta["peak_col"]) - set(peak_cols_present)
    if missing:
        print(f"[warn] {len(missing)} metadata samples have no matching quant column (dropped): "
              f"{sorted(missing)[:5]}{'...' if len(missing) > 5 else ''}", file=sys.stderr)
    meta = meta[meta["peak_col"].isin(peak_cols_present)].copy()
    if len(meta) == 0:
        raise ValueError("no metadata samples matched the quant matrix -- check filename join key")

    feat_ids = hits["row ID"].tolist()
    quant_hits = quant[quant["row ID"].isin(feat_ids)].set_index("row ID")
    missing_ids = set(feat_ids) - set(quant_hits.index)
    if missing_ids:
        print(f"[warn] {len(missing_ids)} annotated row ID(s) absent from quant matrix "
              f"(feature filtered out upstream, e.g. blank/QC removal): {sorted(missing_ids)}", file=sys.stderr)

    long_rows = []
    for row_id, row in quant_hits.iterrows():
        for _, m in meta.iterrows():
            val = row.get(m["peak_col"], np.nan)
            long_rows.append(
                {
                    "row_ID": row_id,
                    "filename": m["filename"],
                    "species": m["ATTRIBUTE_SPECIES"],
                    "strain": m["canonical_strain"],
                    "source": m["ATTRIBUTE_SOURCE"],
                    "peak_area": val,
                }
            )
    long_df = pd.DataFrame(long_rows)
    # ANALYSIS_OK[missingness]: MZmine leaves peak_area blank/NaN when no peak
    # was integrated in that sample -- absence of a peak, not a missing
    # measurement. Treating NaN as 0 for the detected/not-detected call is the
    # correct semantics here; detection_rate below is computed from this flag.
    long_df["detected"] = long_df["peak_area"].fillna(0) > 0

    name_map = hits.set_index("row ID")["sirius_structure_name"].to_dict()
    conf_map = hits.set_index("row ID")["sirius_structure_confidence"].to_dict()
    long_df["compound_name"] = long_df["row_ID"].map(name_map)
    long_df["sirius_confidence"] = long_df["row_ID"].map(conf_map)

    # Species-level summary, split by sample source (cell pellet vs supernatant)
    grp = long_df.groupby(["row_ID", "compound_name", "sirius_confidence", "species", "source"])
    species_summary = grp["peak_area"].agg(
        n_samples="size",
        # ANALYSIS_OK[missingness]: same NaN-as-absent-peak semantics as `detected` above.
        n_detected=lambda s: (s.fillna(0) > 0).sum(),
        mean_peak_area="mean",
        median_peak_area="median",
        max_peak_area="max",
    ).reset_index()
    species_summary["detection_rate"] = species_summary["n_detected"] / species_summary["n_samples"]
    species_summary = species_summary.sort_values(
        ["sirius_confidence", "row_ID", "mean_peak_area"], ascending=[False, True, False]
    )
    species_summary.to_csv(OUT_DIR / "siderophore_abundance_by_species.tsv", sep="\t", index=False)

    # Strain-level summary (finer grain; supernatant + cell pellet pooled per strain,
    # since not every strain has both fractions sampled)
    grp_s = long_df.groupby(["row_ID", "compound_name", "sirius_confidence", "species", "strain"])
    strain_summary = grp_s["peak_area"].agg(
        n_samples="size",
        # ANALYSIS_OK[missingness]: same NaN-as-absent-peak semantics as `detected` above.
        n_detected=lambda s: (s.fillna(0) > 0).sum(),
        mean_peak_area="mean",
        max_peak_area="max",
    ).reset_index()
    strain_summary["detection_rate"] = strain_summary["n_detected"] / strain_summary["n_samples"]
    strain_summary = strain_summary.sort_values(
        ["sirius_confidence", "row_ID", "mean_peak_area"], ascending=[False, True, False]
    )
    strain_summary.to_csv(OUT_DIR / "siderophore_abundance_by_strain.tsv", sep="\t", index=False)

    print(f"\n[done] wrote {OUT_DIR/'siderophore_feature_annotations.tsv'}")
    print(f"[done] wrote {OUT_DIR/'siderophore_abundance_by_species.tsv'} ({len(species_summary)} rows)")
    print(f"[done] wrote {OUT_DIR/'siderophore_abundance_by_strain.tsv'} ({len(strain_summary)} rows)")


if __name__ == "__main__":
    main()
