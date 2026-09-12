#!/usr/bin/env python3
"""
AHL autoinducer / colony morphology hypothesis
(analysis/ahl_autoinducer_search/AHL_AUTOINDUCER_SEARCH.md): build the
strain-level colony-texture (morphology proxy) table from the ingested
`salinity_texture_baseline` dataset (data/raw/salinity_texture_baseline/,
see data/DATA_MANIFEST.md and data/metadata/salinity_texture_baseline/ for
full provenance/caveats).

This is a PROXY for colony morphology (smooth/rough), not a validated
categorical call -- see the ingestion docs. Nothing here should be read
as "morphology" without that caveat attached.

Aggregation: median across replicate colony rows per strain, matching
this project's existing convention for collapsing multi-row phenotype
data to one row per strain (build_strain_phenotype_table.py uses the same
median-collapse pattern for color). Median (not mean) chosen because
Shape_Area shows a clear bimodal split between real colonies and
segmentation-artifact specks (see --min-colony-area below) -- median is
more robust than mean if an artifact slips through the area filter for a
strain with few replicates.

--min-colony-area (default 1000 px): colony-level QC filter applied
BEFORE aggregation. Shape_Area in the raw data is clearly bimodal: 11/997
rows fall below 1000 px (down to 10 px), then the next percentile jumps
straight to ~23,000 px -- the low group is almost certainly segmentation
specks (debris, plate edge artifacts), not small colonies. Rows below the
threshold are dropped entirely (not median-collapsed in), and the count
dropped per strain is logged in the diagnostics file so it's auditable.

Texture summary columns: the 13 `-avg-scale05` (rotation-invariant)
Haralick GLCM metrics only -- NOT the 52 directional (deg000/045/090/135)
columns, which are retained in the raw table for anisotropy diagnostics
but are not this script's summary target. See
data/metadata/salinity_texture_baseline/schema.yaml for metric semantics.

Output: one row per strain with n_colonies_raw, n_colonies_used (post
area-filter), n_colonies_dropped_low_area, and median value for each of
the 13 avg-scale05 metrics plus Shape_Area/Shape_Solidity (kept for QC
and as a known confound axis -- colony size already confounds abundance
comparisons elsewhere in this project, see
.living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md).

Usage:
    python3 analysis/scripts/build_strain_texture_table.py
    python3 analysis/scripts/build_strain_texture_table.py --min-colony-area 500
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
RAW_TEXTURE = REPO / "data" / "raw" / "salinity_texture_baseline" / "salinity0_hours90_texture.csv.gz"
CROSSWALK = (
    REPO / "data" / "metadata" / "EXFAB_UCR-005"
    / "MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz"
)
OUT_DIR = REPO / "analysis" / "ahl_autoinducer_search"
OUT_CSV = OUT_DIR / "strain_texture_table.csv"
DIAG_TXT = OUT_DIR / "strain_texture_table_diagnostics.txt"

TEXTURE_METRICS = [
    "AngularSecondMoment", "Contrast", "Correlation", "HaralickVariance",
    "InverseDifferenceMoment", "SumAverage", "SumVariance", "SumEntropy",
    "Entropy", "DiffVariance", "DiffEntropy", "InfoCorrelation1", "InfoCorrelation2",
]
AVG_COLS = [f"Texture_{m}-avg-scale05" for m in TEXTURE_METRICS]
QC_COLS = ["Shape_Area", "Shape_Solidity"]


def log(msg: str, fh) -> None:
    print(msg)
    fh.write(msg + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-colony-area", type=float, default=1000.0,
                    help="Drop colony rows with Shape_Area below this (px) before aggregation.")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not RAW_TEXTURE.exists():
        raise SystemExit(f"Texture table not found: {RAW_TEXTURE} -- run the ingestion step first.")
    if not CROSSWALK.exists():
        raise SystemExit(f"Strain-ID crosswalk not found: {CROSSWALK}")

    df = pd.read_csv(RAW_TEXTURE)
    n_raw = len(df)
    assert n_raw == 997, f"Expected 997 raw colony rows from the ingested subset, got {n_raw} -- source file changed?"

    missing_cols = [c for c in AVG_COLS + QC_COLS if c not in df.columns]
    if missing_cols:
        raise SystemExit(f"Expected columns missing from raw texture table: {missing_cols}")

    with open(DIAG_TXT, "w") as fh:
        log(f"Raw colony rows: {n_raw}", fh)
        log(f"Raw distinct Strain ID: {df['Strain ID'].nunique()}", fh)

        # crosswalk Strain ID (numeric) -> strain_code
        cross = pd.read_csv(CROSSWALK, sep="\t", usecols=["Strain ID", "Strain"], low_memory=False).dropna()
        cross["Strain ID"] = cross["Strain ID"].astype(float)
        cross = cross.drop_duplicates(subset="Strain ID")
        log(f"Crosswalk rows (Strain ID -> strain_code): {len(cross)}", fh)

        df = df.merge(cross, on="Strain ID", how="left", suffixes=("", "_crosswalk"))
        n_no_code = df["Strain"].isna().sum()
        log(f"Colony rows with no resolvable strain_code after crosswalk join: {n_no_code} (dropped)", fh)
        df = df.dropna(subset=["Strain"]).rename(columns={"Strain": "strain_code"})

        # texture-value completeness check (6 rows expected missing per ingestion docs)
        texture_missing = df[AVG_COLS].isna().any(axis=1)
        log(f"Colony rows missing one or more avg-scale05 texture values: {int(texture_missing.sum())} (dropped)", fh)
        df = df[~texture_missing].copy()

        # area-based QC filter (segmentation-artifact colonies)
        n_before_area_filter = len(df)
        low_area = df["Shape_Area"] < args.min_colony_area
        log(f"Colony rows below --min-colony-area={args.min_colony_area} px: {int(low_area.sum())} (dropped)", fh)
        df = df[~low_area].copy()
        log(f"Colony rows entering aggregation: {len(df)} (of {n_before_area_filter} pre-area-filter, {n_raw} raw)", fh)

        n_colonies_raw = df.groupby("strain_code").size().rename("n_colonies_used")
        agg_cols = AVG_COLS + QC_COLS
        strain_med = df.groupby("strain_code")[agg_cols].median()
        strain_species = df.groupby("strain_code")["Species"].first()

        out = strain_med.join(n_colonies_raw).join(strain_species).reset_index()
        out = out.rename(columns={c: c.replace("Texture_", "").replace("-avg-scale05", "_median") for c in AVG_COLS})
        out = out.rename(columns={"Shape_Area": "shape_area_median", "Shape_Solidity": "shape_solidity_median"})

        assert out["strain_code"].is_unique, "Duplicate strain_code rows in aggregated output -- groupby failed"
        log(f"\nFinal strain-level table: {len(out)} strains", fh)
        log(f"n_colonies_used per strain: min={out['n_colonies_used'].min()}, "
            f"median={out['n_colonies_used'].median()}, max={out['n_colonies_used'].max()}", fh)

        out.to_csv(OUT_CSV, index=False)
        log(f"\nWrote {OUT_CSV}", fh)
        log(f"Wrote {DIAG_TXT}", fh)

    print(out.describe().T[["count", "mean", "std", "min", "max"]].to_string())


if __name__ == "__main__":
    main()
