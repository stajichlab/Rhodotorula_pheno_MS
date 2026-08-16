#!/usr/bin/env python3
"""
Phase 1 (analysis/INTEGRATED_ANALYSIS_STRATEGY.md): build the strain-level
phenotype table that is the response variable for every downstream phase --
one row per strain with L*, a*, b*, C*, h(deg), and a composite orange_score.

Color columns: uses Median_ColorLab_{L*,a*,b*}Mean, the same columns
scripts/pcoa_color_phenotype.py uses for the existing color-phenotype PCoA,
so this table's L*/a*/b* stay numerically comparable to that ordination.
(YPD/color_shape_growth/scripts/00_build_master_table.py used the
Mean_ColorLab_* family instead for a *different*, AA-composition-focused
analysis; that choice isn't disturbed here, this script just doesn't reuse
it, to stay consistent with the canonical color ordination instead.)

A handful of strains (2 as of 2026-08) have more than one YPD2 row
(different plates/replicates) -- collapsed to the mean across rows before
computing derived quantities, matching the convention already used in
YPD/color_shape_growth/scripts/00_build_master_table.py.

Derived quantities:
  - C* = sqrt(a*^2 + b*^2)                          (chroma / color intensity)
  - h_deg = atan2(b*, a*) in degrees, wrapped to [0, 360)   (hue angle)
  - orange_score = z(-L*) + z(a*) + z(b*)            (composite "darker
    orange" score: darker + redder + yellower all push it up equally).
    z-scores are computed once across all IN-SCOPE strains with complete
    L*/a*/b* data (not per-species) -- see --genus-filter below for what
    "in scope" means. The three component z-scores (L_z, a_z, b_z) are
    kept as columns so the composite's weighting is auditable and can be
    redone with different weights without recomputing from raw color values.

--genus-filter (default "Rhodotorula"): the strategy doc's premise is that
R. dairenensis looks darker orange *relative to other Rhodotorula species*
(that's what the existing CIELAB PCoA compared). The raw YPD2 table also
includes a few non-Rhodotorula outgroup strains (Cystobasidium sp.,
Pseudomicrostroma phylloplanum) plus R. araucariae, a distant Rhodotorula
species with a very different baseline color. Computing the z-scores across
ALL 316 phenotyped strains lets those 3 outlier species drag the population
mean/SD and distort the "darker orange within Rhodotorula" ranking -- an
initial full-dataset run of this pipeline had R. dairenensis ranked only
4th of 17 species, which needed exactly this check before being trusted.
Pass --genus-filter "" to disable and use every strain regardless of genus.
Strains outside the filter are still kept in the output table (raw L*/a*/b*/
C*/h_deg computed for everyone) but get NaN for L_z/a_z/b_z/orange_score,
since they were excluded from the population those z-scores are relative to.

Strains with a missing Species (16 as of 2026-08) are kept in the output
(needed for e.g. strain-level within-species tests that don't need species)
but flagged via species_missing=True; they drop out of any species-level
join automatically since they have no species key to join on, and are
excluded from the genus filter (unknown genus).

--source (default "control_90_110"): which phenotype table to read.
  - "control_90_110": the CANONICAL source (chosen 2026-08-15, see
    analysis/examine_phenotype_calling/RECOMMENDATION.md). Control-media
    (Cu=0), latest-timepoint (105/108h) colony color+size table, formally
    ingested into data/raw/control_phenotype_90_110h/ +
    data/metadata/control_phenotype_90_110h/ (schema.yaml, provenance.md,
    summary_stats.md) -- see data/DATA_MANIFEST.md and .living/decisions.md
    for full provenance. 314 rows / 303 distinct strain_code (303 usable
    here; 10 rows have no strain_code and are dropped by the groupby
    below -- see the ingested doc's Known Issues), 2 missing species, no
    missing L/a/b, includes R. evergladensis (absent from ypd2).
  - "ypd2": the original data/metadata/EXFAB_UCR-005/YPD2_phenotypic...csv.gz
    (318 rows, 15-16 strains missing Species as of 2026-08). Legacy;
    kept for cross-validation, not formally re-ingested under the
    data/raw + data/metadata split (predates this repo's mycelium
    adoption -- see the EXFAB_UCR-005 entry in data/DATA_MANIFEST.md).
  - "control_70_80" / "control_80_90": the two earlier-timepoint windows
    control_90_110 was chosen over (see
    analysis/examine_phenotype_calling/RECOMMENDATION.md for the
    comparison). NOT ingested -- still read cross-project from the
    sibling Rhodotorula_phenotypes project, since only the canonical
    choice was formally ingested. Kept here only so the timepoint
    comparison script can still reproduce its result; do not build new
    analysis on these two.

Usage:
    python3 analysis/scripts/build_strain_phenotype_table.py
    python3 analysis/scripts/build_strain_phenotype_table.py --genus-filter ""   # full dataset, no genus scoping
    python3 analysis/scripts/build_strain_phenotype_table.py --source control_90_110
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase1_phenotype"
OUT_CSV = OUT_DIR / "strain_phenotype_table.csv"
DIAG_TXT = OUT_DIR / "strain_phenotype_table_diagnostics.txt"

SOURCES = {
    "ypd2": dict(
        path=REPO / "data" / "metadata" / "EXFAB_UCR-005" / "YPD2_phenotypic.20260702.fixed.csv.gz",
        strain_col="Strain",
        species_col="Species",
        lab_cols={
            "L": "Median_ColorLab_L*Mean", "a": "Median_ColorLab_a*Mean", "b": "Median_ColorLab_b*Mean",
            "area": "Median_Shape_Area",
        },
    ),
    "control_90_110": dict(
        # Ingested 2026-08-15 -- data/raw/control_phenotype_90_110h/, see
        # data/DATA_MANIFEST.md and .living/decisions.md for provenance.
        # This is the canonical source (see module docstring).
        path=REPO / "data" / "raw" / "control_phenotype_90_110h" / "phenotype_control_timepoint_90_110.csv",
        strain_col="strain_code",
        species_col="species",
        lab_cols={"L": "l_median", "a": "a_median", "b": "b_median", "area": "area_median"},
    ),
    **{
        f"control_{window}": dict(
            # NOT ingested into data/raw/ -- read cross-project. Only used
            # for the timepoint-window comparison in
            # analysis/examine_phenotype_calling/; control_90_110 above is
            # the ingested, canonical choice.
            path=Path(
                "/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/analysis/"
                f"control_late_timepoint_phenotype/results/phenotype_control_timepoint_{window}.csv"
            ),
            strain_col="strain_code",
            species_col="species",
            lab_cols={"L": "l_median", "a": "a_median", "b": "b_median"},
        )
        for window in ("70_80", "80_90")
    },
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--genus-filter", default="Rhodotorula",
        help='Only strains whose species starts with this string contribute to the orange_score '
             'z-score population (pass "" to disable and use all strains). Default: Rhodotorula.',
    )
    ap.add_argument(
        "--source", default="control_90_110", choices=sorted(SOURCES),
        help="Which phenotype table to read (see module docstring). Default: control_90_110 -- "
             "see analysis/examine_phenotype_calling/ for why this is the chosen canonical source "
             "(best agreement with YPD2 across every trait, most plateaued growth stage, most "
             "complete species coverage).",
    )
    args = ap.parse_args()
    genus_filter = args.genus_filter or None
    src = SOURCES[args.source]
    strain_col, species_col, lab_cols = src["strain_col"], src["species_col"], src["lab_cols"]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(src["path"])
    n_raw = len(raw)
    n_raw_strains = raw[strain_col].nunique()

    missing_lab = raw[list(lab_cols.values())].isna().any(axis=1)
    dropped_no_lab = int(missing_lab.sum())
    raw = raw.loc[~missing_lab].copy()

    # pandas.groupby drops NaN group keys silently -- a strain_col-NULL row
    # ("unidentified spot" in control_* sources) would otherwise vanish from
    # the output with no warning. Surface it explicitly every run instead.
    missing_strain_id = raw[strain_col].isna()
    dropped_no_strain_id = int(missing_strain_id.sum())
    if dropped_no_strain_id:
        print(
            f"WARNING: {dropped_no_strain_id} row(s) have no {strain_col!r} value and will be "
            f"dropped by groupby (real color/area data, currently unrecoverable here): "
            f"strain_id={raw.loc[missing_strain_id, 'strain_id'].tolist() if 'strain_id' in raw.columns else '(no strain_id column)'}",
            file=sys.stderr,
        )
    raw = raw.loc[~missing_strain_id].copy()

    # collapse duplicate-strain rows (different plates/replicates) to their mean
    agg = {v: "mean" for v in lab_cols.values()}
    agg[species_col] = "first"
    grouped = raw.groupby(strain_col, as_index=False).agg(agg)
    n_replicates = raw.groupby(strain_col).size().rename("n_replicates")
    grouped = grouped.merge(n_replicates, left_on=strain_col, right_index=True)

    grouped = grouped.rename(
        columns={strain_col: "strain_id", species_col: "species", **{v: k for k, v in lab_cols.items()}}
    )
    grouped["species_missing"] = grouped["species"].isna()

    if "area" not in grouped.columns:
        grouped["area"] = np.nan  # sources without an area column (none currently) still get the output column

    L, a, b = grouped["L"].to_numpy(), grouped["a"].to_numpy(), grouped["b"].to_numpy()
    grouped["C"] = np.sqrt(a**2 + b**2)
    grouped["h_deg"] = np.degrees(np.arctan2(b, a)) % 360.0

    in_scope = (
        grouped["species"].notna() & grouped["species"].str.startswith(genus_filter)
        if genus_filter
        else pd.Series(True, index=grouped.index)
    )
    n_in_scope = int(in_scope.sum())

    def zscore_in_scope(x: np.ndarray, mask: np.ndarray) -> pd.Series:
        pop = x[mask]
        z = np.full(x.shape, np.nan)
        z[mask] = (pop - pop.mean()) / pop.std(ddof=1)
        return pd.Series(z, index=grouped.index)

    in_scope_arr = in_scope.to_numpy()
    grouped["L_z"] = zscore_in_scope(-L, in_scope_arr)
    grouped["a_z"] = zscore_in_scope(a, in_scope_arr)
    grouped["b_z"] = zscore_in_scope(b, in_scope_arr)
    grouped["orange_score"] = grouped["L_z"] + grouped["a_z"] + grouped["b_z"]

    grouped = grouped.rename(columns={"L": "L*", "a": "a*", "b": "b*", "C": "C*"})
    out_cols = [
        "strain_id", "species", "species_missing", "n_replicates",
        "L*", "a*", "b*", "C*", "h_deg", "area", "L_z", "a_z", "b_z", "orange_score",
    ]
    grouped = grouped[out_cols].sort_values(["species", "strain_id"], na_position="last")
    grouped.to_csv(OUT_CSV, index=False)

    n_species_missing = int(grouped["species_missing"].sum())
    with DIAG_TXT.open("w") as fh:
        fh.write(f"Source: {args.source} ({src['path']})\n")
        fh.write(f"Input rows: {n_raw + dropped_no_lab} ({n_raw_strains} distinct {strain_col} values)\n")
        fh.write(f"Dropped for missing L*/a*/b*: {dropped_no_lab}\n")
        fh.write(f"Dropped for missing {strain_col!r} (would be silently dropped by groupby otherwise): {dropped_no_strain_id}\n")
        fh.write(f"Collapsed to {len(grouped)} strain rows (mean across replicate rows where >1)\n")
        fh.write(
            f"orange_score genus filter: {genus_filter!r} -- {n_in_scope}/{len(grouped)} strains in scope "
            f"(the rest get NaN L_z/a_z/b_z/orange_score)\n"
        )
        fh.write(f"Strains with missing species: {n_species_missing}\n")
        if n_species_missing:
            fh.write("  " + ", ".join(grouped.loc[grouped["species_missing"], "strain_id"]) + "\n")
        fh.write("\nStrain count by species:\n")
        counts = grouped.loc[~grouped["species_missing"], "species"].value_counts()
        for sp, n in counts.items():
            fh.write(f"  {sp}: {n}\n")

    print(f"Wrote {len(grouped)} strains -> {OUT_CSV}")
    print(f"Diagnostics: {DIAG_TXT}")
    if n_species_missing:
        print(
            f"NOTE: {n_species_missing} strain(s) have no Species assigned in --source {args.source} -- "
            "check ~/projects/Rhodotorula_phenotypes/ for a more complete phenotype table if this repo "
            "hasn't already ingested one (see --source control_90_110), then rerun.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
