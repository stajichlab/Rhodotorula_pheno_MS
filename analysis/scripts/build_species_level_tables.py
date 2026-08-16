#!/usr/bin/env python3
"""
Collapse a strain-level table (one row per strain) to a species-level table
(one row per species), per the "Species-Level Collapse: Procedure" section
of analysis/INTEGRATED_ANALYSIS_STRATEGY.md.

Why: PGLS/PIC/Blomberg's K/Pagel's lambda all model trait evolution along a
tree where each tip is one independent evolutionary lineage. Feeding these
methods multiple strains of the same species as separate tips (or as
repeated rows against one tip) treats non-independent, closely related
strains as independent draws -- pseudoreplication at a finer grain than the
already-flagged R. dairenensis clade-vs-rest confound. Collapsing to one row
per species (matched to one tip per species in the pruned tree produced by
prune_species_tree.R) aligns the analysis unit with what the statistical
model assumes.

Column handling:
  - "binary" columns (declared via --binary-cols): collapsed to species-level
    *prevalence* (fraction of strains in the species with value 1/True), not
    a hard binary collapse. A feature present in 2/9 strains is a different
    claim than 9/9.
  - "count" columns you want proteome-size-normalized (--normalize-cols):
    divided by --normalize-by (e.g. total_protein_count) *per strain first*,
    then averaged. Normalizing after averaging would let strains with
    different proteome sizes contribute unequal weight to the raw count.
  - everything else numeric: species mean/median + SD, per --agg.
  - non-numeric columns (e.g. a free-text species description) are dropped
    from the aggregate other than the grouping key itself.

Every output row also carries n_strains and, for continuous columns, the
within-species SD -- species with n_strains == 1 have no SD estimate; the
diagnostics file flags these explicitly so downstream readers know which
tips have zero within-species variance information, not "genuinely zero
variance."

Usage:
  python build_species_level_tables.py \
      --input analysis/integrated_analysis/phase1_phenotype/strain_phenotype_table.csv \
      --species-col species \
      --id-col strain_id \
      --out analysis/integrated_analysis/phase3_genome_phenotype/species_level/species_phenotype_table.csv

  python build_species_level_tables.py \
      --input analysis/integrated_analysis/phase3_genome_phenotype/genome_feature_matrices/pfam_strain_matrix.csv \
      --species-col species --id-col strain_id \
      --normalize-by total_protein_count --normalize-cols PF00001,PF00002,... \
      --out analysis/integrated_analysis/phase3_genome_phenotype/species_level/species_pfam_matrix.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent


def collapse_to_species(
    df: pd.DataFrame,
    species_col: str,
    id_col: str,
    binary_cols: list[str],
    normalize_cols: list[str],
    normalize_by: str | None,
    agg: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (species_table, diagnostics_table)."""
    df = df.copy()

    if normalize_cols:
        if not normalize_by:
            raise ValueError("--normalize-cols given without --normalize-by")
        if normalize_by not in df.columns:
            raise ValueError(f"--normalize-by column {normalize_by!r} not found in input")
        bad = df[normalize_by].isna() | (df[normalize_by] == 0)
        if bad.any():
            raise ValueError(
                f"{bad.sum()} strain(s) have missing/zero {normalize_by}; "
                "cannot normalize counts by it. Fix upstream or exclude these strains first."
            )
        for col in normalize_cols:
            df[col] = df[col] / df[normalize_by]

    numeric_cols = [
        c for c in df.columns if c not in (species_col, id_col) and pd.api.types.is_numeric_dtype(df[c])
    ]
    binary_set = set(binary_cols)
    unknown_binary = binary_set - set(numeric_cols)
    if unknown_binary:
        raise ValueError(f"--binary-cols not found/not numeric in input: {sorted(unknown_binary)}")
    continuous_cols = [c for c in numeric_cols if c not in binary_set]

    for c in binary_cols:
        bad_vals = set(df[c].dropna().unique()) - {0, 1, 0.0, 1.0}
        if bad_vals:
            raise ValueError(f"--binary-cols column {c!r} has non-0/1 values: {sorted(bad_vals)}")

    grouped = df.groupby(species_col, sort=True)

    out = pd.DataFrame(index=grouped.size().index)
    out["n_strains"] = grouped.size()

    center = grouped[continuous_cols].agg(agg) if continuous_cols else pd.DataFrame(index=out.index)
    spread = grouped[continuous_cols].std() if continuous_cols else pd.DataFrame(index=out.index)
    for c in continuous_cols:
        out[f"{c}_{agg}"] = center[c]
        out[f"{c}_sd"] = spread[c]

    prevalence = grouped[binary_cols].mean() if binary_cols else pd.DataFrame(index=out.index)
    for c in binary_cols:
        out[f"{c}_prevalence"] = prevalence[c]

    out = out.reset_index().rename(columns={species_col: "species"})

    diag = out[["species", "n_strains"]].copy()
    diag["single_strain_species"] = diag["n_strains"] == 1
    n_single = int(diag["single_strain_species"].sum())
    diag.attrs["summary"] = (
        f"{len(diag)} species; {n_single} represented by a single strain "
        "(no within-species SD estimate for those rows)."
    )
    return out, diag


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=Path, help="Strain-level input CSV (one row per strain).")
    ap.add_argument("--species-col", default="species")
    ap.add_argument("--id-col", default="strain_id")
    ap.add_argument(
        "--binary-cols",
        default="",
        help="Comma-separated 0/1 columns to collapse to species-level prevalence instead of mean/SD.",
    )
    ap.add_argument(
        "--normalize-cols",
        default="",
        help="Comma-separated count columns to divide by --normalize-by (per strain) before aggregating.",
    )
    ap.add_argument(
        "--normalize-by",
        default=None,
        help="Column (e.g. total_protein_count) used to normalize --normalize-cols per strain.",
    )
    ap.add_argument("--agg", default="mean", choices=["mean", "median"], help="Central-tendency stat for continuous columns.")
    ap.add_argument("--out", required=True, type=Path, help="Output species-level CSV.")
    args = ap.parse_args()

    binary_cols = [c.strip() for c in args.binary_cols.split(",") if c.strip()]
    normalize_cols = [c.strip() for c in args.normalize_cols.split(",") if c.strip()]

    df = pd.read_csv(args.input)
    if args.species_col not in df.columns:
        raise SystemExit(f"--species-col {args.species_col!r} not found in {args.input}")
    if args.id_col not in df.columns:
        raise SystemExit(f"--id-col {args.id_col!r} not found in {args.input}")

    species_table, diag = collapse_to_species(
        df,
        species_col=args.species_col,
        id_col=args.id_col,
        binary_cols=binary_cols,
        normalize_cols=normalize_cols,
        normalize_by=args.normalize_by,
        agg=args.agg,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    species_table.to_csv(args.out, index=False)

    diag_path = args.out.with_name(args.out.stem + "_diagnostics.csv")
    diag.to_csv(diag_path, index=False)

    print(f"Collapsed {len(df)} strains -> {len(species_table)} species")
    print(diag.attrs["summary"])
    print(f"Wrote {args.out}")
    print(f"Wrote {diag_path}")


if __name__ == "__main__":
    main()
