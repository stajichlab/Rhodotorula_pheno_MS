#!/usr/bin/env python3
"""
Read a completed (merged) SIRIUS run and distill it into one row per
target feature -- top-ranked formula, top-ranked structure (if any), and
the CANOPUS compound-class call -- keyed by row ID, for
build_compound_summary.py to merge into every comparison's
compound_summary.tsv (and the cross-comparison rollup).

Join key: SIRIUS's own "mappingFeatureId" column, present in every
SIRIUS summary table, equals this project's "row ID" -- confirmed
directly against a completed run in the sibling project (
../Rhodotorula_MS2_pheno_explore/analysis/secreted_products/sirius_annotation/):
mappingFeatureId 3760 in formula_identifications.tsv matches FEATURE_ID
3760 in that run's own target MGF, which in turn is that run's row_id
3760. export_sirius_targets_mgf.py in this repo sets FEATURE_ID=<row ID>
for exactly this reason.

Per feature:
  - sirius_formula / sirius_adduct: formulaRank == 1 row from
    formula_identifications.tsv (SIRIUS's best molecular formula guess;
    always present if SIRIUS produced any result for that feature).
  - sirius_structure_name / _smiles / _confidence: structurePerIdRank
    == 1 row from structure_identifications.tsv, if CSI:FingerID found
    any structure candidate at all (not guaranteed -- many features will
    have a formula but no structure hit).
  - sirius_npc_pathway / _class / _classyfire_class: the corresponding
    row (matched by mappingFeatureId) from canopus_structure_summary.tsv
    if present, else canopus_formula_summary.tsv -- CANOPUS's predicted
    compound class, useful even when no exact structure was found.

Usage:
    python3 scripts/import_sirius_annotations.py \
        [--merged-dir analysis/sirius_annotation/sirius_project_container/merged]
    (then re-run build_compound_summary.py to push the annotations into
    every comparison's compound_summary.tsv)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DEFAULT_MERGED_DIR = REPO / "analysis" / "sirius_annotation" / "sirius_project_container" / "merged"
OUT_PATH = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"


def load(merged_dir: Path, fname: str) -> pd.DataFrame | None:
    path = merged_dir / fname
    if not path.exists():
        print(f"skip: {path} not found", file=sys.stderr)
        return None
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"mappingFeatureId": "row ID"})
    df["row ID"] = df["row ID"].astype(int)
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-dir", type=Path, default=DEFAULT_MERGED_DIR)
    args = ap.parse_args()

    formula = load(args.merged_dir, "formula_identifications.tsv")
    structure = load(args.merged_dir, "structure_identifications.tsv")
    canopus_structure = load(args.merged_dir, "canopus_structure_summary.tsv")
    canopus_formula = load(args.merged_dir, "canopus_formula_summary.tsv")

    if formula is None:
        sys.exit(f"no formula_identifications.tsv in {args.merged_dir} -- has SIRIUS finished/been merged?")

    best_formula = (
        formula[formula["formulaRank"] == 1][["row ID", "molecularFormula", "adduct"]]
        .rename(columns={"molecularFormula": "sirius_formula", "adduct": "sirius_adduct"})
    )

    parts = [best_formula]

    if structure is not None:
        best_structure = (
            structure[structure["structurePerIdRank"] == 1]
            [["row ID", "name", "smiles", "ConfidenceScoreExact"]]
            .rename(columns={
                "name": "sirius_structure_name",
                "smiles": "sirius_structure_smiles",
                "ConfidenceScoreExact": "sirius_structure_confidence",
            })
        )
        parts.append(best_structure)

    canopus = canopus_structure if canopus_structure is not None else canopus_formula
    if canopus is not None:
        canopus_cols = {
            "row ID": "row ID",
            "NPC#pathway": "sirius_npc_pathway",
            "NPC#class": "sirius_npc_class",
            "ClassyFire#class": "sirius_classyfire_class",
        }
        available = [c for c in canopus_cols if c in canopus.columns]
        best_canopus = canopus[available].drop_duplicates("row ID").rename(columns=canopus_cols)
        parts.append(best_canopus)

    out = parts[0]
    for p in parts[1:]:
        out = out.merge(p, on="row ID", how="outer")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_PATH, sep="\t", index=False)

    n_formula = out["sirius_formula"].notna().sum()
    n_structure = out.get("sirius_structure_name", pd.Series(dtype=object)).notna().sum()
    print(
        f"{len(out)} features: {n_formula} with a SIRIUS formula, {n_structure} with a structure hit -> {OUT_PATH}",
        file=sys.stderr,
    )
    print("re-run: python3 scripts/build_compound_summary.py   to push these into every compound_summary.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
