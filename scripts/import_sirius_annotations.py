#!/usr/bin/env python3
"""
Distill one or more completed (merged) SIRIUS runs into one row per
feature -- top-ranked formula, top-ranked structure (if any), and the
CANOPUS compound-class call -- keyed by row ID, accumulating into
analysis/sirius_annotation/sirius_annotations.tsv across separate runs
over time (this project's own, plus e.g. relevant runs from the sibling
project) rather than each import clobbering the last.

Join key: SIRIUS's own "mappingFeatureId" column, present in every
SIRIUS summary table, equals this project's "row ID" -- confirmed twice
now against completed runs in the sibling project: mappingFeatureId 3760
in a secreted_products run matches FEATURE_ID 3760 in that run's own
target MGF (this repo's row 3760); mappingFeatureId 4015/26859/6891 in
the pathway_targeted_association "broad" run match this repo's own row
4015 (m/z 275.160), 26859 (m/z 620.302), 6891 (m/z 786.161) exactly on
m/z (RT differs slightly, a minor calibration difference between export
runs -- m/z is the reliable check). Both sibling-project runs used the
same underlying MZmine feature alignment as this project, so their row
numbering lines up directly -- no re-derivation needed, just import.

Per feature, per source run:
  - sirius_formula / sirius_adduct: formulaRank == 1 row from
    formula_identifications.tsv (always present if SIRIUS produced any
    result for that feature).
  - sirius_structure_name / _smiles / _confidence: structurePerIdRank
    == 1 row from structure_identifications.tsv, if CSI:FingerID found
    any structure candidate at all (not guaranteed).
  - sirius_npc_pathway / _class / _classyfire_class: the corresponding
    row from canopus_structure_summary.tsv if present, else
    canopus_formula_summary.tsv.

Accumulation across runs: if a row ID appears in more than one imported
run (this run's own, plus whatever's already in sirius_annotations.tsv
from prior imports), the row with a structure hit wins over one without,
and among structure hits the higher sirius_structure_confidence wins;
ties keep whichever was already accumulated. Every contributing run is
still recorded in source_run (semicolon-joined) even for a row ID whose
predictions came from a different run, so coverage history isn't lost
just because a later run's guess replaced an earlier one's.

Usage:
    python3 scripts/import_sirius_annotations.py [MERGED_DIR ...] [--label LABEL ...] [--fresh]
        (default MERGED_DIR: analysis/sirius_annotation/sirius_results/merged,
        i.e. this repo's own run; --label defaults to each dir's absolute
        path if not given, one --label per MERGED_DIR)
    (then re-run build_compound_summary.py to push the annotations into
    every comparison's compound_summary.tsv)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DEFAULT_MERGED_DIR = REPO / "analysis" / "sirius_annotation" / "sirius_results" / "merged"
OUT_PATH = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"

CHEM_COLS = [
    "sirius_formula", "sirius_adduct", "sirius_structure_name", "sirius_structure_smiles",
    "sirius_structure_confidence", "sirius_npc_pathway", "sirius_npc_class", "sirius_classyfire_class",
]


def load(merged_dir: Path, fname: str) -> pd.DataFrame | None:
    path = merged_dir / fname
    if not path.exists():
        print(f"skip: {path} not found", file=sys.stderr)
        return None
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"mappingFeatureId": "row ID"})
    df["row ID"] = df["row ID"].astype(int)
    return df


def distill_one_run(merged_dir: Path, label: str) -> pd.DataFrame:
    formula = load(merged_dir, "formula_identifications.tsv")
    structure = load(merged_dir, "structure_identifications.tsv")
    canopus_structure = load(merged_dir, "canopus_structure_summary.tsv")
    canopus_formula = load(merged_dir, "canopus_formula_summary.tsv")

    if formula is None:
        sys.exit(f"no formula_identifications.tsv in {merged_dir} -- has SIRIUS finished/been merged?")

    best_formula = (
        formula[formula["formulaRank"] == 1][["row ID", "molecularFormula", "adduct"]]
        .drop_duplicates("row ID")
        .rename(columns={"molecularFormula": "sirius_formula", "adduct": "sirius_adduct"})
    )
    parts = [best_formula]

    if structure is not None:
        best_structure = (
            structure[structure["structurePerIdRank"] == 1]
            [["row ID", "name", "smiles", "ConfidenceScoreExact"]]
            .drop_duplicates("row ID")
            .rename(columns={
                "name": "sirius_structure_name",
                "smiles": "sirius_structure_smiles",
                "ConfidenceScoreExact": "sirius_structure_confidence",
            })
        )
        # SIRIUS reports -Infinity for some rows (a structure candidate
        # exists but confidence was uncomputable) -- treat that the same
        # as "no confidence available" rather than propagating a literal
        # +/-Infinity float downstream (breaks strict JSON consumers,
        # prints oddly, and isn't a meaningful confidence value either way).
        best_structure["sirius_structure_confidence"] = best_structure["sirius_structure_confidence"].replace(
            [float("inf"), float("-inf")], pd.NA
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
    out["source_run"] = label
    return out


def accumulate(existing: pd.DataFrame | None, new_runs: list[pd.DataFrame]) -> pd.DataFrame:
    frames = ([existing] if existing is not None else []) + new_runs
    combined = pd.concat(frames, ignore_index=True)

    # Every contributing source is preserved regardless of which row "wins"
    # below -- semicolon-joined, de-duplicated, order-preserving.
    def join_sources(s: pd.Series) -> str:
        seen = []
        for v in s:
            for part in str(v).split(";"):
                if part and part not in seen:
                    seen.append(part)
        return ";".join(seen)

    all_sources = combined.groupby("row ID")["source_run"].apply(join_sources)

    has_structure = combined["sirius_structure_name"].notna() if "sirius_structure_name" in combined.columns else pd.Series(False, index=combined.index)
    confidence = combined.get("sirius_structure_confidence", pd.Series(dtype=float)).fillna(-1)
    combined = combined.assign(_has_structure=has_structure, _confidence=confidence)
    combined = combined.sort_values(["_has_structure", "_confidence"], ascending=[False, False])
    winners = combined.drop_duplicates("row ID", keep="first").drop(columns=["_has_structure", "_confidence"])

    winners = winners.set_index("row ID")
    winners["source_run"] = all_sources
    return winners.reset_index().sort_values("row ID")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("merged_dirs", nargs="*", type=Path, default=[DEFAULT_MERGED_DIR])
    ap.add_argument("--label", action="append", help="one per MERGED_DIR, in order; defaults to each dir's absolute path")
    ap.add_argument("--fresh", action="store_true", help="ignore any existing sirius_annotations.tsv instead of accumulating into it")
    args = ap.parse_args()

    labels = args.label or [str(d.resolve()) for d in args.merged_dirs]
    if len(labels) != len(args.merged_dirs):
        sys.exit(f"got {len(args.merged_dirs)} MERGED_DIR but {len(labels)} --label -- need one per dir")

    existing = None
    if not args.fresh and OUT_PATH.exists():
        existing = pd.read_csv(OUT_PATH, sep="\t")
        print(f"accumulating onto {len(existing)} previously-imported features from {OUT_PATH}", file=sys.stderr)

    new_runs = []
    for merged_dir, label in zip(args.merged_dirs, labels):
        run_df = distill_one_run(merged_dir, label)
        n_formula = run_df["sirius_formula"].notna().sum()
        n_structure = run_df.get("sirius_structure_name", pd.Series(dtype=object)).notna().sum()
        print(f"[{label}] {len(run_df)} features: {n_formula} formula, {n_structure} structure hits", file=sys.stderr)
        new_runs.append(run_df)

    out = accumulate(existing, new_runs)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_PATH, sep="\t", index=False)

    n_formula = out["sirius_formula"].notna().sum()
    n_structure = out.get("sirius_structure_name", pd.Series(dtype=object)).notna().sum()
    n_runs = out["source_run"].str.split(";").explode().nunique()
    print(
        f"accumulated total: {len(out)} distinct features across {n_runs} source run(s), "
        f"{n_formula} with a formula, {n_structure} with a structure hit -> {OUT_PATH}",
        file=sys.stderr,
    )
    print("re-run: python3 scripts/build_compound_summary.py   to push these into every compound_summary.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
