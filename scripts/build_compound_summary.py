#!/usr/bin/env python3
"""
For each differential-feature comparison, join the plate-blocking-robust
significant features to whatever compound identity is available -- the
Everything-Bagel (EB/GNPS) library search, and (once run) SIRIUS's own
formula/structure/compound-class predictions -- as a "best guess so far".

Called automatically at the end of differential_features_by_species.py
(so every comparison gets compound_summary.tsv for free), and also
runnable standalone here to backfill/refresh existing comparisons -- in
particular, re-run this after a SIRIUS job completes and
scripts/import_sirius_annotations.py has written
analysis/sirius_annotation/sirius_annotations.tsv, to pull those
predictions into every comparison's table without re-running the stats.

Three identity sources, all keyed by MS2 scan number (== this project's
"row ID" everywhere else, confirmed: MGF SCANS=1 has m/z 278.1897 / RT
194.16s, matching row ID 1 in the feature table exactly):
  - feature_library_search: exact EB/GNPS spectral library matches
  - feature_analog_library_search: analog matches (structurally similar,
    not identical -- precursor mass differs from the library hit)
  - sirius_annotations.tsv (optional, see import_sirius_annotations.py):
    SIRIUS's top-ranked formula and, when available, top-ranked
    structure + CANOPUS compound-class prediction, keyed by
    mappingFeatureId == row ID (confirmed directly against a completed
    SIRIUS run in the sibling project: mappingFeatureId 3760 matches
    FEATURE_ID 3760 in that run's own target MGF).
Best-identity precedence (most to least trustworthy): exact library match
> confident SIRIUS structure (ConfidenceScoreExact) > analog library
match > SIRIUS formula only (no confident structure) > unidentified.

Significance filter: q_value_perm_plate < --fdr when that column is
present and non-null for a row (i.e. the plate-block permutation test
actually ran for it -- see differential_features_by_species.py); falls
back to plain Mann-Whitney q_value for rows where it didn't (e.g. one
group too small to block), flagged in the output so a plate-confounded
hit is never silently presented as blocking-confirmed.

Output:
  - <comparison>/compound_summary.tsv, one row per significant feature,
    sorted by (blocking-robust first, then q-value).
  - analysis/differential_features/all_significant_features_summary.tsv:
    every comparison's compound_summary.tsv concatenated with a
    "comparison" column, for a single across-all-pairs view.

Usage:
    python3 scripts/build_compound_summary.py [DIR ...] [--fdr 0.05]
        (default DIRS: every analysis/differential_features/*/ that has a
        differential_features.csv.gz)
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DIFF_ROOT = REPO / "analysis" / "differential_features"
EB_DIR = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9" / "b773ffa18c2b41e5a3484526293a54f9" / "nf_output"
)
LIBRARY_SEARCH = EB_DIR / "feature_library_search" / "merged_feature_library_search_results.tsv"
ANALOG_SEARCH = EB_DIR / "feature_analog_library_search" / "merged_feature_analog_library_search_results.tsv"
SIRIUS_ANNOTATIONS = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"
MASTER_OUT = DIFF_ROOT / "all_significant_features_summary.tsv"

GNPS_KEEP_COLS = ["NAME", "cosine", "matched_peaks", "ADDUCT", "SMILES", "INCHI", "SPECTRUMID", "ORGANISM"]


def load_gnps(path: Path, prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["query_scan"] = df["query_scan"].astype(int)
    df["cosine"] = df["cosine"].astype(float)
    df = df.sort_values("cosine", ascending=False).drop_duplicates("query_scan", keep="first")
    out = df[["query_scan"] + GNPS_KEEP_COLS].rename(columns={c: f"{prefix}_{c}" for c in GNPS_KEEP_COLS})
    return out.rename(columns={"query_scan": "row ID"})


def load_sirius() -> pd.DataFrame | None:
    if not SIRIUS_ANNOTATIONS.exists():
        return None
    df = pd.read_csv(SIRIUS_ANNOTATIONS, sep="\t")
    df["row ID"] = df["row ID"].astype(int)
    return df


def annotate(df: pd.DataFrame, library: pd.DataFrame, analog: pd.DataFrame, sirius: pd.DataFrame | None) -> pd.DataFrame:
    df = df.merge(library, on="row ID", how="left").merge(analog, on="row ID", how="left")
    if sirius is not None:
        df = df.merge(sirius, on="row ID", how="left")

    df["best_identity"] = df["library_NAME"]
    df["best_identity_source"] = df["library_NAME"].notna().map({True: "exact_library_match", False: pd.NA})

    if sirius is not None and "sirius_structure_name" in df.columns:
        confident_structure = df["best_identity"].isna() & df["sirius_structure_name"].notna() & (
            df.get("sirius_structure_confidence", pd.Series(dtype=float)).fillna(0) > 0
        )
        df.loc[confident_structure, "best_identity"] = df.loc[confident_structure, "sirius_structure_name"]
        df.loc[confident_structure, "best_identity_source"] = "sirius_structure"

    still_missing = df["best_identity"].isna() & df["analog_NAME"].notna()
    df.loc[still_missing, "best_identity"] = df.loc[still_missing, "analog_NAME"]
    df.loc[still_missing, "best_identity_source"] = "analog_library_match"

    if sirius is not None and "sirius_formula" in df.columns:
        formula_only = df["best_identity"].isna() & df["sirius_formula"].notna()
        df.loc[formula_only, "best_identity"] = df.loc[formula_only, "sirius_formula"] + " (formula only)"
        df.loc[formula_only, "best_identity_source"] = "sirius_formula_only"

    df.loc[df["best_identity"].isna(), "best_identity_source"] = "unidentified"
    return df


def summarize_one(diff_dir: Path, library: pd.DataFrame, analog: pd.DataFrame, sirius, fdr: float) -> pd.DataFrame:
    diff_path = diff_dir / "differential_features.csv.gz"
    with gzip.open(diff_path, "rt") as fh:
        df = pd.read_csv(fh)

    has_perm = "q_value_perm_plate" in df.columns
    if has_perm:
        df["blocking_robust"] = df["q_value_perm_plate"].notna() & (df["q_value_perm_plate"] < fdr)
        df["q_used"] = df["q_value_perm_plate"].where(df["q_value_perm_plate"].notna(), df["q_value"])
        df["q_used_source"] = df["q_value_perm_plate"].notna().map(
            {True: "q_value_perm_plate", False: "q_value (plate test did not run for this row)"}
        )
    else:
        df["blocking_robust"] = False
        df["q_used"] = df["q_value"]
        df["q_used_source"] = "q_value (no plate-blocking columns in this table)"

    sig = df[df["q_used"] < fdr].copy()
    sig = annotate(sig, library, analog, sirius)
    sig = sig.sort_values(["blocking_robust", "q_used"], ascending=[False, True])

    cols = [
        "row ID", "row m/z", "row retention time", "adduct", "has_ms2",
        "log2FC_a_over_b", "q_value", "q_value_perm_plate", "q_used", "q_used_source", "blocking_robust",
        "best_identity", "best_identity_source",
        "library_NAME", "library_cosine", "library_matched_peaks", "library_SMILES", "library_INCHI", "library_ORGANISM",
        "analog_NAME", "analog_cosine", "analog_matched_peaks", "analog_SMILES", "analog_INCHI", "analog_ORGANISM",
        "sirius_formula", "sirius_adduct", "sirius_structure_name", "sirius_structure_smiles",
        "sirius_structure_confidence", "sirius_npc_pathway", "sirius_npc_class", "sirius_classyfire_class",
    ]
    cols = [c for c in cols if c in sig.columns]
    out_path = diff_dir / "compound_summary.tsv"
    sig[cols].to_csv(out_path, sep="\t", index=False)

    n_robust = int(sig["blocking_robust"].sum())
    n_identified = int(sig["best_identity"].notna().sum())
    print(
        f"{diff_dir.name}: {len(sig)} significant (FDR<{fdr:.0%}, {n_robust} plate-blocking-robust), "
        f"{n_identified}/{len(sig)} identified (library/analog/SIRIUS) -> {out_path}",
        file=sys.stderr,
    )
    return sig[cols].assign(comparison=diff_dir.name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*", help="comparison dirs (default: all with differential_features.csv.gz)")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--no-master", action="store_true", help="skip writing the cross-comparison rollup table")
    args = ap.parse_args()

    if args.dirs:
        diff_dirs = [Path(d) for d in args.dirs]
    else:
        diff_dirs = sorted(p.parent for p in DIFF_ROOT.glob("*/differential_features.csv.gz"))

    library = load_gnps(LIBRARY_SEARCH, "library")
    analog = load_gnps(ANALOG_SEARCH, "analog")
    sirius = load_sirius()
    print(
        f"loaded {len(library)} exact + {len(analog)} analog EB library matches"
        + (f", {len(sirius)} SIRIUS annotations" if sirius is not None else " (no SIRIUS annotations yet)"),
        file=sys.stderr,
    )

    all_sig = [summarize_one(d, library, analog, sirius, args.fdr) for d in diff_dirs]

    if not args.no_master and all_sig:
        master = pd.concat(all_sig, ignore_index=True)
        master.to_csv(MASTER_OUT, sep="\t", index=False)
        print(f"wrote {len(master)} rows across {len(diff_dirs)} comparisons -> {MASTER_OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
