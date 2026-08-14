#!/usr/bin/env python3
"""
For each differential-feature comparison, join the plate-blocking-robust
significant features to whatever compound identity the Everything-Bagel
(EB/GNPS) pipeline already worked out, as a "best guess so far" -- before
spending SIRIUS time on anything.

Two EB tables carry that existing identity, both keyed by MS2 scan number
(== this project's "row ID" everywhere else, confirmed: MGF SCANS=1 has
m/z 278.1897 / RT 194.16s, matching row ID 1 in the feature table exactly):
  - feature_library_search: exact spectral library matches (GNPS-style,
    cosine + matched_peaks against a reference library spectrum)
  - feature_analog_library_search: analog matches (same library, but the
    precursor mass differs -- i.e. "structurally similar, not identical")
Both are already best-hit-per-feature (1 row per query_scan) in this
dataset, so no extra ranking is needed on our side.

Significance filter: q_value_perm_plate < --fdr when that column is
present and non-null for a row (i.e. the plate-block permutation test
actually ran for it -- see differential_features_by_species.py); falls
back to plain Mann-Whitney q_value for rows where it didn't (e.g. one
group too small to block), flagged in the output so a plate-confounded
hit is never silently presented as blocking-confirmed.

Output per comparison folder: compound_summary.tsv, one row per
significant feature, sorted by (blocking-robust first, then q-value).

Usage:
    python3 scripts/build_compound_summary.py [DIR ...] [--fdr 0.05]
        (default DIRS: every analysis/differential_features/*/ that has a
        differential_features.csv.gz)
"""
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

GNPS_KEEP_COLS = ["NAME", "cosine", "matched_peaks", "ADDUCT", "SMILES", "INCHI", "SPECTRUMID", "ORGANISM"]


def load_gnps(path: Path, prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["query_scan"] = df["query_scan"].astype(int)
    df["cosine"] = df["cosine"].astype(float)
    df = df.sort_values("cosine", ascending=False).drop_duplicates("query_scan", keep="first")
    out = df[["query_scan"] + GNPS_KEEP_COLS].rename(columns={c: f"{prefix}_{c}" for c in GNPS_KEEP_COLS})
    return out.rename(columns={"query_scan": "row ID"})


def summarize_one(diff_dir: Path, library: pd.DataFrame, analog: pd.DataFrame, fdr: float) -> int:
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
    sig = sig.merge(library, on="row ID", how="left").merge(analog, on="row ID", how="left")
    sig["best_identity"] = sig["library_NAME"].fillna(sig["analog_NAME"])
    sig["best_identity_source"] = sig["library_NAME"].notna().map({True: "exact_library_match", False: ""})
    no_exact = sig["library_NAME"].isna() & sig["analog_NAME"].notna()
    sig.loc[no_exact, "best_identity_source"] = "analog_library_match"
    sig.loc[sig["best_identity"].isna(), "best_identity_source"] = "no_EB_library_match -- SIRIUS candidate"

    sig = sig.sort_values(["blocking_robust", "q_used"], ascending=[False, True])

    cols = [
        "row ID", "row m/z", "row retention time", "adduct", "has_ms2",
        "log2FC_a_over_b", "q_value", "q_value_perm_plate", "q_used", "q_used_source", "blocking_robust",
        "best_identity", "best_identity_source",
        "library_NAME", "library_cosine", "library_matched_peaks", "library_SMILES", "library_INCHI", "library_ORGANISM",
        "analog_NAME", "analog_cosine", "analog_matched_peaks", "analog_SMILES", "analog_INCHI", "analog_ORGANISM",
    ]
    cols = [c for c in cols if c in sig.columns]
    out_path = diff_dir / "compound_summary.tsv"
    sig[cols].to_csv(out_path, sep="\t", index=False)

    n_robust = int(sig["blocking_robust"].sum())
    n_identified = int(sig["best_identity"].notna().sum())
    print(
        f"{diff_dir.name}: {len(sig)} significant (FDR<{fdr:.0%}, {n_robust} plate-blocking-robust), "
        f"{n_identified}/{len(sig)} have an existing EB library match -> {out_path}",
        file=sys.stderr,
    )
    return len(sig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*", help="comparison dirs (default: all with differential_features.csv.gz)")
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()

    if args.dirs:
        diff_dirs = [Path(d) for d in args.dirs]
    else:
        diff_dirs = sorted(p.parent for p in DIFF_ROOT.glob("*/differential_features.csv.gz"))

    library = load_gnps(LIBRARY_SEARCH, "library")
    analog = load_gnps(ANALOG_SEARCH, "analog")
    print(f"loaded {len(library)} exact + {len(analog)} analog EB library matches", file=sys.stderr)

    for d in diff_dirs:
        summarize_one(d, library, analog, args.fdr)


if __name__ == "__main__":
    main()
