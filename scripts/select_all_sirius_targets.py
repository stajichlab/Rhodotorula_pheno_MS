#!/usr/bin/env python3
"""
Select ALL significant features (with MS2) not already in
analysis/sirius_annotation/sirius_annotations.tsv for a comprehensive
SIRIUS annotation run. Unlike select_sirius_targets.py (which picks only
the top-N unidentified per comparison), this outputs every feature that
still needs SIRIUS -- the full set for sharded, serial array execution.

Output: analysis/sirius_annotation/sirius_targets.csv (overwritten),
one row per target feature with row ID, m/z, RT, and adduct.

Usage:
    python3 scripts/select_all_sirius_targets.py
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SUMMARY = REPO / "analysis" / "differential_features" / "all_significant_features_summary.tsv"
EXISTING = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"
OUT_PATH = REPO / "analysis" / "sirius_annotation" / "sirius_targets.csv"


def main():
    df = pd.read_csv(SUMMARY, sep="\t")

    # SIRIUS requires MS2 spectra
    df = df[df["has_ms2"] == True]  # noqa: E712

    # Skip features already annotated by SIRIUS
    already_done = set()
    if EXISTING.exists():
        existing = pd.read_csv(EXISTING, sep="\t")
        already_done = set(existing["row ID"].astype(int))

    targets = df[~df["row ID"].isin(already_done)]
    targets = targets.drop_duplicates("row ID")

    out_cols = ["row ID", "row m/z", "row retention time", "adduct"]
    out_cols = [c for c in out_cols if c in targets.columns]

    targets[out_cols].to_csv(OUT_PATH, index=False)

    print(
        f"{len(df)} significant features with MS2, "
        f"{len(already_done)} already annotated, "
        f"{len(targets)} targets -> {OUT_PATH}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
