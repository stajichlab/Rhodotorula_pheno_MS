#!/usr/bin/env python3
"""
Pick a bounded, worth-the-compute set of features to run SIRIUS on: the
top --top-n *unidentified* (no existing EB/GNPS exact or analog library
match -- see build_compound_summary.py) plate-blocking-robust significant
features from each given comparison's compound_summary.tsv, deduplicated
by row ID (the same feature often tops more than one comparison, e.g. it
separates a species from several others at once).

Already-identified features are deliberately excluded here -- SIRIUS is
slow and strictly serial (see knowledgebase/sirius.md), so compute goes
to features the existing library search couldn't already name, not ones
it could.

Output: analysis/sirius_annotation/sirius_targets.csv, one row per unique
target feature with its m/z/RT/adduct and which comparison(s)/rank(s) it
came from (provenance, for tracing a SIRIUS hit back to why it was
selected -- same pattern as the prior project's target_provenance.csv).

Usage:
    python3 scripts/select_sirius_targets.py [DIR ...] [--top-n 15]
        (default DIRS: every analysis/differential_features/*/ that has a
        compound_summary.tsv; run build_compound_summary.py first)
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DIFF_ROOT = REPO / "analysis" / "differential_features"
OUT_DIR = REPO / "analysis" / "sirius_annotation"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*")
    ap.add_argument("--top-n", type=int, default=15)
    args = ap.parse_args()

    if args.dirs:
        comp_dirs = [Path(d) for d in args.dirs]
    else:
        comp_dirs = sorted(p.parent for p in DIFF_ROOT.glob("*/compound_summary.tsv"))

    picks = []
    for d in comp_dirs:
        path = d / "compound_summary.tsv"
        if not path.exists():
            print(f"skip {d.name}: no compound_summary.tsv (run build_compound_summary.py first)", file=sys.stderr)
            continue
        df = pd.read_csv(path, sep="\t")
        unident = df[
            df["blocking_robust"] & df["best_identity"].isna() & (df["has_ms2"] == True)  # noqa: E712
        ].sort_values("q_used").head(args.top_n)
        for rank, (_, row) in enumerate(unident.iterrows(), start=1):
            picks.append(
                {
                    "row ID": int(row["row ID"]),
                    "row m/z": row["row m/z"],
                    "row retention time": row["row retention time"],
                    "adduct": row["adduct"],
                    "comparison": d.name,
                    "rank_in_comparison": rank,
                    "q_used": row["q_used"],
                    "log2FC_a_over_b": row["log2FC_a_over_b"],
                }
            )
        print(f"{d.name}: {len(unident)} unidentified top targets picked", file=sys.stderr)

    if not picks:
        sys.exit("no targets picked from any comparison")

    picks_df = pd.DataFrame(picks)
    provenance = (
        picks_df.groupby("row ID")[["comparison", "rank_in_comparison", "q_used"]]
        .apply(lambda g: "; ".join(f"{r.comparison}(rank {r.rank_in_comparison}, q={r.q_used:.2e})" for r in g.itertuples()), include_groups=False)
        .rename("provenance")
        .reset_index()
    )
    unique_targets = picks_df.drop_duplicates("row ID")[
        ["row ID", "row m/z", "row retention time", "adduct"]
    ].merge(provenance, on="row ID").sort_values("row ID")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "sirius_targets.csv"
    unique_targets.to_csv(out_path, index=False)
    print(
        f"{len(picks_df)} picks from {len(comp_dirs)} comparisons -> "
        f"{len(unique_targets)} unique target features -> {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
