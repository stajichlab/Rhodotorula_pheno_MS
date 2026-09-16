#!/usr/bin/env python3
"""
Decoy/permutation null for the AHL and N-acyl amino acid targeted mass
searches (2026-09-16, PI request -- validation step owed to both since
their target lists were built, see AHL_AUTOINDUCER_SEARCH.md and
NACYL_AMINO_ACID_SEARCH.md "what has not been done").

Question: how many raw feature matches would a mass-only search of this
size, at this ppm tolerance, produce in this feature table BY CHANCE,
if the target masses were not real AHL / N-acyl-amino-acid masses?

Decoy design: each permutation independently shifts EVERY target m/z
(already including its adduct) by a random offset, sign randomized,
magnitude drawn from Uniform(15, 60) Da. This range is chosen to move
each decoy target well clear of its own real value and of common
adduct/isotope mass differences (e.g. +1.003 C-13, +21.98 Na-H, +17.03
NH3) that could otherwise let a "decoy" accidentally land back on a
real, chemically meaningful mass difference, while keeping the decoy
targets inside the same broad, densely-populated region of the observed
feature m/z distribution (146-1486 Da, see script docstring below) that
the real target lists (146-493 Da) already sit inside -- so the null
matches the real search's mass region, not an arbitrary/unrepresentative
one. This is a documented, arbitrary-but-motivated choice, not a
statistically unique "correct" decoy design; a different shift range is
a legitimate sensitivity check if wanted.

For each of --n-perm permutations: recompute the same 20 ppm window
search used by the real searches, record (a) total raw matches (a
feature can match >1 shifted target) and (b) distinct matched row IDs.
Empirical p = (permutations with null count >= observed count + 1) /
(n_perm + 1).

Usage:
    python3 analysis/scripts/mass_search_decoy_null.py --which ahl --n-perm 1000
    python3 analysis/scripts/mass_search_decoy_null.py --which nacyl --n-perm 1000
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
FULL_FEATURES_CSV = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9" / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output" / "feature_finding" / "feature_finding_results" / "aligned_features_ms2.csv"
)
OUT_DIR = REPO / "analysis" / "ahl_autoinducer_search"

CONFIGS = {
    "ahl": dict(
        target_list=OUT_DIR / "ahl_target_list.csv",
        observed_raw=103, observed_distinct=103,
        out_csv=OUT_DIR / "ahl_decoy_null.csv",
    ),
    "nacyl": dict(
        target_list=OUT_DIR / "nacyl_amino_acid_target_list.csv",
        observed_raw=918, observed_distinct=541,
        out_csv=OUT_DIR / "nacyl_amino_acid_decoy_null.csv",
    ),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--which", choices=list(CONFIGS), required=True)
    ap.add_argument("--ppm", type=float, default=20.0)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--shift-lo", type=float, default=15.0)
    ap.add_argument("--shift-hi", type=float, default=60.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    cfg = CONFIGS[args.which]
    targets = pd.read_csv(cfg["target_list"])
    target_mz = targets["target_mz"].to_numpy()
    n_targets = len(target_mz)

    feat = pd.read_csv(FULL_FEATURES_CSV, usecols=["row ID", "row m/z"])
    mz_sorted_idx = np.argsort(feat["row m/z"].to_numpy())
    mz_sorted = feat["row m/z"].to_numpy()[mz_sorted_idx]
    row_id_sorted = feat["row ID"].to_numpy()[mz_sorted_idx]

    rng = np.random.default_rng(args.seed)
    n_raw = np.empty(args.n_perm, dtype=np.int64)
    n_distinct = np.empty(args.n_perm, dtype=np.int64)

    for p in range(args.n_perm):
        sign = rng.choice([-1.0, 1.0], size=n_targets)
        magnitude = rng.uniform(args.shift_lo, args.shift_hi, size=n_targets)
        decoy_mz = target_mz + sign * magnitude

        raw_count = 0
        matched_rows = set()
        for dm in decoy_mz:
            tol = dm * args.ppm / 1e6
            lo_i = np.searchsorted(mz_sorted, dm - tol, side="left")
            hi_i = np.searchsorted(mz_sorted, dm + tol, side="right")
            if hi_i > lo_i:
                raw_count += hi_i - lo_i
                matched_rows.update(row_id_sorted[lo_i:hi_i].tolist())
        n_raw[p] = raw_count
        n_distinct[p] = len(matched_rows)

    p_raw = (np.sum(n_raw >= cfg["observed_raw"]) + 1) / (args.n_perm + 1)
    p_distinct = (np.sum(n_distinct >= cfg["observed_distinct"]) + 1) / (args.n_perm + 1)

    out = pd.DataFrame({"perm": np.arange(args.n_perm), "n_raw_matches": n_raw, "n_distinct_rows": n_distinct})
    out.to_csv(cfg["out_csv"], index=False)

    print(f"[{args.which}] n_targets={n_targets}, n_perm={args.n_perm}, ppm={args.ppm}, "
          f"decoy shift=U({args.shift_lo},{args.shift_hi}) Da (signed)")
    print(f"  observed: {cfg['observed_raw']} raw matches, {cfg['observed_distinct']} distinct rows")
    print(f"  null raw matches:      mean={n_raw.mean():.1f} sd={n_raw.std():.1f} "
          f"median={np.median(n_raw):.1f} 95th pct={np.percentile(n_raw, 95):.1f} max={n_raw.max()}")
    print(f"  null distinct rows:    mean={n_distinct.mean():.1f} sd={n_distinct.std():.1f} "
          f"median={np.median(n_distinct):.1f} 95th pct={np.percentile(n_distinct, 95):.1f} max={n_distinct.max()}")
    print(f"  empirical p (raw matches >= observed):      {p_raw:.4f}")
    print(f"  empirical p (distinct rows >= observed):    {p_distinct:.4f}")
    print(f"\nWrote {cfg['out_csv']}")


if __name__ == "__main__":
    main()
