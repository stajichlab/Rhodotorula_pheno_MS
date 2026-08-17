#!/usr/bin/env python3
"""
List every unordered pair of species that has enough samples in a given
fraction to make a Mann-Whitney comparison meaningful, for both fractions.
Feeds the SLURM array job in scripts/run_all_differential_pairs.sbatch --
one line per (fraction, species_a, species_b) triple, no header, so
SLURM_ARRAY_TASK_ID can index it directly with `sed -n "${N}p"`.

Species with fewer than --min-n samples in a fraction are dropped for
that fraction (default 3): below that, a rank-based test and a two-group
median fold-change are barely defined, let alone FDR-worth reporting.

Usage:
    python3 scripts/generate_species_pairs.py [--min-n 3]
"""
import argparse
import itertools
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
LINKED = REPO / "analysis" / "linked_data"
OUT_PATH = REPO / "analysis" / "differential_features" / "pairs_to_run.tsv"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3, help="minimum samples per species per fraction")
    args = ap.parse_args()

    meta = pd.read_csv(LINKED / "sample_metadata.csv.gz")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for fraction in ["cell", "supernatant"]:
        counts = meta.loc[meta["fraction"] == fraction, "Species"].value_counts()
        eligible = sorted(counts[counts >= args.min_n].index)
        for species_a, species_b in itertools.combinations(eligible, 2):
            rows.append((fraction, species_a, species_b))

    with OUT_PATH.open("w") as fh:
        for fraction, species_a, species_b in rows:
            fh.write(f"{fraction}\t{species_a}\t{species_b}\n")

    print(f"{len(rows)} pairs (min_n={args.min_n}) -> {OUT_PATH}")


if __name__ == "__main__":
    main()
