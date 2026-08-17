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
    python3 scripts/generate_species_pairs.py --paired [--min-n 3]
"""
import argparse
import itertools
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
LINKED = REPO / "analysis" / "linked_data"
OUT_PATH = REPO / "analysis" / "differential_features" / "pairs_to_run.tsv"
PAIRED_OUT_PATH = REPO / "analysis" / "differential_features" / "paired_species_to_run.tsv"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3, help="minimum samples per species per fraction")
    ap.add_argument("--paired", action="store_true",
                    help="instead of species-pair rows, list each species with >= --min-n STRAINS that have BOTH a cell and a supernatant sample (for --sup-vs-cell)")
    args = ap.parse_args()

    meta = pd.read_csv(LINKED / "sample_metadata.csv.gz")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    if args.paired:
        cell = meta[(meta["fraction"] == "cell")].drop_duplicates(subset="Strain ID").set_index("Strain ID")
        sup = meta[(meta["fraction"] == "supernatant")].drop_duplicates(subset="Strain ID").set_index("Strain ID")
        both = cell.index.intersection(sup.index)
        species_of = cell.loc[both, "Species"]
        paired_counts = species_of.value_counts()
        eligible = sorted(paired_counts[paired_counts >= args.min_n].index)
        with PAIRED_OUT_PATH.open("w") as fh:
            for sp in eligible:
                fh.write(f"{sp}\n")
        print(f"{len(eligible)} species with >= {args.min_n} paired strains -> {PAIRED_OUT_PATH}")
        return

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
