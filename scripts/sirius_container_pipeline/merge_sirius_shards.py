#!/usr/bin/env python3
"""Concatenate per-shard SIRIUS summary TSVs into one merged file per table.

Each shard's summary tables (formula_identifications.tsv, structure_identifications.tsv,
canopus_formula_summary.tsv, canopus_structure_summary.tsv, denovo_structure_identifications.tsv,
spectral_matches.tsv, spectral_matches_analog.tsv) share one column schema across shards --
their rank columns (e.g. formulaRank, structurePerIdRank) are per-compound, not global, so a
straight concatenation with a single retained header is correct; no renumbering needed.
"""
import argparse
import pathlib
import sys

SUMMARY_FILES = [
    "formula_identifications.tsv",
    "structure_identifications.tsv",
    "denovo_structure_identifications.tsv",
    "canopus_formula_summary.tsv",
    "canopus_structure_summary.tsv",
    "spectral_matches.tsv",
    "spectral_matches_analog.tsv",
]


def merge(shard_root: pathlib.Path, out_dir: pathlib.Path) -> None:
    shard_dirs = sorted(p for p in shard_root.glob("shard_*") if p.is_dir())
    if not shard_dirs:
        raise ValueError(f"No shard_* directories found under {shard_root}")
    out_dir.mkdir(parents=True, exist_ok=True)

    for fname in SUMMARY_FILES:
        header = None
        rows = []
        n_shards_with_file = 0
        for shard_dir in shard_dirs:
            fpath = shard_dir / fname
            if not fpath.exists():
                continue
            n_shards_with_file += 1
            with open(fpath) as fh:
                lines = fh.read().splitlines()
            if not lines:
                continue
            shard_header, *shard_rows = lines
            if header is None:
                header = shard_header
            elif shard_header != header:
                raise ValueError(
                    f"{fpath}: header does not match first shard's header for {fname}\n"
                    f"  expected: {header}\n  got:      {shard_header}"
                )
            rows.extend(shard_rows)

        if header is None:
            print(f"{fname}: not present in any shard, skipping", file=sys.stderr)
            continue

        out_path = out_dir / fname
        with open(out_path, "w") as fh:
            fh.write(header + "\n")
            for row in rows:
                fh.write(row + "\n")
        print(f"{fname}: merged {len(rows)} rows from {n_shards_with_file}/{len(shard_dirs)} shards -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shard-root", required=True, type=pathlib.Path, help="Dir containing shard_NNN/ subdirs (SIRIUS project outputs per shard)")
    ap.add_argument("--out-dir", required=True, type=pathlib.Path, help="Dir to write merged summary TSVs")
    args = ap.parse_args()
    merge(args.shard_root, args.out_dir)
