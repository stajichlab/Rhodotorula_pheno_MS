#!/usr/bin/env python3
"""Split an MGF file into N roughly-equal-spectra-count shard files.

Splitting is safe here only because the downstream SIRIUS chain is
formula -> fingerprint -> canopus -> structures -> write-summaries: every one
of those is a <COMPOUND TOOL> in SIRIUS's own --help text, i.e. computed
per-spectrum. ZODIAC is explicitly a <DATASET TOOL> ("Identify Molecular
formulas of all compounds in a dataset together") and is deliberately never
included in this pipeline -- sharding would silently change its results.
"""
import argparse
import pathlib
import sys


def split_mgf(input_path: pathlib.Path, out_dir: pathlib.Path, n_shards: int) -> list[pathlib.Path]:
    spectra = []
    current = []
    with open(input_path) as fh:
        for line in fh:
            if line.startswith("BEGIN IONS"):
                current = [line]
            elif line.startswith("END IONS"):
                current.append(line)
                spectra.append(current)
                current = []
            else:
                current.append(line)
    if current and any(line.startswith("BEGIN IONS") for line in current):
        raise ValueError(f"{input_path}: trailing spectrum block without END IONS")
    if not spectra:
        raise ValueError(f"{input_path}: no BEGIN IONS/END IONS blocks found")

    n_shards = min(n_shards, len(spectra))
    out_dir.mkdir(parents=True, exist_ok=True)
    shard_paths = []
    for shard_idx in range(n_shards):
        shard_spectra = spectra[shard_idx::n_shards]
        shard_path = out_dir / f"shard_{shard_idx:03d}.mgf"
        with open(shard_path, "w") as fh:
            for block in shard_spectra:
                fh.writelines(block)
        shard_paths.append(shard_path)
        print(f"{shard_path}: {len(shard_spectra)} spectra", file=sys.stderr)

    total_out = sum(len(spectra[i::n_shards]) for i in range(n_shards))
    assert total_out == len(spectra), (
        f"shard spectra count {total_out} != input spectra count {len(spectra)}"
    )
    print(f"Split {len(spectra)} spectra from {input_path} into {n_shards} shards under {out_dir}", file=sys.stderr)
    return shard_paths


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=pathlib.Path, help="Source MGF file")
    ap.add_argument("--out-dir", required=True, type=pathlib.Path, help="Directory to write shard_NNN.mgf files")
    ap.add_argument("--n-shards", required=True, type=int, help="Number of shards to split into")
    args = ap.parse_args()
    split_mgf(args.input, args.out_dir, args.n_shards)
