#!/usr/bin/env python3
"""
Pull the MS2 spectra for analysis/sirius_annotation/sirius_targets.csv's
target features out of the EB pipeline's own gap-filled MGF (SCANS==row
ID, confirmed directly: SCANS=52 has PEPMASS 455.2985, RT 324.27s, exactly
matching row ID 52's row m/z/RT in the feature table) -- no need to
re-extract from raw mzML the way the prior project's pipeline did, since
the EB pipeline already exported one representative MS2 spectrum per
feature.

Usage:
    python3 scripts/export_sirius_targets_mgf.py
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
TARGETS_CSV = REPO / "analysis" / "sirius_annotation" / "sirius_targets.csv"
SOURCE_MGF_ZST = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "Exfab_-_20260130_ExFAB_Rhodo_--b773ffa18c2b41e5a3484526293a54f9-aligned_features_filled.mgf.zst"
)
OUT_MGF = REPO / "analysis" / "sirius_annotation" / "sirius_targets.mgf"


def parse_mgf_blocks(text: str):
    """Yield (scans_id, full_block_text) for every BEGIN IONS...END IONS block."""
    block_lines = []
    scans = None
    for line in text.splitlines(keepends=True):
        if line.startswith("BEGIN IONS"):
            block_lines = [line]
            scans = None
        elif line.startswith("END IONS"):
            block_lines.append(line)
            yield scans, "".join(block_lines)
        else:
            block_lines.append(line)
            if line.startswith("SCANS="):
                scans = int(line.strip().split("=", 1)[1])


def main():
    targets = pd.read_csv(TARGETS_CSV)
    target_ids = set(targets["row ID"].astype(int))
    print(f"looking up {len(target_ids)} target row IDs in {SOURCE_MGF_ZST.name}", file=sys.stderr)

    text = subprocess.run(
        ["zstd", "-dc", str(SOURCE_MGF_ZST)], stdout=subprocess.PIPE, check=True
    ).stdout.decode()

    found = {}
    for scans, block in parse_mgf_blocks(text):
        if scans in target_ids:
            found[scans] = block

    missing = target_ids - found.keys()
    if missing:
        print(f"WARNING: {len(missing)} target row ID(s) not found in the MGF: {sorted(missing)}", file=sys.stderr)

    OUT_MGF.parent.mkdir(parents=True, exist_ok=True)
    with OUT_MGF.open("w") as fh:
        for scans in sorted(found):
            fh.write(found[scans])

    print(f"wrote {len(found)}/{len(target_ids)} spectra -> {OUT_MGF}", file=sys.stderr)


if __name__ == "__main__":
    main()
