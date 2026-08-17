#!/usr/bin/env python3
"""
Sanity-check congruence between the *fixed* phenotype metadata and the
*fixed* EB merged-metadata sample sheet (see scripts/fix_metadata_strain_ids.py,
which derives both and does the actual strain-code repair work).

This script re-derives nothing; it just verifies that every EB sample
row's canonical_strain (written by fix_metadata_strain_ids.py) is present
in the fixed phenotype metadata's Strain column, and reports any that
aren't. Rows with ATTRIBUTE_TYPE in {Blank, QC, SPE_blank} are controls
and are expected to have no phenotype match -- they are excluded up
front. Run fix_metadata_strain_ids.py first if the *.fixed.* files don't
exist yet or are stale.

Usage:
    python3 scripts/check_eb_metadata_congruence.py [--out report.tsv]
"""
import argparse
import csv
import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PHENOTYPE_CSV = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "YPD2_phenotypic.20260702.fixed.csv.gz"
EB_MERGED_METADATA = (
    REPO
    / "data"
    / "processed"
    / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "Exfab_-_20260130_ExFAB_Rhodo_--b773ffa18c2b41e5a3484526293a54f9-merged_metadata.fixed.tsv.gz"
)
CONTROL_TYPES = {"Blank", "QC", "SPE_blank"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out",
        default=str(REPO / "BFD" / "eb_metadata_phenotype_mismatch.tsv"),
        help="path to write the mismatch report TSV",
    )
    args = ap.parse_args()

    if not PHENOTYPE_CSV.exists() or not EB_MERGED_METADATA.exists():
        sys.exit(
            "fixed metadata files not found -- run "
            "scripts/fix_metadata_strain_ids.py first"
        )

    with gzip.open(PHENOTYPE_CSV, "rt", newline="") as fh:
        pheno_strains = {
            row["Strain"].strip() for row in csv.DictReader(fh) if row.get("Strain", "").strip()
        }

    with gzip.open(EB_MERGED_METADATA, "rt", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    n_excluded = 0
    n_matched = 0
    mismatches = []
    for row in rows:
        if row["ATTRIBUTE_TYPE"] in CONTROL_TYPES:
            n_excluded += 1
            continue
        canonical = row.get("canonical_strain", "").strip()
        if canonical and canonical in pheno_strains:
            n_matched += 1
        else:
            mismatches.append(
                (
                    row["filename"],
                    row["ATTRIBUTE_ID_1"],
                    canonical,
                    row.get("canonical_strain_note", ""),
                )
            )

    out_path = Path(args.out)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["filename", "attribute_id_1", "canonical_strain", "note"])
        writer.writerows(mismatches)

    print(f"EB sample rows: {len(rows) - n_excluded}", file=sys.stderr)
    print(f"excluded (Blank/QC/SPE_blank): {n_excluded}", file=sys.stderr)
    print(f"matched to fixed phenotype metadata: {n_matched}", file=sys.stderr)
    print(f"no phenotype match: {len(mismatches)}", file=sys.stderr)
    print(f"report written to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
