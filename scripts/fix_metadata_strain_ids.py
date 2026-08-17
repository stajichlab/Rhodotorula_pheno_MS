#!/usr/bin/env python3
"""
Repair missing/corrupted strain codes in the phenotype metadata and
propagate the fix into the EB merged-metadata sample sheet, writing both
as new "*.fixed.*" files (originals are left untouched).

Background
----------
data/metadata/EXFAB_UCR-005/YPD2_phenotypic.20260702.csv has a "Strain ID"
(integer) column and a free-text "Strain Name" column (first line is the
raw collection code, e.g. "25-327C-3" or "83E1"); a third "Strain" column
is supposed to hold the canonicalized code ("TFCN_25-327C-3"). For ~17
rows "Strain" is blank or wrong -- in a few cases (id 270) because the
value was hand/Excel-corrupted (e.g. "TFCN_7_6_3" for what should be
"TFCN_7-6-2003"), in most cases because it was simply never filled in.

The EB merged-metadata TSV (Exfab_..._merged_metadata.tsv.gz) identifies
samples by filename "C_<id>.mzML" / "SUP_<id>.mzML", where <id> is the
*same* integer as YPD2's "Strain ID" (verified: of 589 sample rows only
one, id 328/JBEI-13807, has no corresponding YPD2 row at all). That makes
"Strain ID" a much more reliable join key than fuzzy string matching on
the corrupted strain text.

Approach
--------
1. Re-derive "Strain" for every YPD2 row from "Strain Name" using the
   same collection-prefix convention already used by the ~300 already-
   correct rows (TFCN_ for bare TFCN-style codes, DBVPG_ for bare
   numbers, NRRL_Y-<n> for bare Y<n> codes; codes that already carry a
   known prefix are normalized, not re-prefixed). Rows where the
   existing "Strain" already agrees are left alone; rows where it's
   blank or disagrees are corrected and flagged in two new audit
   columns (Strain_original, Strain_fix_method).
2. For every EB merged-metadata sample row, look up the fixed YPD2
   Strain by the shared numeric id and write it into a new
   "canonical_strain" column. Control rows (Blank/QC/SPE_blank) get
   canonical_strain = their type, uppercased. Rows with no YPD2 match
   (id 328) are left blank and flagged in a "canonical_strain_note"
   column for manual review -- this script does not guess codes it
   cannot derive from either source file.

Usage:
    python3 scripts/fix_metadata_strain_ids.py
"""
import csv
import gzip
import re
import sys
from pathlib import Path

def resolve_input(*candidates: Path) -> Path:
    """Return the first candidate path that exists. Metadata files in this
    repo get gzip-compressed by an external process from time to time
    (loose .csv today, .csv.gz tomorrow), so scripts resolve at run time
    rather than hardcoding one extension."""
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(f"none of {candidates} exist")


def smart_open(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, newline="") if path.suffix == ".gz" else path.open(
        mode.replace("t", ""), newline=""
    )


REPO = Path(__file__).resolve().parent.parent
YPD2_STEM = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "YPD2_phenotypic.20260702"
YPD2_OUT = YPD2_STEM.with_name(YPD2_STEM.name + ".fixed.csv.gz")
EB_MERGED_IN = (
    REPO
    / "data"
    / "processed"
    / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "Exfab_-_20260130_ExFAB_Rhodo_--b773ffa18c2b41e5a3484526293a54f9-merged_metadata.tsv.gz"
)
EB_MERGED_OUT = (
    REPO
    / "data"
    / "processed"
    / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "Exfab_-_20260130_ExFAB_Rhodo_--b773ffa18c2b41e5a3484526293a54f9-merged_metadata.fixed.tsv.gz"
)
CU_AUC_IN = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "Cu_AUC.20260811.csv.gz"
CU_AUC_OUT = REPO / "data" / "metadata" / "EXFAB_UCR-005" / "Cu_AUC.20260811.fixed.csv.gz"
CONTROL_TYPES = {"Blank", "QC", "SPE_blank"}
NO_MS2_DATA = "No MS2 Data"

KNOWN_PREFIX = re.compile(r"^(DBVPG|NRRL|TFCN|EXF|JES)[_\s-]+(.*)$", re.IGNORECASE)
Y_NUMBER = re.compile(r"^Y[_-]?(\d+)$", re.IGNORECASE)
PURE_DIGITS = re.compile(r"^\d+$")
SAMPLE_FILENAME = re.compile(r"^(C|SUP)_(\d+)\.mzML$")


def derive_strain(strain_name: str):
    """Return (canonical_code, prefix_used, is_clean) from a raw 'Strain
    Name' first line, or (None, None, False) if it's empty.

    is_clean is False when the first line carries more than one
    whitespace-separated token (a secondary BY-collection id riding
    along on the same line, e.g. "102D-1 8130") or non-ASCII characters
    (observed transcription corruption, e.g. "186CL⒖8-2") -- in both
    cases the derived code is still returned for visibility, but callers
    should not use it to silently overwrite an existing value.
    """
    raw = strain_name.split("\n")[0].strip().strip('"').strip()
    if not raw:
        return None, None, False
    is_clean = len(raw.split()) == 1 and raw.isascii()
    m = KNOWN_PREFIX.match(raw)
    if m:
        prefix = m.group(1).upper()
        rest = m.group(2).strip()
        return f"{prefix}_{rest}", prefix, is_clean
    m = Y_NUMBER.match(raw)
    if m:
        return f"NRRL_Y-{m.group(1)}", "NRRL", is_clean
    if PURE_DIGITS.match(raw):
        return f"DBVPG_{raw}", "DBVPG", is_clean
    return f"TFCN_{raw}", "TFCN", is_clean


def existing_prefix(code: str):
    m = re.match(r"^([A-Za-z]+)[_-]", code)
    return m.group(1).upper() if m else None


def normalize(code: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", code).upper()


def fix_ypd2():
    ypd2_in = resolve_input(
        YPD2_STEM.with_name(YPD2_STEM.name + ".csv.gz"),
        YPD2_STEM.with_name(YPD2_STEM.name + ".csv"),
    )
    with smart_open(ypd2_in) as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames) + ["Strain_original", "Strain_fix_method"]
        rows = list(reader)

    n_blank_fixed = 0
    n_blank_unresolved = 0
    n_corrected = 0
    n_unchanged = 0
    review = []
    for row in rows:
        existing = row["Strain"].strip()
        derived, _, is_clean = derive_strain(row["Strain Name"])
        row["Strain_original"] = ""
        row["Strain_fix_method"] = ""

        if not existing:
            if derived is not None and is_clean:
                row["Strain_fix_method"] = "derived_from_blank"
                row["Strain"] = derived
                n_blank_fixed += 1
            else:
                n_blank_unresolved += 1
                review.append(
                    (row["Strain ID"], row["Strain Name"], existing, derived or "", "blank_unresolved")
                )
            continue

        if derived is None or normalize(existing) == normalize(derived):
            n_unchanged += 1
            continue

        # Existing value disagrees with what 'Strain Name' implies. Only
        # auto-correct when the source is unambiguous (single clean token)
        # and the collection-prefix family isn't changing (a family flip,
        # e.g. DBVPG -> NRRL, is exactly the kind of transcription error
        # that needs a human to confirm, not infer).
        if is_clean and existing_prefix(existing) == existing_prefix(derived):
            row["Strain_original"] = existing
            row["Strain_fix_method"] = "corrected_mismatch"
            row["Strain"] = derived
            n_corrected += 1
        else:
            n_unchanged += 1
            review.append(
                (row["Strain ID"], row["Strain Name"], existing, derived, "mismatch_needs_review")
            )

    with gzip.open(YPD2_OUT, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    review_path = YPD2_OUT.with_name("YPD2_phenotypic.20260702.strain_review_needed.tsv.gz")
    with gzip.open(review_path, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["Strain ID", "Strain Name", "existing Strain", "derived guess", "reason"])
        writer.writerows(review)

    print(
        f"YPD2: {len(rows)} rows | unchanged={n_unchanged} "
        f"blank_filled={n_blank_fixed} corrected={n_corrected} "
        f"needs_manual_review={len(review)} (of which blank_unresolved={n_blank_unresolved})",
        file=sys.stderr,
    )
    print(f"wrote {YPD2_OUT}", file=sys.stderr)
    print(f"wrote {review_path}", file=sys.stderr)
    return {row["Strain ID"]: row["Strain"] for row in rows if row["Strain ID"]}


def fix_eb_merged(strain_by_id: dict):
    with gzip.open(EB_MERGED_IN, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames) + [
            "canonical_strain",
            "canonical_strain_note",
        ]
        rows = list(reader)

    n_control = 0
    n_matched = 0
    n_unresolved = 0
    for row in rows:
        row["canonical_strain"] = ""
        row["canonical_strain_note"] = ""
        atype = row["ATTRIBUTE_TYPE"]
        if atype in CONTROL_TYPES:
            row["canonical_strain"] = atype.upper()
            n_control += 1
            continue
        if atype != "sample":
            continue
        m = SAMPLE_FILENAME.match(row["filename"])
        sample_id = m.group(2) if m else None
        strain = strain_by_id.get(sample_id) if sample_id else None
        if strain:
            row["canonical_strain"] = strain
            n_matched += 1
        else:
            row["canonical_strain_note"] = "no_matching_ypd2_strain_id"
            n_unresolved += 1

    with gzip.open(EB_MERGED_OUT, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"EB merged metadata: {len(rows)} rows | control={n_control} "
        f"matched={n_matched} unresolved={n_unresolved}",
        file=sys.stderr,
    )
    print(f"wrote {EB_MERGED_OUT}", file=sys.stderr)


def load_real_sample_filenames():
    """Set of 'sample'-type mzML basenames (no extension) that genuinely
    exist in the EB merged metadata, e.g. {'C_190', 'SUP_190', ...}. Used
    to catch Cu_AUC rows that reference a file that was never generated
    (typically because that well ended up run as a blank/QC control
    instead of the intended strain)."""
    with gzip.open(EB_MERGED_IN, "rt", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {
        row["filename"].removesuffix(".mzML")
        for row in rows
        if row["ATTRIBUTE_TYPE"] == "sample"
    }


def fix_cu_auc(strain_by_id: dict):
    if not CU_AUC_IN.exists():
        return
    real_files = load_real_sample_filenames()

    with gzip.open(CU_AUC_IN, "rt", newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames) + [
            "SAMPLE_NAME_original",
            "MS2_SAMPLE_Cell_original",
            "MS2_SAMPLE_Supernatant_original",
            "fix_note",
        ]
        rows = list(reader)

    n_fixed = 0
    for row in rows:
        row["SAMPLE_NAME_original"] = ""
        row["MS2_SAMPLE_Cell_original"] = ""
        row["MS2_SAMPLE_Supernatant_original"] = ""
        row["fix_note"] = ""
        sid = row["Strain ID"].strip()
        notes = []

        expected_name = strain_by_id.get(sid)
        if expected_name and row["SAMPLE_NAME"].strip() != expected_name:
            row["SAMPLE_NAME_original"] = row["SAMPLE_NAME"]
            row["SAMPLE_NAME"] = expected_name
            notes.append("sample_name_did_not_match_strain_id")

        for field, prefix in (("MS2_SAMPLE_Cell", "C_"), ("MS2_SAMPLE_Supernatant", "SUP_")):
            val = row[field].strip()
            if val == NO_MS2_DATA or not val:
                continue
            expected = f"{prefix}{sid}"
            if val == expected and val in real_files:
                continue
            orig_key = f"{field}_original"
            if expected in real_files:
                # The correct file exists under the strain's own id but this
                # row pointed somewhere else (e.g. another strain's file,
                # or a "*_empty" blank-well variant) -- repoint it.
                row[orig_key] = row[field]
                row[field] = expected
                notes.append(f"{field}_corrected_to_{expected}")
            elif val in real_files:
                # References a real file, but not this strain's own file,
                # and this strain's own file doesn't exist. Ambiguous --
                # don't guess, flag for manual review.
                notes.append(f"{field}_points_to_another_strains_file_needs_review")
            else:
                row[orig_key] = row[field]
                row[field] = NO_MS2_DATA
                notes.append(f"{field}_no_real_file_set_to_no_ms2_data")

        if notes:
            row["fix_note"] = ";".join(notes)
            n_fixed += 1

    with gzip.open(CU_AUC_OUT, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Cu_AUC: {len(rows)} rows | fixed={n_fixed}", file=sys.stderr)
    print(f"wrote {CU_AUC_OUT}", file=sys.stderr)


def main():
    strain_by_id = fix_ypd2()
    fix_eb_merged(strain_by_id)
    fix_cu_auc(strain_by_id)


if __name__ == "__main__":
    main()
