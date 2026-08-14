#!/usr/bin/env python3
"""
Sync genome-annotation symlinks from BFD/input_all/* into BFD/input/*
for strains that have phenotype metadata and/or Everything-Bagel (EB)
MS2 data, and report strain coverage across all three data sources.

Matching problem: BFD/input_all/{cds,dna,gff3,pep}/ holds symlinks named
"<Species>_<strain-code>.<suffix>". Metadata files identify strains by
code alone (e.g. "DBVPG_10075", "TFCN_7W-292-3"), and those codes are
sometimes mangled (underscore/dash swaps, stray splits like "TF_CN_").
Strain codes are matched against the trailing portion of each BFD
basename using a separator-insensitive, boundary-safe comparison
(normalized strain code must match the *end* of the normalized BFD
basename) so "TFCN_2M-1-1" does not spuriously match
"..._TFCN_2M-1-14".

Usage:
    python3 scripts/sync_bfd_symlinks.py [--dry-run] [--out summary.tsv]
"""
import argparse
import csv
import gzip
import io
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BFD_INPUT_ALL = REPO / "BFD" / "input_all"
BFD_INPUT = REPO / "BFD" / "input"
METADATA_DIR = REPO / "data" / "metadata"
EB_FEATURES = (
    REPO
    / "data"
    / "processed"
    / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "aligned_features_ms2.csv.zst"
)
SUBDIRS = ("cds", "dna", "gff3", "pep")


def normalize(code: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", code).upper()


def open_csv_rows(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt", newline="") as fh:
            yield from csv.DictReader(fh)
    else:
        with path.open(newline="") as fh:
            yield from csv.DictReader(fh)


def load_bfd_strains():
    """basename (species_strain, no suffix) -> set of subdirs present in."""
    bases = {}
    suffix_re = {
        "cds": r"\.cds-transcripts\.fa$",
        "dna": r"\.scaffolds\.fa$",
        "gff3": r"\.gff3$",
        "pep": r"\.proteins\.fa$",
    }
    for sub in SUBDIRS:
        d = BFD_INPUT_ALL / sub
        if not d.is_dir():
            continue
        pat = re.compile(suffix_re[sub])
        for entry in d.iterdir():
            if not entry.is_symlink():
                continue
            base = pat.sub("", entry.name)
            bases.setdefault(base, set()).add(sub)
    return bases


def match_bfd(strain_code: str, bfd_norm_index: dict):
    """Return list of BFD basenames whose normalized name ends with the
    normalized strain code (separator/case-insensitive, boundary-safe)."""
    n = normalize(strain_code)
    if not n:
        return []
    return [base for norm_base, base in bfd_norm_index if norm_base.endswith(n)]


def load_phenotype_metadata():
    """Strain codes from every metadata csv(.gz) that has a 'Strain' column,
    excluding files that primarily crosswalk to MS sample IDs (Cu_AUC*)."""
    strains = {}
    for path in sorted(METADATA_DIR.rglob("*.csv*")):
        if "cu_auc" in path.name.lower():
            continue
        try:
            rows = list(open_csv_rows(path))
        except (OSError, csv.Error):
            continue
        if not rows or "Strain" not in rows[0]:
            continue
        for row in rows:
            code = (row.get("Strain") or "").strip()
            if code:
                strains.setdefault(code, set()).add(path.relative_to(REPO).as_posix())
    return strains


def load_eb_ms_strains():
    """Strain codes with MS2 (cell and/or supernatant) data actually present
    as peak-area columns in aligned_features_ms2.csv.zst, via the Cu_AUC
    crosswalk (SAMPLE_NAME <-> MS2_SAMPLE_Cell/Supernatant)."""
    cu_auc_files = sorted(METADATA_DIR.rglob("Cu_AUC*.csv*"))
    if not cu_auc_files or not EB_FEATURES.exists():
        return {}

    header_proc = subprocess.run(
        ["zstd", "-dc", str(EB_FEATURES)],
        stdout=subprocess.PIPE,
        check=True,
    )
    header_line = header_proc.stdout.split(b"\n", 1)[0].decode()
    ms_columns = {c.strip() for c in header_line.split(",")}

    strains = {}
    for cu_auc in cu_auc_files:
        for row in open_csv_rows(cu_auc):
            code = (row.get("SAMPLE_NAME") or "").strip()
            if not code:
                continue
            cell = (row.get("MS2_SAMPLE_Cell") or "").strip()
            sup = (row.get("MS2_SAMPLE_Supernatant") or "").strip()
            present = []
            if cell and f"{cell}.mzML Peak area" in ms_columns:
                present.append(f"cell:{cell}")
            if sup and f"{sup}.mzML Peak area" in ms_columns:
                present.append(f"sup:{sup}")
            if present:
                strains.setdefault(code, set()).update(present)
    return strains


def rsync_symlinks(bfd_base: str, dry_run: bool):
    """rsync the symlinks for one strain (all subdirs it exists in) from
    input_all/<sub>/ into input/<sub>/, preserving symlinks as-is."""
    copied = []
    for sub in SUBDIRS:
        src_dir = BFD_INPUT_ALL / sub
        dst_dir = BFD_INPUT / sub
        matches = sorted(src_dir.glob(f"{bfd_base}.*"))
        if not matches:
            continue
        dst_dir.mkdir(parents=True, exist_ok=True)
        cmd = ["rsync", "-a", "--links"]
        if dry_run:
            cmd.append("--dry-run")
        cmd += [str(m) for m in matches] + [str(dst_dir) + "/"]
        subprocess.run(cmd, check=True)
        copied.extend(m.name for m in matches)
    return copied


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true", help="rsync -n, no changes made")
    ap.add_argument(
        "--out",
        default=str(REPO / "BFD" / "strain_coverage_summary.tsv"),
        help="path to write the strain coverage summary TSV",
    )
    args = ap.parse_args()

    print("Loading BFD/input_all symlink inventory...", file=sys.stderr)
    bfd_bases = load_bfd_strains()
    bfd_norm_index = [(normalize(b), b) for b in bfd_bases]

    print("Loading phenotype metadata...", file=sys.stderr)
    pheno_strains = load_phenotype_metadata()

    print("Loading EB MS2 strain crosswalk...", file=sys.stderr)
    ms_strains = load_eb_ms_strains()

    all_strains = sorted(set(pheno_strains) | set(ms_strains), key=str.upper)
    print(f"{len(all_strains)} distinct strain codes across metadata sources", file=sys.stderr)

    rows = []
    synced_bases = set()
    for strain in all_strains:
        in_pheno = strain in pheno_strains
        in_ms = strain in ms_strains
        matches = match_bfd(strain, bfd_norm_index)
        in_bfd = bool(matches)
        do_sync = in_bfd and (in_pheno or in_ms)
        rows.append(
            {
                "strain": strain,
                "in_phenotype_metadata": "yes" if in_pheno else "no",
                "in_eb_ms_data": "yes" if in_ms else "no",
                "in_bfd_sequences": "yes" if in_bfd else "no",
                "bfd_match": ";".join(matches) if matches else "",
                "rsynced": "yes" if do_sync else "no",
            }
        )
        if do_sync:
            synced_bases.update(matches)

    print(f"Rsyncing symlinks for {len(synced_bases)} matched strains...", file=sys.stderr)
    for base in sorted(synced_bases):
        rsync_symlinks(base, args.dry_run)

    out_path = Path(args.out)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "strain",
                "in_phenotype_metadata",
                "in_eb_ms_data",
                "in_bfd_sequences",
                "bfd_match",
                "rsynced",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    n_pheno = sum(1 for r in rows if r["in_phenotype_metadata"] == "yes")
    n_ms = sum(1 for r in rows if r["in_eb_ms_data"] == "yes")
    n_bfd = sum(1 for r in rows if r["in_bfd_sequences"] == "yes")
    n_synced = sum(1 for r in rows if r["rsynced"] == "yes")
    n_missing_bfd = sum(
        1 for r in rows if r["in_bfd_sequences"] == "no"
    )
    print(f"\nSummary written to {out_path}", file=sys.stderr)
    print(
        f"strains total={len(rows)} in_phenotype_metadata={n_pheno} "
        f"in_eb_ms_data={n_ms} in_bfd_sequences={n_bfd} "
        f"missing_from_bfd={n_missing_bfd} rsynced={n_synced}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
