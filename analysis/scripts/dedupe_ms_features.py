#!/usr/bin/env python3
"""
Group the 16,332 MS2 features in analysis/linked_data/feature_abundance_matrix.csv.gz
into deduplicated compound-level groups, using annotation EverythingBagel
(EB) already computed but that this project's simplified working matrix
doesn't carry: `isotope_source_id`, `adduct_source_id`, `is_default_adduct`,
`is_isf`, `isf_parent_id` in the fuller per-file EB output
(nf_output/feature_finding/feature_finding_results/aligned_features_ms2.csv).

Why this matters (Phase 2, step 0 of analysis/INTEGRATED_ANALYSIS_STRATEGY.md,
per Fable's 2026-08-15 review): 16,332 aligned MS2 "features" are not
16,332 independent chemical entities. The same underlying compound
routinely produces several aligned features -- isotopologues (M+0, M+1,
M+2, ...) and different adducts ([M+H]+, [M+Na]+, [2M+H]+, ...) of the same
molecule, plus in-source fragments (ISF) that are decomposition products of
a co-eluting parent, not independent compounds at all. BH-FDR in Phase 2
assumes roughly independent tests; testing every isotopologue/adduct/ISF of
one color-associated compound as a separate "hit" inflates both the
apparent number of significant features and (via re-identification in
Phase 3) double-counts evidence for whatever SIRIUS NPC class that
compound belongs to.

This is NOT a from-scratch de-dup -- EverythingBagel already grouped
isotopologues/adducts/ISFs during feature finding (aligned_compounds.tsv,
aligned_isf.tsv in the nf_output). This script simply recovers and applies
that existing grouping to the 16,332-row working feature set, rather than
re-deriving it via ad hoc RT/mass-difference heuristics.

Grouping logic, per feature row:
  1. If `is_isf` is True: fold into the group of its parent feature
     (`isf_parent_id` -> that parent row's own group, resolved
     transitively in case a parent is itself an ISF of something else).
  2. Otherwise: group by `adduct_source_id` -- EB's own isotopologue/
     same-adduct-type grouping (a compound's [M+H]+ M+0/M+1/M+2/... peaks
     all share one adduct_source_id).
  3. Cross-adduct-type merging (the same compound detected as BOTH [M+H]+
     and [M+Na]+, ~4.9% of compounds per aligned_compounds.tsv) is NOT
     applied here -- it requires parsing aligned_compounds.tsv's free-text
     `members` field (m/z:isotope/adduct triples, not row IDs) and cross-
     referencing by m/z, which is a separate, noisier join. Left as a
     documented limitation (see module-level NOTE below) rather than
     silently attempted with a fragile parser.

Within each group, the representative feature (used for e.g. Phase 2's
per-group correlation test) is: the `is_default_adduct == True` row if
present, else the row with the highest `total_scans`.

Usage:
    python3 analysis/scripts/dedupe_ms_features.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
EB_FEATURE_FINDING = (
    REPO
    / "data"
    / "processed"
    / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output"
    / "feature_finding"
    / "feature_finding_results"
)
FULL_FEATURES_CSV = EB_FEATURE_FINDING / "aligned_features_ms2.csv"  # has adduct_source_id/isotope_source_id/is_isf
WORKING_MATRIX = REPO / "analysis" / "linked_data" / "feature_abundance_matrix.csv.gz"  # 16,332 rows this project uses
OUT_CSV = REPO / "analysis" / "linked_data" / "ms_feature_dedup_groups.csv"
OUT_SUMMARY = REPO / "analysis" / "linked_data" / "ms_feature_dedup_summary.txt"

# NOTE: cross-adduct-type compound merging (aligned_compounds.tsv,
# feature_finding_results/) is NOT applied here -- see module docstring.
# ~4.9% of EB "compounds" (1,802/36,785, all detected features not just
# the has-ms2 subset used by this project) span >1 adduct type; this
# script's groups are therefore a conservative (slight over-count, not
# under-count) de-dup relative to the true compound count.


def resolve_isf_root(row_id: int, is_isf: pd.Series, isf_parent_id: pd.Series, row_id_index: pd.Index, max_hops: int = 10) -> int:
    """Follow isf_parent_id chains to the non-ISF root feature's row ID."""
    current = row_id
    for _ in range(max_hops):
        if current not in row_id_index or not bool(is_isf.get(current, False)):
            return current
        parent = isf_parent_id.get(current)
        if pd.isna(parent):
            return current
        parent = int(parent)
        if parent == current:
            return current
        current = parent
    return current  # give up after max_hops, avoid infinite loop on a cycle


def main():
    full = pd.read_csv(FULL_FEATURES_CSV)
    working = pd.read_csv(WORKING_MATRIX, usecols=["row ID"])
    working_ids = set(working["row ID"])

    full = full.set_index("row ID", drop=False)
    is_isf = full["is_isf"].fillna(False).astype(bool)
    isf_parent_id = full["isf_parent_id"]

    # resolve each working-set row to its non-ISF root row ID
    roots = {}
    for rid in working_ids:
        if rid not in full.index:
            roots[rid] = rid  # not found in fuller table (shouldn't happen); keep isolated
            continue
        roots[rid] = resolve_isf_root(rid, is_isf, isf_parent_id, full.index)

    # group key: the root's adduct_source_id (falls back to the root row ID
    # itself if adduct_source_id is missing)
    group_key = {}
    for rid, root in roots.items():
        if root in full.index:
            asid = full.at[root, "adduct_source_id"]
            group_key[rid] = f"asid_{int(asid)}" if pd.notna(asid) else f"row_{root}"
        else:
            group_key[rid] = f"row_{root}"

    df = working.copy()
    df["dedup_group_id"] = df["row ID"].map(group_key)
    df["is_isf_member"] = df["row ID"].map(lambda r: bool(is_isf.get(r, False)) if r in full.index else False)
    df["isf_root_row_id"] = df["row ID"].map(roots)
    df["is_default_adduct"] = df["row ID"].map(lambda r: full.at[r, "is_default_adduct"] if r in full.index else None)
    df["total_scans"] = df["row ID"].map(lambda r: full.at[r, "total_scans"] if r in full.index else None)

    group_sizes = df.groupby("dedup_group_id")["row ID"].transform("size")
    df["group_size"] = group_sizes

    # representative: is_default_adduct==True, else highest total_scans, within each group
    df = df.sort_values(["dedup_group_id", "is_default_adduct", "total_scans"], ascending=[True, False, False])
    df["is_group_representative"] = ~df.duplicated("dedup_group_id", keep="first")
    df = df.sort_values("row ID").reset_index(drop=True)

    n_features = len(df)
    n_groups = df["dedup_group_id"].nunique()
    n_isf = int(df["is_isf_member"].sum())
    multi_member_groups = int((df.groupby("dedup_group_id")["row ID"].transform("size") > 1).sum())

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df[["row ID", "dedup_group_id", "group_size", "is_group_representative", "is_isf_member", "isf_root_row_id"]].to_csv(
        OUT_CSV, index=False
    )

    with OUT_SUMMARY.open("w") as fh:
        fh.write("MS2 feature de-duplication summary (adducts/isotopologues/in-source fragments)\n")
        fh.write("=" * 78 + "\n")
        fh.write(f"Source (fuller EB output with grouping annotation): {FULL_FEATURES_CSV}\n")
        fh.write(f"Working matrix (this project's 16,332-feature set): {WORKING_MATRIX}\n\n")
        fh.write(f"Raw features (rows in working matrix): {n_features}\n")
        fh.write(f"Deduplicated groups (isotopologue + same-adduct-type + ISF collapsed): {n_groups}\n")
        fh.write(f"Reduction: {n_features} -> {n_groups} ({100 * (1 - n_groups / n_features):.1f}% collapsed)\n")
        fh.write(f"Features flagged as in-source fragments (folded into their parent's group): {n_isf}\n")
        fh.write(f"Groups with >1 member feature: {multi_member_groups}\n\n")
        fh.write(
            "NOT applied: cross-adduct-type compound merging (e.g. a compound's [M+H]+ "
            "and [M+Na]+ groups sharing one true compound identity) -- see module "
            "docstring in dedupe_ms_features.py. This means the group count above is a "
            "conservative (slightly too high, never too low) de-dup relative to the "
            "true number of independent chemical entities.\n"
        )

    print(f"{n_features} raw features -> {n_groups} deduplicated groups ({n_isf} ISF features folded into parents)")
    print(f"Wrote {OUT_CSV}")
    print(f"Wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
