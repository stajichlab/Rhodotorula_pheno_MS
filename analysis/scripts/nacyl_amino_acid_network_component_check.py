#!/usr/bin/env python3
"""
Molecular-network (GNPS-style) corroboration for the N-acyl amino acid
candidates (2026-09-16, additional validation step -- independent of
everything done so far: this uses neither our own mass-search target
list nor a re-run of SIRIUS, but the raw EB pipeline's own MS2 cosine-
similarity molecular network, built upstream and untouched by this
investigation, analysis/network_components/ already uses the same
network elsewhere in this project).

Question: do rows 4109 (Palmitoyl arginine), 51126 (N-myristoyl-
arginine), and 51152 (N6-Palmitoyl lysine) sit in molecular-network
components with OTHER independently-SIRIUS-annotated members that form a
chemically coherent homologous series (same amino-acid backbone,
different acyl chain length/saturation) -- a third, independent method
(after mass-search + SIRIUS-formula-match, and MS2 fragment check) that
neither of the first two used.

Source: data/processed/.../nf_output/networking/filtered_pairs.tsv
(CLUSTERID1, CLUSTERID2, DeltaMZ, Cosine, ComponentIndex) -- the pipeline's
own cosine-similarity edges between MS2 spectra, thresholded upstream of
this project (not by any script here).

Usage:
    python3 analysis/scripts/nacyl_amino_acid_network_component_check.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
FILTERED_PAIRS = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9" / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output" / "networking" / "filtered_pairs.tsv"
)
FULL_FEATURES_CSV = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9" / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output" / "feature_finding" / "feature_finding_results" / "aligned_features_ms2.csv"
)
SIRIUS_ANNOTATIONS = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"
OUT_DIR = REPO / "analysis" / "ahl_autoinducer_search"
OUT_CSV = OUT_DIR / "nacyl_amino_acid_network_components.csv"

TARGET_ROWS = {
    4109: "N-C16:0-arginine (Palmitoyl arginine)",
    51126: "N-C14:0-arginine (N-myristoyl-arginine)",
    51152: "N-C16:0-lysine (N6-Palmitoyl lysine)",
    24998: "N-C16:0-histidine (supplementary)",
    40740: "N-C16:0-histidine (supplementary)",
}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(FILTERED_PAIRS, sep="\t")
    feat = pd.read_csv(FULL_FEATURES_CSV, usecols=["row ID", "row m/z", "row retention time", "adduct", "is_isf"])
    sirius = pd.read_csv(SIRIUS_ANNOTATIONS, sep="\t")

    rows_out = []
    components_seen = {}
    for row_id, label in TARGET_ROWS.items():
        comp_ids = pairs.loc[(pairs.CLUSTERID1 == row_id) | (pairs.CLUSTERID2 == row_id), "ComponentIndex"].unique()
        if len(comp_ids) == 0:
            print(f"row {row_id} ({label}): singleton, no network edges")
            continue
        for comp_id in comp_ids:
            components_seen[(row_id, label)] = comp_id
            sub = pairs[pairs.ComponentIndex == comp_id]
            members = sorted(set(sub.CLUSTERID1) | set(sub.CLUSTERID2))
            print(f"\nrow {row_id} ({label}) -> ComponentIndex {comp_id}, {len(members)} members")

            mfeat = feat[feat["row ID"].isin(members)].merge(sirius, on="row ID", how="left")
            mfeat = mfeat.sort_values("row m/z")
            for _, r in mfeat.iterrows():
                is_target = r["row ID"] in TARGET_ROWS
                marker = " <-- OUR CANDIDATE" if is_target else ""
                name = r.get("sirius_structure_name")
                cls = r.get("sirius_npc_class")
                print(f"  row {int(r['row ID']):6d}  mz={r['row m/z']:.4f}  rt={r['row retention time']:.2f}  "
                      f"class={cls if pd.notna(cls) else '-':20s} name={name if pd.notna(name) else '-'}{marker}")
                rows_out.append(dict(
                    query_row=row_id, query_label=label, component_index=comp_id,
                    member_row=r["row ID"], member_mz=r["row m/z"], member_rt=r["row retention time"],
                    is_query_row=is_target, sirius_npc_class=cls, sirius_structure_name=name,
                ))

    out = pd.DataFrame(rows_out).drop_duplicates()
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
