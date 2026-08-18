#!/usr/bin/env python3
"""Build a SIRIUS-cross-referenced compound report for the network-component
association analysis (network_component_association.py).

For every molecular-family (GNPS ComponentIndex) that passed BH-FDR<0.05 in
either design (GSEA-style enrichment ``p_es_fdr`` or grouped ``max|rho|
p_maxrho_fdr``) for any of the six within-R. mucilaginosa trait panels --
plus a small curated list of hand-checked components -- this writes:

  reports/component_by_trait.tsv          association stats per (component, trait)
  reports/component_summary.tsv           one row per component (union), with
                                          member count, identified-member count
                                          and per-trait significance flags
  reports/component_feature_identity.tsv/.csv
                                          one row per member feature, joined to
                                          the Everything-Bagel/GNPS library search
                                          (exact + analog) and to SIRIUS
                                          formula/structure/CANOPUS annotations
                                          (same identity schema and precedence as
                                          scripts/build_compound_summary.py)
  reports/component_feature_identity.html
                                          the same table rendered with the shared
                                          sortable/filterable HTML template
                                          (scripts/generate_compound_table_html.py)

Usage:
    python3 analysis/scripts/network_component_sirius_report.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import build_compound_summary as bcs  # noqa: E402
from generate_compound_table_html import build_html  # noqa: E402

OUT = REPO / "analysis" / "network_components" / "outputs"
REPORT = REPO / "analysis" / "network_components" / "reports"
FEAT = REPO / "analysis" / "linked_data" / "feature_abundance_matrix.csv.gz"
TRAITS = ["growth_cell", "growth_supernatant",
          "color_a_cell", "color_a_supernatant",
          "color_C_cell", "color_C_supernatant"]
FDR = 0.05
CURATED = [887, 1642, 1835, 2262, 1604, 1035]  # hand-checked this session

TITLE = "Network-component associations - SIRIUS/GNPS identity table"


def main() -> None:
    REPORT.mkdir(parents=True, exist_ok=True)

    comp = {t: pd.read_csv(OUT / f"{t}_components.tsv", sep="\t") for t in TRAITS}
    rho = {t: pd.read_csv(OUT / f"{t}_rhoperm.tsv", sep="\t") for t in TRAITS}

    # per (component, trait) association stats, traits with >=2 members only
    by_trait = []
    for t in TRAITS:
        r = comp[t].copy()
        r["trait"] = t
        r["sig_es"] = r["p_es_fdr"] < FDR
        r["sig_maxrho"] = r["p_maxrho_fdr"] < FDR
        r["sig_familywise"] = r["p_fdr_max"] < FDR
        by_trait.append(r)
    by_trait = pd.concat(by_trait, ignore_index=True)

    # reported components = significant in any trait, plus curated list
    sig_cids = set(by_trait.loc[by_trait["sig_es"] | by_trait["sig_maxrho"], "ComponentIndex"])
    curated_cids = [c for c in CURATED if c not in sig_cids]
    reported = set(sig_cids) | set(CURATED)

    # one row per (component, trait) among reported components
    bt = by_trait[by_trait["ComponentIndex"].isin(reported)].copy()
    bt.to_csv(REPORT / "component_by_trait.tsv", sep="\t", index=False)

    # member features: union across traits of rows belonging to a reported comp
    member_rows = []
    for t in TRAITS:
        pf = rho[t].copy()
        pf = pf[pf["ComponentIndex"].isin(reported)]
        pf["trait"] = t
        member_rows.append(pf[["row ID", "ComponentIndex", "trait", "rho_abs", "p_feature_perm"]])
    members = pd.concat(member_rows, ignore_index=True)
    # per-feature rho summary across the traits where it was tested
    fm = members.groupby(["row ID", "ComponentIndex"]).apply(
        lambda g: "; ".join(f"{r.trait}={r.rho_abs:.3f}" for _, r in g.iterrows()),
        include_groups=False,
    ).rename("rho_by_trait").reset_index()

    # feature metadata + library + SIRIUS identity (same schema as build_compound_summary)
    feat = pd.read_csv(FEAT, usecols=["row ID", "row m/z", "row retention time", "adduct", "has_ms2"])
    library = bcs.load_gnps(bcs.LIBRARY_SEARCH, "library")
    analog = bcs.load_gnps(bcs.ANALOG_SEARCH, "analog")
    sirius = bcs.load_sirius()
    print(f"loaded {len(library)} exact + {len(analog)} analog library matches, "
          f"{0 if sirius is None else len(sirius)} SIRIUS annotations", file=sys.stderr)

    ann = fm.merge(feat, on="row ID", how="left")
    ann = bcs.annotate(ann, library, analog, sirius)

    # component-level per-trait significance flags as a single cell
    sigmap = bt[bt["sig_es"] | bt["sig_maxrho"]].groupby("ComponentIndex")["trait"].agg(
        lambda s: "; ".join(sorted(set(s)))
    ).rename("component_traits")

    # component summary (union of reported components)
    summ = ann.groupby("ComponentIndex").agg(
        component_n_members=("row ID", "nunique"),
        n_identified=("best_identity", lambda s: int(s.notna().sum())),
        identities=("best_identity", lambda s: "; ".join(
            dict.fromkeys(x for x in s if isinstance(x, str) and x))[:400]),
    ).reset_index()
    summ["n_identified"] = summ["n_identified"].fillna(0)
    summ["identities"] = summ["identities"].fillna("")
    summ = summ.merge(sigmap, on="ComponentIndex", how="left")
    summ["component_traits"] = summ["component_traits"].fillna("")
    summ["curated_only"] = summ["ComponentIndex"].isin(curated_cids)
    summ["report_reason"] = summ["component_traits"].where(
        summ["component_traits"] != "", "curated (not significant)")
    stats_tbl = bt.groupby("ComponentIndex").apply(
        lambda g: "; ".join(
            f"{r.trait} rho={r.max_rho_abs:.3f} es_fdr={r.p_es_fdr:.3f} mx_fdr={r.p_maxrho_fdr:.3f}"
            for _, r in g.iterrows()), include_groups=False,
    ).rename("stats_by_trait").reset_index()
    summ = summ.merge(stats_tbl, on="ComponentIndex", how="left")
    summ = summ.sort_values("ComponentIndex")
    summ.to_csv(REPORT / "component_summary.tsv", sep="\t", index=False)

    # final per-feature compound table
    ann = ann.merge(summ[["ComponentIndex", "component_n_members", "component_traits"]],
                    on="ComponentIndex", how="left")
    out_cols = [c for c in [
        "best_identity", "best_identity_source",
        "row ID", "row m/z", "row retention time", "adduct", "has_ms2",
        "ComponentIndex", "component_n_members", "component_traits", "rho_by_trait",
        "library_NAME", "library_cosine", "library_matched_peaks",
        "library_SMILES", "library_INCHI", "library_ORGANISM",
        "analog_NAME", "analog_cosine", "analog_matched_peaks",
        "analog_SMILES", "analog_INCHI", "analog_ORGANISM",
        "sirius_formula", "sirius_adduct", "sirius_structure_name", "sirius_structure_smiles",
        "sirius_structure_confidence", "sirius_npc_pathway", "sirius_npc_class",
        "sirius_classyfire_class", "source_run",
    ] if c in ann.columns]
    out = ann[out_cols].drop_duplicates("row ID").sort_values(["ComponentIndex", "row ID"])
    out.to_csv(REPORT / "component_feature_identity.tsv", sep="\t", index=False)
    out.to_csv(REPORT / "component_feature_identity.csv", index=False)

    n_ident = int(out["best_identity"].notna().sum())
    print(f"wrote {len(out)} feature rows across {len(reported)} components "
          f"({n_ident} identified; curated-only: {curated_cids})", file=sys.stderr)

    html = build_html(out, TITLE)
    (REPORT / "component_feature_identity.html").write_text(html)
    print(f"wrote {REPORT / 'component_feature_identity.html'}", file=sys.stderr)


if __name__ == "__main__":
    main()
