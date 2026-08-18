#!/usr/bin/env python3
"""Render explanatory figures for the network-component association report.

Reads the component-level association outputs (analysis/network_components/
outputs) and the cross-referenced identity reports (analysis/network_components/
reports) and writes a small set of figures to analysis/network_components/
figures/ (PNG for embedding in REPORT.md, PDF for vector close-up):

  fig1_sig_components_by_trait      sig component counts per trait panel, both
                                    designs, familywise count annotated
  fig2_effect_vs_enrichment_2x3     6 trait panels: max|rho| (effect) vs
                                    -log10(enrichment perm p), significance tinted
  fig3_component887_by_fraction     comp 887's member features in supernatant vs
                                    cell fraction, labelled with identities
  fig4_permutation_calibration      ECDF of per-component permutation p across
                                    traits vs the uniform null (calibration check)

Usage:
    python3 analysis/scripts/network_component_report_figures.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
OUT = REPO / "analysis" / "network_components" / "outputs"
REP = REPO / "analysis" / "network_components" / "reports"
FIG = REPO / "analysis" / "network_components" / "figures"

TRAITS = [
    ("growth_cell", "Growth AUC - cell"),
    ("growth_supernatant", "Growth AUC - supernatant"),
    ("color_a_cell", r"Color $a^*$ - cell"),
    ("color_a_supernatant", r"Color $a^*$ - supernatant"),
    ("color_C_cell", "Color C* - cell"),
    ("color_C_supernatant", "Color C* - supernatant"),
]
COLORS = {
    "both": "#b2182b",
    "enrichment": "#f4a582",
    "maxrho": "#4393c3",
    "none": "#d9d9d9",
}


def sig_state(es, mx):
    if es and mx:
        return "both"
    if es:
        return "enrichment"
    if mx:
        return "maxrho"
    return "none"


def savefig(fig, name: str) -> None:
    fig.savefig(FIG / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}.png/.pdf", flush=True)


def fig1_by_trait() -> None:
    rows = []
    for tname, _ in TRAITS:
        r = pd.read_csv(OUT / f"{tname}_components.tsv", sep="\t")
        rows.append({
            "trait": tname,
            "enrichment": int((r["p_es_fdr"] < 0.05).sum()),
            "maxrho": int((r["p_maxrho_fdr"] < 0.05).sum()),
            "familywise": int((r["p_fdr_max"] < 0.05).sum()),
        })
    df = pd.DataFrame(rows)
    x = np.arange(len(df))
    w = 0.38
    fig, ax = plt.subplots(figsize=(9.5, 5))
    ax.bar(x - w / 2, df["enrichment"], w, label="Enrichment BH-FDR < 0.05", color="#f4a582", edgecolor="black", lw=0.5)
    ax.bar(x + w / 2, df["maxrho"], w, label="max|rho| BH-FDR < 0.05", color="#4393c3", edgecolor="black", lw=0.5)
    for xi, fw in zip(x, df["familywise"]):
        if fw:
            ax.annotate(f"fit.\n{fw}", (xi, df.loc[xi, "maxrho"] + 4), ha="center", va="bottom", fontsize=8, color="#b2182b")
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for _, lbl in TRAITS], rotation=25, ha="right")
    ax.set_ylabel("Components significant (of 1,153 tested)")
    ax.set_title("GNPS component-level association: significant molecular families per trait\n(within R. mucilaginosa, n=1000 permutations)")
    ax.legend(frameon=False)
    ax.set_ylim(0, df[["enrichment", "maxrho"]].max().max() * 1.18)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig(fig, "fig1_sig_components_by_trait")


def fig2_effect_vs_enrichment() -> None:
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5), sharex=True, sharey=True)
    state_labels = {"both": "significant in BOTH designs", "enrichment": "enrichment only",
                    "maxrho": "max|rho| only", "none": "not significant"}
    legend_handles = []
    for ax, (tname, lbl) in zip(axes.ravel(), TRAITS):
        r = pd.read_csv(OUT / f"{tname}_components.tsv", sep="\t")
        es_sig = r["p_es_fdr"] < 0.05
        mx_sig = r["p_maxrho_fdr"] < 0.05
        y = -np.log10(r["p_es_perm"])
        state = [sig_state(a, b) for a, b in zip(es_sig, mx_sig)]
        for st, c in COLORS.items():
            mask = np.array([s == st for s in state])
            if not mask.any():
                continue
            ax.scatter(r["max_rho_abs"][mask], y[mask], s=8, color=c, alpha=0.8, label=state_labels[st])
            if len(legend_handles) < len(COLORS):
                legend_handles.append(ax.collections[-1])
        ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.7)
        ax.set_title(f"{lbl}\n({len(r)} components)", fontsize=10)
        ax.set_xlim(0, 0.38)
        ax.set_ylim(0, 6.2)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.ravel():
        ax.set_xlabel("max $|\\rho|$ over member features (effect)")
        ax.set_ylabel("$-\\log_{10}$(enrichment perm p)")
    fig.legend(handles=legend_handles, loc="lower center", ncol=len(COLORS),
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Component-level association: effect size vs enrichment significance, per trait\n"
                 "dashed grey line = enrichment perm p = 0.05", y=1.0, fontsize=13)
    fig.tight_layout()
    savefig(fig, "fig2_effect_vs_enrichment_2x3")


def fig3_comp887() -> None:
    sup = pd.read_csv(OUT / "growth_supernatant_rhoperm.tsv", sep="\t")
    cel = pd.read_csv(OUT / "growth_cell_rhoperm.tsv", sep="\t")
    ident = pd.read_csv(REP / "component_feature_identity.tsv", sep="\t")[
        ["row ID", "best_identity", "ComponentIndex"]]
    m = sup[sup["ComponentIndex"] == 887][["row ID", "rho_abs", "p_feature_perm"]].merge(
        cel[cel["ComponentIndex"] == 887][["row ID", "rho_abs"]],
        on="row ID", suffixes=("_sup", "_cell"))
    m = m.merge(ident, on="row ID", how="left")
    m = m.sort_values("rho_abs_sup", ascending=False)

    fig, ax = plt.subplots(figsize=(9, 6.2))
    ax.scatter(m["rho_abs_cell"], m["rho_abs_sup"], s=34, color="#4d4d4d",
               edgecolor="white", lw=0.6, zorder=3)
    lim = max(m["rho_abs_sup"].max(), m["rho_abs_cell"].max()) * 1.05
    ax.plot([0, lim], [0, lim], color="#7f7f7f", ls=":", lw=1, label="equal in both fractions")
    for _, r in m.head(10).iterrows():
        lab = (str(r["best_identity"]).split(";")[0].strip() if pd.notna(r["best_identity"]) else "unidentified")
        if len(lab) > 38:
            lab = lab[:36] + "..."
        ax.annotate(f"{r['rho_abs_sup']:.2f} {lab}", (r["rho_abs_cell"], r["rho_abs_sup"]),
                    fontsize=7.5, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("|Spearman rho| vs growth AUC, cell fraction")
    ax.set_ylabel("|Spearman rho| vs growth AUC, supernatant fraction")
    ax.set_title(f"Component 887 ({len(m)} member features) - the F-002 purine-nucleoside family:\n"
                 "growth association lives in the supernatant, not the cell fraction")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig(fig, "fig3_component887_by_fraction")


def fig4_calibration() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    sel = [("growth_supernatant", "GROWTH sup (signal panel)", "#b2182b"),
           ("growth_cell", "GROWTH cell", "#f4a582"),
           ("color_C_supernatant", "Color C* sup (true-null panel)", "#4393c3")]
    for ax, pcol in zip(axes, ["p_es_perm", "p_maxrho_perm"]):
        for tname, lbl, color in sel:
            r = pd.read_csv(OUT / f"{tname}_components.tsv", sep="\t")
            pv = np.sort(r[pcol].to_numpy())
            ecdf = np.arange(1, len(pv) + 1) / len(pv)
            ax.plot(pv, ecdf, color=color, lw=1.8, label=lbl)
        ax.plot([0, 1], [0, 1], color="grey", ls="--", lw=0.8, label="uniform (untested null)")
        ax.set_xlabel("permutation p  (" + ("enrichment" if pcol == "p_es_perm" else "max$|\\rho|$") + ")")
        ax.set_ylabel("cumulative fraction of 1,153 components")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.suptitle("Component permutation p under calibration: excess of small p only in the growth-supernatant panel", y=1.04)
    fig.tight_layout()
    savefig(fig, "fig4_permutation_calibration")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig1_by_trait()
    fig2_effect_vs_enrichment()
    fig3_comp887()
    fig4_calibration()


if __name__ == "__main__":
    main()
