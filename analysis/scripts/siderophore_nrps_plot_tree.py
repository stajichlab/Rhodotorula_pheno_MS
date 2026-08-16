#!/usr/bin/env python3
"""
Render RA_NRPS_candidates.tree.nwk (FastTree, 266 tips) as a labeled,
species-colored figure -- PDF (vector, for close inspection) and PNG (for
quick viewing/embedding in RESULTS.md).

Usage:
    python3 analysis/scripts/siderophore_nrps_plot_tree.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo

REPO = Path(__file__).resolve().parent.parent.parent
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase_siderophore" / "outputs"
GENE_TREE = OUT_DIR / "RA_NRPS_candidates.tree.nwk"
STRAIN_SUMMARY = OUT_DIR / "RA_NRPS_strain_summary.csv"
REF_NAME = "F2DD6D01_006956-T1"
OUTGROUP_TIP = "DBVPG6740_002504-T1"

SPECIES_COLORS = {
    "Rhodotorula mucilaginosa": "#1f77b4", "Rhodotorula paludigena": "#ff7f0e",
    "Rhodotorula toruloides": "#2ca02c", "Rhodotorula sp. clade I": "#d62728",
    "Rhodotorula dairenensis": "#9467bd", "Rhodotorula taiwanensis": "#8c564b",
    "Rhodotorula diobovata": "#e377c2", "Rhodotorula sphaerocarpa": "#7f7f7f",
    "Rhodotorula graminis": "#bcbd22", "Rhodotorula kratochvilovae": "#17becf",
    "Rhodotorula sp. clade XI": "#aec7e8", "Rhodotorula pacifica": "#ffbb78",
    "Rhodotorula evergladensis": "#98df8a", "Rhodotorula glutinis": "#ff9896",
    "Rhodotorula araucariae": "#c5b0d5", "Rhodotorula sp. clade XIII": "#c49c94",
    "Pseudomicrostroma phylloplanum": "#f7b6d2",
}
REF_COLOR = "black"


def main():
    summary = pd.read_csv(STRAIN_SUMMARY)
    protein_to_label = {}
    protein_to_color = {}
    for _, row in summary.iterrows():
        pid = row["sseqid"]
        if pd.isna(pid):
            continue
        sp = row["SPECIES"]
        strain = row["STRAIN"]
        protein_to_label[pid] = f"{sp} {strain}"
        protein_to_color[pid] = SPECIES_COLORS.get(sp, "#333333")
    protein_to_label[REF_NAME] = "REFERENCE R. kratochvilovae Y14"
    protein_to_color[REF_NAME] = REF_COLOR

    tree = Phylo.read(str(GENE_TREE), "newick")

    outgroup = tree.common_ancestor([c for c in tree.get_terminals()
                                     if c.name == OUTGROUP_TIP])
    if tree.root != outgroup:
        tree.root_with_outgroup(outgroup)
    tree.ladderize()

    label_to_color = {}
    for clade in tree.get_terminals():
        color = protein_to_color.get(clade.name, "#999999")
        label = protein_to_label.get(clade.name, clade.name)
        label_to_color[label] = color
        clade.name = label

    n_tips = tree.count_terminals()
    fig_height = max(12, n_tips * 0.11)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    Phylo.draw(tree, do_show=False, axes=ax,
               label_colors=lambda name: label_to_color.get(name, "black"),
               branch_labels=None, show_confidence=False)
    ax.set_title("Rhodotorulic-acid NRPS candidate gene tree (FastTree, protein, 266 tips)\n"
                 "colored by species; branch lengths = substitutions/site", fontsize=10)
    plt.tight_layout()

    png_path = OUT_DIR / "RA_NRPS_candidates.tree.png"
    pdf_path = OUT_DIR / "RA_NRPS_candidates.tree.pdf"
    fig.savefig(png_path, dpi=150)
    fig.savefig(pdf_path)
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
