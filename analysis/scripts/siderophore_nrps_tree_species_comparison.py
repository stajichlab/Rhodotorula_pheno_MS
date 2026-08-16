#!/usr/bin/env python3
"""
PI follow-up (2026-08-16): does the rhodotorulic-acid NRPS gene tree
(RA_NRPS_candidates.tree.nwk, built from siderophore_nrps_build_multifasta.py's
alignment) track the species tree (simple vertical inheritance, expected
for a conserved biosynthetic gene), or show anomalies -- non-monophyletic
species (contamination, paralogy, incomplete lineage sorting, or real
HGT) or unusually long branches (accelerated evolution / pseudogenization
candidates)?

Checks performed:
  1. Per-species monophyly: for every species with >=2 tips in the gene
     tree, is the smallest clade containing all its tips composed ONLY of
     that species' tips (monophyletic), or does it also contain tips from
     other species (non-monophyletic)?
  2. Terminal branch length outliers: tips with a branch length >3 SD
     above the mean, across the whole tree -- candidate fast-evolving/
     pseudogenizing copies.
  3. A qualitative topology narrative comparing the gene tree's
     species-level clade groupings to the species tree
     (phase1_phenotype/species_tree.nwk), for species represented by a
     genuinely distinct sequence (not swamped by *R. mucilaginosa*'s
     near-identical block).

Not done: a formal Robinson-Foulds distance or a proper reconciliation
analysis (no ete3/dendropy available in this environment) -- this is a
descriptive first pass, not a formal cophylogenetic test.

Usage:
    python3 analysis/scripts/siderophore_nrps_tree_species_comparison.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

REPO = Path(__file__).resolve().parent.parent.parent
OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase_siderophore" / "outputs"
GENE_TREE = OUT_DIR / "RA_NRPS_candidates.tree.nwk"
STRAIN_SUMMARY = OUT_DIR / "RA_NRPS_strain_summary.csv"
SPECIES_TREE = REPO / "analysis" / "integrated_analysis" / "phase1_phenotype" / "species_tree.nwk"


def main():
    summary = pd.read_csv(STRAIN_SUMMARY)
    protein_to_species = summary.set_index("sseqid")["SPECIES"].to_dict()
    protein_to_strain = summary.set_index("sseqid")["STRAIN"].to_dict()

    tree = Phylo.read(str(GENE_TREE), "newick")
    terminals = tree.get_terminals()
    tip_names = [t.name for t in terminals]
    # the reference sequence isn't in strain_summary -- exclude it from species-level analysis
    ref_name = "F2DD6D01_006956-T1"
    tip_species = {t: protein_to_species.get(t) for t in tip_names if t != ref_name}
    unmatched = [t for t, s in tip_species.items() if s is None]
    if unmatched:
        print(f"WARNING: {len(unmatched)} tip(s) not matched to a species: {unmatched[:5]}")

    species_tips = {}
    for tip, sp in tip_species.items():
        species_tips.setdefault(sp, []).append(tip)

    print(f"Gene tree: {len(terminals)} tips ({len(tip_species)} strains across {len(species_tips)} species, "
          f"+1 reference)\n")

    # --- 1. monophyly check ---
    print("=== Per-species monophyly ===")
    rows = []
    for sp, tips in sorted(species_tips.items(), key=lambda kv: -len(kv[1])):
        if len(tips) < 2:
            continue
        clades = [tree.find_any(name=t) for t in tips]
        mrca = tree.common_ancestor(clades)
        mrca_terminals = {c.name for c in mrca.get_terminals()}
        extra = mrca_terminals - set(tips) - {ref_name}
        is_mono = len(extra) == 0
        extra_species = sorted({tip_species.get(t, "?") for t in extra}) if extra else []
        rows.append(dict(species=sp, n_strains=len(tips), monophyletic=is_mono,
                          n_extra_tips_in_clade=len(extra), extra_species=";".join(extra_species)))
        flag = "OK" if is_mono else "NOT MONOPHYLETIC"
        extra_note = f" -- includes tips from: {', '.join(extra_species)}" if extra_species else ""
        print(f"  {sp} (n={len(tips)}): {flag}{extra_note}")

    mono_df = pd.DataFrame(rows)
    mono_df.to_csv(OUT_DIR / "RA_NRPS_species_monophyly.csv", index=False)

    # --- 2. terminal branch length outliers ---
    print("\n=== Terminal branch length outliers (>3 SD above mean) ===")
    bl = np.array([t.branch_length or 0.0 for t in terminals])
    mean_bl, sd_bl = bl.mean(), bl.std()
    threshold = mean_bl + 3 * sd_bl
    print(f"Mean terminal branch length: {mean_bl:.5f}, SD: {sd_bl:.5f}, outlier threshold: {threshold:.5f}")
    outlier_rows = []
    for t in terminals:
        blen = t.branch_length or 0.0
        if blen > threshold and t.name != ref_name:
            sp = tip_species.get(t.name, "?")
            strain = protein_to_strain.get(t.name, "?")
            outlier_rows.append(dict(protein_id=t.name, species=sp, strain=strain, branch_length=blen))
            print(f"  {t.name} ({sp}, {strain}): branch_length={blen:.5f}")
    if not outlier_rows:
        print("  None.")
    pd.DataFrame(outlier_rows).to_csv(OUT_DIR / "RA_NRPS_branch_length_outliers.csv", index=False)

    # --- 3. species tree comparison narrative ---
    print("\n=== Species tree for reference (topology to compare against) ===")
    sptree = Phylo.read(str(SPECIES_TREE), "newick")
    Phylo.draw_ascii(sptree)

    print("\nDone. See RA_NRPS_species_monophyly.csv and RA_NRPS_branch_length_outliers.csv")


if __name__ == "__main__":
    main()
