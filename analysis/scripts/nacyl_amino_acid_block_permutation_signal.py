#!/usr/bin/env python3
"""
N-acyl amino acid candidates -- block-permutation phylogenetic signal test
(2026-09-16, PI request: "apply the phase2_metabolome_phenotype [approach]
to test if there is a phylogenetic signal in the production of these
compounds").

Reuses this project's own species-tree block-construction method
(`species_blocks_from_tree`, verbatim from
analysis/scripts/phase2_color_metabolome_association.py) rather than
inventing a new one. That script's own design tests compound abundance
AGAINST AN EXTERNAL PHENOTYPE (color, growth) while controlling for
phylogeny as a nuisance variable -- a different question from "is
production itself phylogenetically structured." This script answers the
literal question asked: does species-tree clade membership (the same
6-clade default split phase2 uses) explain more variance in compound
production than expected by chance?

This is a non-parametric COMPLEMENT to
nacyl_amino_acid_phylogenetic_signal.R's Blomberg's K / Pagel's lambda,
not a replacement -- K/lambda assume a Brownian-motion-like continuous
trait and operate on species-level means (n=16); this test uses raw
STRAIN-level data (n=265, more power) and makes no distributional
assumption, which suits these compounds' zero-inflated distributions
better.

Method: one-way ANOVA F-statistic (log10(peak_area+1) ~ clade) on strain-
level cell-fraction data, clade assignment from the species tree (6
clades, average-linkage on patristic distance, same as phase2's default).
Null: permute clade LABELS across strains (unrestricted -- this tests
whether the REAL tree-derived clade assignment explains more variance
than an arbitrary same-sized grouping, i.e. a direct test of phylogenetic
clustering) --n-perm times, empirical p = fraction of permuted F-stats >=
observed.

Usage:
    python3 analysis/scripts/nacyl_amino_acid_block_permutation_signal.py --n-perm 2000
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import f_oneway

REPO = Path(__file__).resolve().parent.parent.parent
COMPARTMENT_CSV = REPO / "analysis" / "ahl_autoinducer_search" / "nacyl_amino_acid_compartment_species.csv"
SPECIES_TREE = REPO / "analysis" / "integrated_analysis" / "phase1_phenotype" / "species_tree.nwk"
OUT_CSV = REPO / "analysis" / "ahl_autoinducer_search" / "nacyl_amino_acid_block_permutation_signal.csv"

COMPOUND_LABELS = {
    4109: "row_4109_palmitoyl_arginine",
    51126: "row_51126_myristoyl_arginine",
    51152: "row_51152_palmitoyl_lysine",
    24998: "row_24998_histidine_MH",
    40740: "row_40740_histidine_MNa",
}


def species_blocks_from_tree(tree_path: Path, n_clades: int) -> dict[str, str]:
    """Verbatim from phase2_color_metabolome_association.py."""
    from Bio import Phylo

    tree = Phylo.read(str(tree_path), "newick")
    terminals = tree.get_terminals()
    names = [t.name for t in terminals]
    n = len(names)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = tree.distance(terminals[i], terminals[j])
            dist[i, j] = dist[j, i] = d
    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    cluster_ids = fcluster(z, t=n_clades, criterion="maxclust")
    return {name.replace("_", " "): f"clade_{c}" for name, c in zip(names, cluster_ids)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--n-clades", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    d = pd.read_csv(COMPARTMENT_CSV)
    d = d[(d.fraction == "cell") & d.species.notna()].copy()
    d["species"] = d["species"].str.replace(r"^(Rhodotorula) sp_clade_(\w+)$", r"\1 sp. clade \2", regex=True)
    d["log_peak"] = np.log10(d["peak_area"] + 1)

    blocks_by_species = species_blocks_from_tree(SPECIES_TREE, args.n_clades)
    missing = set(d["species"]) - set(blocks_by_species)
    if missing:
        raise SystemExit(f"Species not found in tree-derived blocks: {missing}")
    d["clade"] = d["species"].map(blocks_by_species)

    rng = np.random.default_rng(args.seed)
    rows = []
    for row_id, label in COMPOUND_LABELS.items():
        dd = d[d.row_id == row_id]
        y = dd["log_peak"].to_numpy()
        clade = dd["clade"].to_numpy()
        groups = [y[clade == c] for c in np.unique(clade)]
        f_obs, _ = f_oneway(*groups)

        null_f = np.empty(args.n_perm)
        for p in range(args.n_perm):
            perm_clade = rng.permutation(clade)
            perm_groups = [y[perm_clade == c] for c in np.unique(clade)]
            null_f[p] = f_oneway(*perm_groups)[0]
        emp_p = (np.sum(null_f >= f_obs) + 1) / (args.n_perm + 1)

        rows.append(dict(
            compound=label, n_strains=len(dd), n_clades=len(np.unique(clade)),
            f_observed=f_obs, empirical_p=emp_p,
            null_f_mean=null_f.mean(), null_f_p95=np.percentile(null_f, 95),
        ))
        print(f"{label}: n={len(dd)}, F={f_obs:.3f}, empirical_p={emp_p:.4f} "
              f"(null F mean={null_f.mean():.3f}, 95th pct={np.percentile(null_f, 95):.3f})")

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
