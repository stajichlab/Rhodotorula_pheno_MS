#!/usr/bin/env python3
"""
Block-constrained permutation testing, shared by differential_features_by_species.py
(and anything else that needs it later).

Why: a plain Mann-Whitney U test between two species' strains treats every
strain as an independent draw. Strains of the same species usually are NOT
independent -- they can share a collection batch/MS run (Library Plate),
collection site, or (the big one) recent common ancestry, none of which the
test knows about. A stats-expert review of this project's first differential-
feature run (diobovata n=9 vs mucilaginosa n=209, 45% of tested features
"significant" at FDR<5%) flagged that fraction as implausibly high for
genuinely independent testing -- a plate/relatedness confound inflating
apparent significance was the leading explanation.

The fix implemented here: instead of asking "is species A's median different
from species B's median" against a null that shuffles species labels freely,
ask it against a null that only shuffles labels *within blocks* (same MS
plate, same collection site, or -- once a phylogeny is available -- the same
clade). If two species differ mainly because of which batch/lineage happened
to get sampled, block-permutation p-values reflect that (they can't get more
significant than chance from a confound the null already contains); if the
difference survives block-permutation, it's not explained by that confound.

Two block sources are implemented:
  - derive_blocks_from_metadata(): any categorical column already in
    sample_metadata.csv.gz (e.g. "Library Plate", "Origin"). Works today.
  - derive_blocks_from_tree(): clades cut from a phylogenetic tree (Newick),
    via patristic distance + hierarchical clustering. The strain<->tree-tip
    mapping comes from BFD/strain_coverage_summary.tsv's bfd_match column
    (already built by scripts/sync_bfd_symlinks.py). Written and ready, but
    UNTESTED against a real tree as of this writing, because the phyling
    pipeline building BFD/results/phyling_pep hasn't finished yet -- run a
    quick sanity check (n_clades vs number of species, a couple of known-
    close strain pairs landing in the same clade) the first time a real tree
    is available before trusting its output.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

REPO = Path(__file__).resolve().parent.parent
STRAIN_COVERAGE_TSV = REPO / "BFD" / "strain_coverage_summary.tsv"


def derive_blocks_from_metadata(meta: pd.DataFrame, sample_ids: list[str], column: str) -> dict:
    """sample_id -> block label, read straight from sample_metadata.csv.gz.
    Missing/NaN values become their own block ('__missing__') rather than
    being dropped, since e.g. a handful of strains with unknown Origin are
    still real samples that should be permutable among themselves."""
    lookup = meta.set_index("sample_id")[column]
    blocks = {}
    for sid in sample_ids:
        val = lookup.get(sid)
        blocks[sid] = "__missing__" if pd.isna(val) else str(val)
    return blocks


def _strain_to_tip_map() -> dict:
    """canonical_strain -> tree tip name, from BFD/strain_coverage_summary.tsv
    (bfd_match column; first match if a strain resolved to more than one
    BFD basename -- see sync_bfd_symlinks.py)."""
    cov = pd.read_csv(STRAIN_COVERAGE_TSV, sep="\t", dtype=str)
    cov = cov[cov["bfd_match"].fillna("") != ""]
    return {row["strain"]: row["bfd_match"].split(";")[0] for _, row in cov.iterrows()}


def derive_blocks_from_tree(
    tree_path: Path,
    sample_to_strain: dict,
    n_clades: int | None = None,
    clade_height: float | None = None,
) -> dict:
    """sample_id -> 'clade_<k>', cut from a Newick tree's patristic distances
    via average-linkage hierarchical clustering. Exactly one of n_clades
    (target number of clusters, scipy criterion='maxclust') or clade_height
    (cut at this patristic distance, criterion='distance') must be given.

    The tree does NOT need to be rooted: patristic (tip-to-tip, summed
    branch-length) distance is a property of the tree's topology and branch
    lengths alone, not of where (or whether) a root is placed, so an
    unrooted phyling/IQ-TREE-style output works as-is.

    Samples whose strain has no sequence in BFD (so no tree tip) are simply
    absent from the returned dict -- callers should drop them from the
    blocked test (with a printed count) rather than error, since incomplete
    BFD coverage is expected (see BFD/strain_coverage_summary.tsv:
    missing_from_bfd).
    """
    from Bio import Phylo  # deferred import: only needed on this code path

    if (n_clades is None) == (clade_height is None):
        raise ValueError("give exactly one of n_clades or clade_height")

    strain_to_tip = _strain_to_tip_map()
    tip_by_sample = {}
    for sid, strain in sample_to_strain.items():
        tip = strain_to_tip.get(strain)
        if tip is not None:
            tip_by_sample[sid] = tip
    n_missing = len(sample_to_strain) - len(tip_by_sample)
    if n_missing:
        print(
            f"derive_blocks_from_tree: {n_missing}/{len(sample_to_strain)} samples "
            "have no BFD sequence -> no tree tip -> dropped from blocked test",
            file=sys.stderr,
        )

    tree = Phylo.read(str(tree_path), "newick")
    tip_names_in_tree = {t.name for t in tree.get_terminals()}
    usable = {sid: tip for sid, tip in tip_by_sample.items() if tip in tip_names_in_tree}
    n_not_in_tree = len(tip_by_sample) - len(usable)
    if n_not_in_tree:
        print(
            f"derive_blocks_from_tree: {n_not_in_tree} more sample(s) have a BFD "
            "match not present as a tip in this tree -> also dropped",
            file=sys.stderr,
        )
    if len(usable) < 3:
        raise ValueError(f"only {len(usable)} samples have a usable tree tip -- too few to cluster")

    sample_ids = list(usable.keys())
    tips = [usable[sid] for sid in sample_ids]
    terminal_by_name = {t.name: t for t in tree.get_terminals()}
    n = len(tips)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = tree.distance(terminal_by_name[tips[i]], terminal_by_name[tips[j]])
            dist[i, j] = dist[j, i] = d

    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    if n_clades is not None:
        cluster_ids = fcluster(z, t=n_clades, criterion="maxclust")
    else:
        cluster_ids = fcluster(z, t=clade_height, criterion="distance")

    return {sid: f"clade_{c}" for sid, c in zip(sample_ids, cluster_ids)}


def block_permutation_pvalues(
    mat_a: np.ndarray,
    mat_b: np.ndarray,
    blocks_a: list,
    blocks_b: list,
    n_perm: int = 2000,
    seed: int = 0,
) -> np.ndarray:
    """Two-sided permutation p-values (one per feature/column) for the
    difference in group medians, permuting the group A/B label only within
    each block (so every permuted dataset has the exact same per-block A/B
    counts as the real data -- a block that is 100% one group, e.g. an MS
    plate diobovata was never run on, correctly contributes no information
    either way, rather than being silently dropped or treated as if it did).
    """
    mat_all = np.concatenate([mat_a, mat_b], axis=0)
    labels = np.concatenate([np.zeros(mat_a.shape[0], dtype=bool), np.ones(mat_b.shape[0], dtype=bool)])
    blocks = np.array(list(blocks_a) + list(blocks_b))

    observed = np.abs(np.median(mat_a, axis=0) - np.median(mat_b, axis=0))
    n_features = mat_all.shape[1]
    exceed_count = np.zeros(n_features, dtype=np.int64)

    rng = np.random.default_rng(seed)
    block_indices = [np.where(blocks == b)[0] for b in np.unique(blocks)]

    perm_labels = labels.copy()
    for _ in range(n_perm):
        for idx in block_indices:
            perm_labels[idx] = rng.permutation(labels[idx])
        stat = np.abs(
            np.median(mat_all[~perm_labels], axis=0) - np.median(mat_all[perm_labels], axis=0)
        )
        exceed_count += stat >= observed

    return (exceed_count + 1) / (n_perm + 1)
