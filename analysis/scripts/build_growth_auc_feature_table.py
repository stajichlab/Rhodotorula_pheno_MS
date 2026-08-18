#!/usr/bin/env python3
"""
Grow-rate (Cu-AUC) <-> supernatant metabolome: focused table for the 10
within-species features.

Recomputes the analysis in
.living/findings/abundance-axis-growth-rate-relationships.md (F-002/F-003):
- whole-cell-extract ('cell') and 'supernatant' fractions analyzed separately
- TSS normalization per sample, replicates collapsed by canonical_strain
- adduct/isotopologue de-duplication via ms_feature_dedup_groups.csv
- Spearman rank correlation of each feature's relative abundance vs
  mean_auc_rate (liquid-growth AUC in Copper media), restricted to
  R. mucilaginosa (the only species with adequate n)
- keeps the 10 supernatant features with |rho| > 0.3 (significantly above
  permutation null; cell fraction has 0 -- entirely between-species)

Outputs (focusing on what can be investigated):
  outputs/within_mucilaginosa_auc_features.tsv.gz  - all features w/ rho vs AUC
  outputs/compound_summary.tsv   - the 10 features annotated (same identity
                                   pipeline + column schema as the pairwise
                                   differential comparisons)
  outputs/compound_summary.html  - sortable/filterable table (same renderer as
                                   the differential comparisons)
"""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import build_compound_summary as bcs  # noqa: E402
import generate_compound_table_html as gcth  # noqa: E402

META = REPO / "analysis/linked_data/sample_metadata.csv.gz"
FEAT = REPO / "analysis/linked_data/feature_abundance_matrix.csv.gz"
DEDUP = REPO / "analysis/linked_data/ms_feature_dedup_groups.csv"
OUT = REPO / "analysis/growth_rate_auc/outputs"

RHO_THRESH = 0.3
N_PERM = 1000


def build(frac: str, restrict: set | None = None):
    meta = pd.read_csv(META)
    feat = pd.read_csv(FEAT)
    dedup = pd.read_csv(DEDUP)
    rep_ids = set(dedup[dedup["is_group_representative"]]["row ID"])
    feat = feat[feat["row ID"].isin(rep_ids)]

    m = meta[meta["fraction"] == frac]
    cols = [c for c in feat.columns if c in set(m["sample_id"])]
    mat = feat[cols].to_numpy(dtype=float)
    cs = mat.sum(axis=0)
    cs[cs == 0] = 1.0
    mat = mat / cs
    smap = m.set_index("sample_id")["canonical_strain"]
    strains = [smap[c] for c in cols]
    sdf = pd.DataFrame(mat.T, index=strains).groupby(level=0).mean()
    span = sdf.max() - sdf.min()
    sdf = sdf.loc[:, span > 0]
    fids = feat["row ID"].to_numpy()[np.where(span > 0)[0]].astype(int)
    dic = m.drop_duplicates("canonical_strain").set_index("canonical_strain")
    auc = dic["mean_auc_rate"]
    sprc = dic["Species"]
    if restrict is not None:
        sdf = sdf[sdf.index.isin(restrict)]
    common = [s for s in sdf.index.intersection(auc.index) if pd.notna(auc[s])]
    return sdf, fids, auc.loc[common], sprc, meta


def rho_signed(sdf, y):
    """Vectorized Spearman rho per feature vs phenotype on the shared strains."""
    Xr = sdf.loc[sdf.index.intersection(y.index)].rank(axis=0)
    yr = y.rank()
    yc = (yr - yr.mean()) / np.sqrt(((yr - yr.mean()) ** 2).sum())
    Xc = Xr - Xr.mean(axis=0)
    Xc = Xc / np.sqrt((Xc**2).sum(axis=0) + 1e-300)
    return (Xc.T @ yc).to_numpy(), Xr.index


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)

    # --- within-mucilaginosa: compute rho for every feature in each fraction ---
    spfull = pd.read_csv(META).drop_duplicates("canonical_strain").set_index("canonical_strain")["Species"]
    muc = spfull[spfull == "Rhodotorula mucilaginosa"].index

    for frac in ["cell", "supernatant"]:
        sdf, fids, auc, sprc, _ = build(frac, restrict=muc)
        rho, st = rho_signed(sdf, auc)
        rho = pd.Series(rho, index=fids).rename("rho_auc")
        # within-muc permutation null (presserve no structure: single species)
        hit = np.abs(rho.to_numpy()) > RHO_THRESH
        n_perm = N_PERM
        cnt = []
        for _ in range(n_perm):
            yp = pd.Series(rng.permutation(auc.loc[st].to_numpy()), index=st)
            a, _ = rho_signed(sdf.loc[st], yp)
            cnt.append((np.abs(a) > RHO_THRESH).sum())
        cnt = np.array(cnt)
        obs = int(hit.sum())
        emp_p = float((cnt >= obs).mean())
        out = pd.DataFrame({"row ID": fids, "rho_auc": rho.to_numpy()})
        gz = gzip.open(OUT / f"within_mucilaginosa_{frac}_auc_features.tsv.gz", "wt")
        out.to_csv(gz, sep="\t", index=False)
        gz.close()
        print(f"{frac}: n features={len(fids)}, |rho|>{RHO_THRESH} obs={obs}, "
              f"perm null mean={cnt.mean():.2f} sd={cnt.std():.2f}, emp_p={emp_p:.4f}", file=sys.stderr)

    # --- the 10 features: annotation + focused outputs ---
    sdf, fids, auc, sprc, meta = build("supernatant", restrict=muc)
    rho, st = rho_signed(sdf, auc)
    rho = pd.Series(rho, index=fids).rename("rho_auc")
    hits = rho[np.abs(rho) > RHO_THRESH].sort_values()
    print(f"\n{len(hits)} feature hits", file=sys.stderr)

    # empirical p per hit via the shared within-muc null (single null, same count)
    # -> per-feature p via permutation of AUC (free swap within muc)
    n_perm = N_PERM
    rng = np.random.default_rng(7)
    per_feat = np.zeros(len(hits))
    Xr = sdf.loc[st].rank(axis=0).to_numpy()
    Xc = Xr - Xr.mean(axis=0)
    Xc = Xc / np.sqrt((Xc**2).sum(axis=0) + 1e-300)
    hit_col = [np.where(fids == r)[0][0] for r in hits.index]
    for _ in range(n_perm):
        yr = pd.Series(rng.permutation(auc.loc[st].to_numpy()), index=st).rank()
        yc = (yr - yr.mean()) / np.sqrt(((yr - yr.mean()) ** 2).sum())
        a = np.abs(Xc.T @ yc)
        per_feat += (a[hit_col] >= np.abs(hits.to_numpy())).astype(float)
    emp_p_per = (per_feat + 1) / (n_perm + 1)

    rowmeta = meta[meta["fraction"] == "supernatant"].drop_duplicates("canonical_strain")
    an = pd.DataFrame({
        "row ID": hits.index,
        "rho_auc": hits.round(4).to_numpy(),
        "emp_p_perm_within_muc": np.round(emp_p_per, 6),
    })
    print(an.to_string(index=False), file=sys.stderr)

    # --- annotate with the same identity pipeline as differential comparisons ---
    library = bcs.load_gnps(bcs.LIBRARY_SEARCH, "library")
    analog = bcs.load_gnps(bcs.ANALOG_SEARCH, "analog")
    sirius = bcs.load_sirius()
    sig = an.merge(
        pd.read_csv(FEAT)[["row ID", "row m/z", "row retention time", "adduct", "has_ms2"]],
        on="row ID", how="left",
    )
    sig = bcs.annotate(sig, library, analog, sirius)

    # align with compound_summary column schema used by the differential comparisons
    cols = [
        "row ID", "row m/z", "row retention time", "adduct", "has_ms2",
        "rho_auc", "emp_p_perm_within_muc",
        "best_identity", "best_identity_source",
        "library_NAME", "library_cosine", "library_matched_peaks", "library_SMILES", "library_INCHI", "library_ORGANISM",
        "analog_NAME", "analog_cosine", "analog_matched_peaks", "analog_SMILES", "analog_INCHI", "analog_ORGANISM",
        "sirius_formula", "sirius_adduct", "sirius_structure_name", "sirius_structure_smiles",
        "sirius_structure_confidence", "sirius_npc_pathway", "sirius_npc_class", "sirius_classyfire_class",
    ]
    sig = sig[[c for c in cols if c in sig.columns]].sort_values("rho_auc", ascending=False)
    tsv_out = OUT / "compound_summary.tsv"
    sig.to_csv(tsv_out, sep="\t", index=False)
    print(f"wrote {tsv_out} ({len(sig)} rows)", file=sys.stderr)

    # --- focused HTML (same renderer; primary cols swapped for rho/emp-p) ---
    gcth.PRIMARY_COLS = [
        ("best_identity_source", "ID", "glyph"),
        ("row ID", "row ID", "int"),
        ("row m/z", "m/z", "float"),
        ("row retention time", "RT (min)", "float"),
        ("adduct", "adduct", "text"),
        ("rho_auc", "rho (vs AUC)", "float"),
        ("emp_p_perm_within_muc", "emp. perm p", "sci"),
        ("best_identity", "best identity", "text"),
    ]
    gcth.SECONDARY_LABELS = {**gcth.SECONDARY_LABELS,
        "rho_auc": "Within-species Spearman rho (R. mucilaginosa, Cu AUC)",
        "emp_p_perm_within_muc": "Empirical permutation p (within R. mucilaginosa, 1000 perms)",
    }

    html_out = OUT / "compound_summary.html"
    df_html = pd.read_csv(tsv_out, sep="\t")
    title = "Growth-rate (Cu-AUC) associated features — R. mucilaginosa supernatant (n=10)"
    html_out.write_text(gcth.build_html(df_html, title))
    print(f"wrote {html_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
