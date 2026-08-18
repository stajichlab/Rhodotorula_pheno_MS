#!/usr/bin/env python3
"""Test network (GNPS/MS2 molecular-family) components for phenotype association.

Groups features by their connected component in the MS2 molecular network
(``filtered_pairs.tsv`` ComponentIndex) instead of testing ~10.9 k features
independently.  Two component-level designs are computed and compared:

  1) GSEA-style enrichment:
       per-feature permutation-derived p (normal approx. on Spearman rank rho)
       -> component score = mean(-log10 p) over members.
     Interpretation: are a component's members collectively more associated
     with the phenotype than the panel-wide null?

  2) max|rho|:
       component statistic = max |Spearman rho| over members; null calibrated
       by permuting the phenotype within preserved component structure.
     Interpretation: does the best member of a molecular family exceed the
     expectation for the largest association in any component?

Both are reported per component with permutation p-values; max|rho| also gets a
component-level FDR estimated from the whole-panel max null.  Components are
retained only if they have >= MIN_MEMBERS dedup-representative members in the
analysis set; network singletons are reported in the per-feature rho file.

Targets (within R. mucilaginosa):
  * growth AUC (``mean_auc_rate``) per fraction (cell / supernatant)
  * color ``a*`` and colourfulness ``C*`` per fraction

Outputs: analysis/network_components/outputs/<trait>_{components,rhoperm}.tsv
and a comparison summary table.

USAGE:  analysis/scripts/network_component_association.py [--nperm N]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import build_compound_summary as bcs  # noqa: E402

EB = bcs.EB_DIR
META = REPO / "analysis/linked_data/sample_metadata.csv.gz"
FEAT = REPO / "analysis/linked_data/feature_abundance_matrix.csv.gz"
DEDUP = REPO / "analysis/linked_data/ms_feature_dedup_groups.csv"
NET = EB / "networking" / "filtered_pairs.tsv"
OUT = REPO / "analysis" / "network_components" / "outputs"

MUC = "Rhodotorula mucilaginosa"
MIN_MEMBERS = 2

# trait -> (phenotype column, fraction, within-muc only)
TRAITS = {
    "growth_cell": ("mean_auc_rate", "cell", True),
    "growth_supernatant": ("mean_auc_rate", "supernatant", True),
    "color_a_cell": ("Mean_ColorLab_a*Mean", "cell", True),
    "color_a_supernatant": ("Mean_ColorLab_a*Mean", "supernatant", True),
    "color_C_cell": ("chroma", "cell", True),
    "color_C_supernatant": ("chroma", "supernatant", True),
}


def zrank(s: pd.Series) -> np.ndarray:
    r = s.rank().to_numpy(dtype=float)
    r = r - r.mean()
    n = np.sqrt((r**2).sum())
    return r / n if n > 0 else r


def load_network_components() -> pd.DataFrame:
    ed = pd.read_csv(NET, sep="\t", usecols=["CLUSTERID1", "CLUSTERID2", "ComponentIndex"])
    lo = pd.DataFrame({"row ID": ed.CLUSTERID1, "ComponentIndex": ed.ComponentIndex})
    hi = pd.DataFrame({"row ID": ed.CLUSTERID2, "ComponentIndex": ed.ComponentIndex})
    m = pd.concat([lo, hi], ignore_index=True).drop_duplicates("row ID")
    return m


def collate_per_component(
    fids: np.ndarray, fc: np.ndarray, rho: np.ndarray, nlp: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sort by component; return component ids (sorted), sizes, start + end
    indices into the *sorted* arrays, and the stable sort order.

    Members' values live at ``rho[order[s:e]]`` etc.; use ``order`` on the
    original arrays so observed and permutation data are sliced identically."""
    order = np.argsort(fc, kind="stable")
    fc_s = fc[order]
    edges = np.where(fc_s[1:] != fc_s[:-1])[0] + 1
    starts = np.concatenate([[0], edges])
    ends = np.concatenate([edges, [len(fc_s)]])
    cids = fc_s[starts]
    sizes = ends - starts
    return cids, sizes, starts, ends, order


def run_trait(
    tname: str,
    pheno_col: str,
    frac: str,
    in_muc: bool,
    feat_row: pd.DataFrame,
    meta: pd.DataFrame,
    comp_df: pd.DataFrame,
    fids_analysis: np.ndarray,
    comp_of: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    m = meta[meta["fraction"] == frac].copy()
    if in_muc:
        m = m[m["Species"] == MUC]
    if pheno_col == "chroma":
        m["pheno"] = np.sqrt(m["Mean_ColorLab_a*Mean"] ** 2 + m["Mean_ColorLab_b*Mean"] ** 2)
    else:
        m["pheno"] = m[pheno_col]
    # strain-level phenotype
    phen = m.groupby("canonical_strain")["pheno"].mean().dropna()
    if len(phen) < 5:
        print(f"[{tname}] skipping: only {len(phen)} strains", flush=True)
        return None
    samp2strain = m.set_index("sample_id")["canonical_strain"]
    cols = [c for c in feat_row.columns if c in samp2strain.index]
    mat = feat_row[cols].to_numpy(dtype=float)
    cs = mat.sum(axis=0)
    cs[cs == 0] = 1.0
    mat = mat / cs
    fid_col = feat_row["row ID"].to_numpy()
    sdf = pd.DataFrame(mat.T, index=samp2strain.loc[cols].to_numpy(),
                       columns=fid_col).groupby(level=0).mean()
    common = sdf.index.intersection(phen.index)
    X = sdf.loc[common].rank(axis=0)
    y = phen.loc[common]
    Xc = X - X.mean(axis=0)
    Xc = Xc / np.sqrt((Xc**2).sum(axis=0) + 1e-300)
    T = Xc.to_numpy()            # n_strains x n_feat
    yc = zrank(y)                # n_strains
    n = len(yc)
    fcol = X.columns.to_numpy()
    # keep only analysis features actually present in this matrix
    fids = np.array([f for f in fids_analysis if f in X.columns])
    fidx = np.array([np.where(fcol == f)[0][0] for f in fids])
    T = T[:, fidx]
    fc = comp_of[np.searchsorted(fids_analysis, fids)]
    mapped = fc > 0  # in a connected network component (size may be 1)

    # --- per-feature stats over ALL features in this trait matrix -----------
    rho_full = T.T @ yc
    with np.errstate(divide="ignore"):
        t_ = rho_full * np.sqrt((n - 2) / np.maximum(1 - rho_full**2, 1e-12))
        p_full = 2 * stats.t.sf(np.abs(t_), n - 2)
        nlp_full = -np.log10(p_full + 1e-300)
    pf = pd.DataFrame({
        "row ID": fids,
        "ComponentIndex": fc,          # 0 => no network mapping
        "rho_abs": np.abs(rho_full),
        "es_member": nlp_full,
    })

    # --- observed (component testing uses only network-mapped features) -----
    Tm = T[:, mapped]
    fids_m = fids[mapped]
    fc_m = fc[mapped]
    rho_obs = Tm.T @ yc
    with np.errstate(divide="ignore"):
        t_ = rho_obs * np.sqrt((n - 2) / np.maximum(1 - rho_obs**2, 1e-12))
        p_obs = 2 * stats.t.sf(np.abs(t_), n - 2)
        nlp_obs = -np.log10(p_obs + 1e-300)

    cids, sizes, starts, ends, order = collate_per_component(fids_m, fc_m, rho_obs, nlp_obs)
    obs_maxrho = np.array([np.abs(rho_obs[order[s:e]]).max() for s, e in zip(starts, ends)])
    obs_es = np.array([nlp_obs[order[s:e]].mean() for s, e in zip(starts, ends)])

    # --- permutation null ---------------------------------------------------
    null_maxrho = np.empty((n_perm, len(cids)))
    null_es = np.empty((n_perm, len(cids)))
    perm_top = np.empty(n_perm)  # whole-panel max of per-comp max|rho|
    pn = np.empty((n_perm, len(fids_m)))
    pfn = np.empty((n_perm, len(fids)))  # full per-feature |rho| for perm p
    for i in range(n_perm):
        yp = zrank(pd.Series(rng.permutation(yc), index=y.index))
        rp = Tm.T @ yp
        rp_sorted = rp[order]
        perm_top[i] = np.max(np.abs(rp))
        pn[i] = np.abs(rp)
        pfn[i] = np.abs(T.T @ yp)
        nlp_p = -np.log10(2 * stats.t.sf(np.abs(rp_sorted) * np.sqrt(
            (n - 2) / np.maximum(1 - rp_sorted**2, 1e-12)), n - 2) + 1e-300)
        for j, (s, e) in enumerate(zip(starts, ends)):
            g = rp_sorted[s:e]
            null_maxrho[i, j] = g[np.abs(g).argmax()]
            null_es[i, j] = nlp_p[s:e].mean()

    keep = sizes >= MIN_MEMBERS
    rdf = pd.DataFrame({
        "ComponentIndex": cids,
        "n_members": sizes,
        "max_rho_abs": obs_maxrho,
        "es": obs_es,
    })
    rdf = rdf[keep].copy()
    if len(rdf) == 0:
        print(f"[{tname}] no components with >= {MIN_MEMBERS} members", flush=True)
        return None
    rdf["p_maxrho_perm"] = [
        min(1.0, float((null_maxrho[:, j] >= obs_maxrho[j]).mean()) + 1e-6)
        for j in rdf.index
    ]
    rdf["p_es_perm"] = [
        min(1.0, float((null_es[:, j] >= obs_es[j]).mean()) + 1e-6)
        for j in rdf.index
    ]
    # component-level FDR for max|rho| via whole-panel max null
    rdf["p_fdr_max"] = [min(1.0, float((perm_top >= obs_m).mean()) + 1e-6) for obs_m in rdf["max_rho_abs"]]
    rdf["p_es_fdr"] = stats.false_discovery_control(rdf["p_es_perm"].to_numpy())
    rdf["p_maxrho_fdr"] = stats.false_discovery_control(rdf["p_maxrho_perm"].to_numpy())

    # per-feature permutation p / BH-FDR (all features incl. singletons)
    pfn_abs = np.abs(pfn)
    pf["p_feature_perm"] = [
        min(1.0, float((pfn_abs[:, i] >= v).mean()) + 1e-6)
        for i, v in enumerate(pf["rho_abs"])
    ]
    pf["p_feature_fdr"] = stats.false_discovery_control(pf["p_feature_perm"].to_numpy())
    return rdf, pf


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nperm", type=int, default=1000)
    args = ap.parse_args()
    rng = np.random.default_rng(42)

    print("loading data ...", flush=True)
    meta = pd.read_csv(META)
    feat = pd.read_csv(FEAT)
    dedup = pd.read_csv(DEDUP)
    reps = dedup[dedup["is_group_representative"]]
    rep_ids = set(reps["row ID"])
    comp_df = load_network_components()
    comp_of = np.zeros(feat["row ID"].max() + 1, dtype=np.int64)
    comp_of[comp_df["row ID"].to_numpy()] = comp_df["ComponentIndex"].to_numpy()
    fids_analysis = np.array(sorted(rep_ids), dtype=np.int64)
    n_net = (comp_of[fids_analysis] > 0).sum()
    print(f"  features {len(feat)}   analysis reps {len(fids_analysis)}   "
          f"in-network reps {n_net}", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    summary = []
    for tname, (pheno_col, frac, in_muc) in TRAITS.items():
        print(f"[{tname}] ...", flush=True)
        res = run_trait(tname, pheno_col, frac, in_muc, feat, meta, comp_df,
                        fids_analysis, comp_of, args.nperm, rng)
        if res is None:
            continue
        rdf, pf = res
        rdf.to_csv(OUT / f"{tname}_components.tsv", sep="\t", index=False)
        pf.to_csv(OUT / f"{tname}_rhoperm.tsv", sep="\t", index=False)
        sub = rdf.sort_values("p_es_fdr").head(8)
        summary.append(f"\n[{tname}] top components (by FDR-adjusted enrichment):")
        for _, r in sub.iterrows():
            summary.append(
                f"  comp {int(r.ComponentIndex):>5}  n={int(r.n_members):>3}  "
                f"max|rho|={r.max_rho_abs:.3f} (p={r.p_maxrho_perm:.3f}, "
                f"BH={r.p_maxrho_fdr:.3f})  es={r.es:.2f} "
                f"(p={r.p_es_perm:.3f}, BH={r.p_es_fdr:.3f})")
        n_sig = int((rdf["p_es_fdr"] < 0.05).sum())
        n_sig_m = int((rdf["p_maxrho_fdr"] < 0.05).sum())
        summary.append(
            f"  {len(rdf)} components >= {MIN_MEMBERS} members; "
            f"enrichment FDR<0.05: {n_sig};  max|rho| FDR<0.05: {n_sig_m}")

    txt = "\n".join(summary)
    print(txt, flush=True)
    (OUT / "README.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
