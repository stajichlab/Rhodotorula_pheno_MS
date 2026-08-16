"""Idea 3: parse hmmsearch --domtblout output from the PI's 28-profile
pigment_protein_hmms panel into a per-protein best-hit table and a
per-strain family presence/copy-number matrix.

Usage:
    python3 analysis/integrated_analysis/phase5_genome_linkage/idea3_pigment_hmm_search/scripts/parse_hmm_hits.py
"""
from pathlib import Path

import duckdb
import pandas as pd

OUT_DIR = Path("analysis/integrated_analysis/phase5_genome_linkage/idea3_pigment_hmm_search/outputs")
DOMTBLOUT = OUT_DIR / "pigment_hmm_hits.domtblout"
BFD_DUCKDB = Path("BFD/db/BFD.duckdb")

COLS = ["target_name", "target_acc", "tlen", "query_name", "query_acc", "qlen",
        "full_evalue", "full_score", "full_bias", "dom_num", "dom_of",
        "c_evalue", "i_evalue", "dom_score", "dom_bias",
        "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"]


def load_domtblout(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split(None, len(COLS) - 1)
            rows.append(parts[:len(COLS)])
    df = pd.DataFrame(rows, columns=COLS)
    for c in ["tlen", "qlen", "dom_num", "dom_of", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to"]:
        df[c] = df[c].astype(int)
    for c in ["full_evalue", "full_score", "full_bias", "c_evalue", "i_evalue", "dom_score", "dom_bias"]:
        df[c] = df[c].astype(float)
    return df


def main():
    df = load_domtblout(DOMTBLOUT)
    print(f"Parsed {len(df)} domain hit rows, {df['query_name'].nunique()} distinct HMM queries hit, "
          f"{df['target_name'].nunique()} distinct proteins hit")

    # best domain per (protein, HMM family) -- lowest i-evalue
    best = df.sort_values("i_evalue").drop_duplicates(subset=["target_name", "query_name"])
    best["species_prefix"] = best["target_name"].str.split("_").str[0]

    con = duckdb.connect(str(BFD_DUCKDB), read_only=True)
    species = con.execute("SELECT LOCUSTAG, GENUS, SPECIES, STRAIN FROM species").df()
    con.close()
    best = best.merge(species, left_on="species_prefix", right_on="LOCUSTAG", how="left")
    unmatched = best["LOCUSTAG"].isna().sum()
    if unmatched:
        print(f"WARNING: {unmatched} hit rows did not match a species_prefix -> strain mapping "
              f"(check protein-ID naming convention for those strains)")

    out_cols = ["query_name", "target_name", "species_prefix", "GENUS", "SPECIES", "STRAIN",
                "full_evalue", "full_score", "i_evalue", "dom_score", "ali_from", "ali_to", "qlen", "tlen"]
    best_out = best[out_cols].sort_values(["query_name", "GENUS", "SPECIES", "i_evalue"])
    best_out.to_csv(OUT_DIR / "pigment_hmm_hits_by_protein.csv", index=False)
    print(f"Wrote {OUT_DIR / 'pigment_hmm_hits_by_protein.csv'} ({len(best_out)} rows)")

    summary = (
        best.groupby(["GENUS", "SPECIES", "STRAIN", "query_name"])
        .size().reset_index(name="copy_number")
        .pivot_table(index=["GENUS", "SPECIES", "STRAIN"], columns="query_name", values="copy_number", fill_value=0)
    )
    summary = species.set_index(["GENUS", "SPECIES", "STRAIN"])[[]].join(summary, how="left").fillna(0).astype(int)
    summary.to_csv(OUT_DIR / "pigment_hmm_strain_summary.csv")
    print(f"Wrote {OUT_DIR / 'pigment_hmm_strain_summary.csv'} ({summary.shape[0]} strains x {summary.shape[1]} families)")

    print("\nHits per HMM family (proteins, all strains, i-Evalue<1e-5):")
    print(best["query_name"].value_counts().to_string())

    print("\nStrains with >=1 hit per family (out of 278):")
    presence = (summary > 0).sum(axis=0)
    print(presence.sort_values(ascending=False).to_string())

    print("\nPer-species mean copy number:")
    per_species = summary.reset_index().groupby(["GENUS", "SPECIES"])[summary.columns.tolist()].mean()
    print(per_species.round(2).to_string())


if __name__ == "__main__":
    main()
