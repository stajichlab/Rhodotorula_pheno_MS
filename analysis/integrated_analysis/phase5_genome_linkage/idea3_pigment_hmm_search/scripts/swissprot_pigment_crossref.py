"""Idea 3: quick keyword cross-reference of BFD's per-protein SwissProt
diamond/blast hits against pigment-pathway-relevant SwissProt entries.

Independent of the custom HMM panel search -- uses only what's already in
BFD/tables/swissprot.parquet (per-protein best hits) joined against
swissprot_annot.parquet (SwissProt entry metadata: protein_name/gene_name)
and BFD.duckdb's species table (strain/species mapping). Not a replacement
for the HMM panel (SwissProt is a generic reference DB, not fungal-
pigment-specific), but a fast independent sanity check of what's already
in BFD.

Usage:
    python3 analysis/integrated_analysis/phase5_genome_linkage/idea3_pigment_hmm_search/scripts/swissprot_pigment_crossref.py
"""
from pathlib import Path

import duckdb
import pandas as pd

BFD_TABLES = Path("BFD/tables")
BFD_DUCKDB = Path("BFD/db/BFD.duckdb")
OUT_DIR = Path("analysis/integrated_analysis/phase5_genome_linkage/idea3_pigment_hmm_search/outputs")

# Keywords matched against SwissProt protein_name / gene_name (case-insensitive substring).
# Mirrors the pigment_protein_hmms panel families where a SwissProt analog is plausible.
KEYWORDS = [
    "carotenoid", "phytoene", "lycopene", "carotene desaturase",
    "crtb", "crti", "crte", "crtr", "crto", "crtp", "crtq", "crty",
    "tyrosinase", "laccase", "melanin", "dhn-melanin", "scytalone",
    "tetrahydroxynaphthalene", "1,3,6,8-tetrahydroxynaphthalene",
    "trihydroxynaphthalene reductase", "hydroxynaphthalene reductase",
    "homogentisate", "hydroxyphenylpyruvate dioxygenase", "hppd",
    "polyketide synthase", "ayg1", "ochre", "melanogenesis",
]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    annot = pd.read_parquet(BFD_TABLES / "swissprot_annot.parquet")
    pat = "|".join(KEYWORDS)
    hay = (annot["protein_name"].fillna("") + " " + annot["gene_name"].fillna("")).str.lower()
    mask = hay.str.contains(pat, regex=True, na=False)
    hit_annot = annot[mask].copy()
    print(f"SwissProt entries matching pigment keywords: {len(hit_annot)} / {len(annot)}")
    hit_annot.to_csv(OUT_DIR / "swissprot_pigment_keyword_entries.csv", index=False)

    hits = pd.read_parquet(BFD_TABLES / "swissprot.parquet")
    print(f"Total BFD protein->SwissProt hit rows: {len(hits)}")

    merged = hits.merge(
        hit_annot[["accession", "entry_name", "protein_name", "gene_name", "organism"]],
        left_on="swissprot_acc", right_on="accession", how="inner",
    )
    print(f"BFD protein hits landing on a pigment-keyword SwissProt entry: {len(merged)}")

    # keep best (lowest evalue) hit per query protein per matched SwissProt entry family
    merged = merged.sort_values("evalue").drop_duplicates(subset=["protein_id", "swissprot_acc"])

    con = duckdb.connect(str(BFD_DUCKDB), read_only=True)
    species = con.execute("SELECT LOCUSTAG, GENUS, SPECIES, STRAIN FROM species").df()
    con.close()
    merged["species_prefix_clean"] = merged["species_prefix"]
    merged = merged.merge(species, left_on="species_prefix_clean", right_on="LOCUSTAG", how="left")

    cols = ["protein_id", "species_prefix", "GENUS", "SPECIES", "STRAIN",
            "swissprot_acc", "entry_name", "protein_name", "gene_name",
            "pident", "qcovhsp", "evalue", "bitscore"]
    merged = merged[cols].sort_values(["gene_name", "GENUS", "SPECIES", "evalue"])
    merged.to_csv(OUT_DIR / "swissprot_pigment_hits_by_protein.csv", index=False)
    print(f"Wrote {OUT_DIR / 'swissprot_pigment_hits_by_protein.csv'} ({len(merged)} rows)")

    # summary: strain x matched-gene-family presence/count
    summary = (
        merged.groupby(["GENUS", "SPECIES", "STRAIN", "gene_name"])
        .size().reset_index(name="n_hits")
        .pivot_table(index=["GENUS", "SPECIES", "STRAIN"], columns="gene_name", values="n_hits", fill_value=0)
    )
    summary.to_csv(OUT_DIR / "swissprot_pigment_strain_summary.csv")
    print(f"Wrote {OUT_DIR / 'swissprot_pigment_strain_summary.csv'} ({summary.shape[0]} strains x {summary.shape[1]} gene families)")
    print("\nGene-family hit counts (proteins, across all strains):")
    print(merged["gene_name"].value_counts().to_string())


if __name__ == "__main__":
    main()
