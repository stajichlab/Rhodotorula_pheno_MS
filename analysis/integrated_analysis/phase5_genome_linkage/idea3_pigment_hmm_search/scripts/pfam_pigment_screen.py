"""Idea 3: coarse Pfam-domain screen for pigment-pathway-relevant families
in BFD.duckdb, independent of the PI's custom HMM panel.

Runs the pre-filter step originally scaffolded (never executed against
real data) in phase5_candidate_gene_genotyping.py's CANDIDATE_GENES dict,
plus a few additional pigment-relevant Pfam families (laccase multicopper
oxidase domains, tyrosinase) matching the HMM panel's melanin-pathway
coverage. This is a coarse presence screen (domain-level, not
ortholog-confirmed) -- a starting triage, not a substitute for the HMM
panel or a real ortholog-confirmation step (MSA + tree).

Usage:
    python3 analysis/integrated_analysis/phase5_genome_linkage/idea3_pigment_hmm_search/scripts/pfam_pigment_screen.py
"""
from pathlib import Path

import duckdb
import pandas as pd

BFD_DUCKDB = Path("BFD/db/BFD.duckdb")
OUT_DIR = Path("analysis/integrated_analysis/phase5_genome_linkage/idea3_pigment_hmm_search/outputs")

# pfam_acc in BFD carries a version suffix (e.g. PF00494.21) -- match by prefix.
CANDIDATE_PFAM = {
    "crtYB": ["PF00494"],                    # Squalene/phytoene synthase
    "crtI": ["PF01593"],                     # Phytoene dehydrogenase / FAD-dependent oxidoreductase
    "crtS/crtR-hydroxylase": ["PF00067"],    # Cytochrome P450
    "crtR-CBR-type": ["PF00970", "PF00258"], # Ferredoxin reductase / flavodoxin-like
    "HMGR": ["PF00368"],                     # HMG-CoA reductase
    "GGPPS": ["PF00348"],                    # Polyprenyl synthetase
    "laccase/multicopper-oxidase": ["PF00394", "PF07731", "PF07732"],
    "tyrosinase": ["PF00264"],
}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(BFD_DUCKDB), read_only=True)

    species = con.execute("SELECT LOCUSTAG, GENUS, SPECIES, STRAIN FROM species").df()

    rows = []
    for gene, pfam_ids in CANDIDATE_PFAM.items():
        like_clause = " OR ".join(f"pfam_acc LIKE '{p}.%'" for p in pfam_ids)
        df = con.execute(
            f"""
            SELECT species_prefix, protein_id, pfam_acc, full_seq_e_value, full_seq_score
            FROM pfam
            WHERE {like_clause}
            """
        ).fetchdf()
        df["gene"] = gene
        rows.append(df)
        n_strains = df["species_prefix"].nunique()
        print(f"{gene} ({','.join(pfam_ids)}): {len(df)} domain hits, {n_strains}/278 strains")

    hits = pd.concat(rows, ignore_index=True)
    hits = hits.merge(species, left_on="species_prefix", right_on="LOCUSTAG", how="left")
    hits = hits[["gene", "pfam_acc", "protein_id", "species_prefix", "GENUS", "SPECIES", "STRAIN",
                 "full_seq_e_value", "full_seq_score"]]
    hits.to_csv(OUT_DIR / "pfam_pigment_hits_by_protein.csv", index=False)
    print(f"\nWrote {OUT_DIR / 'pfam_pigment_hits_by_protein.csv'} ({len(hits)} rows)")

    # per-strain copy number matrix
    summary = (
        hits.groupby(["GENUS", "SPECIES", "STRAIN", "gene"])
        .size().reset_index(name="copy_number")
        .pivot_table(index=["GENUS", "SPECIES", "STRAIN"], columns="gene", values="copy_number", fill_value=0)
    )
    # ensure all 278 strains appear even with 0 hits everywhere
    summary = species.set_index(["GENUS", "SPECIES", "STRAIN"])[[]].join(summary, how="left").fillna(0).astype(int)
    summary.to_csv(OUT_DIR / "pfam_pigment_strain_summary.csv")
    print(f"Wrote {OUT_DIR / 'pfam_pigment_strain_summary.csv'} ({summary.shape[0]} strains x {summary.shape[1]} families)")

    print("\nPer-species mean copy number:")
    per_species = summary.reset_index().groupby(["GENUS", "SPECIES"])[list(CANDIDATE_PFAM.keys())].mean()
    print(per_species.round(2).to_string())

    con.close()


if __name__ == "__main__":
    main()
