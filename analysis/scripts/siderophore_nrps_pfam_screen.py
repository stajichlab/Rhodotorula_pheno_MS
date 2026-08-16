#!/usr/bin/env python3
"""
PI request (2026-08-16): coarse genome-level screen for the NRPS/
ornithine-N5-monooxygenase machinery needed to make rhodotorulic acid,
as a stand-in until a specific characterized reference sequence for the
actual rhodotorulic-acid-synthesizing NRPS is available (PI: "we know the
NRPS gene that makes Rhodotorulic acid" -- once given an accession/
sequence, this should be superseded by a real ortholog search, same as
Idea 3's custom HMM panel superseded its own initial Pfam pre-filter).

Rhodotorulic acid biosynthesis (from characterized homologs in other
basidiomycetes, e.g. Ustilago maydis Sid1/Sid2) requires, in order:
  1. L-ornithine N5-monooxygenase (hydroxylase) -- a flavin-dependent
     monooxygenase, Pfam PF00743 "FMO-like".
  2. An NRPS that activates and condenses 2 units of the resulting
     N5-acetyl-N5-hydroxyornithine into the cyclic dihydroxamate
     diketopiperazine -- detected via the co-occurrence, on the SAME
     protein, of adenylation (PF00501, AMP-binding), condensation
     (PF00668), and phosphopantetheine-attachment/PCP (PF00550) domains.
     Rhodotorulic acid is a DIMER of a single amino acid building block,
     so its synthetase is expected to be a small, ~2-module NRPS -- i.e.
     a protein with ~2 adenylation domains, in contrast to larger
     multi-module NRPS (e.g. the 6-module ferrichrome synthetase, or
     unrelated secondary-metabolite NRPS elsewhere in the genome). This
     script reports the adenylation-domain-count distribution per
     candidate NRPS protein specifically to distinguish these cases.

CAVEAT (explicit, not a footnote): this is a coarse, generic domain-level
screen. PF00501/PF00668/PF00550 are large superfamilies (any NRPS or even
some unrelated multi-domain enzymes can carry AMP-binding-like domains);
PF00743 covers FMO-like monooxygenases broadly, not specifically ornithine
hydroxylases. This CANNOT confirm true orthology -- it can only narrow
the search space (candidate proteins/strains) for the real ortholog
search once a reference sequence is available. Do not treat "N adenylation
domains, PF00743 also present" as proof this is the rhodotorulic acid
synthetase.

Usage:
    python3 analysis/scripts/siderophore_nrps_pfam_screen.py
"""
from __future__ import annotations

from pathlib import Path

import duckdb
import pandas as pd

BFD_DUCKDB = Path("BFD/db/BFD.duckdb")
OUT_DIR = Path("analysis/integrated_analysis/phase_siderophore")

NRPS_DOMAIN_PFAM = {
    "adenylation": "PF00501",
    "condensation": "PF00668",
    "PCP_thiolation": "PF00550",
}
HYDROXYLASE_PFAM = "PF00743"  # FMO-like


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(BFD_DUCKDB), read_only=True)
    species = con.execute("SELECT LOCUSTAG, GENUS, SPECIES, STRAIN FROM species").df()

    # ornithine-hydroxylase-like presence per strain
    hydrox = con.execute(
        f"SELECT species_prefix, protein_id, pfam_acc, full_seq_e_value FROM pfam WHERE pfam_acc LIKE '{HYDROXYLASE_PFAM}.%'"
    ).fetchdf()
    hydrox = hydrox.merge(species, left_on="species_prefix", right_on="LOCUSTAG", how="left")
    hydrox.to_csv(OUT_DIR / "siderophore_hydroxylase_hits_by_protein.csv", index=False)
    print(f"Ornithine-hydroxylase-like ({HYDROXYLASE_PFAM}) hits: {len(hydrox)} proteins, "
          f"{hydrox['species_prefix'].nunique()}/278 strains")

    # NRPS domain co-occurrence: pull all 3 domain types, find proteins carrying all 3
    all_doms = []
    for name, acc in NRPS_DOMAIN_PFAM.items():
        d = con.execute(
            f"SELECT species_prefix, protein_id, pfam_acc, hmm_from, hmm_to, ali_from, ali_to, domain_score "
            f"FROM pfam WHERE pfam_acc LIKE '{acc}.%'"
        ).fetchdf()
        d["domain_type"] = name
        all_doms.append(d)
        print(f"{name} ({acc}): {len(d)} domain hits, {d['species_prefix'].nunique()}/278 strains")
    doms = pd.concat(all_doms, ignore_index=True)
    doms.to_csv(OUT_DIR / "siderophore_nrps_domain_hits.csv", index=False)

    # candidate NRPS proteins: carry adenylation + condensation + PCP domains together
    domain_types_per_protein = doms.groupby("protein_id")["domain_type"].apply(set)
    candidate_proteins = domain_types_per_protein[domain_types_per_protein.apply(lambda s: len(s) == 3)].index
    print(f"\nCandidate NRPS proteins (all 3 domain types present): {len(candidate_proteins)}")

    a_domain_counts = (
        doms[(doms["protein_id"].isin(candidate_proteins)) & (doms["domain_type"] == "adenylation")]
        .groupby("protein_id").size()
    )
    print("\nAdenylation (A-)domain count per candidate NRPS protein (module number proxy):")
    print(a_domain_counts.value_counts().sort_index().to_string())

    two_module = a_domain_counts[a_domain_counts == 2].index
    print(f"\nCandidate 2-module NRPS proteins (best structural match to rhodotorulic acid's "
          f"dimeric architecture): {len(two_module)}")

    cand = doms[doms["protein_id"].isin(candidate_proteins)][["species_prefix", "protein_id"]].drop_duplicates()
    cand["n_adenylation_domains"] = cand["protein_id"].map(a_domain_counts)
    cand["is_two_module"] = cand["protein_id"].isin(two_module)
    cand = cand.merge(species, left_on="species_prefix", right_on="LOCUSTAG", how="left")
    cand.to_csv(OUT_DIR / "siderophore_candidate_nrps_proteins.csv", index=False)

    # per-strain summary: has hydroxylase AND has >=1 candidate 2-module NRPS
    strains_with_hydrox = set(hydrox["species_prefix"])
    strains_with_2mod_nrps = set(cand.loc[cand["is_two_module"], "species_prefix"])
    summary = species.copy()
    summary["has_hydroxylase"] = summary["LOCUSTAG"].isin(strains_with_hydrox)
    summary["has_2module_nrps_candidate"] = summary["LOCUSTAG"].isin(strains_with_2mod_nrps)
    summary["has_both"] = summary["has_hydroxylase"] & summary["has_2module_nrps_candidate"]
    summary.to_csv(OUT_DIR / "siderophore_pathway_strain_summary.csv", index=False)
    print(f"\nStrains with hydroxylase-like hit: {summary['has_hydroxylase'].sum()}/278")
    print(f"Strains with a candidate 2-module NRPS: {summary['has_2module_nrps_candidate'].sum()}/278")
    print(f"Strains with BOTH: {summary['has_both'].sum()}/278")
    print(f"\nWrote outputs to {OUT_DIR}")

    con.close()


if __name__ == "__main__":
    main()
