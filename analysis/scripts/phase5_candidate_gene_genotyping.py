#!/usr/bin/env python3
"""
Idea 3 (analysis/ideas/2026-08-15-color-metabolome-genome-null-brainstorm/
DEVELOPMENT_PLAN.md, Part B): candidate carotenoid-pathway gene genotyping
against color, across all BFD strains. The first genome<->color test in
this project (strategy doc's Phase 5, "Genome Linkage" -- moved forward
per PI request 2026-08-15, ahead of Phase 3's metabolome-side enrichment).

Two tiers, per DEVELOPMENT_PLAN.md Part A (does this need SNP calling?
No -- see that doc for the full reasoning):
  Tier 1a: gene copy-number / presence-absence from existing BFD Pfam
           annotations, tested against color via the existing PGLS
           infrastructure. No new bioinformatics beyond ortholog curation.
  Tier 1b: loss-of-function / catalytic-residue variant calls from a
           multiple-sequence alignment of each candidate gene's ALREADY-
           ANNOTATED protein sequences across strains (assemblies are
           already called; no raw-read remapping or variant caller needed).

BLOCKING DEPENDENCY: BFD/db/BFD.duckdb is currently being rebuilt by a
separate nextflow pipeline (PI, 2026-08-15) -- every DB-dependent step
below is marked "BLOCKED: needs BFD.duckdb" and will raise/exit rather
than silently querying a stale or partial copy. Do not point this script
at any other duckdb path without explicit confirmation the rebuild is
finished and that copy is the intended one.

Candidate genes (priority order, see DEVELOPMENT_PLAN.md Part B Step 1
for the full table and rationale):
  crtYB  bifunctional phytoene synthase/lycopene cyclase   Pfam PF00494 (coarse)
  crtI   phytoene desaturase                                Pfam PF01593 (coarse)
  crtS   torulene -> torularhodin cytochrome P450            Pfam PF00067 (coarse)
  crtR   crtS's cognate P450 reductase                       Pfam PF00970+PF00258 (coarse)
  HMGR   HMG-CoA reductase (precursor supply)                 Pfam PF00368
  GGPPS  geranylgeranyl PP synthase (precursor supply)        Pfam PF00348

These Pfam IDs are deliberately coarse (documented false-positive risk --
e.g. PF00494 also hits unrelated squalene synthases, PF01593 hits many
unrelated FAD-oxidoreductases) and MUST be filtered by reciprocal-best-hit
against curated reference proteins (Step 2) before being trusted as
"this strain has crtYB" -- Pfam domain match alone is not ortholog
confirmation.

STATUS (2026-08-15): reference sequences/catalytic-residue positions for
the candidate genes have NOT yet been sourced (literature/PI task, tracked
in DEVELOPMENT_PLAN.md Part D). REFERENCE_FASTA below is a placeholder --
fill in before running Step 2.

Tools required (module system, checked available 2026-08-15):
    module load diamond/2.1.24
    module load mafft/7.505

Usage (once BFD.duckdb and reference sequences are both ready):
    python3 analysis/scripts/phase5_candidate_gene_genotyping.py --step ortholog
    python3 analysis/scripts/phase5_candidate_gene_genotyping.py --step copy-number
    python3 analysis/scripts/phase5_candidate_gene_genotyping.py --step msa-lof
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
BFD_DUCKDB = REPO / "BFD" / "db" / "BFD.duckdb"
BFD_PEP_DIR = REPO / "BFD" / "input_all" / "pep"  # per-strain *.proteins.fa, one per strain

OUT_DIR = REPO / "analysis" / "integrated_analysis" / "phase5_genome_linkage"
REFERENCE_FASTA = REPO / "analysis" / "scripts" / "data" / "carotenoid_pathway_references.fasta"  # NOT YET POPULATED

# Coarse Pfam pre-filter per candidate gene -- see module docstring for the
# false-positive caveat. Values are Pfam accessions (PFxxxxx), not names.
CANDIDATE_GENES = {
    "crtYB": ["PF00494"],
    "crtI": ["PF01593"],
    "crtS": ["PF00067"],
    "crtR": ["PF00970", "PF00258"],
    "HMGR": ["PF00368"],
    "GGPPS": ["PF00348"],
}


def require_bfd_duckdb() -> Path:
    if not BFD_DUCKDB.exists():
        sys.exit(
            f"BLOCKED: {BFD_DUCKDB} does not exist yet -- the BFD rebuild pipeline is still running "
            "(PI, 2026-08-15). Do not point this script at any other duckdb copy without explicit "
            "confirmation the rebuild is finished. Re-run this script once it lands at that exact path."
        )
    return BFD_DUCKDB


def require_reference_fasta() -> Path:
    if not REFERENCE_FASTA.exists():
        sys.exit(
            f"BLOCKED: {REFERENCE_FASTA} not found. Reference protein sequences for "
            "crtYB/crtI/crtS/crtR (e.g. from Xanthophyllomyces dendrorhous) and catalytic-residue "
            "positions have not been sourced yet -- see DEVELOPMENT_PLAN.md Part B Step 1 / Part D. "
            "Populate this FASTA (header format: >GENE_NAME|accession|species) before running --step ortholog."
        )
    return REFERENCE_FASTA


def step_pfam_prefilter(con) -> "pd.DataFrame":
    """Coarse Pfam-based candidate list per strain -- fast filter before the
    expensive ortholog-confirmation step. Returns strain x candidate-gene hit counts."""
    import pandas as pd

    all_pfam_ids = [pid for ids in CANDIDATE_GENES.values() for pid in ids]
    placeholders = ",".join(f"'{p}'" for p in all_pfam_ids)
    df = con.execute(
        f"""
        SELECT species_prefix, pfam_acc, protein_id
        FROM pfam
        WHERE pfam_acc IN ({placeholders})
        """
    ).fetchdf()
    print(f"Pfam pre-filter: {len(df)} raw hits across {df['species_prefix'].nunique()} strains", file=sys.stderr)
    return df


def step_ortholog_confirmation(prefilter_hits, reference_fasta: Path):
    """Reciprocal-best-hit diamond blastp of each strain's Pfam-pre-filtered
    candidate proteins against the curated reference set. Requires:
        module load diamond/2.1.24
    Not yet implemented beyond the pipeline skeleton below -- fill in once
    reference_fasta is populated and the strain protein FASTAs' exact path
    convention is confirmed against the rebuilt BFD.duckdb's ASMID scheme."""
    diamond_check = subprocess.run(["which", "diamond"], capture_output=True, text=True)
    if diamond_check.returncode != 0:
        sys.exit("`diamond` not on PATH -- run `module load diamond/2.1.24` first.")
    raise NotImplementedError(
        "Ortholog confirmation pipeline not yet built -- needs: (1) reference_fasta populated, "
        "(2) per-strain candidate protein sequences extracted (from BFD_PEP_DIR or gene_proteins "
        "table, matched via prefilter_hits' protein_id), (3) `diamond makedb` on the reference set, "
        "(4) `diamond blastp` each strain's candidates against it, (5) reciprocal best hit confirmed "
        "the other direction. See DEVELOPMENT_PLAN.md Part B Step 2."
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step", choices=["ortholog", "copy-number", "msa-lof"], required=True)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.step == "ortholog":
        import duckdb

        db_path = require_bfd_duckdb()
        ref_fasta = require_reference_fasta()
        con = duckdb.connect(str(db_path), read_only=True)
        prefilter_hits = step_pfam_prefilter(con)
        step_ortholog_confirmation(prefilter_hits, ref_fasta)

    elif args.step == "copy-number":
        sys.exit(
            "BLOCKED: --step copy-number depends on --step ortholog's confirmed gene calls "
            "(analysis/integrated_analysis/phase5_genome_linkage/ortholog_calls.csv, not yet produced). "
            "Run --step ortholog first once its NotImplementedError is resolved."
        )

    elif args.step == "msa-lof":
        sys.exit(
            "BLOCKED: --step msa-lof depends on --step ortholog's confirmed gene calls, plus "
            "catalytic-residue position annotations (not yet sourced -- see module docstring) and "
            "`module load mafft/7.505`. Not yet implemented."
        )


if __name__ == "__main__":
    main()
