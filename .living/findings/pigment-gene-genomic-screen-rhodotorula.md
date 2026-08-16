---
topic: pigment-gene-genomic-screen-rhodotorula
description: Genome-level presence/copy-number screening of carotenoid and melanin pathway genes across the 278-strain BFD panel, via Pfam domains, SwissProt best-hits, and a PI-curated custom HMM panel.
created: 2026-08-16
last_updated: 2026-08-16
status: active
---

# Pigment-pathway gene screening across the Rhodotorula BFD genome panel

## F-001: Three independent screens agree the melanin pathway route in this genus is laccase-based, not tyrosinase-based
**Status:** preliminary
**Claim:** A coarse Pfam-domain screen (multicopper oxidase domains
PF00394/PF07731/PF07732 vs. tyrosinase PF00264), a SwissProt best-hit
keyword cross-reference (LAC*/LCC* families vs. TYR/tyr1/melC2), and the
PI's custom 28-profile pigment HMM panel (`laccase` vs. `tyrosinase`
profiles) **all independently agree**: laccase/multicopper-oxidase hits
are abundant and present in all 278 genomes (Pfam: 3,437 hits/278 strains;
SwissProt: LAC*/LCC* families >10,000 hits; HMM panel: 1,220 hits/278
strains, ~2-5.5 copies/genome), while tyrosinase hits are essentially
absent in all three (Pfam: 2 hits/1 strain; SwissProt: 7 hits total across
3 gene-name variants; HMM panel: 3 hits/1 strain).
**Implications:** Convergent agreement across 3 independent, methodologically
different screens is much stronger evidence than any one method alone.
This genus appears to melanize via the laccase route rather than the
tyrosinase route -- relevant context for interpreting any future
melanin-pathway <-> color association test, and for prioritizing which
gene families are worth deeper ortholog-confirmation effort.
**Tags:** melanin-pathway, laccase, tyrosinase, pfam, swissprot, hmm-panel, convergent-evidence

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Pfam screen (pfam_pigment_screen.py) + SwissProt crossref (swissprot_pigment_crossref.py) + custom HMM panel (hmmsearch + parse_hmm_hits.py) | BFD.duckdb pfam/swissprot tables + BFD/input/pep proteomes (278 strains, 2,188,032 proteins) | Rhodotorula_pheno_MS | Laccase present in ~278/278 strains across all 3 methods; tyrosinase near-absent (0-1 strains) across all 3 methods | refines |

## F-002: The custom HMM panel's zero-hit families likely reflect lineage-specific paralog discrimination, not panel failure -- but 4 gene names need a PI sanity check
**Status:** preliminary
**Claim:** Of the PI's 28 pigment HMM profiles, 10 found zero hits genome-wide
at E<1e-5: `ayg1`, `crtB`, `hppd`, `scd`, `scyB/C/D/E/F`, `mysB`. Within
each multi-paralog family (scyA-F, mysA-E), only one member consistently
hits (scyA at ~1 copy/genome; mysA/C/D/E at ~1, 0, 0, ~15 copies/genome
respectively) while the others are silent -- consistent with the panel
correctly discriminating fungal-lineage-specific paralogs rather than a
threshold or search failure. `crtB` (phytoene synthase, core carotenogenesis)
being silent while `crtP`/`crt_fungal_psy` hit broadly is plausibly explained
by many basidiomycetes using a fused bifunctional crtYB
(phytoene-synthase + lycopene-cyclase) instead of a separate bacterial-type
crtB -- a textbook fungal carotenogenesis distinction -- but this is an
inference, not confirmed.
**Implications:** Do not treat the 10 zero-hit families as "gene absent
from the genus" without PI confirmation; `ayg1`/`hppd`/`scd`/the missing
scy/mys paralogs specifically need a sanity check on whether they were
expected to be present in *Rhodotorula* at all, since the panel may have
been built from a different fungal lineage's pigment gene set.
**Tags:** hmm-panel, paralog-discrimination, needs-pi-confirmation, crtB, ayg1

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Custom HMM panel search (parse_hmm_hits.py) | Same as F-001 | Rhodotorula_pheno_MS | 10/28 profiles zero hits genome-wide; within-family paralog selectivity observed | refines |

### Open Questions
- Are `ayg1`, `hppd`, `scd`, `scyB-F` (except A), `mysB` genuinely expected
  to be absent from *Rhodotorula*, or were these HMMs built from a
  different fungal lineage's sequences and simply too divergent to hit at
  this threshold? (PI has the seed-alignment provenance for these HMMs.)

## F-003: Several HMM-panel and Pfam families show implausibly high copy numbers, suggesting superfamily-level (not ortholog-level) domain matches
**Status:** preliminary
**Claim:** `t3hnr` (~20-45 copies/genome), `t4hnr` (~35-64), `pks_melanin`
(~18-27), `mysE` (~12-18), `crtP` (~6-14) from the HMM panel, and
`crtS/crtR-hydroxylase` (PF00067, cytochrome P450, ~9-35/genome) and
`crtI` (PF01593, FAD-oxidoreductase, ~7-11/genome) from the Pfam screen,
are far too high to represent single-gene ortholog copy number. These
Pfam domains in particular (P450, FAD-oxidoreductase) are large
superfamilies genome-wide; the HMM-panel families likely score a shared
catalytic domain (SDR reductase fold for t3hnr/t4hnr, ketosynthase domain
for pks_melanin) shared across a broader gene family.
**Implications:** None of these counts should be used as "gene copy
number" in any downstream genome<->phenotype test without first tightening
orthology confirmation (stricter threshold, reciprocal-best-hit, or a
gene tree) -- using the raw counts would silently conflate a real
pathway-gene signal with generic superfamily size variation.
**Tags:** hmm-panel, pfam, superfamily-domain, ortholog-confirmation-needed, false-copy-number-risk

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | HMM panel + Pfam screen (same runs as F-001/F-002) | Same as F-001 | Rhodotorula_pheno_MS | 5 HMM-panel families and 2 Pfam families show 6-64 copies/genome, implausible for single orthologs | refines |

### Open Questions
- Which of `crtP`/`laccase`/`mysE` (the copy-number-variable-but-plausible
  families, as opposed to the clearly-superfamily ones) are worth
  promoting to a real ortholog-confirmation step (diamond blastp +
  mafft MSA + gene tree, as originally scaffolded in
  `phase5_candidate_gene_genotyping.py`'s `step_ortholog_confirmation()`)?
- Since most hit families are present in ~all 278 strains (no presence/
  absence variance), the only immediately testable genome<->color
  candidates are the sparse-presence families (`crtR` 52/278, `crtQ`
  48/278, `hgd` 28/278) -- has anyone checked whether that presence
  pattern tracks species or color yet? (Not yet done as of this entry.)
