# Species re-assignment record (2026-08-16)

Whole-genome fastANI confirms 4 strains are mis-specified. Screen method:
whole-genome fungi_odb10 protein tree monophyly flagging -> targeted fastANI
(1.33, fragment len 3000) against labeled vs. closest species.

ANI criterion: >95% = same species; ~78% = sister-species background.

## Confirmed re-assignments

| Strain ID | Strain            | Old species             | New species             | ANI to new species (best) | ANI to old species (best) |
|-----------|-------------------|-------------------------|-------------------------|---------------------------|---------------------------|
| 28        | TFCN_1A-1-3       | R. pacifica             | R. mucilaginosa         | 99.9986 (TFCN_98C-8)     | n/a (no R. pacifica genome in set) |
| 80        | TFCN_1B-1-2       | R. pacifica             | R. mucilaginosa         | 99.9936 (TFCN_152A-12)   | n/a                        |
| 157       | TFCN_1A-1-2       | R. paludigena           | R. mucilaginosa         | 99.9990 (TFCN_98A-10)    | 77.98 (TFCN_43A-4)        |
| 105       | TFCN_152C-6       | R. toruloides           | R. taiwanensis          | 99.9994 (TFCN_357-1)     | 77.68 (DBVPG_6121)        |

Notes:
- The two R. pacifica-labeled strains are mucilaginosa; no *true* R. pacifica
  genome is present in BFD/input/dna, so no within-species comparison.
- These four were caught because their RA-NRPS gene copies fell inside the
  R. mucilaginosa (or R. taiwanensis) clade of the NRPS gene tree.

## Unresolved / needs data

| Strain ID | Strain      | Labeled species | Status                                              |
|-----------|-------------|-----------------|-----------------------------------------------------|
| 288       | TFCN_1A_1_5 | R. pacifica     | No genome in BFD/input/dna; ANI unverifiable. Likely R. mucilaginosa given the other 2 pacifica strains are, but needs a genome to confirm. |

## R. paludigena sanity check

The 14 remaining R. paludigena strains all have ~69-71% pident RA-NRPS and
were NOT flagged by the whole-genome tree. TFCN_1A-1-2 was the only paludigena
that clustered with mucilaginosa. Screen shows R. paludigena species clade is
intact after removing TFCN_1A-1-2 (n=15 -> 14).

## Records updated

- `BFD/samples.csv` (genome/ASMID record: SPECIES_IN, SPECIES, NCBI_TAXONID)
- `data/metadata/EXFAB_UCR-005/YPD2_phenotypic.20260702.fixed.csv.gz` (phenotype Species col)
- `data/metadata/EXFAB_UCR-005/Cu_AUC.20260811.fixed.csv.gz` (SPECIES col)
- `data/metadata/EXFAB_UCR-005/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz`
  (Verified Species / db_species_full / Species cols)

## ANI evidence

- `ANI_check/ANI.out` (3 mucilaginosa-confirmed queries vs 217 refs)
- `ANI_check/t152/ANI.out` (TFCN_152C-6 vs 15 toruloides+taiwanensis refs)
