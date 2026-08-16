# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-08-15] pandas groupby silently drops rows with a NULL group key

**Category**: gotcha

**What happened**: `build_strain_phenotype_table.py` groups the raw
phenotype CSV by `strain_code` (or `Strain` for the ypd2 source) to
collapse replicate rows. The `control_90_110` source has 10/314 rows with
`strain_code = NaN` ("unidentified spots" in the upstream pipeline's own
terms). `df.groupby(strain_col)` drops NaN keys by default — these 10
rows vanish from the output with no warning, no error, no row-count
mismatch that would be obvious without deliberately checking
`len(df)` vs. `df[strain_col].notna().sum()` before and after.

**Why it matters**: silent data loss. These 10 rows have real color/area
measurements; if any of them happen to be strains of analytical interest
(e.g. one more *R. dairenensis* replicate, or a species with few other
strains), they're invisibly absent from every downstream table without a
flag anywhere. Found only by manually auditing `data/raw/control_phenotype_90_110h/`
during formal ingestion — would not have been caught by just checking the
output row count "looked about right."

**Resolution**: documented in `data/metadata/control_phenotype_90_110h/schema.yaml`
(strain_code column notes) and `provenance.md` (Known Issues), AND fixed
in `build_strain_phenotype_table.py`: it now explicitly drops
`strain_col`-NULL rows *before* the groupby (rather than relying on
groupby's implicit NaN-key drop) and prints a `WARNING` with the affected
`strain_id`s to stderr every run, plus records the count in
`strain_phenotype_table_diagnostics.txt`. The 10 rows are still excluded
from the output (no strain_code to key them on), but the exclusion is now
loud, not silent.

**Tags**: pandas, groupby, data-loss, phenotype-ingestion, silent-failure

**mitigation_type**: structural

**structural_mitigation_candidate**: (shipped, see Resolution) —
`build_strain_phenotype_table.py` prints a WARNING + strain_id list for
any strain_col-NULL row and logs the count to the diagnostics file.

### [2026-08-15] NaN correlation values silently produce spuriously significant permutation p-values

**Category**: gotcha

**What happened**: In `phase2_color_metabolome_association.py`'s
block-permutation test, constant-abundance features (zero variance across
strains after TSS normalization -- 533-785 of 10,949 depending on
fraction) produce a Spearman rho of NaN (0/0 in the correlation formula).
The exceed-count comparison `np.abs(perm_rho) >= np.abs(observed_rho)`
silently evaluates to `False` for any comparison involving NaN in numpy
(no error, no warning propagated to the accumulator) -- so these features
never "exceeded" across any permutation, giving them the minimum possible
empirical p-value `1/(n_perm+1)` regardless of the fact that their
correlation is literally undefined. First real run put several NaN-rho
features at the very top of the "most significant" list.

**Why it matters**: this is a general trap for any permutation-test
implementation using `>=`/`<=` comparisons against a vectorized statistic
that can be NaN for degenerate inputs (constant columns, all-zero rows,
etc.) -- NaN doesn't raise, it silently participates in comparisons as
"never true," which for a *count of times exceeded* means NaN inputs look
like the *strongest* possible hits, exactly backwards from correct
(undefined, not significant).

**Resolution**: explicitly compute a `valid = ~np.isnan(observed_rho)`
mask before the permutation loop, restrict the permutation comparison and
all output columns (empirical_p, empirical_fdr, asymptotic_p,
asymptotic_fdr) to valid features only, and set invalid features to NaN in
every output column rather than letting them flow through the exceed-count
logic at all.

**Tags**: permutation-test, nan-handling, numpy, silent-failure, statistics

**mitigation_type**: structural

**structural_mitigation_candidate**: (shipped, see Resolution) — explicit
`valid` mask computed and applied before any permutation comparison in
`phase2_color_metabolome_association.py`.

### [2026-08-15] `genome_strain_species_busco_map.csv` has no generating script — species-level analyses can't pick up new BFD taxa on their own

**Category**: gotcha

**What happened**: While rebuilding phylogeny-dependent outputs after the
PHYling tree grew from 276 to 278 taxa, `prune_species_tree.R` silently
(well, with a warning) dropped both new taxa from the rebuilt
`species_tree.nwk` because
`analysis/integrated_analysis/phase1_phenotype/genome_strain_species_busco_map.csv`
doesn't list them. Searched the repo for whatever produced that CSV — no
script, no `run.sh`, nothing in `.living/decisions.md` beyond a note that
it was "built from BFD directly." It's a hand-built artifact with no
reproduction path.

**Why it matters**: any future BFD taxon-count bump (this project already
went through one at 276->278) will look "rebuilt" at the strain-tree
level while silently staying stale at every species-level analysis
(`species_tree.nwk`, `phylogenetic_signal.R`, `convergent_color_test.R`,
and anything built on `species_phenotype_table.csv`) until someone
manually notices the warning and regenerates the map by hand.

**Resolution**: not fixed this session (out of scope for a tree-pointer
rebuild) — flagged as a todo instead. A real fix would be a small script
(e.g. `analysis/scripts/build_busco_map.R` or `.py`) that queries BFD
directly for `species_prefix`/`complete_pct` and joins to a
strain-name-normalization step, so this map can be regenerated with one
command instead of whatever ad hoc process built it originally.

**Tags**: phylogeny, reproducibility, missing-script, species-level-analysis, busco-map

**mitigation_type**: none

**structural_mitigation_candidate**: write `build_busco_map.{R,py}` to
make `genome_strain_species_busco_map.csv` regeneration a one-command
step (see linked todo).

### [2026-08-15] "Hard gate" negative-control checks need to verify statistical adequacy, not just file freshness

**Category**: failure

**What happened**: The Phase 2 within-species scripts' hard gate (refuse
to run a real predictor without a fresh negative-control decoy run)
checked only that the decoy output file existed and was newer than the
input data. A `--n-perm 20` smoke test (run to verify the code worked,
not as a real analysis) satisfied both conditions and got silently
treated as a valid negative control. The resulting "0 hits" from that
under-powered run was written into `WITHIN_SPECIES_MUCILAGINOSA.md` and a
`.living/findings/` entry as "the area confound doesn't reproduce within-
species" — a properly-powered rerun (`--n-perm 200`) later showed 2,025
hits, contradicting that claim. Caught only because a newly-written,
independent test (the ANOVA/pattern-group script, whose decoy WAS run at
full power from the start) gave a wildly different decoy result for the
same predictor/species/fraction, which triggered a direct investigation
rather than being accepted.

**Why it matters**: "the negative control ran" and "the negative control
was adequate" are different claims, and a hard gate that only checks the
former gives false confidence — exactly the failure mode the gate was
built to prevent, just moved one level up. A smoke test's output living
at the same file path as the real deliverable is an easy trap: nothing
about the file's existence or timestamp distinguishes "I verified the
code runs" from "I ran the actual negative control."

**Resolution**: `phase2_within_species_association.py` and
`phase2_anova_pattern_association.py` now write an `n_perm` column into
every output row and the hard gate additionally requires
`decoy_df["n_perm"].min() >= 100` before permitting a real-predictor run.
Not yet backported to `phase2_color_metabolome_association.py` or
`phase2_multivariate_association.py` (see `todo/TODOLIST.md`).

**Tags**: hard-gate, negative-control, smoke-test, permutation-test, statistical-power, silent-failure

**mitigation_type**: structural

**structural_mitigation_candidate**: (partially shipped, see Resolution)
— `n_perm` recording + minimum-permutation-count gate check, in 2 of 4
Phase 2 scripts so far.
