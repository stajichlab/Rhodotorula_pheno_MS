"""Shared helpers for joining BFD genome-derived amino-acid composition
to strain-level phenotype metadata via the phyling protein tree.

Strain identifiers differ in format across the three source files:
  - tree tip labels:  Rhodotorula_mucilaginosa_DBVPG_3776.proteins
  - aa_freq filename:  DBVPG3776.aa_freq.csv.gz  (species_prefix column matches)
  - phenotype metadata SAMPLE_NAME / Strain: DBVPG_3776

We match by a normalized key (strip underscores/dashes/spaces, uppercase)
and require the normalized aa_freq/metadata key to be a *substring* of the
normalized tree tip label (tree tips carry extra species-name prefix text).
"""
import csv
import glob
import gzip
import os
import re

REPO = "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS"
TREEFILE = os.path.join(
    REPO,
    "BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/"
    "protein-Rhodotorula-taxa_276.fungi_odb10.fasttree.support.treefile",
)
AA_FREQ_GLOB = os.path.join(REPO, "BFD/results/genome_stats/aa_freq/*/*.aa_freq.csv.gz")
CU_AUC_CSV = os.path.join(REPO, "data/metadata/EXFAB_UCR-005/Cu_AUC.20260811.fixed.csv.gz")
YPD2_CSV = os.path.join(REPO, "data/metadata/EXFAB_UCR-005/YPD2_phenotypic.20260702.fixed.csv.gz")

AA20 = list("ACDEFGHIKLMNPQRSTVWY")


def norm(s):
    return re.sub(r"[_\-\s]", "", s).upper()


def load_tree_tips(treefile=TREEFILE):
    with open(treefile) as fh:
        s = fh.read()
    tips = re.findall(r"[(,]([A-Za-z0-9_.\-]+):", s)
    return sorted(set(tips))


def match_to_tree(key, tip_norm_map):
    """key: raw strain id string. Returns tree tip label or None.

    Tree tip labels are `<species words>_<strain code>.proteins`; the strain
    code is always the trailing segment before `.proteins`. We therefore
    require the normalized key to be a *suffix* of the tip's normalized
    label (after stripping the trailing "PROTEINS" marker), not merely a
    substring anywhere in it. A plain substring test is too permissive:
    e.g. normalized "TFCN48D1" is a substring of "...TFCN48D10.PROTEINS"
    (the "10" strain), so a naive `in` check wrongly matches both
    "TFCN_48D-1" and "TFCN_48D-10" to the same "...48D10" tip.
    """
    nk = norm(key)
    hits = []
    for tnorm, tip in tip_norm_map.items():
        stripped = re.sub(r"\.?PROTEINS$", "", tnorm)
        if stripped.endswith(nk):
            hits.append(tip)
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        return sorted(hits, key=len)[0]
    return None


def load_aa_freq_table():
    """Returns dict: aa_freq_prefix -> {AA: frequency}"""
    table = {}
    # glob enumerates the full, fixed set of per-strain aa_freq outputs from
    # one BFD run (one file per strain, no latest/backup ambiguity to
    # resolve). The hashed subdirectories under aa_freq/ could in principle
    # produce a duplicate strain prefix, which would silently overwrite a
    # dict entry; the explicit check below turns that into a hard failure
    # (not an `assert`, so it survives python -O) instead of a silent bad join.
    # ANALYSIS_OK[file-selection]: fixed BFD run output dir, duplicate-prefix checked below
    files = glob.glob(AA_FREQ_GLOB)
    prefixes = [os.path.basename(f).split(".")[0] for f in files]
    if len(prefixes) != len(set(prefixes)):
        raise ValueError("duplicate aa_freq strain prefix across hashed subdirectories")
    for f in files:
        prefix = os.path.basename(f).split(".")[0]
        freqs = {}
        with gzip.open(f, "rt") as fh:
            r = csv.DictReader(fh)
            for row in r:
                freqs[row["amino_acid"]] = float(row["frequency"])
        table[prefix] = freqs
    return table


def load_cu_auc():
    rows = []
    with gzip.open(CU_AUC_CSV, "rt") as fh:
        r = csv.DictReader(fh)
        for row in r:
            rows.append(row)
    return rows


def load_ypd2():
    rows = []
    with gzip.open(YPD2_CSV, "rt") as fh:
        r = csv.DictReader(fh)
        for row in r:
            rows.append(row)
    return rows
