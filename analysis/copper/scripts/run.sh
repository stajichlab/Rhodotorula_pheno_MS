#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

python3 00_build_master_table.py
python3 01_naive_correlation.py
Rscript 02_pgls_analysis.R
python3 03_compare_naive_vs_pgls.py
Rscript 04_sensitivity_species_balance.R
