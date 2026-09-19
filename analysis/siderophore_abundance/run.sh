#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${0}")/../.."
python3 analysis/siderophore_abundance/scripts/siderophore_abundance_table.py
