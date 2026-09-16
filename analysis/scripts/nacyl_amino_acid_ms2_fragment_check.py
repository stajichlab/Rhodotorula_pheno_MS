#!/usr/bin/env python3
"""
MS2 fragment confirmation for the 3 structurally-corroborated N-acyl
amino acid candidates (2026-09-16, PI-requested validation step 1):
  - row 4109  -- N-C16:0-arginine  (SIRIUS: "Palmitoyl arginine")
  - row 51126 -- N-C14:0-arginine  (SIRIUS: "N2-(1-Oxotetradecyl)-L-arginine")
  - row 51152 -- N-C16:0-lysine    (SIRIUS: "N6-Palmitoyl lysine")

Pulls each row's real MS2 spectrum from the raw EB pipeline's
aligned_features.mgf (matched by FEATURE_ID, not re-derived or
simulated), then checks observed peaks against the diagnostic fragment
ions expected for a free arginine or lysine backbone -- the classic
"loss of the fatty acyl chain, revealing the underlying amino acid and
its own characteristic fragmentation" pattern used to confirm N-acyl
amino acids in the literature, analogous to how an intact AHL is
confirmed by its homoserine-lactone ring-loss fragment.

Diagnostic fragment set (monoisotopic, positive mode), computed from
first principles, not looked up as magic numbers:
  Arginine (C6H14N4O2, [M+H]+=175.11896):
    175.11896  [Arg+H]+                          (full acyl-chain loss)
    158.09241  [Arg+H-NH3]+                       (175.11896 - NH3)
    116.07061  [Arg+H-CH5N3]+                     (175.11896 - guanidine, CH5N3=59.04801)
     70.06511  C4H8N+                             (Arg side-chain pyrrolinium marker)
  Lysine (C6H14N2O2, [M+H]+=147.11281):
    147.11281  [Lys+H]+                           (full acyl-chain loss)
    130.08626  [Lys+H-NH3]+                       (147.11281 - NH3)
    129.10152  [Lys+H-H2O]+                       (147.11281 - H2O, cyclic lactam-type ion)
     84.08078  C5H10N+                            (secondary lysine marker, loss of NH3+CO from 129)

A hit is called at 15 mDa absolute tolerance (loose enough for this
Orbitrap-class low-mass-fragment data without being loose enough to match
arbitrary noise peaks; observed ppm error is reported for every hit so
the reader can judge fit quality directly rather than trust the
tolerance choice blindly).

Usage:
    python3 analysis/scripts/nacyl_amino_acid_ms2_fragment_check.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
MGF = (
    REPO / "data" / "processed" / "EB_20260130_ExFAB_Rhodo_Sup_and_Pellet"
    / "b773ffa18c2b41e5a3484526293a54f9" / "b773ffa18c2b41e5a3484526293a54f9"
    / "nf_output" / "feature_finding" / "feature_finding_results" / "aligned_features.mgf"
)
OUT_DIR = REPO / "analysis" / "ahl_autoinducer_search"
OUT_CSV = OUT_DIR / "nacyl_amino_acid_ms2_fragment_check.csv"

TARGET_ROWS = {
    4109: ("arginine", "N-C16:0-arginine (Palmitoyl arginine)"),
    51126: ("arginine", "N-C14:0-arginine (N-myristoyl-arginine)"),
    51152: ("lysine", "N-C16:0-lysine (N6-Palmitoyl lysine)"),
}

DIAGNOSTIC_FRAGMENTS = {
    "arginine": {
        175.11896: "[Arg+H]+ (full acyl-chain loss)",
        158.09241: "[Arg+H-NH3]+",
        116.07061: "[Arg+H-CH5N3]+ (guanidine loss)",
        70.06511: "C4H8N+ (Arg side-chain pyrrolinium)",
    },
    "lysine": {
        147.11281: "[Lys+H]+ (full acyl-chain loss)",
        130.08626: "[Lys+H-NH3]+",
        129.10152: "[Lys+H-H2O]+ (cyclic lactam-type ion)",
        84.08078: "C5H10N+ (secondary lysine marker)",
    },
}

TOLERANCE_DA = 0.015


def extract_spectra(mgf_path: Path, feature_ids: set[int]) -> dict[int, dict]:
    spectra = {}
    cur_lines: list[str] = []
    cur_id = None
    cur_meta = {}
    with open(mgf_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "BEGIN IONS":
                cur_lines = []
                cur_id = None
                cur_meta = {}
            elif line == "END IONS":
                if cur_id in feature_ids:
                    peaks = []
                    for pl in cur_lines:
                        parts = pl.split()
                        if len(parts) == 2:
                            try:
                                peaks.append((float(parts[0]), float(parts[1])))
                            except ValueError:
                                pass
                    spectra[cur_id] = {**cur_meta, "peaks": peaks}
            elif "=" in line and not line[0].isdigit():
                k, _, v = line.partition("=")
                if k == "FEATURE_ID":
                    cur_id = int(v)
                cur_meta[k] = v
            else:
                cur_lines.append(line)
    return spectra


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    spectra = extract_spectra(MGF, set(TARGET_ROWS))
    missing = set(TARGET_ROWS) - set(spectra)
    if missing:
        raise SystemExit(f"Spectra not found in MGF for row IDs: {missing}")

    rows = []
    for row_id, (aa, label) in TARGET_ROWS.items():
        sp = spectra[row_id]
        peaks = sorted(sp["peaks"], key=lambda p: -p[1])  # sort by intensity desc
        max_intensity = max(p[1] for p in peaks) if peaks else 1.0
        print(f"\n{'='*70}\nrow {row_id} -- {label}")
        print(f"PEPMASS={sp.get('PEPMASS')} SOURCE_FILE={sp.get('SOURCE_FILE')} n_peaks={len(peaks)}")
        diag = DIAGNOSTIC_FRAGMENTS[aa]
        n_hit = 0
        for expected_mz, fragment_label in diag.items():
            best = None
            for obs_mz, obs_intensity in peaks:
                if abs(obs_mz - expected_mz) <= TOLERANCE_DA:
                    if best is None or obs_intensity > best[1]:
                        best = (obs_mz, obs_intensity)
            if best is not None:
                obs_mz, obs_intensity = best
                ppm_error = (obs_mz - expected_mz) / expected_mz * 1e6
                rel_intensity = obs_intensity / max_intensity * 100
                n_hit += 1
                print(f"  HIT  {fragment_label:45s} expected={expected_mz:.5f} observed={obs_mz:.4f} "
                      f"ppm_error={ppm_error:+.1f} rel_intensity={rel_intensity:.1f}%")
                rows.append(dict(row_id=row_id, compound=label, amino_acid=aa, fragment=fragment_label,
                                  expected_mz=expected_mz, observed_mz=obs_mz, ppm_error=ppm_error,
                                  observed_intensity=obs_intensity, rel_intensity_pct=rel_intensity, hit=True))
            else:
                print(f"  miss {fragment_label:45s} expected={expected_mz:.5f} -- no peak within {TOLERANCE_DA*1000:.0f} mDa")
                rows.append(dict(row_id=row_id, compound=label, amino_acid=aa, fragment=fragment_label,
                                  expected_mz=expected_mz, observed_mz=None, ppm_error=None,
                                  observed_intensity=None, rel_intensity_pct=None, hit=False))
        print(f"  {n_hit}/{len(diag)} diagnostic {aa} fragments matched")

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
