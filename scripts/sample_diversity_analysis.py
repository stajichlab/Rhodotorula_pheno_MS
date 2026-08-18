#!/usr/bin/env python3
"""
Analyze diversity of non-zero peak area values within each R. mucilaginosa sample.
Focus on what the non-zero values look like - are they all targeted to the same value?
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
LINKED = REPO / "analysis" / "linked_data"
OUT_DIR = REPO / "analysis" / "sample_diversity"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Load metadata and filter to R. mucilaginosa
meta = pd.read_csv(LINKED / "sample_metadata.csv.gz")
mucil_meta = meta[meta["Species"] == "Rhodotorula mucilaginosa"].copy()
mucil_meta["sample_type"] = mucil_meta["sample_id"].apply(
    lambda x: "SUP" if x.startswith("SUP") else "C"
)

# Load feature matrix
feature_df = pd.read_csv(LINKED / "feature_abundance_matrix.csv.gz", compression="gzip")
meta_cols = ["row ID", "row m/z", "row retention time", "adduct", 
             "is_default_adduct", "has_ms2", "detection_count", "detection_rate"]
sample_cols = [c for c in feature_df.columns if c not in meta_cols]
mucil_sample_ids = mucil_meta["sample_id"].tolist()

def compute_nonzero_diversity(sample_id):
    """Compute diversity metrics for non-zero values in a sample."""
    vals = feature_df[sample_id].values.astype(float)
    nonzero = vals[vals > 0]
    total = len(vals)
    n_nonzero = len(nonzero)
    
    result = {
        "sample_id": sample_id,
        "total_features": total,
        "n_nonzero": n_nonzero,
        "prop_nonzero": n_nonzero / total if total > 0 else 0,
    }
    
    if n_nonzero < 10:
        result.update({
            "nonzero_mode_prop": np.nan,
            "nonzero_n_unique": 0,
            "nonzero_cv": np.nan,
            "nonzero_median": np.nan,
            "nonzero_mean": np.nan,
        })
        return result
    
    # For non-zero values: what fraction share the most common value?
    unique_nz, counts_nz = np.unique(nonzero, return_counts=True)
    mode_prop_nz = counts_nz.max() / n_nonzero
    
    # Coefficient of variation of non-zero values
    cv = nonzero.std() / nonzero.mean() if nonzero.mean() > 0 else np.nan
    
    result.update({
        "nonzero_mode_prop": mode_prop_nz,
        "nonzero_n_unique": len(unique_nz),
        "nonzero_cv": cv,
        "nonzero_median": np.median(nonzero),
        "nonzero_mean": np.mean(nonzero),
    })
    return result

print("Computing non-zero diversity metrics...")
results = [compute_nonzero_diversity(sid) for sid in mucil_sample_ids]
metrics_df = pd.DataFrame(results)

# Merge with metadata
metrics_df = metrics_df.merge(mucil_meta[["sample_id", "sample_type"]], on="sample_id")

# Save
metrics_df.to_csv(OUT_DIR / "nonzero_diversity_metrics.csv", index=False)

# Separate by sample type
sup_data = metrics_df[metrics_df["sample_type"] == "SUP"]
c_data = metrics_df[metrics_df["sample_type"] == "C"]

print(f"\nSUP (n={len(sup_data)}): nonzero median mode proportion = {sup_data['nonzero_mode_prop'].median():.4f}")
print(f"C (n={len(c_data)}): nonzero median mode proportion = {c_data['nonzero_mode_prop'].median():.4f}")

# Create plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Plot 1: Mode proportion of non-zero values
ax = axes[0, 0]
ax.hist(sup_data["nonzero_mode_prop"].dropna(), bins=30, alpha=0.7, label="SUP", color="orange")
ax.hist(c_data["nonzero_mode_prop"].dropna(), bins=30, alpha=0.7, label="C", color="blue")
ax.set_xlabel("Mode Proportion (Non-Zero Values)")
ax.set_ylabel("Count")
ax.set_title("Non-Zero Value Uniformity")
ax.legend()
ax.axvline(sup_data["nonzero_mode_prop"].median(), color="orange", linestyle="--")
ax.axvline(c_data["nonzero_mode_prop"].median(), color="blue", linestyle="--")

# Plot 2: Number of unique non-zero values
ax = axes[0, 1]
ax.hist(sup_data["nonzero_n_unique"], bins=30, alpha=0.7, label="SUP", color="orange")
ax.hist(c_data["nonzero_n_unique"], bins=30, alpha=0.7, label="C", color="blue")
ax.set_xlabel("Number of Unique Non-Zero Values")
ax.set_ylabel("Count")
ax.set_title("Non-Zero Value Diversity")
ax.legend()

# Plot 3: Coefficient of variation
ax = axes[0, 2]
ax.hist(sup_data["nonzero_cv"].dropna(), bins=30, alpha=0.7, label="SUP", color="orange")
ax.hist(c_data["nonzero_cv"].dropna(), bins=30, alpha=0.7, label="C", color="blue")
ax.set_xlabel("Coefficient of Variation")
ax.set_ylabel("Count")
ax.set_title("Non-Zero Value Spread")
ax.legend()

# Plot 4: Proportion non-zero
ax = axes[1, 0]
ax.hist(sup_data["prop_nonzero"], bins=30, alpha=0.7, label="SUP", color="orange")
ax.hist(c_data["prop_nonzero"], bins=30, alpha=0.7, label="C", color="blue")
ax.set_xlabel("Proportion Non-Zero")
ax.set_ylabel("Count")
ax.set_title("Feature Detection Rate")
ax.legend()

# Plot 5: Box plots comparing mode proportion
ax = axes[1, 1]
data_to_plot = [sup_data["nonzero_mode_prop"].dropna(), c_data["nonzero_mode_prop"].dropna()]
bp = ax.boxplot(data_to_plot, labels=["SUP", "C"], patch_artist=True)
bp["boxes"][0].set_facecolor("orange")
bp["boxes"][1].set_facecolor("blue")
ax.set_ylabel("Mode Proportion (Non-Zero)")
ax.set_title("Non-Zero Value Uniformity Comparison")

# Plot 6: Scatter - mode prop vs detection rate
ax = axes[1, 2]
ax.scatter(sup_data["prop_nonzero"], sup_data["nonzero_mode_prop"], 
           alpha=0.5, label="SUP", color="orange", s=30)
ax.scatter(c_data["prop_nonzero"], c_data["nonzero_mode_prop"], 
           alpha=0.5, label="C", color="blue", s=30)
ax.set_xlabel("Proportion Non-Zero")
ax.set_ylabel("Mode Proportion (Non-Zero)")
ax.set_title("Uniformity vs Detection Rate")
ax.legend()

plt.suptitle("Non-Zero Value Diversity: R. mucilaginosa SUP vs Cell", fontsize=14)
plt.tight_layout()
fig.savefig(OUT_DIR / "nonzero_diversity_histograms.png", dpi=150, bbox_inches="tight")
fig.savefig(OUT_DIR / "nonzero_diversity_histograms.pdf", bbox_inches="tight")

# Statistical test
t_stat, p_val = stats.ttest_ind(
    sup_data["nonzero_mode_prop"].dropna(),
    c_data["nonzero_mode_prop"].dropna()
)
print(f"\nT-test for mode proportion: t={t_stat:.3f}, p={p_val:.2e}")

print(f"\nSaved to {OUT_DIR}")
