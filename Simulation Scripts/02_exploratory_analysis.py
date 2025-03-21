#!/usr/bin/env python3
# Exploratory Analysis for Simulated GWAS data (Optional)
# This script will not be included in the final pipeline.
# ==== Initially done in Jupyter Notebook === #

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- Handle command-line arguments ---
if len(sys.argv) != 4:
    print("Usage: python3 02_exploratory_analysis.py <exposure_gwas.csv> <outcome_gwas.csv> <output_dir>")
    sys.exit(1)

exposure_path = sys.argv[1]
outcome_path = sys.argv[2]
output_dir = sys.argv[3]
os.makedirs(output_dir, exist_ok=True)

# --- Helper to parse GWAS CSV --- #
def parse_gwas(path):
    return pd.read_csv(os.path.expanduser(path))

# --- Load data ---
exposure = parse_gwas(exposure_path)
outcome = parse_gwas(outcome_path)

print(f"exposure shape: {exposure.shape}, outcome shape: {outcome.shape}")
print("\nExposure dtypes:\n", exposure.dtypes)
print("\nOutcome dtypes:\n", outcome.dtypes)

# --- Check for missing values --- #
print("\nMissing values in exposure:")
print(exposure.isna().sum())
print("\nMissing values in outcome:")
print(outcome.isna().sum())

# ================================ #
# --- Manhattan Plot: Exposure --- #
# ================================ #
df = exposure[["BP", "PVALUE", "CHR"]].copy()
df["-log10(PVALUE)"] = -np.log10(df["PVALUE"])

sns.set_style("whitegrid")
plt.figure(figsize=(16, 8), dpi=300)

chromosomes = sorted(df["CHR"].unique())
colors = ["#1f77b4", "#d62728"] * (len(chromosomes) // 2 + 1)

x_labels = []
x_ticks = []
x_offset = 0

for i, chrom in enumerate(chromosomes):
    subset = df[df["CHR"] == chrom]
    plt.scatter(
        subset["BP"] + x_offset,
        subset["-log10(PVALUE)"],
        color=colors[i % 2],
        s=10,
        alpha=0.75,
        edgecolors="none"
    )
    x_labels.append(chrom)
    x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
    x_offset += subset["BP"].max() - subset["BP"].min() + 1

plt.axhline(y=-np.log10(5e-8), color="black", linestyle="dashed", linewidth=1.5,
            label="Genome-wide significance (5e-8)")

plt.xticks(x_ticks, x_labels, rotation=0, fontsize=12, fontweight="bold")
plt.yticks(fontsize=12, fontweight="bold")
plt.xlabel("Chromosome", fontsize=14, fontweight="bold", labelpad=10)
plt.ylabel("-log10(P-value)", fontsize=14, fontweight="bold", labelpad=10)
plt.title("Manhattan Plot of Exposure GWAS", fontsize=18, fontweight="bold", pad=15)

sns.despine()
plt.legend(frameon=True, fontsize=12, loc="upper right", edgecolor="black")
plt.tight_layout()
exposure_plot_path = os.path.join(output_dir, "exposure_manhattan.png")
plt.savefig(exposure_plot_path, dpi=300)
plt.show()
print(f"[INFO] Saved: {exposure_plot_path}")

# =============================== #
# --- Manhattan Plot: Outcome --- #
# =============================== #

df2 = outcome[["BP", "PVALUE", "CHR"]].copy()
df2["-log10(PVALUE)"] = -np.log10(df2["PVALUE"])

sns.set_style("whitegrid")
plt.figure(figsize=(16, 8), dpi=300)

chromosomes = sorted(df2["CHR"].unique())
colors = ["#1f77b4", "#d62728"] * (len(chromosomes) // 2 + 1)

x_labels = []
x_ticks = []
x_offset = 0

for i, chrom in enumerate(chromosomes):
    subset = df2[df2["CHR"] == chrom]
    plt.scatter(
        subset["BP"] + x_offset,
        subset["-log10(PVALUE)"],
        color=colors[i % 2],
        s=10,
        alpha=0.75,
        edgecolors="none"
    )
    x_labels.append(chrom)
    x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
    x_offset += subset["BP"].max() - subset["BP"].min() + 1

plt.axhline(y=-np.log10(5e-8), color="black", linestyle="dashed", linewidth=1.5,
            label="Genome-wide significance (5e-8)")

plt.xticks(x_ticks, x_labels, rotation=0, fontsize=12, fontweight="bold")
plt.yticks(fontsize=12, fontweight="bold")
plt.xlabel("Chromosome", fontsize=14, fontweight="bold", labelpad=10)
plt.ylabel("-log10(P-value)", fontsize=14, fontweight="bold", labelpad=10)
plt.title("Manhattan Plot of Outcome GWAS", fontsize=18, fontweight="bold", pad=15)

sns.despine()
plt.legend(frameon=True, fontsize=12, loc="upper right", edgecolor="black")
plt.tight_layout()
outcome_plot_path = os.path.join(output_dir, "outcome_manhattan.png")
plt.savefig(outcome_plot_path, dpi=300)
plt.show()
print(f"[INFO] Saved: {outcome_plot_path}")

print("\n[INFO] Exploratory analyses complete!")
