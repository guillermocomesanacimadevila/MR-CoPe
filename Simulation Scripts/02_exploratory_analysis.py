#!/usr/bin/env python3
# Exploratory Analysis for Simulated GWAS data (Optional)
# This script will not be inlcluded in the final pipeline.
# ==== Initially done in Jupyter Notebook === #

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Handle command-line arguments ---
if len(sys.argv) != 4:
    print("Usage: python3 02_exploratory_analysis.py <exposure_gwas.csv> <outcome_gwas.csv> <output_dir>")
    sys.exit(1)

exposure_path = sys.argv[1]
outcome_path = sys.argv[2]
output_dir = sys.argv[3]
# Note: output_dir kept for compatibility, but not used in file paths

# --- Helper to parse GWAS CSV ---
def parse_gwas(location):
    with open(os.path.expanduser(location), "r") as file:
        return pd.read_csv(file, sep=",", index_col=0)

# --- Load data ---
exposure = parse_gwas(exposure_path)
outcome = parse_gwas(outcome_path)

print(f"exposure shape: {exposure.shape}, outcome shape: {outcome.shape}")
print("\nExposure dtypes:")
print(exposure.dtypes)
print("\nOutcome dtypes:")
print(outcome.dtypes)

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

plt.xticks(x_ticks, x_labels, rotation=90, fontsize=12, fontweight="bold")
plt.yticks(fontsize=12, fontweight="bold")
plt.xlabel("Chromosome", fontsize=14, fontweight="bold", labelpad=10)
plt.ylabel("-log10(P-value)", fontsize=14, fontweight="bold", labelpad=10)
plt.title("Manhattan Plot of Exposure GWAS", fontsize=18, fontweight="bold", pad=15)

plt.legend(frameon=True, fontsize=12, loc="upper right", edgecolor="black")
plt.tight_layout()
plt.savefig("exposure_manhattan.png", dpi=300)
plt.show()
print("[INFO] Saved: exposure_manhattan.png")

# =============================== #
# --- Manhattan Plot: Outcome --- #
# =============================== #
df2 = outcome[["BP", "PVALUE", "CHR"]].copy()
df2["-log10(PVALUE)"] = -np.log10(df2["PVALUE"])

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

plt.xticks(x_ticks, x_labels, rotation=90, fontsize=12, fontweight="bold")
plt.yticks(fontsize=12, fontweight="bold")
plt.xlabel("Chromosome", fontsize=14, fontweight="bold", labelpad=10)
plt.ylabel("-log10(P-value)", fontsize=14, fontweight="bold", labelpad=10)
plt.title("Manhattan Plot of Outcome GWAS", fontsize=18, fontweight="bold", pad=15)

plt.legend(frameon=True, fontsize=12, loc="upper right", edgecolor="black")
plt.tight_layout()
plt.savefig("outcome_manhattan.png", dpi=300)
plt.show()
print("[INFO] Saved: outcome_manhattan.png")

print("\n[INFO] Exploratory analyses complete!")
