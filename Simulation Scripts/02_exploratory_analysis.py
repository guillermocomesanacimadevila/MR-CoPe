# Exploratory Analysis for Simulated GWAS data (Optiomal)
# This script will not be inlcluded in the final pipeline.
# ==== Initially done in Jupyter Notebook === #

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_gwas(location):
    with open(os.path.expanduser(location), "r") as file:
        return pd.read_csv(file, sep=",", index_col=0)

# ! pwd
# ! ls -lh

ldl_c = parse_gwas("/home/jovyan/MR_simulation/Scripts/ldl_gwas.csv")
ad = parse_gwas("/home/jovyan/MR_simulation/Scripts/ad_gwas.csv")

print(f"ldl_c shape: {ldl_c.shape}, ad shape: {ad.shape}")
print(ldl_c.dtypes)
print(ad.dtypes)

# Check for missing values
print(ldl_c.isna().sum())
print(ad.isna().sum())

# Plot Manhattan / Scatter
df = ldl_c[["BP", "PVALUE", "CHR"]]
df2 = ad[["BP", "PVALUE", "CHR"]]

# Transform pvalues to f(p) = -log10(p)
df["-log10(PVALUE)"] = -np.log10(df["PVALUE"])
df2["-log10(PVALUE)"] = -np.log10(df2["PVALUE"])

# Set up the Manhattan plot -> Exposure
plt.figure(figsize=(14, 6))

# Define colors for alternating chromosomes
chromosomes = df["CHR"].unique()
colors = ["blue", "red"] * (len(chromosomes) // 2 + 1)

x_labels = []
x_ticks = []
x_offset = 0

# Sort chromosomes for consistent plotting
for i, chrom in enumerate(sorted(chromosomes)):
    subset = df[df["CHR"] == chrom]
    plt.scatter(subset["BP"] + x_offset, subset["-log10(PVALUE)"], color=colors[i % 2], s=10)
    x_labels.append(chrom)
    x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
    x_offset += subset["BP"].max() - subset["BP"].min() + 1

# Add genome-wide significance threshold line (-log10(5e-8) = 7.3)
plt.axhline(y=-np.log10(5e-8), color="black", linestyle="--", label="Genome-wide significance (5e-8)")

# Formatting
plt.xticks(x_ticks, x_labels, rotation=90)
plt.xlabel("Chromosome")
plt.ylabel("-log10(P-value)")
plt.title("Manhattan Plot of GWAS Results")
plt.tight_layout()
plt.legend()

# Show plot
plt.show()

# --- Outcome --- #

# Set up the Manhattan plot -> Outcome
plt.figure(figsize=(14, 6))

# Define colors for alternating chromosomes
chromosomes = df["CHR"].unique()
colors = ["blue", "red"] * (len(chromosomes) // 2 + 1)

x_labels = []
x_ticks = []
x_offset = 0

# Sort chromosomes for consistent plotting
for i, chrom in enumerate(sorted(chromosomes)):
    subset = df2[df2["CHR"] == chrom]
    plt.scatter(subset["BP"] + x_offset, subset["-log10(PVALUE)"], color=colors[i % 2], s=10)
    x_labels.append(chrom)
    x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
    x_offset += subset["BP"].max() - subset["BP"].min() + 1

# Add genome-wide significance threshold line (-log10(5e-8) = 7.3)
plt.axhline(y=-np.log10(5e-8), color="black", linestyle="--", label="Genome-wide significance (5e-8)")

# Formatting
plt.xticks(x_ticks, x_labels, rotation=90)
plt.xlabel("Chromosome")
plt.ylabel("-log10(P-value)")
plt.title("Manhattan Plot of GWAS Results")
plt.tight_layout()
plt.legend()

# Show plot
plt.show()
