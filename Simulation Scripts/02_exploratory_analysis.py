# Exploratory Analysis for Simulated GWAS data (Optiomal)
# This script will not be inlcluded in the final pipeline.
# ==== Initially done in Jupyter Notebook === #

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

# --- Exposure --- #
# Set up the Manhattan plot -> Exposure
sns.set_style("whitegrid")
plt.figure(figsize=(16, 8), dpi=300)

# Define colors for alternating chromosomes
chromosomes = sorted(df["CHR"].unique())
colors = ["#1f77b4", "#d62728"] * (len(chromosomes) // 2 + 1)

x_labels = []
x_ticks = []
x_offset = 0

# Plot each chromosome separately 
for i, chrom in enumerate(chromosomes):
    subset = df[df["CHR"] == chrom]
    plt.scatter(
        subset["BP"] + x_offset, subset["-log10(PVALUE)"], 
        color=colors[i % 2], s=10, alpha=0.75, edgecolors="none"
    )
    x_labels.append(chrom)
    x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
    x_offset += subset["BP"].max() - subset["BP"].min() + 1

# Add genome-wide significance threshold line (-log10(5e-8) = 7.3)
threshold = -np.log10(5e-8)
plt.axhline(y=threshold, color="black", linestyle="dashed", linewidth=1.5, label="Genome-wide significance (5e-8)")

# Formatting 
plt.xticks(x_ticks, x_labels, rotation=90, fontsize=12, fontweight="bold")
plt.yticks(fontsize=12, fontweight="bold")
plt.xlabel("Chromosome", fontsize=14, fontweight="bold", labelpad=10)
plt.ylabel("-log10(P-value)", fontsize=14, fontweight="bold", labelpad=10)
plt.title("Manhattan Plot of GWAS Results", fontsize=18, fontweight="bold", pad=15)

# Remove top and right spines 
sns.despine()

# Add a legend
plt.legend(frameon=True, fontsize=12, loc="upper right", edgecolor="black")

# Optimise layout
plt.tight_layout()

# Save Figure
plt.savefig("exposure_mahattan.png", dpi=300)

# Show plot
plt.show()

# --- Outcome --- #
# Set up the Manhattan plot -> Outcome
sns.set_style("whitegrid")
plt.figure(figsize=(16, 8), dpi=300)

# Define colors for alternating chromosomes
chromosomes = sorted(df2["CHR"].unique())
colors = ["#1f77b4", "#d62728"] * (len(chromosomes) // 2 + 1)

x_labels = []
x_ticks = []
x_offset = 0

# Plot each chromosome separately 
for i, chrom in enumerate(chromosomes):
    subset = df2[df2["CHR"] == chrom]
    plt.scatter(
        subset["BP"] + x_offset, subset["-log10(PVALUE)"], 
        color=colors[i % 2], s=10, alpha=0.75, edgecolors="none"
    )
    x_labels.append(chrom)
    x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
    x_offset += subset["BP"].max() - subset["BP"].min() + 1

# Add genome-wide significance threshold line (-log10(5e-8) = 7.3)
threshold = -np.log10(5e-8)
plt.axhline(y=threshold, color="black", linestyle="dashed", linewidth=1.5, label="Genome-wide significance (5e-8)")

# Formatting 
plt.xticks(x_ticks, x_labels, rotation=90, fontsize=12, fontweight="bold")
plt.yticks(fontsize=12, fontweight="bold")
plt.xlabel("Chromosome", fontsize=14, fontweight="bold", labelpad=10)
plt.ylabel("-log10(P-value)", fontsize=14, fontweight="bold", labelpad=10)
plt.title("Manhattan Plot of GWAS Results", fontsize=18, fontweight="bold", pad=15)

# Remove top and right spines 
sns.despine()

# Add a legend
plt.legend(frameon=True, fontsize=12, loc="upper right", edgecolor="black")

# Optimise layout
plt.tight_layout()

# Save Figure
plt.savefig("outcome_mahattan.png", dpi=300)

# Show plot
plt.show()

print("Exploratory analyses complete!")
