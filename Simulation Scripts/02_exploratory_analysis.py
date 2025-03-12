# Exploratory Analysis for Simulated GWAS data (Optiomal)
# This script will not be inlcluded in the final pipeline.
# ==== Initially done in Jupyter Notebook === #

import os
import pandas as pd
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

# Check for pval < 0.05
plt.figure(figsize=(12, 6))

# Exposure GWAS
# Scatter plot with chromosome-based coloring
chromosomes = df["CHR"].unique()
colors = ["blue", "red"] * (len(chromosomes) // 2 + 1)

x_labels = []
x_ticks = []
x_offset = 0

for i, chrom in enumerate(sorted(chromosomes)):
    subset = df[df["CHR"] == chrom]
    plt.scatter(subset["BP"] + x_offset, subset["PVALUE"], color=colors[i % 2], s=10)
    x_labels.append(chrom)
    x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
    x_offset += subset["BP"].max() - subset["BP"].min() + 1

plt.xticks(x_ticks, x_labels, rotation=90)
plt.axhline(y=0.05, color="black", linestyle="--", label="P = 0.05")
plt.title("Manhattan Plot (P-values)")
plt.xlabel("Chromosome")
plt.ylabel("P-value")
plt.tight_layout()

plt.show()

# Outcome GWAS
plt.figure(figsize=(12, 6))

chromosomes = df2["CHR"].unique()
colors = ["blue", "red"] * (len(chromosomes) // 2 + 1)

x_labels = []
x_ticks = []
x_offset = 0

for i, chrom in enumerate(sorted(chromosomes)):
    subset = df2[df2["CHR"] == chrom]
    plt.scatter(subset["BP"] + x_offset, subset["PVALUE"], color=colors[i % 2], s=10)
    x_labels.append(chrom)
    x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
    x_offset += subset["BP"].max() - subset["BP"].min() + 1

plt.xticks(x_ticks, x_labels, rotation=90)
plt.axhline(y=0.05, color="black", linestyle="--", label="P = 0.05")
plt.title("Manhattan Plot (P-values)")
plt.xlabel("Chromosome")
plt.ylabel("P-value")
plt.tight_layout()

plt.show()
