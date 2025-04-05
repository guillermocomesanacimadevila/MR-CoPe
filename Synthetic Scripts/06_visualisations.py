#!/usr/bin/env python3
# Visualisation of MR Results
# This script will be used in the main pipeline.
# ==== Initially done in Jupyter Notebook === #

import sys
import pandas as pd
import matplotlib.pyplot as plt

# --- Handle command-line arguments ---
if len(sys.argv) != 4:
    print("Usage: python3 05_visualisations.py <MR_Formatted_Results.csv> <MR_IVW_OR_Per_SNP.csv> <output_dir>")
    sys.exit(1)

results_file = sys.argv[1]
snp_file = sys.argv[2]
output_dir = sys.argv[3]  # Not used for saving anymore (but kept for structure)

# ============================= #
# --- Load MR summary table --- #
# ============================= #
df = pd.read_csv(results_file)
print(f"[INFO] Loaded MR summary: {df.shape}")
print(df.head())

# --- Extract method-level values ---
IVW_OR = df["IVW_OR"].values[0]
WM_OR = df["WM_OR"].values[0]
Egger_OR = df["Egger_OR"].values[0]

IVW_lower = df["IVW_Lower_95"].values[0]
WM_lower = df["WM_Lower_95"].values[0]
Egger_lower = df["Egger_Lower_95"].values[0]

IVW_upper = df["IVW_Upper_95"].values[0]
WM_upper = df["WM_Upper_95"].values[0]
Egger_upper = df["Egger_Upper_95"].values[0]

# Build dataframe for plotting
methods = ["IVW", "Weighted Median", "Egger"]
ORs = [IVW_OR, WM_OR, Egger_OR]
Lower = [IVW_lower, WM_lower, Egger_lower]
Upper = [IVW_upper, WM_upper, Egger_upper]

# ================================ #
# --- Plot method-level summary --- #
# ================================ #
plt.figure(figsize=(10, 6), dpi=300)
x = range(len(methods))
yerr = [[OR - low for OR, low in zip(ORs, Lower)],
        [up - OR for OR, up in zip(ORs, Upper)]]

plt.errorbar(x, ORs, yerr=yerr, fmt='o', color='black', capsize=5, markersize=8, linewidth=2)
plt.axhline(y=1, linestyle='--', color='gray', linewidth=1.5)

plt.xticks(x, methods, fontsize=12, fontweight='bold')
plt.ylabel("Odds Ratio (95% CI)", fontsize=13)
plt.title("MR Effect Estimates (Odds Ratios)", fontsize=16, fontweight='bold', pad=15)
plt.tight_layout()

plt.savefig("mr_summary_estimates.png", dpi=300)
plt.close()
print("[INFO] Saved: mr_summary_estimates.png")

# =============================== #
# --- Load per-SNP IVW results --- #
# =============================== #
snp_df = pd.read_csv(snp_file)
print(f"[INFO] Loaded SNP-level IVW results: {snp_df.shape[0]} SNPs")

# Save full table for inspection
snp_df.to_csv("ivw_all_snp_ORs.csv", index=False)

# ========================== #
# --- Plot top 30 forest --- #
# ========================== #
snp_df["Lower_CI"] = snp_df["IVW_Lower_95"]
snp_df["Upper_CI"] = snp_df["IVW_Upper_95"]
top_snp_df = snp_df.nlargest(30, "IVW_OR").copy()
top_snp_df = top_snp_df.sort_values("IVW_OR", ascending=False).reset_index(drop=True)

plt.figure(figsize=(8, 0.35 * len(top_snp_df) + 2), dpi=300)
y = range(len(top_snp_df))

plt.errorbar(
    top_snp_df["IVW_OR"], y,
    xerr=[
        top_snp_df["IVW_OR"] - top_snp_df["Lower_CI"],
        top_snp_df["Upper_CI"] - top_snp_df["IVW_OR"]
    ],
    fmt='o', color='black', ecolor='gray', elinewidth=1.5, capsize=3
)

plt.axvline(x=1, linestyle="--", color="gray")
plt.yticks(y, top_snp_df["SNP"], fontsize=8)
plt.xlabel("IVW Odds Ratio (95% CI)", fontsize=12)
plt.title("Top 30 SNPs by IVW OR", fontsize=14, fontweight="bold")
plt.tight_layout()

plt.savefig("ivw_per_snp_forest_plot.png", dpi=300)
plt.close()
print("[INFO] Saved: ivw_per_snp_forest_plot.png")

print("\n[INFO] Both visualisations generated successfully.")