#!/usr/bin/env python3
# Visualisation of MR Results
# This script will be used in the main pipeline.
# ==== Initially done in Jupyter Notebook === #

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# --- Handle command-line arguments --- # 
if len(sys.argv) != 4:
    print("Usage: python3 05_visualisations.py <MR_Formatted_Results.csv> <MR_IVW_OR_Per_SNP.csv> <output_dir>")
    sys.exit(1)

summary_file = sys.argv[1]
snp_file = sys.argv[2]
output_dir = sys.argv[3]
os.makedirs(output_dir, exist_ok=True)

# =============================== #
# --- Plot 1: MR Summary Box ---- #
# =============================== #

# --- Load MR summary results ---
df = pd.read_csv(summary_file)
print(f"[INFO] Loaded MR summary: {df.shape}")

# --- Extract estimates and CIs --- #
methods = ["IVW", "Weighted Median", "Egger"]
ORs = [df["IVW_OR"][0], df["WM_OR"][0], df["Egger_OR"][0]]
Lower = [df["IVW_Lower_95"][0], df["WM_Lower_95"][0], df["Egger_Lower_95"][0]]
Upper = [df["IVW_Upper_95"][0], df["WM_Upper_95"][0], df["Egger_Upper_95"][0]]

fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
x = range(len(methods))
yerr = [[OR - low for OR, low in zip(ORs, Lower)],
        [up - OR for OR, up in zip(ORs, Upper)]]

ax.errorbar(x, ORs, yerr=yerr, fmt='o', color='black', capsize=5, markersize=8, linewidth=2)
ax.axhline(y=1, color='gray', linestyle='--', linewidth=1.5)

ax.set_xticks(x)
ax.set_xticklabels(methods, fontsize=12, fontweight='bold')
ax.set_ylabel("Odds Ratio (95% CI)", fontsize=13)
ax.set_title("MR Effect Estimates (Odds Ratios)", fontsize=16, fontweight='bold', pad=15)
ax.tick_params(axis='y', labelsize=11)
fig.tight_layout()

summary_plot_path = os.path.join(output_dir, "mr_summary_estimates.png")
plt.savefig(summary_plot_path, dpi=300)
plt.close()
print(f"[INFO] Saved: {summary_plot_path}")

# =============================== #
# --- Plot 2: IVW OR per SNP ---- #
# =============================== #

# --- Load per-SNP results --- #
snp_df = pd.read_csv(snp_file)
print(f"[INFO] Loaded SNP-level IVW results: {snp_df.shape[0]} SNPs")

# Save full table for reference
snp_df.to_csv(os.path.join(output_dir, "ivw_all_snp_ORs.csv"), index=False)

# --- Clean and subset top 30 SNPs --- #
snp_df["Lower_CI"] = snp_df["IVW_Lower_95"]
snp_df["Upper_CI"] = snp_df["IVW_Upper_95"]
top_snp_df = snp_df.nlargest(30, "IVW_OR").copy()
top_snp_df = top_snp_df.sort_values("IVW_OR", ascending=False).reset_index(drop=True)

# --- Forest plot of top 30 SNPs --- #
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

snp_plot_path = os.path.join(output_dir, "ivw_per_snp_forest_plot.png")
plt.savefig(snp_plot_path, dpi=300)
plt.close()
print(f"[INFO] Saved: {snp_plot_path}")

print("\n[INFO] Both visualisations generated successfully.")
