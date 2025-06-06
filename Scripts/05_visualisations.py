#!/usr/bin/env python3
"""
MR-CoPe | Visualisation of Mendelian Randomisation Results
-----------------------------------------------------------
Author: Guillermo Comesaña & Christian Pepler
Date: 2025

Usage:
    python3 05_visualisations.py <MR_Formatted_Results.csv> <MR_IVW_OR_Per_SNP.csv> <output_dir>

Description:
- Generates method-level MR results plot (IVW, WM, Egger).
- Generates SNP-level forest plot of top 30 SNPs by IVW OR.
- Outputs all plots and full SNP table.
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt


def validate_inputs(paths, labels):
    for path, label in zip(paths, labels):
        if not os.path.isfile(path):
            print(f"❌ ERROR: {label} file not found at: {path}")
            sys.exit(1)


def safe_load_csv(path, label):
    try:
        df = pd.read_csv(path)
        if df.empty:
            print(f"⚠️ WARNING: {label} file is empty. Skipping visualisation.")
            sys.exit(0)
        return df
    except Exception as e:
        print(f"❌ ERROR loading {label}: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    results_file, snp_file, output_dir = sys.argv[1], sys.argv[2], sys.argv[3]

    print("\n==============================================================")
    print("MR-CoPe | Visualisation of MR Results")
    print("==============================================================\n")

    validate_inputs([results_file, snp_file], ["MR summary results", "SNP-level results"])

    os.makedirs(output_dir, exist_ok=True)

    # --- Load MR Summary Data --- #
    print("📥 Loading MR summary results...")
    df = safe_load_csv(results_file, "MR summary results")
    print(f"✅ Loaded: {df.shape}\n")

    required_cols = ["IVW_OR", "WM_OR", "Egger_OR"]
    if not all(col in df.columns for col in required_cols):
        print(f"❌ ERROR: Required columns missing from MR summary results. Found columns: {df.columns.tolist()}")
        sys.exit(1)

    methods = ["IVW", "Weighted Median", "Egger"]
    ORs = [df["IVW_OR"].iloc[0], df["WM_OR"].iloc[0], df["Egger_OR"].iloc[0]]
    Lower = [df["IVW_CI_Lower"].iloc[0], df["WM_CI_Lower"].iloc[0], df["Egger_CI_Lower"].iloc[0]]
    Upper = [df["IVW_CI_Upper"].iloc[0], df["WM_CI_Upper"].iloc[0], df["Egger_CI_Upper"].iloc[0]]

    # --- Plot Method-level Summary --- #
    print("📊 Creating method-level summary plot...")

    plt.figure(figsize=(10, 6), dpi=300)
    x = range(len(methods))
    yerr = [[OR - low for OR, low in zip(ORs, Lower)],
            [up - OR for OR, up in zip(ORs, Upper)]]

    plt.errorbar(x, ORs, yerr=yerr, fmt='o', color='black', capsize=5, markersize=8, linewidth=2)
    plt.axhline(y=1, linestyle='--', color='gray', linewidth=1.5)
    plt.xticks(x, methods, fontsize=12, fontweight='bold')
    plt.ylabel("Odds Ratio (95% CI)", fontsize=13)
    plt.title("MR Effect Estimates", fontsize=16, fontweight='bold', pad=15)
    plt.tight_layout()

    output_path1 = os.path.join(output_dir, "mr_summary_estimates.png")
    plt.savefig(output_path1, dpi=300)
    plt.close()
    print(f"✅ Saved: {output_path1}\n")

    # --- Load SNP-level IVW Results --- #
    print("📥 Loading SNP-level IVW results...")
    snp_df = safe_load_csv(snp_file, "SNP-level results")
    print(f"✅ Loaded SNPs: {snp_df.shape[0]}\n")

    snp_df["Lower_CI"] = snp_df["IVW_Lower_95"]
    snp_df["Upper_CI"] = snp_df["IVW_Upper_95"]

    snp_output = os.path.join(output_dir, "ivw_all_snp_ORs.csv")
    snp_df.to_csv(snp_output, index=False)
    print(f"✅ Full SNP results saved: {snp_output}\n")

    # --- Forest Plot of Top 30 SNPs --- #
    print("📊 Creating forest plot of top 30 SNPs...")

    top_snp_df = snp_df.nlargest(30, "IVW_OR").copy().sort_values("IVW_OR", ascending=False).reset_index(drop=True)

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

    output_path2 = os.path.join(output_dir, "ivw_per_snp_forest_plot.png")
    plt.savefig(output_path2, dpi=300)
    plt.close()
    print(f"✅ Saved: {output_path2}\n")

    print("🎉 Visualisation complete!")
    print(f"All outputs saved in: {output_dir}")
    print("==============================================================\n")


if __name__ == "__main__":
    main()
