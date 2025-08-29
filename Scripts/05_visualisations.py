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
import matplotlib as mpl
import matplotlib.pyplot as plt

# =========================
# Global visual aesthetics
# =========================
mpl.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 11,
    "font.family": "DejaVu Sans",
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "axes.linewidth": 1.0,
    "xtick.labelsize": 10,
    "ytick.labelsize": 9,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.color": "0.85",
    "grid.linestyle": "-",
    "grid.linewidth": 0.8,
})

# Subtle, colorblind-friendly palette
COL_POINT = "#2B6CB0"   # blue for method points
COL_CI    = "#4A5568"   # slate/gray for CIs
COL_REF   = "#A0AEC0"   # reference line
COL_POS   = "#2F855A"   # green (OR>1)
COL_NEG   = "#C53030"   # red (OR<1)
COL_BAND  = "#F7FAFC"   # banding in forest plot


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

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.grid(axis="y")

    x = range(len(methods))
    yerr = [
        [OR - low for OR, low in zip(ORs, Lower)],
        [up - OR for OR, up in zip(ORs, Upper)]
    ]

    # Draw CIs and then points (cleaner layering)
    ax.errorbar(x, ORs, yerr=yerr, fmt="none", ecolor=COL_CI, elinewidth=2, capsize=6, zorder=1, label="95% CI")
    ax.scatter(x, ORs, s=60, color=COL_POINT, zorder=2, label="Point estimate")

    # Reference line at null
    ax.axhline(y=1, linestyle='--', color=COL_REF, linewidth=1.2, zorder=0)

    ax.set_xticks(list(x), methods)
    ax.set_ylabel("Odds Ratio (95% CI)")
    ax.set_title("MR Effect Estimates")

    # Gentle y padding
    ymin = min(Lower) - 0.05
    ymax = max(Upper) + 0.05
    ax.set_ylim(ymin, ymax)

    # Minimal legend (since we keep labels off the figure)
    ax.legend(loc="upper left", ncols=2, frameon=False)

    plt.tight_layout()
    output_path1 = os.path.join(output_dir, "mr_summary_estimates.png")
    plt.savefig(output_path1, dpi=300, bbox_inches="tight")
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

    fig, ax = plt.subplots(figsize=(8, 0.35 * len(top_snp_df) + 2), dpi=300)

    y = list(range(len(top_snp_df)))

    # Alternating bands to guide the eye
    for i in y:
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=COL_BAND, zorder=0)

    # CI lines (horizontal) and points colored by direction vs null
    ax.errorbar(
        top_snp_df["IVW_OR"], y,
        xerr=[
            top_snp_df["IVW_OR"] - top_snp_df["Lower_CI"],
            top_snp_df["Upper_CI"] - top_snp_df["IVW_OR"]
        ],
        fmt='none', ecolor=COL_CI, elinewidth=2, capsize=3, zorder=1
    )

    colors = (top_snp_df["IVW_OR"] >= 1).map({True: COL_POS, False: COL_NEG}).values
    ax.scatter(top_snp_df["IVW_OR"], y, s=40, color=colors, zorder=2)

    # Reference line
    ax.axvline(x=1, linestyle="--", color=COL_REF, linewidth=1.2, zorder=0)

    ax.set_yticks(y, top_snp_df["SNP"])
    ax.set_xlabel("IVW Odds Ratio (95% CI)")
    ax.set_title("Top 30 SNPs by IVW OR")

    # Tidy x-limits with small padding
    xmin = float((top_snp_df["IVW_OR"] - (top_snp_df["IVW_OR"] - top_snp_df["Lower_CI"])).min())
    xmax = float((top_snp_df["IVW_OR"] + (top_snp_df["Upper_CI"] - top_snp_df["IVW_OR"])).max())
    xmin = min(top_snp_df["Lower_CI"].min(), xmin)
    xmax = max(top_snp_df["Upper_CI"].max(), xmax)
    xpad = (xmax - xmin) * 0.05 if xmax > xmin else 0.1
    ax.set_xlim(left=max(0.5, xmin - xpad), right=xmax + xpad)

    ax.tick_params(axis="y", pad=4)
    ax.grid(axis="x")

    plt.tight_layout()
    output_path2 = os.path.join(output_dir, "ivw_per_snp_forest_plot.png")
    plt.savefig(output_path2, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved: {output_path2}\n")

    print("🎉 Visualisation complete!")
    print(f"All outputs saved in: {output_dir}")
    print("==============================================================\n")


if __name__ == "__main__":
    main()
