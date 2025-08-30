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
- ALSO generates F-stat histogram + density (from ld_pruned_SNPs.csv)
  and a funnel plot (from harmonised_data.csv).
- Outputs all plots and full SNP table.
"""

import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Optional, for KDE and nicer hist aesthetics
try:
    import seaborn as sns
    _HAS_SNS = True
except Exception:
    _HAS_SNS = False

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
COL_POINT = "#2B6CB0"   # blue for points / lines
COL_CI    = "#4A5568"   # slate/gray for CIs
COL_REF   = "#A0AEC0"   # reference line
COL_POS   = "#2F855A"   # green (OR>1)
COL_NEG   = "#C53030"   # red (OR<1)
COL_BAND  = "#F7FAFC"   # banding in forest plot

# Extra palette for F-stats/funnel
COL_BLUE   = "#4C72B0"  # hist, IVW line, points
COL_GREEN  = "#55A868"  # density fill
COL_RED    = "#D62728"  # threshold/null line
GRID_COLOR = "#E6E8EC"


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


def save_all(fig, outdir: Path, stem: str):
    """Save a figure as PNG, PDF, and SVG with white background."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.png", facecolor="white", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", facecolor="white", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.svg", facecolor="white", bbox_inches="tight")


def method_summary_plot(df: pd.DataFrame, output_dir: Path):
    required_cols = ["IVW_OR", "WM_OR", "Egger_OR",
                     "IVW_CI_Lower", "IVW_CI_Upper",
                     "WM_CI_Lower", "WM_CI_Upper",
                     "Egger_CI_Lower", "Egger_CI_Upper"]
    if not all(col in df.columns for col in required_cols):
        print(f"❌ ERROR: Required columns missing from MR summary results. Found columns: {df.columns.tolist()}")
        sys.exit(1)

    methods = ["IVW", "Weighted Median", "Egger"]
    ORs     = [df["IVW_OR"].iloc[0], df["WM_OR"].iloc[0], df["Egger_OR"].iloc[0]]
    Lower   = [df["IVW_CI_Lower"].iloc[0], df["WM_CI_Lower"].iloc[0], df["Egger_CI_Lower"].iloc[0]]
    Upper   = [df["IVW_CI_Upper"].iloc[0], df["WM_CI_Upper"].iloc[0], df["Egger_CI_Upper"].iloc[0]]

    print("📊 Creating method-level summary plot...")
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.grid(axis="y")

    x = range(len(methods))
    yerr = [
        [OR - low for OR, low in zip(ORs, Lower)],
        [up - OR for OR, up in zip(ORs, Upper)]
    ]

    ax.errorbar(x, ORs, yerr=yerr, fmt="none", ecolor=COL_CI, elinewidth=2, capsize=6, zorder=1, label="95% CI")
    ax.scatter(x, ORs, s=60, color=COL_POINT, zorder=2, label="Point estimate")
    ax.axhline(y=1, linestyle='--', color=COL_REF, linewidth=1.2, zorder=0)

    ax.set_xticks(list(x), methods)
    ax.set_ylabel("Odds Ratio (95% CI)")
    ax.set_title("MR Effect Estimates")

    ymin = min(Lower) - 0.05
    ymax = max(Upper) + 0.05
    ax.set_ylim(ymin, ymax)

    ax.legend(loc="upper left", ncols=2, frameon=False)
    plt.tight_layout()
    save_all(fig, output_dir, "mr_summary_estimates")
    plt.close(fig)
    print(f"✅ Saved: {output_dir/'mr_summary_estimates.png'}\n")


def snp_forest_plot(snp_df: pd.DataFrame, output_dir: Path):
    print("📊 Creating forest plot of top 30 SNPs...")

    snp_df = snp_df.copy()
    snp_df["Lower_CI"] = snp_df["IVW_Lower_95"]
    snp_df["Upper_CI"] = snp_df["IVW_Upper_95"]

    (output_dir / "ivw_all_snp_ORs.csv").write_text(snp_df.to_csv(index=False))
    print(f"✅ Full SNP results saved: {output_dir/'ivw_all_snp_ORs.csv'}\n")

    top_snp_df = snp_df.nlargest(30, "IVW_OR").copy().sort_values("IVW_OR", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 0.35 * len(top_snp_df) + 2), dpi=300)

    y = list(range(len(top_snp_df)))

    for i in y:
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=COL_BAND, zorder=0)

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

    ax.axvline(x=1, linestyle="--", color=COL_REF, linewidth=1.2, zorder=0)
    ax.set_yticks(y, top_snp_df["SNP"])
    ax.set_xlabel("IVW Odds Ratio (95% CI)")
    ax.set_title("Top 30 SNPs by IVW OR")

    xmin = float(min(top_snp_df["Lower_CI"]))
    xmax = float(max(top_snp_df["Upper_CI"]))
    xpad = (xmax - xmin) * 0.05 if xmax > xmin else 0.1
    ax.set_xlim(left=max(0.5, xmin - xpad), right=xmax + xpad)

    ax.tick_params(axis="y", pad=4)
    ax.grid(axis="x")

    plt.tight_layout()
    save_all(fig, output_dir, "ivw_per_snp_forest_plot")
    plt.close(fig)
    print(f"✅ Saved: {output_dir/'ivw_per_snp_forest_plot.png'}\n")


def fstat_plots(output_dir: Path):
    """
    Read ld_pruned_SNPs.csv from output_dir and make:
      - Fstat_Histogram.(png|pdf|svg)
      - Fstat_Density.(png|pdf|svg)
    """
    ld_path = output_dir / "ld_pruned_SNPs.csv"
    if not ld_path.exists():
        print("⚠️ ld_pruned_SNPs.csv not found — skipping F-stat plots.")
        return

    ld = pd.read_csv(ld_path)
    if "F_STAT" not in ld.columns:
        print("⚠️ 'F_STAT' column missing in ld_pruned_SNPs.csv — skipping F-stat plots.")
        return

    if _HAS_SNS:
        sns.set_theme(context="talk", style="whitegrid")
    print("\nSummary of F-statistics:\n", ld["F_STAT"].describe())

    # X-axis cap to avoid long tails squashing the bulk
    fmax = ld["F_STAT"].quantile(0.995)
    xmax = float(np.ceil(max(10, fmax) / 5) * 5)

    # Histogram
    bins = np.arange(0, max(xmax, ld["F_STAT"].max()) + 5, 5)
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    if _HAS_SNS:
        import seaborn as sns  # local
        sns.histplot(ld["F_STAT"], bins=bins, color=COL_BLUE, edgecolor="white", alpha=0.8, ax=ax)
    else:
        ax.hist(ld["F_STAT"], bins=bins, color=COL_BLUE, edgecolor="white", alpha=0.8)
    ax.axvline(10, color=COL_RED, linestyle="--", linewidth=1.2)
    ax.set_title("Distribution of F-statistics for Instruments")
    ax.set_xlabel("F-statistic")
    ax.set_ylabel("Count")
    ax.set_xlim(0, xmax)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    plt.tight_layout()
    save_all(fig, output_dir, "Fstat_Histogram")
    plt.close(fig)

    # Density
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    if _HAS_SNS:
        sns.kdeplot(ld["F_STAT"], fill=True, color=COL_GREEN, alpha=0.7, linewidth=1.5, ax=ax)
    else:
        # Fallback to a simple KDE via numpy (rough)
        from scipy.stats import gaussian_kde
        xs = np.linspace(0, xmax, 512)
        kde = gaussian_kde(ld["F_STAT"].dropna())
        ax.fill_between(xs, kde(xs), color=COL_GREEN, alpha=0.7)
        ax.plot(xs, kde(xs), color="black", linewidth=1.0)
    ax.axvline(10, color=COL_RED, linestyle="--", linewidth=1.2)
    ax.set_title("Density of F-statistics for Instruments")
    ax.set_xlabel("F-statistic")
    ax.set_ylabel("Density")
    ax.set_xlim(0, xmax)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    plt.tight_layout()
    save_all(fig, output_dir, "Fstat_Density")
    plt.close(fig)

    print(f"✅ Saved: {output_dir/'Fstat_Histogram.png'} and {output_dir/'Fstat_Density.png'}")


def funnel_plot(output_dir: Path):
    """
    Read harmonised_data.csv from output_dir and make:
      - Funnel_Plot.(png|pdf|svg)
    """
    harm_path = output_dir / "harmonised_data.csv"
    if not harm_path.exists():
        print("⚠️ harmonised_data.csv not found — skipping funnel plot.")
        return

    harm = pd.read_csv(harm_path)

    # Flexible column detection (TwoSampleMR/MR-CoPe styles)
    beta_exp_col = "BETA_EXP" if "BETA_EXP" in harm.columns else "beta.exposure"
    beta_out_col = "BETA_OUT" if "BETA_OUT" in harm.columns else "beta.outcome"
    se_out_col   = "SE_OUT"   if "SE_OUT"   in harm.columns else "se.outcome"

    if not all(col in harm.columns for col in [beta_exp_col, beta_out_col, se_out_col]):
        print("⚠️ harmonised_data.csv missing required columns — skipping funnel plot.")
        return

    harm = harm.dropna(subset=[beta_exp_col, beta_out_col, se_out_col])
    harm = harm[harm[beta_exp_col] != 0]

    if harm.empty:
        print("⚠️ No valid rows for funnel plot after filtering — skipping.")
        return

    # Wald ratio per SNP & corresponding SE on ratio scale
    harm["wald_ratio"] = harm[beta_out_col] / harm[beta_exp_col]
    harm["se_ratio"]   = harm[se_out_col] / harm[beta_exp_col].abs()

    # IVW estimate on raw beta scale
    w  = 1.0 / (harm[se_out_col] ** 2)
    bx = harm[beta_exp_col].values
    by = harm[beta_out_col].values
    beta_ivw = float((w * bx * by).sum() / (w * (bx ** 2)).sum())

    # Build funnel cones: se = |wald - beta_ivw| / z
    z = 1.96
    x_min = float(harm["wald_ratio"].quantile(0.005))
    x_max = float(harm["wald_ratio"].quantile(0.995))
    x_vals = np.linspace(x_min, x_max, 400)
    cone = np.abs(x_vals - beta_ivw) / z

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    # Cone shading (95% pseudo-CI around IVW line)
    ax.fill_between(x_vals, 0, cone, color="#B0B7C3", alpha=0.25, linewidth=0)

    # Scatter
        # Scatter
    ax.scatter(
        harm["wald_ratio"], harm["se_ratio"],
        s=40,
        color=COL_BLUE,
        edgecolor="black",
        linewidth=0.4,
        alpha=0.8
    )

    # Reference & IVW lines
    ax.axvline(0,        color=COL_RED,  linestyle="--", linewidth=1.2, label="Null (0)")
    ax.axvline(beta_ivw, color=COL_BLUE, linestyle="-",  linewidth=1.4, label=f"IVW = {beta_ivw:.3f}")

    # Invert y so higher precision (smaller SE) is at the top
    ax.invert_yaxis()

    ax.set_title("Funnel Plot of SNP Effects")
    ax.set_xlabel("Wald Ratio Estimate (βY / βX)")
    ax.set_ylabel("Standard Error")
    ax.set_xlim(x_min, x_max)

    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.legend(loc="upper right")

    plt.tight_layout()
    save_all(fig, output_dir, "Funnel_Plot")
    plt.close(fig)

    print(f"✅ Saved: {output_dir/'Funnel_Plot.png'}")
    return

def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    results_file, snp_file, out_dir_str = sys.argv[1], sys.argv[2], sys.argv[3]
    output_dir = Path(out_dir_str)

    print("\n==============================================================")
    print("MR-CoPe | Visualisation of MR Results")
    print("==============================================================\n")

    # Validate inputs for the core plots
    validate_inputs([results_file, snp_file], ["MR summary results", "SNP-level results"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load and plot method-level summary --- #
    print("📥 Loading MR summary results...")
    df = safe_load_csv(results_file, "MR summary results")
    print(f"✅ Loaded: {df.shape}\n")
    method_summary_plot(df, output_dir)

    # --- Load and plot SNP-level forest --- #
    print("📥 Loading SNP-level IVW results...")
    snp_df = safe_load_csv(snp_file, "SNP-level results")
    print(f"✅ Loaded SNPs: {snp_df.shape[0]}\n")
    snp_forest_plot(snp_df, output_dir)

    # --- Additional visualisations from results dir --- #
    # F-statistics (ld_pruned_SNPs.csv in the same results folder)
    fstat_plots(output_dir)

    # Funnel plot (harmonised_data.csv in the same results folder)
    funnel_plot(output_dir)

    print("🎉 Visualisation complete!")
    print(f"All outputs saved in: {output_dir}")
    print("==============================================================\n")


if __name__ == "__main__":
    main()
