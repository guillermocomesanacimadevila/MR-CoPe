#!/usr/bin/env python3

"""
MR-CoPe | Exploratory Analysis of GWAS Summary Statistics
---------------------------------------------------------
Author: Guillermo Comesaña & Christian Pepler
Date: 2025

Description:
Performs QC checks and generates Manhattan plots for exposure and outcome
GWAS summary statistics provided by the user.

Usage:
    python3 01_exploratory_analysis.py <exposure_gwas.csv> <outcome_gwas.csv> <output_dir> <log10_flag>

Inputs:
    - exposure_gwas.csv : Exposure GWAS summary statistics
    - outcome_gwas.csv  : Outcome GWAS summary statistics
    - output_dir        : Directory to save output files
    - log10_flag        : "y" if p-values are already -log10(p), else "n"
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy import stats

# ----------------- Auto Column Mapping ----------------- #

AUTO_COLUMN_MAP = {
    'SNP': ['snp', 'rsid', 'marker', 'rs_number', 'rsids'],
    'CHR': ['chr', 'chromosome', 'chrom'],
    'BP': ['bp', 'position', 'pos'],
    'A1': ['a1', 'effect_allele', 'ea', 'alt'],
    'A2': ['a2', 'other_allele', 'oa', 'ref'],
    'BETA': ['beta', 'effect_size', 'b'],
    'SE': ['se', 'stderr', 'standard_error'],
    'PVALUE': ['pval', 'p_value', 'p'],
    'EAF': ['eaf', 'effect_allele_freq', 'freq', 'riskfrequency', 'maf']
}

def auto_map_columns(df, label):
    mapping = {}
    df_cols_lower = [c.lower() for c in df.columns]

    print(f"🔎 Auto-mapping columns for {label} GWAS...")

    for standard_name, possible_names in AUTO_COLUMN_MAP.items():
        for possible in possible_names:
            if possible in df_cols_lower:
                matched_col = df.columns[df_cols_lower.index(possible)]
                mapping[standard_name] = matched_col
                break
        if standard_name not in mapping:
            print(f"⚠️ WARNING: Column for {standard_name} not found in {label} GWAS")

    print(f"✅ Column mapping for {label}:")
    for k, v in mapping.items():
        print(f"  {k} --> {v}")
    return mapping

def parse_custom_gwas(df, label):
    print(f"⚙️ Parsing custom GWAS structure for {label}...")
    df["SNP"] = df["riskAllele"].str.split("-").str[0]
    df["CHR"] = df["locations"].str.split(":").str[0]
    raw_bp = df["locations"].str.split(":").str[1]

    if raw_bp.str.contains(r"\D", regex=True).any():
        print("🧹 Detected non-numeric BP values — extracting digits...")
        df["BP"] = raw_bp.str.extract(r"(\d+)")[0]
    else:
        df["BP"] = raw_bp

    df["BP"] = pd.to_numeric(df["BP"], errors='coerce')
    df["PVALUE"] = pd.to_numeric(df["pValue"], errors='coerce')
    df.dropna(subset=["SNP", "CHR", "BP", "PVALUE"], inplace=True)
    df["CHR"] = df["CHR"].astype(str)
    df["BP"] = df["BP"].astype(int)
    print(f"✅ Custom parsing complete for {label}.\n")
    return df

def load_gwas(path, label):
    print(f"📅 Loading {label} GWAS: {path}")
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = df.columns.str.strip()
    print(f"✅ Loaded {label} GWAS | Shape: {df.shape}\n")

    mapping = auto_map_columns(df, label)

    if not all(k in mapping for k in ["SNP", "CHR", "BP", "PVALUE"]):
        if {"riskAllele", "locations", "pValue"}.issubset(df.columns):
            df = parse_custom_gwas(df, label)
        else:
            print(f"❌ ERROR: Missing critical columns in {label} GWAS and no fallback possible.")
            sys.exit(1)
    else:
        for std_col, actual_col in mapping.items():
            if std_col != actual_col:
                df.rename(columns={actual_col: std_col}, inplace=True)

    if "BETA" in df.columns and "SE" in df.columns and "PVALUE" not in df.columns:
        print("🧠 Computing PVALUE from BETA and SE...")
        df["Z"] = df["BETA"] / df["SE"]
        df["PVALUE"] = 2 * norm.sf(np.abs(df["Z"]))

    return df

def qc_report(df, label):
    print(f"🔎 QC Summary: {label} GWAS")
    print("-" * 40)
    print(df.info())
    print("\nMissing values per column:")
    print(df.isna().sum())
    print("\n")

def manhattan_plot(df, output_path, title, is_log10_input=False):
    if "PVALUE" not in df.columns:
        print(f"❌ Cannot plot Manhattan — PVALUE column missing in {title}.")
        return

    df = df.dropna(subset=["CHR", "BP", "PVALUE"])
    df = df[df["PVALUE"] > 0]

    if not is_log10_input:
        df["-log10(PVALUE)"] = -np.log10(df["PVALUE"])
    else:
        df["-log10(PVALUE)"] = df["PVALUE"]

    df["CHR"] = df["CHR"].astype(str)
    df = df.sort_values(["CHR", "BP"])
    df["ind"] = range(len(df))
    df_grouped = df.groupby("CHR")

    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)
    colors = ["#4C72B0", "#55A868"]
    x_labels, x_labels_pos = [], []

    for i, (chrom, group) in enumerate(df_grouped):
        ax.scatter(
            group["ind"], group["-log10(PVALUE)"],
            color=colors[i % 2],
            s=6, alpha=0.8, edgecolor='none'
        )
        mid_pos = (group["ind"].min() + group["ind"].max()) / 2
        x_labels.append(chrom)
        x_labels_pos.append(mid_pos)

    ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', linewidth=1)
    ax.axhline(y=-np.log10(1e-5), color='orange', linestyle='--', linewidth=1)

    ymax = df["-log10(PVALUE)"].max()
    ax.set_ylim([0, ymax + 0.1 * ymax])
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, fontsize=6)
    ax.set_xlim([0, len(df)])
    ax.set_xlabel("Chromosome", fontsize=10)
    ax.set_ylabel("-log10(p)", fontsize=10)
    ax.set_title(title, fontsize=12, weight='bold', pad=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='y', labelsize=7)
    ax.tick_params(axis='x', labelsize=6)
    ax.legend().set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"✅ Manhattan plot saved: {output_path}\n")

def qq_plot(df, output_path, title):
    if "PVALUE" not in df.columns:
        print(f"❌ Cannot plot Q-Q — PVALUE column missing in {title}.")
        return

    pvals = df["PVALUE"].dropna()
    pvals = np.clip(pvals, 1e-300, 1.0)
    n = len(pvals)
    expected = -np.log10(np.linspace(1 / (n + 1), 1, n))
    observed = -np.log10(np.sort(pvals))

    quantiles = np.arange(1, n + 1) / (n + 1)
    ci_low = -np.log10(stats.beta.ppf(0.025, quantiles * n, (1 - quantiles) * n))
    ci_high = -np.log10(stats.beta.ppf(0.975, quantiles * n, (1 - quantiles) * n))

    chisq = stats.chi2.isf(pvals, df=1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)

    cmap = plt.cm.viridis_r
    colors = cmap(np.linspace(0.1, 0.85, n))

    plt.figure(figsize=(7.5, 7.5), dpi=400)
    plt.fill_between(expected, ci_low, ci_high, color="#EAEAF2", alpha=0.6, label="95% CI")
    plt.scatter(expected, observed, c=colors, s=10, edgecolor='none', alpha=0.8)
    plt.plot([0, max(expected)], [0, max(expected)], linestyle="--", color="black", linewidth=1.3)

    plt.xlabel(r"Expected $-\log_{10}$(P)", fontsize=13, weight="bold")
    plt.ylabel(r"Observed $-\log_{10}$(P)", fontsize=13, weight="bold")
    plt.title(title, fontsize=16, weight="bold", pad=10)

    plt.text(0.95, 0.05, f"$\\lambda_{{GC}}$ = {lambda_gc:.3f}",
             transform=plt.gca().transAxes, fontsize=11,
             ha="right", va="bottom",
             bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'))

    plt.grid(True, linestyle=':', linewidth=0.5, alpha=0.6)
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.legend(frameon=False, fontsize=10, loc="upper left")

    plt.tight_layout()
    plt.savefig(output_path, dpi=400)
    plt.close()
    print(f"✅ Q-Q plot saved: {output_path}\n")

def validate_inputs(exposure_path, outcome_path, output_dir):
    if not os.path.isfile(exposure_path):
        print(f"❌ ERROR: Exposure file not found: {exposure_path}")
        sys.exit(1)
    if not os.path.isfile(outcome_path):
        print(f"❌ ERROR: Outcome file not found: {outcome_path}")
        sys.exit(1)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)

    exposure_path, outcome_path, output_dir, log10_flag_input = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    is_log10_input = log10_flag_input.lower().startswith("y")

    print("\n" + "=" * 60)
    print("MR-CoPe | Exploratory Analysis of GWAS Summary Statistics")
    print("=" * 60 + "\n")

    validate_inputs(exposure_path, outcome_path, output_dir)

    exposure = load_gwas(exposure_path, "Exposure")
    outcome = load_gwas(outcome_path, "Outcome")

    qc_report(exposure, "Exposure")
    qc_report(outcome, "Outcome")

    print("📊 Generating Manhattan and Q-Q plots...\n")
    manhattan_plot(exposure, os.path.join(output_dir, "exposure_manhattan.png"), "Exposure GWAS Manhattan Plot", is_log10_input)
    manhattan_plot(outcome, os.path.join(output_dir, "outcome_manhattan.png"), "Outcome GWAS Manhattan Plot", is_log10_input)

    qq_plot(exposure, os.path.join(output_dir, "exposure_qq.png"), "Exposure GWAS Q-Q Plot")
    qq_plot(outcome, os.path.join(output_dir, "outcome_qq.png"), "Outcome GWAS Q-Q Plot")

    print("🎉 Exploratory analysis completed successfully!")
    print(f"All outputs saved in: {output_dir}\n")
    print("=" * 60)

if __name__ == "__main__":
    main()
