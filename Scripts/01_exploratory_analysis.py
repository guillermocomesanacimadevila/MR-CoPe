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
    python3 01_exploratory_analysis.py <exposure_gwas.csv> <outcome_gwas.csv> <output_dir>

Inputs:
    - exposure_gwas.csv : Exposure GWAS summary statistics
    - outcome_gwas.csv  : Outcome GWAS summary statistics
    - output_dir        : Directory to save output files
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

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

def manhattan_plot(df, output_path, title):
    if "PVALUE" not in df.columns:
        print(f"❌ Cannot plot Manhattan — PVALUE column missing in {title}.")
        return

    df = df.dropna(subset=["CHR", "BP", "PVALUE"])
    df["-log10(PVALUE)"] = -np.log10(df["PVALUE"])
    df["CHR"] = df["CHR"].astype(str)

    df = df.sort_values(["CHR", "BP"])
    df["ind"] = range(len(df))
    df_grouped = df.groupby("CHR")

    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)
    colors = ["#4C72B0", "#55A868"]
    x_labels = []
    x_labels_pos = []

    for i, (chrom, group) in enumerate(df_grouped):
        ax.scatter(
            group["ind"], group["-log10(PVALUE)"],
            color=colors[i % 2],
            s=6, alpha=0.8, edgecolor='none'
        )
        mid_pos = (group["ind"].min() + group["ind"].max()) / 2
        x_labels.append(chrom)
        x_labels_pos.append(mid_pos)

    # Significance thresholds
    genomewide_sig = 5e-8
    suggestive_sig = 1e-5
    ax.axhline(y=-np.log10(genomewide_sig), color='red', linestyle='--', linewidth=1)
    ax.axhline(y=-np.log10(suggestive_sig), color='orange', linestyle='--', linewidth=1)

    ax.text(df["ind"].max() * 0.99, -np.log10(genomewide_sig) + 0.2,
            "Genome-wide sig (5e-8)", ha='right', fontsize=7, color='red')
    ax.text(df["ind"].max() * 0.99, -np.log10(suggestive_sig) + 0.2,
            "Suggestive (1e-5)", ha='right', fontsize=7, color='orange')

    # Axes formatting
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, fontsize=6)
    ax.set_xlim([0, len(df)])

    ymax = df["-log10(PVALUE)"].max()
    ax.set_ylim([0, ymax + 0.1 * ymax])

    ax.set_xlabel("Chromosome", fontsize=10)
    ax.set_ylabel("-log10(p)", fontsize=10)
    ax.set_title(title, fontsize=12, weight='bold', pad=15)

    # Aesthetic cleanup
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
    n = len(pvals)
    expected = -np.log10(np.linspace(1 / (n + 1), 1, n))
    observed = -np.log10(np.sort(pvals))

    plt.figure(figsize=(8, 8), dpi=300)
    plt.plot(expected, observed, marker='o', linestyle='none', markersize=3, alpha=0.6)
    plt.plot([0, max(expected)], [0, max(expected)], linestyle='--', color='red', linewidth=1.5)

    plt.xlabel("Expected -log10(P)", fontsize=14, fontweight="bold")
    plt.ylabel("Observed -log10(P)", fontsize=14, fontweight="bold")
    plt.title(title, fontsize=16, fontweight="bold")
    plt.grid(True, linestyle=':', linewidth=0.7)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
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
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    exposure_path, outcome_path, output_dir = sys.argv[1], sys.argv[2], sys.argv[3]

    print("\n" + "=" * 60)
    print("MR-CoPe | Exploratory Analysis of GWAS Summary Statistics")
    print("=" * 60 + "\n")

    validate_inputs(exposure_path, outcome_path, output_dir)

    exposure = load_gwas(exposure_path, "Exposure")
    outcome = load_gwas(outcome_path, "Outcome")

    qc_report(exposure, "Exposure")
    qc_report(outcome, "Outcome")

    print("📊 Generating Manhattan and Q-Q plots...\n")
    manhattan_plot(exposure, os.path.join(output_dir, "exposure_manhattan.png"), "Exposure GWAS Manhattan Plot")
    manhattan_plot(outcome, os.path.join(output_dir, "outcome_manhattan.png"), "Outcome GWAS Manhattan Plot")

    qq_plot(exposure, os.path.join(output_dir, "exposure_qq.png"), "Exposure GWAS Q-Q Plot")
    qq_plot(outcome, os.path.join(output_dir, "outcome_qq.png"), "Outcome GWAS Q-Q Plot")

    print("🎉 Exploratory analysis completed successfully!")
    print(f"All outputs saved in: {output_dir}\n")
    print("=" * 60)

if __name__ == "__main__":
    main()
