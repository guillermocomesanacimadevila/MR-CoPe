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

# ----------------- Auto Column Mapping ----------------- #

AUTO_COLUMN_MAP = {
    'SNP': ['snp', 'rsid', 'marker', 'rs_number'],
    'CHR': ['chr', 'chromosome'],
    'BP': ['bp', 'position', 'pos'],
    'A1': ['a1', 'effect_allele', 'ea'],
    'A2': ['a2', 'other_allele', 'oa'],
    'BETA': ['beta', 'effect_size', 'b'],
    'SE': ['se', 'stderr', 'standard_error'],
    'PVALUE': ['pval', 'p_value', 'p'],
    'EAF': ['eaf', 'effect_allele_freq', 'freq', 'riskfrequency']
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
    print(f"📥 Loading {label} GWAS: {path}")
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = df.columns.str.strip()
    print(f"✅ Loaded {label} GWAS | Shape: {df.shape}\n")

    mapping = auto_map_columns(df, label)

    # If key columns missing → custom parsing
    required = ["SNP", "CHR", "BP", "PVALUE"]
    if not all(k in mapping for k in required):
        if {"riskAllele", "locations", "pValue"}.issubset(df.columns):
            df = parse_custom_gwas(df, label)
        else:
            print(f"❌ ERROR: Missing critical columns in {label} GWAS and no fallback possible.")
            sys.exit(1)
    else:
        df.rename(columns=mapping, inplace=True)

    return df

def qc_report(df, label):
    print(f"🔎 QC Summary: {label} GWAS")
    print("-" * 40)
    print(df.info())
    print("\nMissing values per column:")
    print(df.isna().sum())
    print("\n")

def manhattan_plot(df, output_path, title):
    df["-log10(PVALUE)"] = -np.log10(df["PVALUE"])
    chromosomes = sorted(df["CHR"].unique())
    colors = ["#1f77b4", "#d62728"] * (len(chromosomes) // 2 + 1)

    plt.figure(figsize=(16, 8), dpi=300)
    x_labels, x_ticks, x_offset = [], [], 0

    for i, chrom in enumerate(chromosomes):
        subset = df[df["CHR"] == chrom]
        plt.scatter(subset["BP"] + x_offset, subset["-log10(PVALUE)"],
                    color=colors[i % 2], s=10, alpha=0.75, edgecolors="none")
        x_labels.append(chrom)
        x_ticks.append(x_offset + (subset["BP"].max() - subset["BP"].min()) / 2)
        x_offset += subset["BP"].max() - subset["BP"].min() + 1

    plt.axhline(y=-np.log10(5e-8), color="black", linestyle="dashed", linewidth=1.5,
                label="Genome-wide significance (5e-8)")

    plt.xticks(x_ticks, x_labels, rotation=90, fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Chromosome", fontsize=14, fontweight="bold")
    plt.ylabel("-log10(P-value)", fontsize=14, fontweight="bold")
    plt.title(title, fontsize=18, fontweight="bold")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"✅ Manhattan plot saved: {output_path}\n")

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

    print("📊 Generating Manhattan plots...\n")
    manhattan_plot(exposure, os.path.join(output_dir, "exposure_manhattan.png"), "Exposure GWAS Manhattan Plot")
    manhattan_plot(outcome, os.path.join(output_dir, "outcome_manhattan.png"), "Outcome GWAS Manhattan Plot")

    print("🎉 Exploratory analysis completed successfully!")
    print(f"All outputs saved in: {output_dir}\n")
    print("=" * 60)

if __name__ == "__main__":
    main()
