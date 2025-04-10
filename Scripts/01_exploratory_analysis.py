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


def print_header():
    print("\n" + "="*60)
    print("MR-CoPe | GWAS Exploratory Analysis")
    print("="*60 + "\n")


def validate_inputs(exposure_path, outcome_path, output_dir):
    for path, label in zip([exposure_path, outcome_path], ["Exposure", "Outcome"]):
        if not os.path.isfile(path):
            print(f"❌ ERROR: {label} file not found at: {path}")
            sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    print(f"📁 Output directory created (if not existing): {output_dir}\n")


def load_gwas(path, label):
    print(f"📥 Loading {label} GWAS data...")
    df = pd.read_csv(path)
    print(f"✅ {label} GWAS loaded successfully | Shape: {df.shape}\n")
    return df


def qc_report(df, label):
    print(f"🔎 QC Summary: {label} GWAS")
    print("-" * 40)
    print(df.info())
    print("\nMissing values per column:")
    print(df.isna().sum())
    print("\n")


def manhattan_plot(df, output_path, title):
    if not {'CHR', 'BP', 'PVALUE'}.issubset(df.columns):
        print(f"❌ ERROR: Missing required columns in GWAS data: {output_path}")
        sys.exit(1)

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

    plt.axhline(y=-np.log10(5e-8), color="black", linestyle="dashed",
                linewidth=1.5, label="Genome-wide significance (5e-8)")

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


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    exposure_path, outcome_path, output_dir = sys.argv[1], sys.argv[2], sys.argv[3]

    print_header()
    validate_inputs(exposure_path, outcome_path, output_dir)

    exposure = load_gwas(exposure_path, "Exposure")
    outcome  = load_gwas(outcome_path, "Outcome")

    qc_report(exposure, "Exposure")
    qc_report(outcome, "Outcome")

    print("📊 Generating Manhattan plots...\n")
    manhattan_plot(exposure, os.path.join(output_dir, "exposure_manhattan.png"), "Exposure GWAS Manhattan Plot")
    manhattan_plot(outcome, os.path.join(output_dir, "outcome_manhattan.png"), "Outcome GWAS Manhattan Plot")

    print("🎉 Exploratory analysis completed successfully!")
    print(f"All outputs saved in: {output_dir}\n")
    print("="*60)


if __name__ == "__main__":
    main()