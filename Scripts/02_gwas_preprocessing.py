#!/usr/bin/env python3

"""
MR-CoPe | GWAS Harmonisation & Filtering Script
------------------------------------------------
Author: Guillermo Comesaña & Christian Pepler
Date: 2025

Description:
- Filters & harmonises exposure and outcome GWAS summary statistics.
- Keeps only SNPs present in both datasets.
- Drops missing values.
- Filters out INDELs (keeps only SNPs with alleles: A, T, C, G).
- Outputs a cleaned merged dataset ready for LD pruning and MR analysis.

Usage:
    python3 02_gwas_preprocessing.py <exposure_gwas.csv> <outcome_gwas.csv> <output_filtered.csv>
"""

import sys
import os
import pandas as pd

def print_header():
    print("\n" + "=" * 60)
    print("MR-CoPe | GWAS Harmonisation & F-Statistic Filtering")
    print("=" * 60 + "\n")


def validate_inputs(paths, labels):
    for path, label in zip(paths, labels):
        if not os.path.isfile(path):
            print(f"❌ ERROR: {label} file not found at: {path}")
            sys.exit(1)


def is_snp_allele(a):
    return a in {"A", "T", "C", "G"}


def load_gwas(path, label):
    print(f"📥 Loading {label} GWAS...")
    df = pd.read_csv(path)
    print(f"✅ {label} GWAS loaded | Shape: {df.shape}\n")
    return df


def calculate_f_statistics(df):
    df["F_stat"] = (df["BETA"] ** 2) / (df["SE"] ** 2)
    return df


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    exposure_path, outcome_path, output_path = sys.argv[1], sys.argv[2], sys.argv[3]

    print_header()
    validate_inputs([exposure_path, outcome_path], ["Exposure", "Outcome"])

    exposure = load_gwas(exposure_path, "Exposure")
    outcome  = load_gwas(outcome_path, "Outcome")

    # Identify common SNPs
    common_snps = set(exposure["SNP"]).intersection(set(outcome["SNP"]))
    print(f"🔗 SNPs shared between exposure & outcome: {len(common_snps)}\n")

    exposure = exposure[exposure["SNP"].isin(common_snps)].dropna()
    outcome = outcome[outcome["SNP"].isin(common_snps)].dropna()

    # Filter INDELs
    exposure = exposure[exposure["A1"].apply(is_snp_allele) & exposure["A2"].apply(is_snp_allele)]
    outcome = outcome[outcome["A1"].apply(is_snp_allele) & outcome["A2"].apply(is_snp_allele)]

    print(f"🧹 Exposure after INDEL removal: {exposure.shape}")
    print(f"🧹 Outcome after INDEL removal: {outcome.shape}\n")

    # Calculate F-statistic for Exposure SNPs
    print("🧮 Calculating F-statistics for exposure SNPs...")
    exposure = calculate_f_statistics(exposure)

    weak_snps = exposure[exposure["F_stat"] < 10]["SNP"]
    print(f"🚫 Removing weak SNPs with F < 10: {len(weak_snps)} SNPs\n")

    # Filter out weak SNPs from both exposure & outcome
    exposure = exposure[exposure["F_stat"] >= 10]
    outcome = outcome[outcome["SNP"].isin(exposure["SNP"])]

    print(f"📏 Exposure after F-stat filtering: {exposure.shape}")
    print(f"📏 Outcome after F-stat filtering: {outcome.shape}\n")

    # Merge exposure & outcome
    merged = pd.merge(exposure, outcome, on="SNP", suffixes=("_exp", "_out"))
    print(f"🔀 Final merged dataset shape: {merged.shape}\n")

    # Save output
    merged.to_csv(output_path, index=False)
    print(f"✅ Filtered & harmonised GWAS saved to: {output_path}")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()