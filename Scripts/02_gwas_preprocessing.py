#!/usr/bin/env python3

"""
MR-CoPe | GWAS Harmonisation & Filtering Script
------------------------------------------------
Author: Guillermo Comesaña & Christian Pepler
Date: 2025

Description:
- Automatically detects GWAS column names
- Filters & harmonises exposure and outcome GWAS summary statistics
- Filters for SNPs in both datasets
- Removes INDELs (keeps only SNPs with alleles A/T/C/G)
- Filters weak SNPs (F-stat < 10)
- Outputs a merged dataset ready for MR analysis

Usage:
    python3 02_gwas_preprocessing.py <exposure_gwas.csv> <outcome_gwas.csv> <output_filtered.csv>
"""

#!/usr/bin/env python3

"""
MR-CoPe | GWAS Harmonisation & Filtering Script
------------------------------------------------
Author: Guillermo Comesaña & Christian Pepler
Date: 2025

Description:
- Automatically detects GWAS column names
- Filters & harmonises exposure and outcome GWAS summary statistics
- Filters for SNPs in both datasets
- Removes INDELs (keeps only SNPs with alleles A/T/C/G)
- Filters weak SNPs (F-stat < 10)
- Outputs a merged dataset ready for MR analysis

Usage:
    python3 02_gwas_preprocessing.py <exposure_gwas.csv> <outcome_gwas.csv> <output_filtered.csv>
"""

import sys
import os
import pandas as pd
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

def print_header():
    print("\n" + "=" * 60)
    print("MR-CoPe | GWAS Harmonisation & F-Statistic Filtering")
    print("=" * 60 + "\n")

def validate_inputs(paths, labels):
    for path, label in zip(paths, labels):
        if not os.path.isfile(path):
            print(f"❌ ERROR: {label} file not found at: {path}")
            sys.exit(1)

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

def is_snp_allele(a):
    return a in {"A", "T", "C", "G"}

def infer_se_if_missing(df):
    if "SE" not in df.columns and "BETA" in df.columns and "PVALUE" in df.columns:
        print("🧠 Inferring SE from BETA and PVALUE...")
        df["SE"] = abs(df["BETA"] / norm.ppf(df["PVALUE"] / 2))
        df["SE"] = pd.to_numeric(df["SE"], errors="coerce")
    return df

def load_gwas(path, label):
    print(f"📥 Loading {label} GWAS: {path}")
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
        for std_col, original_col in mapping.items():
            if std_col != original_col and original_col in df.columns:
                df.rename(columns={original_col: std_col}, inplace=True)

    df = infer_se_if_missing(df)
    return df

def calculate_f_statistics(df):
    if "SE" in df.columns:
        df["F_stat"] = (df["BETA"] ** 2) / (df["SE"] ** 2)
    else:
        print("⚠️ WARNING: SE column not found — Skipping F-statistic filtering")
        df["F_stat"] = 9999
    return df

def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    exposure_path, outcome_path, output_path = sys.argv[1], sys.argv[2], sys.argv[3]

    print_header()
    validate_inputs([exposure_path, outcome_path], ["Exposure", "Outcome"])

    exposure = load_gwas(exposure_path, "Exposure")
    outcome = load_gwas(outcome_path, "Outcome")

    common_snps = set(exposure["SNP"]).intersection(set(outcome["SNP"]))
    print(f"🔗 Shared SNPs between exposure & outcome: {len(common_snps)}\n")

    exposure = exposure[exposure["SNP"].isin(common_snps)].dropna()
    outcome = outcome[outcome["SNP"].isin(common_snps)].dropna()

    if "A1" in exposure.columns and "A2" in exposure.columns:
        exposure = exposure[exposure["A1"].apply(is_snp_allele) & exposure["A2"].apply(is_snp_allele)]
        outcome = outcome[outcome["A1"].apply(is_snp_allele) & outcome["A2"].apply(is_snp_allele)]
        print(f"🧹 Exposure after INDEL removal: {exposure.shape}")
        print(f"🧹 Outcome after INDEL removal: {outcome.shape}\n")
    else:
        print("⚠️ WARNING: No A1/A2 columns found — Skipping INDEL filtering\n")

    print("🧮 Calculating F-statistics for exposure SNPs...")
    exposure = calculate_f_statistics(exposure)

    weak_snps = exposure[exposure["F_stat"] < 10]["SNP"]
    print(f"🚫 Removing weak SNPs with F < 10: {len(weak_snps)} SNPs\n")

    exposure = exposure[exposure["F_stat"] >= 10]
    outcome = outcome[outcome["SNP"].isin(exposure["SNP"])]

    print(f"📏 Exposure after F-stat filtering: {exposure.shape}")
    print(f"📏 Outcome after F-stat filtering: {outcome.shape}\n")

    merged = pd.merge(exposure, outcome, on="SNP", suffixes=("_exp", "_out"))

    if "riskFrequency_exp" in merged.columns:
        merged.rename(columns={"riskFrequency_exp": "EAF_exp"}, inplace=True)
    if "riskFrequency_out" in merged.columns:
        merged.rename(columns={"riskFrequency_out": "EAF_out"}, inplace=True)

    print(f"🔀 Final merged dataset shape: {merged.shape}\n")
    merged.to_csv(output_path, index=False)
    print(f"✅ Filtered & harmonised GWAS saved to: {output_path}")
    print("=" * 60 + "\n")

if __name__ == "__main__":
    main()
