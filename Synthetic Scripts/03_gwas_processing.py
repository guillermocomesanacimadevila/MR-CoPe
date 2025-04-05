#!/usr/bin/env python3
# GWAS Harmonisation and Filtering Script
# Originally done in Jupyter Notebook
# Now pipeline-ready!

import sys
import os
import pandas as pd

# --- Handle command-line arguments --- #
if len(sys.argv) != 4:
    print("Usage: python3 03_gwas_processing.py <exposure_gwas.csv> <outcome_gwas.csv> <output_filtered.csv>")
    sys.exit(1)

exposure_path = sys.argv[1]
outcome_path = sys.argv[2]
output_path = sys.argv[3]

# --- Load exposure and outcome datasets --- #
exposure = pd.read_csv(os.path.expanduser(exposure_path))
outcome = pd.read_csv(os.path.expanduser(outcome_path))

print(f"[INFO] Loaded exposure: {exposure.shape}, outcome: {outcome.shape}")
print(f"[INFO] Common SNPs before merge: {len(set(exposure['SNP']).intersection(set(outcome['SNP'])))}")

# --- Keep only SNPs present in both datasets --- #
common_snps = set(exposure["SNP"]).intersection(set(outcome["SNP"]))
exposure = exposure[exposure["SNP"].isin(common_snps)]
outcome = outcome[outcome["SNP"].isin(common_snps)]

print(f"[INFO] SNPs in both datasets: {len(common_snps)}")

# --- Remove empty entries --- #
exposure = exposure.dropna()
outcome = outcome.dropna()

print(f"[INFO] Exposure after dropna: {exposure.shape}")
print(f"[INFO] Outcome after dropna: {outcome.shape}")

# --- Remove possible INDELs --- #
def is_valid_allele(a):
    return a in ["A", "T", "C", "G"]

exposure = exposure[
    exposure["A1"].apply(is_valid_allele) &
    exposure["A2"].apply(is_valid_allele)
]

outcome = outcome[
    outcome["A1"].apply(is_valid_allele) &
    outcome["A2"].apply(is_valid_allele)
]

print(f"[INFO] Exposure after indel filtering: {exposure.shape}")
print(f"[INFO] Outcome after indel filtering: {outcome.shape}")

# --- Merge exposure and outcome on SNP --- #
merged = pd.merge(
    exposure, outcome,
    on="SNP",
    suffixes=("_exp", "_out")
)

print(f"[INFO] Merged dataset shape: {merged.shape}")

# --- Rename relevant columns for downstream MR analysis --- #
merged = merged.rename(columns={
    "A1_exp": "A1_exp",
    "A2_exp": "A2_exp",
    "EAF_exp": "EAF_exp",
    "BETA_exp": "BETA_exp",
    "SE_exp": "SE_exp",
    "PVALUE_exp": "PVALUE_exp",
    "A1_out": "A1_out",
    "A2_out": "A2_out",
    "EAF_out": "EAF_out",
    "BETA_out": "BETA_out",
    "SE_out": "SE_out",
    "PVALUE_out": "PVALUE_out"
})

# --- Save merged and cleaned dataset --- #
merged.to_csv(output_path, index=False)
print(f"[INFO] Filtered SNP dataset saved to: {output_path}")