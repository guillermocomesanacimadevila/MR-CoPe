#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import csv

# --- Handle arguments --- #
if len(sys.argv) != 3:
    print("Usage: python3 01_generate_simulated_data.py <exposure_output.csv> <outcome_output.csv>")
    sys.exit(1)

exposure_output_path = sys.argv[1]
outcome_output_path = sys.argv[2]

# --- Config --- #
np.random.seed(42)
num_snps = 10_000
num_strong = 8000       # p < 5e-8 in exposure
num_indep = 1200        # SNPs with independent genotypes
num_individuals = 500   # for LD simulation

# --- SNP metadata --- #
snp_ids = [f"rs{1_000_000 + i}" for i in range(num_snps)]
chr_values = np.random.randint(1, 23, size=num_snps)
bp_values = np.random.randint(1_000_000, 50_000_000, size=num_snps)
nucleotides = ["A", "T", "C", "G"]

# Ensure A1 ≠ A2 to prevent weird MR behavior
a1 = np.random.choice(nucleotides, size=num_snps)
a2 = np.array([
    np.random.choice([n for n in nucleotides if n != ref])
    for ref in a1
])

eafs = np.random.uniform(0.05, 0.5, size=num_snps)

# --- Simulate genotypes to enforce LD diversity (not saved) --- #
def simulate_genotypes(p, n, seed=None):
    if seed is not None:
        np.random.seed(seed)
    probs = [(1 - p)**2, 2 * p * (1 - p), p**2]
    return np.random.choice([0, 1, 2], size=n, replace=True, p=probs)

genotypes = []
for i in range(num_snps):
    seed = i + 1000 if i < num_indep else (i % num_indep) + 1000
    genotypes.append(simulate_genotypes(eafs[i], num_individuals, seed=seed))

# --- Simulate exposure GWAS stats --- #
beta_exp = np.zeros(num_snps)
beta_exp[:num_strong] = np.random.normal(0.08, 0.01, num_strong)
beta_exp[num_strong:] = np.random.normal(0.001, 0.01, num_snps - num_strong)
se_exp = np.random.uniform(0.01, 0.03, num_snps)
pval_exp = np.ones(num_snps)
pval_exp[:num_strong] = np.random.uniform(1e-12, 5e-8, num_strong)

# --- Simulate outcome GWAS stats --- #
beta_out = np.random.normal(0.001, 0.01, num_snps)
se_out = np.random.uniform(0.01, 0.03, num_snps)
pval_out = np.ones(num_snps)
pval_out[:num_strong] = np.random.uniform(0.1, 1.0, num_strong)

# --- Build exposure DataFrame --- #
exposure_df = pd.DataFrame({
    "SNP": snp_ids,
    "CHR": chr_values,
    "BP": bp_values,
    "A1": a1,
    "A2": a2,
    "EAF": eafs,
    "BETA": beta_exp,
    "SE": se_exp,
    "PVALUE": pval_exp
})

# --- Build outcome DataFrame --- #
outcome_df = pd.DataFrame({
    "SNP": snp_ids,
    "CHR": chr_values,
    "BP": bp_values,
    "A1": a1,
    "A2": a2,
    "EAF": eafs,
    "BETA": beta_out,
    "SE": se_out,
    "PVALUE": pval_out
})

# --- Enforce allele column types (extra safety) --- #
for df in [exposure_df, outcome_df]:
    df["A1"] = df["A1"].astype(str)
    df["A2"] = df["A2"].astype(str)

# --- Save to CSV --- #
exposure_df.to_csv(exposure_output_path, index=False, quoting=csv.QUOTE_NONNUMERIC)
outcome_df.to_csv(outcome_output_path, index=False, quoting=csv.QUOTE_NONNUMERIC)

# --- Logs --- #
print(f"[✅] Simulated {num_snps:,} SNPs")
print(f"[🔗] {num_indep:,} SNPs are independent (genotype-seeded)")
print(f"[📉] {np.sum(pval_exp < 5e-8):,} SNPs with p < 5e-8 in exposure")
print(f"[📈] {np.sum(pval_out > 5e-8):,} SNPs with p > 5e-8 in outcome")
print(f"[💾] Exposure saved to: {exposure_output_path}")
print(f"[💾] Outcome saved to:  {outcome_output_path}")