#!/usr/bin/env python3

# Import libraries
import sys
import pandas as pd
import numpy as np

# --- Handle arguments --- #
if len(sys.argv) != 3:
    print("Usage: python3 01_generate_simulated_data.py <exposure_output.csv> <outcome_output.csv>")
    sys.exit(1)

exposure_output_path = sys.argv[1]
outcome_output_path = sys.argv[2]

# --- Set seed for reproducibility --- #
np.random.seed(42)

# Number of SNPs - 10K per GWAS
num_snps = 10000

# Generate random SNP ids
snp_ids = [f"rs{1000000 + i}" for i in range(num_snps)]

# Generate random GWAS summary statistics - LDL-c
exposure_df = pd.DataFrame({
    "SNP": snp_ids,
    "CHR": np.random.randint(1, 23, size=num_snps),
    "BP": np.random.randint(1e6, 5e7, size=num_snps),
    "A1": np.random.choice(["A", "T", "C", "G"], size=num_snps),
    "A2": np.random.choice(["A", "T", "C", "G"], size=num_snps),
    "EAF": np.random.uniform(0.01, 0.99, size=num_snps),
    "BETA": np.random.normal(0, 0.1, size=num_snps),
    "SE": np.random.uniform(0.02, 0.05, size=num_snps),  # Adjusted SE values to avoid extreme F-stat
})

# Ensure that at least 50% of SNPs have p-value < 0.05 in LDL GWAS
p_values = np.random.uniform(0, 1, size=num_snps)
p_values[: num_snps // 2] = np.random.uniform(0, 0.05, size=num_snps // 2)  # Ensure 50% are < 0.05
p_values[:500] = np.random.uniform(1e-10, 5e-08, 500)  # Ensure at least 500 SNPs have p < 5e-08
exposure_df["PVALUE"] = p_values

# Generate random GWAS summary statistics - AD, ensuring P-values are never < 0.05
outcome_df = pd.DataFrame({
    "SNP": snp_ids,
    "CHR": np.random.randint(1, 23, size=num_snps),
    "BP": np.random.randint(1e6, 5e7, size=num_snps),
    "A1": np.random.choice(["A", "T", "C", "G"], size=num_snps),
    "A2": np.random.choice(["A", "T", "C", "G"], size=num_snps),
    "EAF": np.random.uniform(0.01, 0.99, size=num_snps),
    "BETA": np.random.normal(0, 0.1, size=num_snps),
    "SE": np.random.uniform(0.01, 0.05, size=num_snps),
    "PVALUE": np.random.uniform(0.05, 1, size=num_snps)  # Ensuring P-values ≥ 0.05
})

# Save to CSV
exposure_df.to_csv(exposure_output_path, index=False)
outcome_df.to_csv(outcome_output_path, index=False)

# --- Logging --- #
print(f"[INFO] Simulated exposure GWAS saved to: {exposure_output_path}")
print(f"[INFO] Simulated outcome GWAS saved to:  {outcome_output_path}")
print("[INFO] Exposure GWAS: ≥50% SNPs with p < 0.05, ≥500 SNPs with p < 5e-08.")
print("[INFO] Outcome GWAS:  No SNPs with p < 0.05.")
