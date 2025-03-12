# Import libraries
import os
import pandas as pd
import numpy as np

# First set seed for reproducibility
np.random.seed(42)

# Number of SNPs - 10K per GWAS
num_snps = 10000

# Generate random SNP ids
snp_ids = [f"rs{1000000 + i}" for i in range(num_snps)]

# Generate random GWAS summary statistics - LDL-c
ldl_gwas = pd.DataFrame({
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
ldl_gwas["PVALUE"] = p_values

# Generate random GWAS summary statistics - AD, ensuring P-values are never < 0.05
ad_gwas = pd.DataFrame({
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
ldl_gwas.to_csv("ldl_gwas.csv", index=False, sep=",")
ad_gwas.to_csv("ad_gwas.csv", index=False, sep=",")

print("Simulated GWAS statistics generated! At least 50% of SNPs in LDL GWAS have P-value < 0.05.")
print("At least 500 SNPs have P-value < 5e-08 in LDL GWAS.")