import os
import pandas as pd

# Load GWAS data
def parse_csv(location):
    with open(os.path.expanduser(location), "r") as file:
        df = pd.read_csv(file, sep=",")
        df.columns = df.columns.str.strip()  # Ensure column names are clean
        return df

# Remove INDELs (Ensure A1 and A2 are single bases)
def clean_gwas(df):
    return df[(df["A1"].str.len() == 1) & (df["A2"].str.len() == 1)]

# Load GWAS Data
ldl_c = parse_csv("/home/jovyan/MR_simulation/Scripts/ldl_gwas.csv")
ad = parse_csv("/home/jovyan/MR_simulation/Scripts/ad_gwas.csv")

# Print column names before processing
print("LDL-C GWAS Columns:", ldl_c.columns)
print("AD GWAS Columns:", ad.columns)

# Clean Data
ldl_c = clean_gwas(ldl_c)
ad = clean_gwas(ad)

# Ensure SNP column is explicitly retained
if "SNP" not in ldl_c.columns or "SNP" not in ad.columns:
    raise KeyError("SNP column is missing from the dataset!")

# Check shape
print(f"LDL GWAS shape: {ldl_c.shape}")
print(f"AD GWAS shape: {ad.shape}")

# Match SNPs from both GWAS
matched_snps = pd.merge(ldl_c, ad, on="SNP", suffixes=("_exp", "_out"))

# Ensure SNP column is still present
if "SNP" not in matched_snps.columns:
    raise KeyError("SNP column was lost after merging!")

# Check matched SNPs
print(f"Common SNPs: {matched_snps.shape[0]}")

# Compute F-statistic
matched_snps["F_stat"] = (matched_snps["BETA_exp"] ** 2) / (matched_snps["SE_exp"] ** 2)

# **Debugging: Count SNPs failing each filter**
print(f"SNPs failing F-stat > 10: {(matched_snps['F_stat'] <= 10).sum()}")
print(f"SNPs failing EAF >= 0.01: {(matched_snps['EAF_exp'] < 0.01).sum()}")
print(f"SNPs failing P-value < 5e-08: {(matched_snps['PVALUE_exp'] >= 5e-08).sum()}")
print(f"SNPs failing Outcome P-value > 5e-08: {(matched_snps['PVALUE_out'] < 5e-08).sum()}")

# Apply filtering
filtered_snps = matched_snps[
    (matched_snps["F_stat"] > 10) &  # Instrument strength check
    (matched_snps["EAF_exp"] >= 0.01) &  # Ensure MAF > 0.01
    (matched_snps["PVALUE_exp"] < 5e-08) &  # Only significant exposure SNPs
    (matched_snps["PVALUE_out"] >= 5e-08)  # Avoid outcome-associated SNPs
]

# Ensure SNP column is still present
if "SNP" not in filtered_snps.columns:
    raise KeyError("SNP column was lost after filtering!")

# Check filtered SNPs
print(f"Filtered SNPs: {filtered_snps.shape}")

# Save the Filtered SNPs for MR
filtered_snps.to_csv("filtered_SNPs_for_MR.csv", index=False)
print("Filtered SNPs saved successfully.")
