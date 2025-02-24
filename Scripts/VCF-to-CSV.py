import os
import csv
import pandas as pd
import numpy as np

# ==== Function accumulator ==== #

def read_file(location):
    with open(os.path.expanduser(location), "r") as file:
        return pd.read_csv(file, sep=",")

# Function to Convert VCF to CSV
def vcf_to_csv(vcf_file_path, csv_file_path):
    data = []  # Initialize data to avoid unbound local variable error

    with open(vcf_file_path, "r") as vcf_file:
        with open(csv_file_path, "w", newline="") as csv_file:
            csv_writer = csv.writer(csv_file)
            for line in vcf_file:
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    header = line.strip("#").strip().split("\t")
                    csv_writer.writerow(header)
                else:
                    data = line.strip().split("\t")
                    csv_writer.writerow(data)

    return data if data else None  # Ensure it returns something, even if empty

# Correct file paths
exposure_vcf_path = "~/cpep_MR/Data/exposure.vcf"
exposure_csv_path = "~/cpep_MR/Data/exposure.csv"

outcome_vcf_path = "~/cpep_MR/Data/outcome.vcf"
outcome_csv_path = "~/cpep_MR/Data/outcome.csv"

# Expand user paths and run conversion
vcf_to_csv(os.path.expanduser(exposure_vcf_path), os.path.expanduser(exposure_csv_path))
vcf_to_csv(os.path.expanduser(outcome_vcf_path), os.path.expanduser(outcome_csv_path))

# load up cleaned VCFs as a pd df
exposure = read_file("~/cpep_MR/Cleaned Data/cleaned_exposure.csv")
outcome = read_file("~/cpep_MR/Cleaned Data/cleaned_outcome.csv")

# Separate final columns to get (Beta, SE, etc...)
column_name = "ieu-b-5067"

# ==== Exposure Formatting ==== #
df1 = exposure[column_name].str.split(":", expand=True)
df1.columns = ["Beta", "SE", "LP", "AF", "ID"]

# Convert to numeric
df1["Beta"] = pd.to_numeric(df1["Beta"], errors="coerce")
df1["SE"] = pd.to_numeric(df1["SE"], errors="coerce")
df1["LP"] = pd.to_numeric(df1["LP"], errors="coerce")

# Compute p-value
df1["P_VALUE"] = 10 ** (-df1["LP"])

# Add SNP id and Position
df1["SNP"] = exposure["ID"]
df1["POS"] = exposure["POS"]  

print(df1.head(n=5), f"Shape: {df1.shape}")

# ==== Outcome Formatting ==== #
df2 = outcome[column_name].str.split(":", expand=True)
df2.columns = ["Beta", "SE", "LP", "AF", "ID"]

# Convert to numeric
df2["Beta"] = pd.to_numeric(df2["Beta"], errors="coerce")
df2["SE"] = pd.to_numeric(df2["SE"], errors="coerce")
df2["LP"] = pd.to_numeric(df2["LP"], errors="coerce")

# Compute p-value
df2["P_VALUE"] = 10 ** (-df2["LP"])

# Add SNP id and Position
df2["SNP"] = outcome["ID"]
df2["POS"] = outcome["POS"]

print(df2.head(n=5), f"Shape: {df2.shape}")

# Save updated DataFrames
df1.to_csv("~/cpep_MR/Cleaned Data/exp_stats.csv", index=False)
df2.to_csv("~/cpep_MR/Cleaned Data/outcome_stats.csv", index=False)

# Check column names of original VCF
print(exposure.columns)
