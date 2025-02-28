# https://www.ebi.ac.uk/gwas/search?query=cholesterol

import os
import pandas as pd

def read_tsv(location):
    with open(os.path.expanduser(location), "r") as file:
        return pd.read_csv(location, sep="\t")

def parse_csv(location):
    with open(os.path.expanduser(location), "r") as file:
        return pd.read_csv(location, sep=",")

def common_snps(exp, out):
    common_snps = set(exp["riskAllele"]).intersection(set(out["riskAllele"]))
    return list(common_snps)

# Import exposure and outcome GWAS as TSV

exposure = read_tsv("/Users/guillermocomesanacimadevila/Desktop/T_Chol_GWAS.tsv")
outcome = read_tsv("/Users/guillermocomesanacimadevila/Desktop/AD_GWAS.tsv")

# import the dataframes as a CSV

exposure.to_csv("exposure.csv", header=True, sep=",")
outcome.to_csv("outcome.csv", header=True, sep=",")

# Now test the common SNPs

df_exposure = parse_csv("/Users/guillermocomesanacimadevila/Desktop/exposure.csv")
df_outcome = parse_csv("/Users/guillermocomesanacimadevila/Desktop/outcome.csv")

print(df_exposure.head(n=5))
print(common_snps(exposure, df_exposure))
print(f"Number of common SNPs: {len((common_snps(exposure, df_exposure)))}") # 5780
print(f"Total Number of SNPs Exposure: {df_exposure.shape[0]}, Total Number of SNPs Outcome: {df_outcome.shape[0]}") # 9640 vs 3169
