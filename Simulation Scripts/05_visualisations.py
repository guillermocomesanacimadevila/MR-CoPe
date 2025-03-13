import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load MR summary results
summary_file = "/Users/guillermocomesanacimadevila/Desktop/CRyPTIC_cleaning/MR R/MR_Formatted_Results.csv"
mr_results = pd.read_csv(summary_file)

# Load per-SNP IVW OR results
snp_file = "/Users/guillermocomesanacimadevila/Desktop/CRyPTIC_cleaning/MR R/MR_IVW_OR_Per_SNP.csv"
snp_results = pd.read_csv(snp_file)

# Ensure numeric conversion
snp_results[["IVW_OR", "IVW_Lower_95", "IVW_Upper_95"]] = snp_results[["IVW_OR", "IVW_Lower_95", "IVW_Upper_95"]].apply(pd.to_numeric, errors="coerce")

# Sort by deviation from OR = 1 (most extreme effects first)
snp_results["Deviation"] = abs(snp_results["IVW_OR"] - 1)
snp_results = snp_results.sort_values(by="Deviation", ascending=False)

# Select top 20 SNPs
top_snp_results = snp_results.head(20)

# Create vertical box plot
plt.figure(figsize=(6, 10))
sns.boxplot(data=top_snp_results, y="SNP", x="IVW_OR", color="lightblue", orient="h")

# Add individual points (scatter) for better visualization
sns.stripplot(data=top_snp_results, y="SNP", x="IVW_OR", color="black", alpha=0.7, jitter=True, orient="h")

# Add confidence intervals as error bars
for i, row in enumerate(top_snp_results.itertuples()):
    plt.plot([row.IVW_Lower_95, row.IVW_Upper_95], [i, i], color="black", linestyle="-", linewidth=1)

# Add reference line at OR = 1
plt.axvline(x=1, linestyle="--", color="red", label="Null Effect (OR = 1)")

# Formatting
plt.xlabel("IVW OR (95% CI)")
plt.ylabel("SNP")
plt.title("Top 20 SNPs by IVW OR Deviation from 1")
plt.legend()
plt.tight_layout()
plt.show()

print("Visualisation completed!")