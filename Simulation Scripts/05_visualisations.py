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

# Set professional styling
sns.set_style("whitegrid")
plt.figure(figsize=(7, 10), dpi=300)

# Plot boxplot
sns.boxplot(data=top_snp_results, y="SNP", x="IVW_OR", color="skyblue", orient="h", width=0.6)

# Overlay individual points (strip plot)
sns.stripplot(data=top_snp_results, y="SNP", x="IVW_OR", color="black", alpha=0.7, jitter=True, orient="h", size=5)

# Add confidence intervals as error bars
for i, row in enumerate(top_snp_results.itertuples()):
    plt.plot([row.IVW_Lower_95, row.IVW_Upper_95], [i, i], color="black", linestyle="-", linewidth=1.5)

# Add reference line at OR = 1
plt.axvline(x=1, linestyle="--", color="red", linewidth=1.2, label="Null Effect (OR = 1)")

# Formatting
plt.xlabel("IVW Odds Ratio (95% CI)", fontsize=14, fontweight="bold")
plt.ylabel("SNP", fontsize=14, fontweight="bold")
plt.xticks(fontsize=12)
plt.yticks(fontsize=10)
plt.title("Top 20 SNPs by IVW OR Deviation", fontsize=16, fontweight="bold")

# Remove top and right spines for a cleaner look
sns.despine()

# Add legend
plt.legend(fontsize=12)

# Adjust layout for better spacing
plt.tight_layout()

# Save as high-quality image
plt.savefig("Top20_SNPs_IVW_OR_Boxplot.png", dpi=300, bbox_inches="tight")

# Show plot
plt.show()

print("Visualisations completed!")
