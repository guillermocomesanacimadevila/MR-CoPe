import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load MR results
file_path = "/Users/guillermocomesanacimadevila/Desktop/CRyPTIC_cleaning/MR R/MR_Formatted_Results.csv"
mr_results = pd.read_csv(file_path)

# Convert relevant columns to numeric
numeric_cols = ["IVW_OR", "IVW_Lower_95%", "IVW_Upper_95%", "IVW_Pval",
                "WM_OR", "WM_Lower_95%", "WM_Upper_95%", "WM_Pval",
                "Egger_OR", "Egger_Lower_95%", "Egger_Upper_95%", "Egger_Pval",
                "I2_Statistic", "Mean_F_Statistic"]

mr_results[numeric_cols] = mr_results[numeric_cols].apply(pd.to_numeric, errors="coerce")

# Function to add confidence intervals to plots
def plot_forest(method, or_col, lower_col, upper_col, pval_col):
    plt.figure(figsize=(8, 5))
    sns.pointplot(x=mr_results[or_col], y=[method], join=False, color="blue")
    plt.errorbar(mr_results[or_col], [method], xerr=[
        mr_results[or_col] - mr_results[lower_col],
        mr_results[upper_col] - mr_results[or_col]
    ], fmt='o', color='blue', label=f"{method} OR (95% CI)")

    plt.axvline(x=1, linestyle="--", color="red", label="Null Effect (OR = 1)")
    plt.xlabel("Odds Ratio (OR)")
    plt.ylabel("Method")
    plt.title(f"Forest Plot: {method}")
    plt.legend()
    plt.show()

# Forest plots for each method
plot_forest("IVW", "IVW_OR", "IVW_Lower_95%", "IVW_Upper_95%", "IVW_Pval")
plot_forest("Weighted Median", "WM_OR", "WM_Lower_95%", "WM_Upper_95%", "WM_Pval")
plot_forest("Egger", "Egger_OR", "Egger_Lower_95%", "Egger_Upper_95%", "Egger_Pval")

# Funnel Plot
plt.figure(figsize=(6, 6))
sns.scatterplot(x=mr_results["IVW_OR"], y=-np.log10(mr_results["IVW_Pval"]), color="blue", label="IVW")
sns.scatterplot(x=mr_results["WM_OR"], y=-np.log10(mr_results["WM_Pval"]), color="green", label="WM")
sns.scatterplot(x=mr_results["Egger_OR"], y=-np.log10(mr_results["Egger_Pval"]), color="red", label="Egger")
plt.axvline(x=1, linestyle="--", color="black")
plt.xlabel("Odds Ratio")
plt.ylabel("-log10(P-value)")
plt.title("Funnel Plot for MR Methods")
plt.legend()
plt.show()

# Scatter Plot for IVW & Egger Regression
plt.figure(figsize=(7, 5))
sns.regplot(x=mr_results["IVW_OR"], y=mr_results["Egger_OR"], scatter_kws={"s": 40}, line_kws={"color": "red"})
plt.xlabel("IVW OR")
plt.ylabel("Egger OR")
plt.title("Scatter Plot: IVW vs Egger Regression")
plt.axhline(y=1, linestyle="--", color="gray")
plt.axvline(x=1, linestyle="--", color="gray")
plt.show()

# Leave-One-Out Analysis (if multiple SNPs)
if "N_SNPs" in mr_results.columns and mr_results["N_SNPs"].values[0] > 1:
    sns.boxplot(y=mr_results["IVW_OR"])
    plt.title("Leave-One-Out Sensitivity Analysis")
    plt.show()

print("Visualization completed!")