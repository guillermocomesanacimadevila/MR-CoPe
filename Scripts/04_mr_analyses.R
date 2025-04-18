#!/usr/bin/env Rscript

# ================================================================
# MR-CoPe | Mendelian Randomisation Analysis using TwoSampleMR
# ================================================================
# Author: Guillermo Comesaña & Christian Pepler
# Date: 2025
#
# Usage:
# Rscript 04_mr_analysis.R <input_ld_pruned_snps.csv>
#
# Description:
# Performs MR analyses (IVW, Egger, Weighted Median) on LD-pruned,
# harmonised GWAS summary statistics.
#
# Outputs:
# - exposure_dat.csv
# - outcome_dat.csv
# - harmonised_data.csv
# - MR_Formatted_Results.csv
# - MR_IVW_OR_Per_SNP.csv
# ================================================================

# ---- Required Libraries ----
required_packages <- c("TwoSampleMR", "optparse", "tidyverse")

missing <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing) > 0) {
  cat("📦 Installing missing R packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org/")
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
})


# ---- Handle Arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript 04_mr_analysis.R <input_ld_pruned_snps.csv>")
}

input_file <- args[1]


cat("\n==============================================================\n")
cat("MR-CoPe | Mendelian Randomisation Analysis\n")
cat("==============================================================\n\n")

cat("📥 Loading input file:", input_file, "\n")
filtered_snps <- read.csv(input_file)
cat("📊 SNPs in analysis:", nrow(filtered_snps), "\n\n")


# ---- Check & Clean Column Names ----
colnames(filtered_snps) <- gsub("\\s+", "", colnames(filtered_snps))


# ---- Rename Columns for TwoSampleMR ----
rename_cols <- c(
  A1_exp = "effect_allele",
  A2_exp = "other_allele",
  EAF_exp = "eaf",
  A1_out = "effect_allele.outcome",
  A2_out = "other_allele.outcome",
  EAF_out = "eaf.outcome"
)

filtered_snps <- filtered_snps %>% rename(any_of(rename_cols))


# ---- Exposure Dataset ----
exposure_dat <- filtered_snps %>%
  select(SNP, beta = BETA_exp, se = SE_exp, pval = PVALUE_exp,
         eaf, effect_allele, other_allele) %>%
  drop_na()

write.csv(exposure_dat, "exposure_dat.csv", row.names = FALSE)


# ---- Outcome Dataset ----
outcome_dat <- filtered_snps %>%
  select(SNP, beta = BETA_out, se = SE_out, pval = PVALUE_out,
         eaf = eaf.outcome, effect_allele = effect_allele.outcome, other_allele = other_allele.outcome) %>%
  drop_na()

write.csv(outcome_dat, "outcome_dat.csv", row.names = FALSE)


# ---- Format Data for MR ----
exposure_dat <- format_data(exposure_dat, type = "exposure")
outcome_dat <- format_data(outcome_dat, type = "outcome")


# ---- Harmonise Datasets ----
harmonised_data <- harmonise_data(exposure_dat, outcome_dat, action = 2)

write.csv(harmonised_data, "harmonised_data.csv", row.names = FALSE)

cat("📊 SNPs after harmonisation:", nrow(harmonised_data), "\n\n")


# ---- Run MR Analyses ----
cat("⚙️ Running MR methods...\n")
methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
results_df <- mr(harmonised_data, method_list = methods)


# ---- Sensitivity Analyses ----
heterogeneity <- mr_heterogeneity(harmonised_data)
pleiotropy <- mr_pleiotropy_test(harmonised_data)


# ---- IVW OR per SNP ----
IVW_SNP_results <- data.frame(
  SNP = harmonised_data$SNP,
  IVW_OR = exp(harmonised_data$beta.exposure / harmonised_data$se.exposure),
  IVW_Lower_95 = exp((harmonised_data$beta.exposure - 1.96 * harmonised_data$se.exposure) / harmonised_data$se.exposure),
  IVW_Upper_95 = exp((harmonised_data$beta.exposure + 1.96 * harmonised_data$se.exposure) / harmonised_data$se.exposure)
)

write.csv(IVW_SNP_results, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)


# ---- MR Summary Table ----
results <- tibble(
  N_SNPs = nrow(harmonised_data),

  IVW_OR = results_df %>% filter(method == "Inverse variance weighted") %>% pull(b) %>% exp(),
  IVW_CI_Lower = results_df %>% filter(method == "Inverse variance weighted") %>% mutate(lower = b + qnorm(.025) * se) %>% pull(lower) %>% exp(),
  IVW_CI_Upper = results_df %>% filter(method == "Inverse variance weighted") %>% mutate(upper = b + qnorm(.975) * se) %>% pull(upper) %>% exp(),
  IVW_Pval = results_df %>% filter(method == "Inverse variance weighted") %>% pull(pval),

  WM_OR = results_df %>% filter(method == "Weighted median") %>% pull(b) %>% exp(),
  WM_CI_Lower = results_df %>% filter(method == "Weighted median") %>% mutate(lower = b + qnorm(.025) * se) %>% pull(lower) %>% exp(),
  WM_CI_Upper = results_df %>% filter(method == "Weighted median") %>% mutate(upper = b + qnorm(.975) * se) %>% pull(upper) %>% exp(),
  WM_Pval = results_df %>% filter(method == "Weighted median") %>% pull(pval),

  Egger_OR = results_df %>% filter(method == "MR Egger") %>% pull(b) %>% exp(),
  Egger_CI_Lower = results_df %>% filter(method == "MR Egger") %>% mutate(lower = b + qnorm(.025) * se) %>% pull(lower) %>% exp(),
  Egger_CI_Upper = results_df %>% filter(method == "MR Egger") %>% mutate(upper = b + qnorm(.975) * se) %>% pull(upper) %>% exp(),
  Egger_Pval = results_df %>% filter(method == "MR Egger") %>% pull(pval),

  Egger_Intercept_Pval = pleiotropy$pval,
  IVW_Het_Q_Pval = heterogeneity$Q_pval[heterogeneity$method == "Inverse variance weighted"],
  I2_Statistic = heterogeneity$I2[heterogeneity$method == "Inverse variance weighted"]
)

write.csv(results, "MR_Formatted_Results.csv", row.names = FALSE)


cat("✅ MR analysis completed successfully!\n")
cat("📝 Outputs generated:\n")
cat("- exposure_dat.csv\n")
cat("- outcome_dat.csv\n")
cat("- harmonised_data.csv\n")
cat("- MR_Formatted_Results.csv\n")
cat("- MR_IVW_OR_Per_SNP.csv\n")
cat("==============================================================\n\n")
