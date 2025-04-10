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
  stop("Usage: Rscript 05_mr_analysis.R <input_ld_pruned_snps.csv>")
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


# ---- Manual Harmonisation ----
harmonised_data <- merge(exposure_dat, outcome_dat, by = "SNP")
harmonised_data$mr_keep <- TRUE  # Keep all SNPs (already harmonised upstream)

write.csv(harmonised_data, "harmonised_data.csv", row.names = FALSE)

cat("📊 SNPs after harmonisation:", nrow(harmonised_data), "\n\n")


# ---- Run MR Analyses ----
cat("⚙️ Running MR methods...\n")
res_ivw <- mr(harmonised_data, method = "mr_ivw")
res_egger <- mr(harmonised_data, method = "mr_egger_regression")
res_wmedian <- mr(harmonised_data, method = "mr_weighted_median")


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
results <- data.frame(
  N_SNPs = nrow(harmonised_data),

  IVW_OR = if (!is.null(res_ivw$b)) exp(res_ivw$b) else NA,
  IVW_CI_Lower = if (!is.null(res_ivw$b)) exp(res_ivw$b + qnorm(.025) * res_ivw$se) else NA,
  IVW_CI_Upper = if (!is.null(res_ivw$b)) exp(res_ivw$b + qnorm(.975) * res_ivw$se) else NA,
  IVW_Pval = res_ivw$pval,

  WM_OR = if (!is.null(res_wmedian$b)) exp(res_wmedian$b) else NA,
  WM_CI_Lower = if (!is.null(res_wmedian$b)) exp(res_wmedian$b + qnorm(.025) * res_wmedian$se) else NA,
  WM_CI_Upper = if (!is.null(res_wmedian$b)) exp(res_wmedian$b + qnorm(.975) * res_wmedian$se) else NA,
  WM_Pval = res_wmedian$pval,

  Egger_OR = if (!is.null(res_egger$b)) exp(res_egger$b) else NA,
  Egger_CI_Lower = if (!is.null(res_egger$b)) exp(res_egger$b + qnorm(.025) * res_egger$se) else NA,
  Egger_CI_Upper = if (!is.null(res_egger$b)) exp(res_egger$b + qnorm(.975) * res_egger$se) else NA,
  Egger_Pval = res_egger$pval,

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