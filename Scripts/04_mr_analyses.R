#!/usr/bin/env Rscript

# ================================================================
# MR-CoPe | Mendelian Randomisation Analysis using TwoSampleMR
# ================================================================
# Author: Guillermo Comesaña & Christian Pepler
#
# Usage:
#   Rscript 04_mr_analyses.R <input_ld_pruned_snps.csv>
#
# Description:
#   Performs MR analyses (IVW, Egger, Weighted Median) on LD-pruned,
#   harmonised GWAS summary statistics.
#
# Outputs:
# - exposure_dat.csv
# - outcome_dat.csv
# - harmonised_data.csv
# - MR_Formatted_Results.csv
# - MR_IVW_OR_Per_SNP.csv
# ================================================================

# ---- Required Libraries (with auto-install) ----
install_if_missing <- function(pkg, repo = "https://cloud.r-project.org") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = repo)
  }
}

install_if_missing("nloptr")
install_if_missing("lme4")
install_if_missing("meta")
install_if_missing("remotes")

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never", lib = .libPaths()[1])
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(optparse)
})

# ---- Handle Arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript 04_mr_analyses.R <input_ld_pruned_snps.csv>")
}

input_file <- args[1]

cat("\n==============================================================\n")
cat("MR-CoPe | Mendelian Randomisation Analysis\n")
cat("==============================================================\n\n")

if (!file.exists(input_file) || file.info(input_file)$size == 0) {
  cat("❌ ERROR: Input file is missing or empty — no SNPs passed LD pruning.\n")
  cat("🛑 Skipping MR analysis step.\n")

  empty_df <- data.frame()
  write.csv(empty_df, "MR_Formatted_Results.csv", row.names = FALSE)
  write.csv(empty_df, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)
  write.csv(empty_df, "exposure_dat.csv", row.names = FALSE)
  write.csv(empty_df, "outcome_dat.csv", row.names = FALSE)
  write.csv(empty_df, "harmonised_data.csv", row.names = FALSE)

  quit(status = 0)
}

cat("📥 Loading input file:", input_file, "\n")
filtered_snps <- read.csv(input_file)
cat("📊 SNPs in analysis:", nrow(filtered_snps), "\n\n")

colnames(filtered_snps) <- gsub("\\s+", "", colnames(filtered_snps))

cat("🧾 Columns available in merged dataset:\n")
print(colnames(filtered_snps))

# ---- Fallback Renaming for EAF ----
if (!"eaf" %in% names(filtered_snps)) {
  eaf_candidates <- c("EAF_exp", "effect_allele_freq", "eaf_exposure")
  existing <- intersect(eaf_candidates, names(filtered_snps))
  if (length(existing) > 0) {
    names(filtered_snps)[names(filtered_snps) == existing[1]] <- "eaf"
  } else {
    stop("❌ Could not find a column to map as 'eaf'.")
  }
}

if (!"eaf.outcome" %in% names(filtered_snps)) {
  eaf_out_candidates <- c("EAF_out", "effect_allele_freq.outcome", "eaf_outcome")
  existing_out <- intersect(eaf_out_candidates, names(filtered_snps))
  if (length(existing_out) > 0) {
    names(filtered_snps)[names(filtered_snps) == existing_out[1]] <- "eaf.outcome"
  } else {
    stop("❌ Could not find a column to map as 'eaf.outcome'.")
  }
}

# ---- Rename for harmonisation ----
rename_cols <- c(
  A1_exp = "effect_allele",
  A2_exp = "other_allele",
  A1_out = "effect_allele.outcome",
  A2_out = "other_allele.outcome"
)
filtered_snps <- filtered_snps %>% rename(any_of(rename_cols))

cat("🔄 Renamed: a1_exp → effect_allele \n")
cat("🔄 Renamed: a2_exp → other_allele \n")
cat("🔄 Renamed: a1_out → effect_allele.outcome \n")
cat("🔄 Renamed: a2_out → other_allele.outcome \n\n")

# ---- Exposure Dataset ----
exposure_dat <- filtered_snps %>%
  select(SNP, beta = BETA_exp, se = SE_exp, pval = PVALUE_exp,
         eaf, effect_allele, other_allele) %>%
  drop_na()

write.csv(exposure_dat, "exposure_dat.csv", row.names = FALSE)

# ---- Outcome Dataset ----
outcome_dat <- filtered_snps %>%
  select(SNP, beta = BETA_out, se = SE_out, pval = PVALUE_out,
         eaf = eaf.outcome, effect_allele = effect_allele.outcome,
         other_allele = other_allele.outcome) %>%
  drop_na()

write.csv(outcome_dat, "outcome_dat.csv", row.names = FALSE)

# ---- Format & Harmonise ----
exposure_dat <- format_data(exposure_dat, type = "exposure")
outcome_dat <- format_data(outcome_dat, type = "outcome")

harmonised_data <- harmonise_data(exposure_dat, outcome_dat, action = 2)
write.csv(harmonised_data, "harmonised_data.csv", row.names = FALSE)

cat("📊 SNPs after harmonisation:", nrow(harmonised_data), "\n\n")

if (nrow(harmonised_data) == 0) {
  cat("❌ No data available for MR analysis after filtering. Check SNP overlap and F-stat filtering.\n")
  empty_df <- data.frame()
  write.csv(empty_df, "MR_Formatted_Results.csv", row.names = FALSE)
  write.csv(empty_df, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)
  quit(status = 0)
}

# ---- MR Analysis ----
cat("⚙️ Running MR methods...\n")
methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
results_df <- mr(harmonised_data, method_list = methods)

# ---- Sensitivity ----
heterogeneity <- mr_heterogeneity(harmonised_data)
pleiotropy <- mr_pleiotropy_test(harmonised_data)

# ---- Per-SNP OR ----
IVW_SNP_results <- data.frame(
  SNP = harmonised_data$SNP,
  IVW_OR = exp(harmonised_data$beta.exposure / harmonised_data$se.exposure),
  IVW_Lower_95 = exp((harmonised_data$beta.exposure - 1.96 * harmonised_data$se.exposure) / harmonised_data$se.exposure),
  IVW_Upper_95 = exp((harmonised_data$beta.exposure + 1.96 * harmonised_data$se.exposure) / harmonised_data$se.exposure)
)
write.csv(IVW_SNP_results, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)

# ---- Summary ----
results <- tibble(
  N_SNPs = nrow(harmonised_data),

  IVW_OR = results_df %>% filter(method == "Inverse variance weighted") %>% pull(b) %>% exp(),
  IVW_CI_Lower = results_df %>% filter(method == "Inverse variance weighted") %>%
    mutate(lower = b + qnorm(.025) * se) %>% pull(lower) %>% exp(),
  IVW_CI_Upper = results_df %>% filter(method == "Inverse variance weighted") %>%
    mutate(upper = b + qnorm(.975) * se) %>% pull(upper) %>% exp(),
  IVW_Pval = results_df %>% filter(method == "Inverse variance weighted") %>% pull(pval),

  WM_OR = results_df %>% filter(method == "Weighted median") %>% pull(b) %>% exp(),
  WM_CI_Lower = results_df %>% filter(method == "Weighted median") %>%
    mutate(lower = b + qnorm(.025) * se) %>% pull(lower) %>% exp(),
  WM_CI_Upper = results_df %>% filter(method == "Weighted median") %>%
    mutate(upper = b + qnorm(.975) * se) %>% pull(upper) %>% exp(),
  WM_Pval = results_df %>% filter(method == "Weighted median") %>% pull(pval),

  Egger_OR = results_df %>% filter(method == "MR Egger") %>% pull(b) %>% exp(),
  Egger_CI_Lower = results_df %>% filter(method == "MR Egger") %>%
    mutate(lower = b + qnorm(.025) * se) %>% pull(lower) %>% exp(),
  Egger_CI_Upper = results_df %>% filter(method == "MR Egger") %>%
    mutate(upper = b + qnorm(.975) * se) %>% pull(upper) %>% exp(),
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
