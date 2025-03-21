#!/usr/bin/env Rscript
# MR Analysis using TwoSampleMR
# This script will be used in the main pipeline.

# --- Load Required Package ---
suppressPackageStartupMessages(library(TwoSampleMR))

# --- Handle Command-Line Arguments --- #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 04_mr_analyses.R <filtered_SNPs_for_MR.csv>")
}
input_file <- args[1]

# --- Load CSV --- #
filtered_snps <- read.csv(input_file, header = TRUE)
filtered_snps

# --- Check it has been imported correctly --- #
dim(filtered_snps)
print(colnames(filtered_snps))

# --- Ensure column names are clean --- #
colnames(filtered_snps) <- gsub("\\s+", "", colnames(filtered_snps))
print(colnames(filtered_snps))

# --- Rename columns to match TwoSampleMR requirements --- #
colnames(filtered_snps)[colnames(filtered_snps) == "A1_exp"] <- "effect_allele"
colnames(filtered_snps)[colnames(filtered_snps) == "A2_exp"] <- "other_allele"
colnames(filtered_snps)[colnames(filtered_snps) == "EAF_exp"] <- "eaf"

colnames(filtered_snps)[colnames(filtered_snps) == "A1_out"] <- "effect_allele.outcome"
colnames(filtered_snps)[colnames(filtered_snps) == "A2_out"] <- "other_allele.outcome"
colnames(filtered_snps)[colnames(filtered_snps) == "EAF_out"] <- "eaf.outcome"

# --- Ensure required columns exist --- #
required_cols <- c("SNP", "BETA_exp", "SE_exp", "PVALUE_exp", "effect_allele", "other_allele", "eaf")
missing_cols <- setdiff(required_cols, colnames(filtered_snps))

if (length(missing_cols) > 0) {
  stop(paste("ERROR: Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# --- Define Exposure and Outcome datasets --- #
exposure_dat <- na.omit(
  data.frame(
    SNP = filtered_snps$SNP,
    beta = as.numeric(filtered_snps$BETA_exp),
    se = as.numeric(filtered_snps$SE_exp),
    pval = as.numeric(filtered_snps$PVALUE_exp),
    eaf = as.numeric(filtered_snps$eaf),
    effect_allele = as.character(filtered_snps$effect_allele),
    other_allele = as.character(filtered_snps$other_allele)
  )
)

outcome_dat <- na.omit(
  data.frame(
    SNP = filtered_snps$SNP,
    beta = as.numeric(filtered_snps$BETA_out),
    se = as.numeric(filtered_snps$SE_out),
    pval = as.numeric(filtered_snps$PVALUE_out),
    eaf = as.numeric(filtered_snps$eaf.outcome),
    effect_allele = as.character(filtered_snps$effect_allele.outcome),
    other_allele = as.character(filtered_snps$other_allele.outcome)
  )
)

# --- Save formatted exposure and outcome data --- #
write.csv(exposure_dat, "exposure_dat.csv", row.names = FALSE)
write.csv(outcome_dat, "outcome_dat.csv", row.names = FALSE)

# --- Format data for MR analysis --- #
exposure_dat <- format_data(exposure_dat, type = "exposure")
outcome_dat <- format_data(outcome_dat, type = "outcome")

# --- Merge datasets manually (skip harmonise_data) --- #
harmonised_data <- merge(exposure_dat, outcome_dat, by = "SNP")
harmonised_data$mr_keep <- TRUE 

# --- Run MR methods --- #
res_ivw <- mr(harmonised_data, method = "mr_ivw")  
res_egger <- mr(harmonised_data, method = "mr_egger_regression")  
res_wmedian <- mr(harmonised_data, method = "mr_weighted_median")

# --- Sensitivity analyses --- #
heterogeneity <- mr_heterogeneity(harmonised_data)
pleiotropy <- mr_pleiotropy_test(harmonised_data)

# --- Compute IVW OR per SNP --- #
IVW_SNP_results <- data.frame(
  SNP = harmonised_data$SNP,
  IVW_OR = exp(harmonised_data$beta.exposure / harmonised_data$se.exposure),
  IVW_Lower_95 = exp((harmonised_data$beta.exposure - 1.96 * harmonised_data$se.exposure) / harmonised_data$se.exposure),
  IVW_Upper_95 = exp((harmonised_data$beta.exposure + 1.96 * harmonised_data$se.exposure) / harmonised_data$se.exposure)
)

# --- Create summary results table --- #
results <- data.frame(
  Exposure = "LDL-C",
  Outcome = "Alzheimer's Disease",
  N_SNPs = nrow(harmonised_data),
  
  IVW_OR = if (!is.null(res_ivw$b)) exp(res_ivw$b) else NA,
  IVW_Lower_95 = if (!is.null(res_ivw$b)) exp(res_ivw$b + qnorm(.025) * res_ivw$se) else NA,
  IVW_Upper_95 = if (!is.null(res_ivw$b)) exp(res_ivw$b + qnorm(.975) * res_ivw$se) else NA,
  IVW_Pval = if (!is.null(res_ivw$pval)) res_ivw$pval else NA,
  
  WM_OR = if (!is.null(res_wmedian$b)) exp(res_wmedian$b) else NA,
  WM_Lower_95 = if (!is.null(res_wmedian$b)) exp(res_wmedian$b + qnorm(.025) * res_wmedian$se) else NA,
  WM_Upper_95 = if (!is.null(res_wmedian$b)) exp(res_wmedian$b + qnorm(.975) * res_wmedian$se) else NA,
  WM_Pval = if (!is.null(res_wmedian$pval)) res_wmedian$pval else NA,
  
  Egger_OR = if (!is.null(res_egger$b)) exp(res_egger$b) else NA,
  Egger_Lower_95 = if (!is.null(res_egger$b)) exp(res_egger$b + qnorm(.025) * res_egger$se) else NA,
  Egger_Upper_95 = if (!is.null(res_egger$b)) exp(res_egger$b + qnorm(.975) * res_egger$se) else NA,
  Egger_Pval = if (!is.null(res_egger$pval)) res_egger$pval else NA,
  
  Egger_Intercept_Pval = pleiotropy$pval,
  IVW_Q_Statistic_P = ifelse("Inverse variance weighted" %in% heterogeneity$method, heterogeneity$Q_pval[heterogeneity$method == "Inverse variance weighted"], NA),
  I2_Statistic = if ("Inverse variance weighted" %in% heterogeneity$method && !is.null(heterogeneity$I2[heterogeneity$method == "Inverse variance weighted"])) heterogeneity$I2[heterogeneity$method == "Inverse variance weighted"] else NA,
  Mean_F_Statistic = mean(harmonised_data$F_stat, na.rm = TRUE)
)

# --- Save results --- #
write.csv(results, "MR_Formatted_Results.csv", row.names = FALSE)
write.csv(IVW_SNP_results, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)

print("MR Analysis completed successfully. Results saved to 'MR_Formatted_Results.csv' and 'MR_IVW_OR_Per_SNP.csv'.")
