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
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = repo)
}
install_if_missing("nloptr"); install_if_missing("lme4"); install_if_missing("meta")
install_if_missing("remotes"); install_if_missing("dplyr"); install_if_missing("readr")
install_if_missing("tibble");  install_if_missing("tidyr")

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never", lib = .libPaths()[1])
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript 04_mr_analyses.R <input_ld_pruned_snps.csv>")
input_file <- args[1]

cat("\n==============================================================\n")
cat("MR-CoPe | Mendelian Randomisation Analysis\n")
cat("==============================================================\n\n")

# ---- guard empty ----
if (!file.exists(input_file) || file.info(input_file)$size == 0) {
  cat("❌ Input missing/empty — skipping MR.\n")
  empty_df <- data.frame()
  write.csv(empty_df, "MR_Formatted_Results.csv", row.names = FALSE)
  write.csv(empty_df, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)
  write.csv(empty_df, "exposure_dat.csv", row.names = FALSE)
  write.csv(empty_df, "outcome_dat.csv", row.names = FALSE)
  write.csv(empty_df, "harmonised_data.csv", row.names = FALSE)
  quit(status = 0)
}

cat("📥 Loading input file:", input_file, "\n")
dat <- suppressMessages(read_csv(input_file, show_col_types = FALSE))
names(dat) <- toupper(gsub("\\s+", "", names(dat)))
cat("📊 SNPs in analysis:", nrow(dat), "\n\n")
cat("🧾 Columns:\n"); print(names(dat))

# ---- required columns ----
req <- c("SNP","BETA_EXP","SE_EXP","PVALUE_EXP","BETA_OUT","SE_OUT","PVALUE_OUT")
miss <- setdiff(req, names(dat))
if (length(miss)) stop(paste0("❌ Missing required columns: ", paste(miss, collapse = ", ")))

# ---- EAF mapping (be tolerant if absent) ----
if (!"EAF" %in% names(dat))          dat[["EAF"]] <- if ("EAF_EXP" %in% names(dat)) dat[["EAF_EXP"]] else NA_real_
if (!"EAF.OUTCOME" %in% names(dat))  dat[["EAF.OUTCOME"]] <- if ("EAF_OUT" %in% names(dat)) dat[["EAF_OUT"]] else NA_real_

# ---- allele names for harmonisation ----
if (all(c("A1_EXP","A2_EXP","A1_OUT","A2_OUT") %in% names(dat))) {
  dat <- dat %>%
    rename(effect_allele = A1_EXP,
           other_allele = A2_EXP,
           `effect_allele.outcome` = A1_OUT,
           `other_allele.outcome`  = A2_OUT)
  cat("🔄 Renamed allele columns for harmonisation.\n\n")
} else {
  cat("⚠️  A1/A2 columns not fully present; proceeding anyway.\n\n")
}

# ---- build exposure/outcome tables ----
exposure_dat <- dat %>%
  select(SNP,
         beta = BETA_EXP, se = SE_EXP, pval = PVALUE_EXP,
         eaf = EAF, effect_allele, other_allele) %>%
  drop_na(SNP, beta, se, pval)

outcome_dat <- dat %>%
  select(SNP,
         beta = BETA_OUT, se = SE_OUT, pval = PVALUE_OUT,
         eaf = `EAF.OUTCOME`,
         effect_allele = `effect_allele.outcome`,
         other_allele  = `other_allele.outcome`) %>%
  drop_na(SNP, beta, se, pval)

write.csv(exposure_dat, "exposure_dat.csv", row.names = FALSE)
write.csv(outcome_dat,  "outcome_dat.csv",  row.names = FALSE)

# ---- harmonise ----
exposure_dat  <- format_data(exposure_dat, type = "exposure")
outcome_dat   <- format_data(outcome_dat,  type = "outcome")
harm <- harmonise_data(exposure_dat, outcome_dat, action = 2)
write.csv(harm, "harmonised_data.csv", row.names = FALSE)

cat("📊 SNPs after harmonisation:", nrow(harm), "\n\n")
if (nrow(harm) == 0) {
  cat("❌ No data after harmonisation — stopping.\n")
  empty_df <- data.frame()
  write.csv(empty_df, "MR_Formatted_Results.csv", row.names = FALSE)
  write.csv(empty_df, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)
  quit(status = 0)
}

# ---- MR main methods ----
cat("⚙️ Running MR methods (IVW, Egger, Weighted Median)...\n")
methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
res_df <- mr(harm, method_list = methods)
het    <- mr_heterogeneity(harm)
pleio  <- mr_pleiotropy_test(harm)

# ---- Per-SNP ORs (robust block) ----
per_snp_or <- NULL
ok <- TRUE
snp_try <- try(mr_singlesnp(harm), silent = TRUE)
if (!inherits(snp_try, "try-error") && is.data.frame(snp_try) && "SNP" %in% names(snp_try)) {
  if ("b" %in% names(snp_try) && "se" %in% names(snp_try)) {
    keep_idx <- rep(TRUE, nrow(snp_try))
    if ("method" %in% names(snp_try)) keep_idx <- grepl("Wald", snp_try$method, ignore.case = TRUE)
    snp_sub <- snp_try[keep_idx, , drop = FALSE]
    per_snp_or <- tibble(
      SNP          = snp_sub$SNP,
      IVW_OR       = exp(snp_sub$b),
      IVW_Lower_95 = exp(snp_sub$b + qnorm(0.025) * snp_sub$se),
      IVW_Upper_95 = exp(snp_sub$b + qnorm(0.975) * snp_sub$se)
    )
  } else ok <- FALSE
} else ok <- FALSE

if (!ok) {
  # Manual Wald ratio per SNP (delta-method SE)
  eps <- 1e-12
  hd <- harm %>%
    transmute(
      SNP,
      be = beta.exposure,
      se_be = se.exposure,
      bo = beta.outcome,
      se_bo = se.outcome
    ) %>%
    filter(is.finite(be), abs(be) > eps,
           is.finite(se_be), is.finite(bo), is.finite(se_bo))

  b  <- hd$bo / hd$be
  # FIXED LINE (balanced parentheses):
  se <- sqrt((hd$se_bo^2 / (hd$be^2)) + ((hd$bo^2) * (hd$se_be^2) / (hd$be^4)))

  per_snp_or <- tibble(
    SNP          = hd$SNP,
    IVW_OR       = exp(b),
    IVW_Lower_95 = exp(b + qnorm(0.025) * se),
    IVW_Upper_95 = exp(b + qnorm(0.975) * se)
  )
}
write.csv(per_snp_or, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)

# ---- Summary table ----
get_row <- function(d, m) d %>% filter(method == m)
safe_exp <- function(x) if (!length(x)) NA_real_ else exp(x)
safe_ci  <- function(b, se) if (!length(b) || !length(se)) c(NA_real_, NA_real_) else c(b + qnorm(.025)*se, b + qnorm(.975)*se)

ivw   <- get_row(res_df, "Inverse variance weighted")
wm    <- get_row(res_df, "Weighted median")
egger <- get_row(res_df, "MR Egger")

ivw_ci   <- safe_ci(ivw$b, ivw$se)
wm_ci    <- safe_ci(wm$b,  wm$se)
egger_ci <- safe_ci(egger$b, egger$se)

results <- tibble(
  N_SNPs = nrow(harm),

  IVW_OR       = safe_exp(ivw$b),
  IVW_CI_Lower = safe_exp(ivw_ci[1]),
  IVW_CI_Upper = safe_exp(ivw_ci[2]),
  IVW_Pval     = if (nrow(ivw)) ivw$pval else NA_real_,

  WM_OR        = safe_exp(wm$b),
  WM_CI_Lower  = safe_exp(wm_ci[1]),
  WM_CI_Upper  = safe_exp(wm_ci[2]),
  WM_Pval      = if (nrow(wm)) wm$pval else NA_real_,

  Egger_OR       = safe_exp(egger$b),
  Egger_CI_Lower = safe_exp(egger_ci[1]),
  Egger_CI_Upper = safe_exp(egger_ci[2]),
  Egger_Pval     = if (nrow(egger)) egger$pval else NA_real_,

  Egger_Intercept_Pval = if (!is.null(pleio$pval)) pleio$pval else NA_real_,
  IVW_Het_Q_Pval = {
    hv <- het %>% filter(method == "Inverse variance weighted")
    if (nrow(hv)) hv$Q_pval else NA_real_
  },
  I2_Statistic = {
    hv <- het %>% filter(method == "Inverse variance weighted")
    if (nrow(hv)) hv$I2 else NA_real_
  }
)

write.csv(results, "MR_Formatted_Results.csv", row.names = FALSE)

cat("✅ MR analysis completed successfully!\n")
cat("📝 Outputs:\n")
cat("- exposure_dat.csv\n- outcome_dat.csv\n- harmonised_data.csv\n- MR_Formatted_Results.csv\n- MR_IVW_OR_Per_SNP.csv\n")
cat("==============================================================\n\n")
