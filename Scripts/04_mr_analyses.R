#!/usr/bin/env Rscript

# ================================================================
# MR-CoPe | Mendelian Randomisation Analysis using TwoSampleMR
# + MR-PRESSO integration (https://github.com/rondolab/MR-PRESSO)
# ================================================================
# Usage:
#   Rscript 04_mr_analyses.R <input_ld_pruned_snps.csv>
#
# Inputs:
#   - <input_ld_pruned_snps.csv>: output of LD step (or no-LD copy)
#
# Outputs:
#   - exposure_dat.csv
#   - outcome_dat.csv
#   - harmonised_data.csv
#   - MR_Formatted_Results.csv          (now includes PRESSO fields)
#   - MR_IVW_OR_Per_SNP.csv
#   - MR_PRESSO_Summary.csv             (Global p, #outliers, Distortion p)
#   - MR_PRESSO_Outliers.csv            (list of outlier SNPs; empty if none)
#   - (optional) MR_PRESSO_Outlier_Plot.png (placeholder plot if outliers)
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
install_if_missing("dplyr")
install_if_missing("readr")
install_if_missing("tibble")
install_if_missing("tidyr")

# TwoSampleMR (use CRAN if available, else GitHub)
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never", lib = .libPaths()[1])
}

# MR-PRESSO (GitHub)
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  remotes::install_github("rondolab/MR-PRESSO", upgrade = "never", lib = .libPaths()[1])
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MRPRESSO)
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
cat("MR-CoPe | Mendelian Randomisation Analysis (+ MR-PRESSO)\n")
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
  write.csv(tibble(GlobalTest_P=NA_real_, NbOutliers=NA_integer_, Distortion_P=NA_real_),
            "MR_PRESSO_Summary.csv", row.names = FALSE)
  write.csv(data.frame(), "MR_PRESSO_Outliers.csv", row.names = FALSE)
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
  write.csv(tibble(GlobalTest_P=NA_real_, NbOutliers=NA_integer_, Distortion_P=NA_real_),
            "MR_PRESSO_Summary.csv", row.names = FALSE)
  write.csv(data.frame(), "MR_PRESSO_Outliers.csv", row.names = FALSE)
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
  se <- sqrt((hd$se_bo^2 / (hd$be^2)) + ((hd$bo^2) * (hd$se_be^2) / (hd$be^4)))

  per_snp_or <- tibble(
    SNP          = hd$SNP,
    IVW_OR       = exp(b),
    IVW_Lower_95 = exp(b + qnorm(0.025) * se),
    IVW_Upper_95 = exp(b + qnorm(0.975) * se)
  )
}
write.csv(per_snp_or, "MR_IVW_OR_Per_SNP.csv", row.names = FALSE)

# ---- Build core summary (IVW / WM / Egger) ----
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

# ---- MR-PRESSO (pleiotropy / outlier correction) ----
presso_summary <- tibble(
  GlobalTest_P   = NA_real_,
  NbOutliers     = NA_integer_,
  Distortion_P   = NA_real_
)

presso_adj <- tibble(
  IVW_OR_PressoAdj    = NA_real_,
  IVW_CI_L_PressoAdj  = NA_real_,
  IVW_CI_U_PressoAdj  = NA_real_,
  IVW_P_PressoAdj     = NA_real_
)

try({
  # run only if enough SNPs and required columns for PRESSO exist
  if (nrow(harm) >= 4 &&
      all(c("beta.outcome","beta.exposure","se.outcome","se.exposure") %in% names(harm))) {
    set.seed(1234)
    pr <- mr_presso(
      BetaOutcome     = "beta.outcome",
      BetaExposure    = "beta.exposure",
      SdOutcome       = "se.outcome",
      SdExposure      = "se.exposure",
      OUTLIERtest     = TRUE,
      DISTORTIONtest  = TRUE,
      NbDistribution  = 1000,             # consider 5000+ for final reports
      SignifThreshold = 0.05,
      data = harm
    )

    # Global test p
    gt <- tryCatch(pr$`Main MR results`$`Global Test`$Pvalue, error = function(e) NA_real_)
    # Outliers
    ot <- tryCatch(pr$`Outlier Test`, error = function(e) NULL)
    n_out <- if (!is.null(ot) && is.data.frame(ot)) nrow(ot) else 0L
    # Distortion test p
    dp <- tryCatch(pr$`Distortion Test`$Pvalue, error = function(e) NA_real_)

    presso_summary <- tibble(
      GlobalTest_P = suppressWarnings(as.numeric(gt)),
      NbOutliers   = as.integer(n_out),
      Distortion_P = suppressWarnings(as.numeric(dp))
    )

    # PRESSO-adjusted causal estimate (beta → OR)
    est <- tryCatch(pr$`Main MR results`$`Causal Estimate`, error = function(e) NA_real_)
    se  <- tryCatch(pr$`Main MR results`$Sd,                error = function(e) NA_real_)
    pv  <- tryCatch(pr$`Main MR results`$Pvalue,            error = function(e) NA_real_)

    if (is.finite(est)[1] && is.finite(se)[1] && is.finite(pv)[1]) {
      b <- as.numeric(est)[1]
      s <- as.numeric(se)[1]
      p <- as.numeric(pv)[1]
      presso_adj <- tibble(
        IVW_OR_PressoAdj    = exp(b),
        IVW_CI_L_PressoAdj  = exp(b + qnorm(0.025)*s),
        IVW_CI_U_PressoAdj  = exp(b + qnorm(0.975)*s),
        IVW_P_PressoAdj     = p
      )
    }

    # Persist PRESSO CSVs for the HTML report
    write.csv(presso_summary, "MR_PRESSO_Summary.csv", row.names = FALSE)
    if (!is.null(ot) && is.data.frame(ot)) {
      write.csv(ot, "MR_PRESSO_Outliers.csv", row.names = FALSE)
      # Optional tiny outlier plot (placeholder)
      try({
        if (nrow(ot) > 0 && "SNP" %in% names(ot)) {
          png("MR_PRESSO_Outlier_Plot.png", width = 1100, height = 700, res = 150)
          par(mar = c(5,4,2,1))
          plot(seq_len(nrow(ot)), rep(0, nrow(ot)), pch = 19,
               xlab = "Outlier index", ylab = "Residual (placeholder)",
               main = "MR-PRESSO: Outlier SNPs")
          text(seq_len(nrow(ot)), rep(0, nrow(ot)), labels = ot$SNP, pos = 3, cex = 0.7)
          abline(h = 0, lty = 2, col = "gray50")
          dev.off()
        }
      }, silent = TRUE)
    } else {
      write.csv(data.frame(), "MR_PRESSO_Outliers.csv", row.names = FALSE)
    }
  } else {
    write.csv(presso_summary, "MR_PRESSO_Summary.csv", row.names = FALSE)
    write.csv(data.frame(),    "MR_PRESSO_Outliers.csv", row.names = FALSE)
  }
}, silent = TRUE)

# if PRESSO summary wasn't written (e.g., error), ensure files exist
if (!file.exists("MR_PRESSO_Summary.csv")) {
  write.csv(presso_summary, "MR_PRESSO_Summary.csv", row.names = FALSE)
}
if (!file.exists("MR_PRESSO_Outliers.csv")) {
  write.csv(data.frame(), "MR_PRESSO_Outliers.csv", row.names = FALSE)
}

# ---- Merge PRESSO into main summary (safe if NA) ----
results <- results %>%
  mutate(
    PRESSO_Global_P    = presso_summary$GlobalTest_P[1],
    PRESSO_N_Outliers  = presso_summary$NbOutliers[1],
    PRESSO_Distort_P   = presso_summary$Distortion_P[1],
    IVW_OR_PressoAdj   = presso_adj$IVW_OR_PressoAdj[1],
    IVW_CI_L_PressoAdj = presso_adj$IVW_CI_L_PressoAdj[1],
    IVW_CI_U_PressoAdj = presso_adj$IVW_CI_U_PressoAdj[1],
    IVW_P_PressoAdj    = presso_adj$IVW_P_PressoAdj[1],
    PRESSO_Distort_Pct = {
      adj <- presso_adj$IVW_OR_PressoAdj[1]
      base <- IVW_OR[1]
      if (is.na(adj) || is.na(base) || base == 0) NA_real_ else 100 * (adj - base) / base
    }
  )

# ---- Write outputs ----
write.csv(results, "MR_Formatted_Results.csv", row.names = FALSE)

cat("✅ MR analysis completed successfully!\n")
cat("📝 Outputs:\n")
cat("- exposure_dat.csv\n- outcome_dat.csv\n- harmonised_data.csv\n")
cat("- MR_Formatted_Results.csv (includes PRESSO)\n- MR_IVW_OR_Per_SNP.csv\n")
cat("- MR_PRESSO_Summary.csv\n- MR_PRESSO_Outliers.csv\n")
cat("==============================================================\n\n")
