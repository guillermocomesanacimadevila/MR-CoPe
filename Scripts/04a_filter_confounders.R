#!/usr/bin/env Rscript

# =============================================================================
# MR-CoPe | Automated SNP Filtering Using PhenoScanner
# -----------------------------------------------------------------------------
# Description:
#   This script queries the PhenoScanner V2 API to remove SNPs that are
#   significantly associated (P < 5e-8) with traits unrelated to a given
#   target exposure. This helps reduce horizontal pleiotropy in MR analyses.
#
# Usage:
#   Rscript 04a_filter_confounders.R <input_csv> <output_csv> "<target_trait>"
#
# Example:
#   Rscript 04a_filter_confounders.R ld_pruned_SNPs.csv filtered_snps.csv "CRP"
#
# Requirements:
#   - Input CSV must include a 'SNP' column (e.g. from LD clumping output)
#   - Output is a filtered version of the same file
#
# Author: MR-CoPe Team (2025)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(httr)
})

# ----------------------- Argument Parsing -----------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("❌ ERROR: Invalid number of arguments.\n")
  cat("Usage: Rscript 04a_filter_confounders.R <input_csv> <output_csv> '<target_trait>'\n")
  quit(status = 1)
}

input_file  <- args[1]
output_file <- args[2]
target_trait <- tolower(args[3])

if (!file.exists(input_file)) {
  cat("❌ ERROR: Input file not found:", input_file, "\n")
  quit(status = 1)
}

# ----------------------- Load Input -----------------------
cat("📥 Loading input SNP data from:", input_file, "\n")
gwas <- read_csv(input_file, show_col_types = FALSE)

if (!"SNP" %in% colnames(gwas)) {
  cat("❌ ERROR: Input file must contain a 'SNP' column.\n")
  quit(status = 1)
}

snps <- unique(gwas$SNP)
cat("🔬 SNPs to evaluate:", length(snps), "\n")
cat("🎯 Target trait to retain:", target_trait, "\n\n")

# ----------------------- Query PhenoScanner -----------------------
query_phenoscanner <- function(snp) {
  tryCatch({
    res <- GET("http://www.phenoscanner.medschl.cam.ac.uk/api",
               query = list(
                 query = snp,
                 catalogue = "GWAS",
                 build = 37,
                 p = 5e-8,
                 proxies = "no"
               ),
               timeout(15))
    if (status_code(res) != 200 || !grepl("^SNP", content(res, "text"))) return(NULL)
    read_tsv(content(res, "text"), show_col_types = FALSE)
  }, error = function(e) {
    cat("⚠️ Warning: Failed to query SNP", snp, "|", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ----------------------- Filter Logic -----------------------
snps_to_remove <- c()

for (i in seq_along(snps)) {
  snp <- snps[i]
  cat(sprintf("[%03d/%03d] 🔎 Querying SNP: %s\n", i, length(snps), snp))

  dat <- query_phenoscanner(snp)
  Sys.sleep(0.5)  # API-friendly rate

  if (!is.null(dat) && "trait" %in% names(dat)) {
    traits <- tolower(dat$trait)
    if (!any(grepl(target_trait, traits))) {
      cat("🚫 Removing SNP associated with unrelated traits:", snp, "\n")
      snps_to_remove <- c(snps_to_remove, snp)
    } else {
      cat("✅ SNP retained (matches target trait):", snp, "\n")
    }
  }
}

# ----------------------- Output Filtered File -----------------------
filtered_gwas <- gwas %>% filter(!SNP %in% snps_to_remove)
write_csv(filtered_gwas, output_file)

cat("\n========================================\n")
cat("🧹 SNPs removed due to off-target traits:", length(snps_to_remove), "\n")
cat("✅ SNPs retained for MR analysis       :", nrow(filtered_gwas), "\n")
cat("💾 Output written to:", output_file, "\n")
cat("========================================\n")
