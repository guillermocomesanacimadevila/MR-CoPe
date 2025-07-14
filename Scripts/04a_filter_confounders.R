#!/usr/bin/env Rscript

# =============================================================================
# MR-CoPe | PhenoScanner-Based Confounder Filtering (Strict Mode)
# -----------------------------------------------------------------------------
# Description:
#   Removes SNPs associated with any trait that does not match the specified
#   exposure keyword(s). Keeps only "clean" instruments with no known pleiotropy.
#
# Usage:
#   Rscript 04a_filter_confounders.R <input_ld_pruned.csv> <output_filtered.csv> "<target_keyword>"
#
# Example:
#   Rscript 04a_filter_confounders.R ld_pruned_SNPs.csv filtered_SNPs.csv "LDL|cholesterol"
#
# Author: MR-CoPe Team (updated July 2025)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(httr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("❌ ERROR: Invalid number of arguments.\n")
  cat("Usage: Rscript 04a_filter_confounders.R <input_csv> <output_csv> '<target_keyword>'\n")
  quit(status = 1)
}

input_file  <- args[1]
output_file <- args[2]
target_keyword <- tolower(args[3])

if (!file.exists(input_file)) {
  cat("❌ ERROR: Input file not found:", input_file, "\n")
  quit(status = 1)
}

# -----------------------
# Load input SNPs
# -----------------------
cat("📥 Loading:", input_file, "\n")
gwas <- read_csv(input_file, show_col_types = FALSE)

# Robustly ensure 'SNP' column exists
if (!"SNP" %in% colnames(gwas)) {
  if ("rsid" %in% colnames(gwas)) {
    cat("🛠️  Renaming 'rsid' to 'SNP'...\n")
    gwas <- gwas %>% rename(SNP = rsid)
  } else {
    cat("❌ ERROR: Input file must contain a 'SNP' (or 'rsid') column.\n")
    quit(status = 1)
  }
}

if (nrow(gwas) == 0) {
  cat("⚠️  Input file has 0 rows. Writing empty output and exiting.\n")
  file.create(output_file)
  quit(status = 0)
}

snps <- unique(gwas$SNP)
cat("🔬 SNPs to evaluate:", length(snps), "\n")
cat("🎯 Target trait pattern:", target_keyword, "\n\n")

if (length(snps) == 0) {
  cat("⚠️  No SNPs to query. Writing empty output and exiting.\n")
  file.create(output_file)
  quit(status = 0)
}

# -----------------------
# Query PhenoScanner API
# -----------------------
query_phenoscanner <- function(snp) {
  tryCatch({
    res <- GET("http://www.phenoscanner.medschl.cam.ac.uk/api",
               query = list(query = snp, catalogue = "GWAS", p = 5e-8, build = 37, proxies = "no"),
               timeout(15))
    if (status_code(res) == 200 && grepl("^SNP", content(res, "text"))) {
      read_tsv(content(res, "text"), show_col_types = FALSE)
    } else {
      NULL
    }
  }, error = function(e) {
    cat("⚠️  Warning: Failed to query SNP", snp, "|", conditionMessage(e), "\n")
    NULL
  })
}

# -----------------------
# Filter SNPs strictly
# -----------------------
snps_to_remove <- c()

for (i in seq_along(snps)) {
  snp <- snps[i]
  cat(sprintf("[%03d/%03d] 🔎 Checking SNP: %s\n", i, length(snps), snp))

  dat <- query_phenoscanner(snp)
  Sys.sleep(0.5)  # Respect API limits

  if (!is.null(dat) && "trait" %in% colnames(dat)) {
    traits <- tolower(dat$trait)

    # STRICT mode: all traits must match the exposure keyword
    if (!all(grepl(target_keyword, traits))) {
      cat("🚫 Removing due to off-target traits:", snp, "\n")
      snps_to_remove <- c(snps_to_remove, snp)
    } else {
      cat("✅ Retained (only associated with exposure):", snp, "\n")
    }
  }
}

# -----------------------
# Save filtered results
# -----------------------
gwas_filtered <- gwas %>% filter(!SNP %in% snps_to_remove)

# Ensure the output file always has a 'SNP' column
if (!"SNP" %in% colnames(gwas_filtered)) {
  if ("rsid" %in% colnames(gwas_filtered)) {
    gwas_filtered <- gwas_filtered %>% rename(SNP = rsid)
  }
}

write_csv(gwas_filtered, output_file)

cat("\n🧹 SNPs removed:", length(snps_to_remove), "\n")
cat("✅ SNPs retained:", nrow(gwas_filtered), "\n")
cat("💾 Output written to:", output_file, "\n")
cat("========================================================================\n")
