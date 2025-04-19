#!/usr/bin/env Rscript

# ========================================================================
# MR-CoPe | Confounder Filtering using PhenoScanner
# ========================================================================
# Author: Guillermo Comesaña & Christian Pepler
# Date: 2025
#
# Description:
# Identifies SNPs associated with potential confounding traits using PhenoScanner.
# Flags or removes SNPs with associations to specified phenotypes.
#
# Usage:
# Rscript 07_filter_confounders.R <input_snps.csv> <output_filtered.csv> <output_report.csv>
#
# Inputs:
# - input_snps.csv         : SNP dataset after LD pruning (must include 'SNP' column)
# - output_filtered.csv    : SNPs not linked to confounding traits
# - output_report.csv      : Full PhenoScanner result with flags
# ========================================================================

suppressPackageStartupMessages({
  library(phenoscanner)
  library(dplyr)
  library(stringr)
  library(readr)
})

# ---------------------------- Args ----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript 07_filter_confounders.R <input_snps.csv> <output_filtered.csv> <output_report.csv>")
}

input_file  <- args[1]
output_file <- args[2]
report_file <- args[3]

cat("\n===============================================================\n")
cat("MR-CoPe | Confounder Filtering with PhenoScanner\n")
cat("===============================================================\n\n")

# ---------------------------- Load SNPs ----------------------------
if (!file.exists(input_file)) {
  stop(paste("❌ ERROR: Input SNP file not found:", input_file))
}

snps <- read_csv(input_file, show_col_types = FALSE)

if (!"SNP" %in% names(snps)) {
  stop("❌ ERROR: Input file must contain a 'SNP' column.")
}

cat("📥 Loaded SNPs:", nrow(snps), "\n\n")

# ---------------------------- Confounder List ----------------------------
confounders <- c(
  "bmi", "body mass index", "obesity", "smoking", "alcohol",
  "physical activity", "education", "socioeconomic",
  "lipids", "cholesterol", "hdl", "ldl", "triglycerides",
  "crp", "inflammation", "il6", "white blood cell count"
)

# ---------------------------- PhenoScanner Query ----------------------------
cat("🔎 Querying PhenoScanner API...\n")
results <- phenoscanner(snps = snps$SNP, catalogue = "gwas", pvalue = 5e-5)

if (is.null(results$associations)) {
  stop("❌ ERROR: No results returned from PhenoScanner.")
}

associations <- results$associations %>%
  mutate(confounder_hit = str_detect(tolower(trait), str_c(confounders, collapse = "|")).
                                & p <= 5e-5)

flagged_snps <- unique(associations$SNP[associations$confounder_hit])

cat("🚫 Confounder-associated SNPs detected:", length(flagged_snps), "\n\n")

# ---------------------------- Filtering ----------------------------
snps_filtered <- snps %>%
  filter(!SNP %in% flagged_snps)

write_csv(snps_filtered, output_file)
write_csv(associations, report_file)

cat("✅ Filtered SNPs saved to:", output_file, "\n")
cat("📝 Full report saved to:", report_file, "\n")
cat("===============================================================\n\n")