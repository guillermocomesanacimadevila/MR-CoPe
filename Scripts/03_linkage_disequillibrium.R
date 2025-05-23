#!/usr/bin/env Rscript

# ===================================================================
# MR-CoPe | Linkage Disequilibrium Pruning (Real LD Clumping)
# ===================================================================
# Authors: Guillermo Comesaña & Christian Pepler
#
# Description:
# Performs LD pruning using real linkage disequilibrium structure via
# the TwoSampleMR::ld_clump() function, using the 1000 Genomes EUR panel.
#
# Usage:
#   Rscript 03_linkage_disequillibrium.R <input_filtered_SNPs.csv> <output_ld_pruned.csv> <clump_kb> <clump_r2>
#
# Notes:
# - Requires 'SNP' and 'PVALUE_exp' columns in the input file
# ===================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ieugwasr)
  library(readr)
})

# ------------------ Parse Arguments ------------------ #

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  cat("\n❌ ERROR: Incorrect number of arguments.\n")
  cat("Usage: Rscript 03_linkage_disequillibrium.R <input_filtered_SNPs.csv> <output_ld_pruned.csv> <clump_kb> <clump_r2>\n\n")
  quit(status = 1)
}

input_file  <- args[1]
output_file <- args[2]
clump_kb    <- as.numeric(args[3])
clump_r2    <- as.numeric(args[4])

# ------------------ Banner ------------------ #

cat("\n==================================================================\n")
cat("MR-CoPe | Linkage Disequilibrium Pruning (LD Clumping via TwoSampleMR)\n")
cat("==================================================================\n\n")

# ------------------ Load Data ------------------ #

cat("📥 Reading input summary statistics from:\n  →", input_file, "\n")

if (!file.exists(input_file)) {
  stop(paste0("❌ Input file not found: ", input_file))
}

gwas <- read_csv(input_file, show_col_types = FALSE)
cat("📊 SNPs loaded:", nrow(gwas), "\n\n")

# ------------------ Validate Required Columns ------------------ #

required_cols <- c("SNP", "PVALUE_exp")

missing <- setdiff(required_cols, colnames(gwas))
if (length(missing) > 0) {
  stop(paste0("❌ Missing required columns in input: ", paste(missing, collapse = ", ")))
}

# ------------------ Prepare Input for Clumping ------------------ #

cat("🧬 Preparing input for LD clumping...\n")
clump_input <- gwas %>%
  transmute(
    rsid = SNP,
    pval = PVALUE_exp,
    id = "exposure",   # Required by ld_clump
    chr = NA,
    pos = NA
  )

# ------------------ Perform Clumping ------------------ #

cat("🔗 Performing LD clumping with parameters:\n")
cat("   ➤ Window     :", clump_kb, "kb\n")
cat("   ➤ R² cutoff  :", clump_r2, "\n")
cat("   ➤ P-value    : 1.0 (include all SNPs)\n\n")

clumped <- tryCatch({
  ld_clump(
    clump_input,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = 1.0
  )
}, error = function(e) {
  cat("❌ ERROR during clumping:\n", e$message, "\n")
  quit(status = 1)
})

# ------------------ Filter and Save Results ------------------ #

if (nrow(clumped) == 0) {
  cat("⚠️ No SNPs passed LD clumping.\n")
  cat("💾 Writing empty output file so pipeline can gracefully skip MR.\n")
  file.create(output_file)
  quit(status = 0)
}

gwas_pruned <- gwas %>% filter(SNP %in% clumped$rsid)

cat("✅ SNPs retained after LD clumping:", nrow(gwas_pruned), "\n")
cat("🚫 SNPs removed due to LD:", nrow(gwas) - nrow(gwas_pruned), "\n")

cat("\n💾 Writing pruned dataset to:\n  →", output_file, "\n")
write_csv(gwas_pruned, output_file)

cat("\n✅ LD pruning complete.\n")
cat("==================================================================\n\n")
