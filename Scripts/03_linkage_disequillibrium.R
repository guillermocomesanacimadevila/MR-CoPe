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
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  cat("\n❌ ERROR: Incorrect number of arguments.\n")
  cat("Usage: Rscript 03_linkage_disequillibrium.R <input> <output> <clump_kb> <clump_r2>\n\n")
  quit(status = 1)
}

input_file  <- args[1]
output_file <- args[2]
clump_kb    <- as.numeric(args[3])
clump_r2    <- as.numeric(args[4])

cat("\n==============================================================\n")
cat("MR-CoPe | Linkage Disequilibrium Pruning\n")
cat("==============================================================\n")

# ------------------ Load Data ------------------

if (!file.exists(input_file)) {
  stop(paste0("❌ Input file not found: ", input_file))
}

cat("📥 Reading input from:\n  →", input_file, "\n")
gwas <- read_csv(input_file, show_col_types = FALSE)

cat("📊 SNPs loaded:", nrow(gwas), "\n\n")

# ------------------ Check Required Columns ------------------

if (!"SNP" %in% names(gwas) || !"PVALUE_exp" %in% names(gwas)) {
  stop("❌ Input must contain columns: 'SNP' and 'PVALUE_exp'.")
}

# ------------------ Build clump_input ------------------

has_chr_pos <- all(c("CHR_exp", "BP_exp") %in% names(gwas))

if (!has_chr_pos) {
  cat("⚠️ No CHR/BP info found. Will rely on default LD panel via API lookup.\n")
}

clump_input <- gwas %>%
  transmute(
    rsid = SNP,
    pval = PVALUE_exp,
    id   = "exposure",
    chr  = if (has_chr_pos) CHR_exp else NA,
    pos  = if (has_chr_pos) BP_exp else NA
  )

cat("🧪 Sample SNPs going into ld_clump():\n")
print(head(clump_input))

# ------------------ Run LD Clumping ------------------

cat("🔗 Clumping with:\n")
cat("   ➤ Window :", clump_kb, "kb\n")
cat("   ➤ R²     :", clump_r2, "\n")
cat("   ➤ P      : 1.0 (keep all SNPs)\n\n")

clumped <- tryCatch({
  ld_clump(clump_input, clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = 1.0)
}, error = function(e) {
  cat("❌ Clumping failed:\n", e$message, "\n")
  quit(status = 1)
})

# ------------------ Filter Original ------------------

if (nrow(clumped) == 0) {
  cat("⚠️ No SNPs passed LD clumping.\n")
  cat("💾 Writing empty file so pipeline can skip MR.\n")
  file.create(output_file)
  quit(status = 0)
}

gwas_pruned <- gwas %>% filter(SNP %in% clumped$rsid)

cat("✅ SNPs after clumping:", nrow(gwas_pruned), "\n")
cat("🚫 SNPs removed due to LD:", nrow(gwas) - nrow(gwas_pruned), "\n")

cat("💾 Saving to:\n  →", output_file, "\n")
write_csv(gwas_pruned, output_file)

cat("✅ LD pruning complete.\n")
cat("==============================================================\n\n")
