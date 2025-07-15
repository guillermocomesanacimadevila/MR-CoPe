#!/usr/bin/env Rscript

# ===================================================================
# MR-CoPe | Linkage Disequilibrium Pruning (Real LD Clumping)
# ===================================================================
# Authors: Guillermo Comesaña & Christian Pepler, Updated by ChatGPT
#
# Description:
# Performs LD pruning using real linkage disequilibrium structure via
# the TwoSampleMR::clump_data() function, using the user-supplied 1000 Genomes pop panel.
#
# Usage:
#   Rscript 03_linkage_disequillibrium.R <input_filtered_SNPs.csv> <output_ld_pruned.csv> <clump_kb> <clump_r2> <ld_pop>
#
# Notes:
# - Requires 'SNP' and 'PVALUE_exp' columns in the input file (or 'rsid', which will be renamed).
# ===================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  cat("\n❌ ERROR: Incorrect number of arguments.\n")
  cat("Usage: Rscript 03_linkage_disequillibrium.R <input> <output> <clump_kb> <clump_r2> <ld_pop>\n\n")
  quit(status = 1)
}

input_file  <- args[1]
output_file <- args[2]
clump_kb    <- as.numeric(args[3])
clump_r2    <- as.numeric(args[4])
ld_pop      <- args[5]

cat("\n==============================================================\n")
cat("MR-CoPe | Linkage Disequilibrium Pruning\n")
cat("==============================================================\n")

# ------------------ Load Data ------------------

if (!file.exists(input_file)) {
  stop(paste0("❌ Input file not found: ", input_file))
}

cat("📥 Reading input from:\n  →", input_file, "\n")
gwas <- read_csv(input_file, show_col_types = FALSE)

# Handle empty input
if (nrow(gwas) == 0) {
  cat("⚠️  WARNING: Input file is empty. Writing empty output and exiting.\n")
  file.create(output_file)
  quit(status = 0)
}

# Robust column naming for SNPs
if (!"SNP" %in% names(gwas)) {
  if ("rsid" %in% names(gwas)) {
    cat("🛠️  Renaming 'rsid' to 'SNP'...\n")
    gwas <- gwas %>% rename(SNP = rsid)
  } else {
    stop("❌ Input must contain a 'SNP' (or 'rsid') column.")
  }
}

if (!"PVALUE_exp" %in% names(gwas)) {
  stop("❌ Input must contain column: 'PVALUE_exp'.")
}

cat("📊 SNPs loaded:", nrow(gwas), "\n\n")

# ------------------ Build clump_input ------------------

has_chr_pos <- all(c("CHR_exp", "BP_exp") %in% names(gwas))

if (!has_chr_pos) {
  cat("⚠️ No CHR/BP info found. Will rely on default LD panel via API lookup.\n")
}

clump_input <- gwas %>%
  transmute(
    SNP  = SNP,
    chr  = if (has_chr_pos) CHR_exp else NA,
    pos  = if (has_chr_pos) BP_exp else NA,
    pval = PVALUE_exp,
    id   = "exposure"
  )

cat("🧪 Sample SNPs going into clump_data():\n")
print(head(clump_input))

# ------------------ Run LD Clumping ------------------

cat("🔗 Clumping with:\n")
cat("   ➤ Window :", clump_kb, "kb\n")
cat("   ➤ R²     :", clump_r2, "\n")
cat("   ➤ Pop    :", ld_pop, "\n")
cat("   ➤ P      : 1.0 (keep all SNPs)\n\n")

clumped <- tryCatch({
  clump_data(
    dat      = clump_input,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p1 = 1.0,
    clump_p2 = 1.0,
    pop      = ld_pop
  )
}, error = function(e) {
  cat("❌ Clumping failed:\n", e$message, "\n")
  file.create(output_file)
  quit(status = 1)
})

# ------------------ Filter Original ------------------

if (nrow(clumped) == 0) {
  cat("⚠️ No SNPs passed LD clumping.\n")
  cat("💾 Writing empty file so pipeline can skip MR.\n")
  file.create(output_file)
  quit(status = 0)
}

gwas_pruned <- gwas %>% filter(SNP %in% clumped$SNP)

# =========== ENSURE COLUMN IS CALLED 'SNP' ==========
if (!"SNP" %in% colnames(gwas_pruned)) {
  if ("rsid" %in% colnames(gwas_pruned)) {
    gwas_pruned <- gwas_pruned %>% rename(SNP = rsid)
  } else {
    stop("❌ ERROR: Output file must contain a 'SNP' column (even after clumping)!")
  }
}
# =========== END FIX ==========

cat("✅ SNPs after clumping:", nrow(gwas_pruned), "\n")
cat("🚫 SNPs removed due to LD:", nrow(gwas) - nrow(gwas_pruned), "\n")

cat("💾 Saving to:\n  →", output_file, "\n")
write_csv(gwas_pruned, output_file)

cat("✅ LD pruning complete.\n")
cat("==============================================================\n\n")
