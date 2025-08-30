#!/usr/bin/env Rscript

# ===================================================================
# MR-CoPe | Linkage Disequilibrium Pruning (Real LD Clumping)
# ===================================================================
# Authors: Guillermo Comesaña & Christian Pepler
#
# Description:
# Performs LD pruning using real linkage disequilibrium structure via
# the TwoSampleMR::clump_data() function, using the user-supplied
# 1000 Genomes population panel.
#
# Usage:
#   Rscript 03_linkage_disequillibrium.R <input_filtered_SNPs.csv> <output_ld_pruned.csv> <clump_kb> <clump_r2> <ld_pop>
#
# Notes:
# - Requires 'SNP' and 'PVALUE_EXP' columns in the input file (or 'rsid', which will be renamed).
# ===================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
})

# ------------------ Parse Arguments ------------------
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

cat("📥 Reading input from:\n  → ", input_file, "\n")
gwas <- read_csv(input_file, show_col_types = FALSE)

# Normalize all column names to uppercase for consistency
names(gwas) <- toupper(names(gwas))

# Handle empty input
if (nrow(gwas) == 0) {
  cat("⚠️  WARNING: Input file is empty. Writing empty output and exiting.\n")
  file.create(output_file)
  quit(status = 0)
}

# ------------------ Ensure Required Columns ------------------
if (!"SNP" %in% names(gwas)) {
  if ("RSID" %in% names(gwas)) {
    cat("🛠️  Renaming 'RSID' to 'SNP'...\n")
    gwas <- gwas %>% rename(SNP = RSID)
  } else {
    stop("❌ Input must contain a 'SNP' (or 'RSID') column.")
  }
}

if (!"PVALUE_EXP" %in% names(gwas)) {
  stop("❌ Input must contain column: 'PVALUE_EXP'.")
}

cat("📊 SNPs loaded: ", nrow(gwas), "\n\n")

# ------------------ Build clump_input ------------------
has_chr_pos <- all(c("CHR_EXP", "BP_EXP") %in% names(gwas))

if (!has_chr_pos) {
  cat("⚠️ No CHR/BP info found. Will rely on default LD panel via API lookup.\n")
}

clump_input <- gwas %>%
  transmute(
    SNP  = SNP,
    chr  = if (has_chr_pos) CHR_EXP else NA,
    pos  = if (has_chr_pos) BP_EXP else NA,
    pval = PVALUE_EXP,
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

# ------------------ Filter Original Data ------------------
if (nrow(clumped) == 0) {
  cat("⚠️ No SNPs passed LD clumping.\n")
  cat("💾 Writing empty file so pipeline can skip MR.\n")
  file.create(output_file)
  quit(status = 0)
}

gwas_pruned <- gwas %>% filter(SNP %in% clumped$SNP)

# Ensure SNP column exists
if (!"SNP" %in% colnames(gwas_pruned)) {
  if ("RSID" %in% colnames(gwas_pruned)) {
    gwas_pruned <- gwas_pruned %>% rename(SNP = RSID)
  } else {
    stop("❌ ERROR: Output file must contain a 'SNP' column (even after clumping)!")
  }
}

cat("✅ SNPs after clumping: ", nrow(gwas_pruned), "\n")
cat("🚫 SNPs removed due to LD: ", nrow(gwas) - nrow(gwas_pruned), "\n")

# ------------------ Save Output ------------------
cat("💾 Saving to:\n  → ", output_file, "\n")
write_csv(gwas_pruned, output_file)

cat("✅ LD pruning complete.\n")
cat("==============================================================\n\n")
