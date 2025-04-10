#!/usr/bin/env Rscript

# ======================================================
# MR-CoPe | LD Pruning Script
# ======================================================
# Author: Guillermo Comesaña & Christian Pepler
# Date: 2025
#
# Description:
# Removes SNPs in linkage disequilibrium (LD) using simulated genotypes
# based on EAF, to avoid artificial LD in MR analyses.
#
# Usage:
# Rscript 03_ld_pruning.R <input_filtered_SNPs.csv> <output_ld_pruned.csv>
# ======================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript 04_ld_pruning.R <input_filtered_SNPs.csv> <output_ld_pruned.csv>")
}

input_file <- args[1]
output_file <- args[2]

cat("\n==============================================================\n")
cat("MR-CoPe | Linkage Disequilibrium Pruning\n")
cat("==============================================================\n\n")

cat("📥 Loading filtered SNPs file:", input_file, "\n")
gwas <- read.csv(input_file)
cat("📊 Input SNPs:", nrow(gwas), "\n\n")


# Ensure alleles are characters
allele_cols <- c("A1_exp", "A2_exp", "A1_out", "A2_out")
gwas[allele_cols] <- lapply(gwas[allele_cols], as.character)


# --- Simulate genotypes from EAF --- #
set.seed(42)

simulate_genotypes <- function(eaf, n = 100000) {
  probs <- c((1 - eaf)^2, 2 * eaf * (1 - eaf), eaf^2)
  sample(0:2, size = n, replace = TRUE, prob = probs)
}

cat("🧬 Simulating genotypes for LD calculation (n = 100,000)...\n")
geno_matrix <- sapply(gwas$EAF_exp, simulate_genotypes)
colnames(geno_matrix) <- gwas$SNP
geno_matrix <- t(geno_matrix)


# --- Compute pairwise LD (r²) matrix --- #
cat("🔗 Calculating pairwise LD matrix (r²)...\n")
ld_r2 <- cor(t(geno_matrix))^2
ld_r2[upper.tri(ld_r2, diag = TRUE)] <- NA  # Lower triangle only


# --- Prune SNPs with r² >= 0.001 --- #
cat("🧹 Pruning SNPs with any r² ≥ 0.001...\n")
keep_snps_idx <- which(apply(ld_r2, 1, function(row) all(is.na(row) | row < 0.001)))
gwas_pruned <- gwas[keep_snps_idx, ]

cat("✅ SNPs retained after LD pruning:", nrow(gwas_pruned), "\n")
cat("🚫 SNPs removed:", nrow(gwas) - nrow(gwas_pruned), "\n\n")

# --- Save output --- #
write.csv(gwas_pruned, output_file, row.names = FALSE)
cat("💾 LD-pruned SNP dataset saved to:", output_file, "\n")
cat("==============================================================\n\n")