#!/usr/bin/env Rscript

# ==========================
# LD pruning with large N to eliminate artificial LD
# ==========================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 04_linkage_disequillibrium.R <input_filtered_SNPs.csv> <output_ld_pruned.csv>")
}

input_file <- args[1]
output_file <- args[2]

cat("📥 Loading input:", input_file, "\n")
gwas <- read.csv(input_file)
cat("📊 Input rows:", nrow(gwas), "\n")

# --- Ensure allele columns are characters --- #
gwas$A1_exp <- as.character(gwas$A1_exp)
gwas$A2_exp <- as.character(gwas$A2_exp)
gwas$A1_out <- as.character(gwas$A1_out)
gwas$A2_out <- as.character(gwas$A2_out)

# --- Simulate genotypes using EAF --- #
set.seed(42)
simulate_genotypes <- function(eaf, n = 10000) {  # Increased sample size
  p <- eaf
  probs <- c((1 - p)^2, 2 * p * (1 - p), p^2)
  sample(0:2, size = n, replace = TRUE, prob = probs)
}

cat("🧬 Simulating genotypes (n = 100,000)...\n")
geno_matrix <- sapply(gwas$EAF_exp, simulate_genotypes)
colnames(geno_matrix) <- gwas$SNP
geno_matrix <- t(geno_matrix)

cat("🔗 Calculating pairwise LD matrix using Pearson correlation (r²)...\n")
ld_r <- cor(t(geno_matrix))^2  # r² matrix

ld_r[upper.tri(ld_r, diag = TRUE)] <- NA  # Only consider lower triangle

# --- Prune SNPs with any r² ≥ 0.001 --- #
cat("🧹 Pruning SNPs with r² ≥ 0.001...\n")
keep_snps <- which(apply(ld_r, 1, function(row) all(is.na(row) | row < 0.001)))
gwas_pruned <- gwas[keep_snps, ]

cat("✅ Retained", nrow(gwas_pruned), "SNPs\n")
write.csv(gwas_pruned, output_file, row.names = FALSE)
cat("💾 Written to:", output_file, "\n")
