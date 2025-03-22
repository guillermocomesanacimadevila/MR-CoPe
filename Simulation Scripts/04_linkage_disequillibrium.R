#!/usr/bin/env Rscript
# Allocate - 0, 1, 2 depending on genotype to each SNP 
# Loop through A1 and A2 - if df["A1"] == df["A2"] - label [0-2]

# LD Filtering Script using simulated genotypes based on EAF
# Input: filtered SNPs after gwas_processing.py
# Output: LD-pruned SNP list (r² < 0.001)

# --- Load Required Packages --- #
required_packages <- c("LDcorSV")
missing <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org/")
}
library(LDcorSV)

# --- Handle Args --- #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 04_linkage_disequillibrium.r <input_filtered_SNPs.csv> <output_ld_pruned.csv>")
}
input_file <- args[1]
output_file <- args[2]

cat("📥 Loading input:", input_file, "\n")
gwas <- read.csv(input_file)
cat("📊 Input rows:", nrow(gwas), "\n")

# --- Simulate genotype dosages using EAF --- #
set.seed(42)
simulate_genotypes <- function(eaf, n = 500) {
  p <- eaf
  geno_probs <- c((1 - p)^2, 2 * p * (1 - p), p^2)
  sample(0:2, size = n, replace = TRUE, prob = geno_probs)
}

cat("🧬 Simulating genotypes based on EAF...\n")
geno_matrix <- sapply(gwas$EAF_exp, simulate_genotypes)
colnames(geno_matrix) <- gwas$SNP
geno_matrix <- t(geno_matrix)

# --- Compute LD Matrix --- #
cat("🔗 Computing LD matrix...\n")
ld_matrix <- LDcorSV:::LD.cor(geno_matrix, snp.data = TRUE)  # ← fix here
ld_matrix[upper.tri(ld_matrix, diag = TRUE)] <- 0

# --- Prune SNPs with r² >= 0.001 --- #
keep_snps <- which(apply(ld_matrix < 0.001, 2, all))
gwas_pruned <- gwas[keep_snps, ]

cat("✅ Retained", nrow(gwas_pruned), "SNPs after LD pruning\n")
write.csv(gwas_pruned, output_file, row.names = FALSE)
cat("💾 Output written to:", output_file, "\n")
