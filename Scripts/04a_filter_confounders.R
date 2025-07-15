#!/usr/bin/env Rscript

# =============================================================================
# MR-CoPe | PhenoScanner-Based Confounder Filtering (Strict Mode, Batch, Robust)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(phenoscanner)
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
filter_log_file <- sub("\\.csv$", "_filter_log.csv", output_file)

cat(sprintf("\n=== MR-CoPe PhenoScanner Strict Mode Filter ===\n"))
cat(sprintf("Input:      %s\n", input_file))
cat(sprintf("Output:     %s\n", output_file))
cat(sprintf("Keyword(s): %s\n", target_keyword))
cat(sprintf("Filter log: %s\n\n", filter_log_file))

if (!file.exists(input_file)) {
  cat("❌ ERROR: Input file not found:", input_file, "\n")
  quit(status = 1)
}

gwas <- read_csv(input_file, show_col_types = FALSE)

# Ensure SNP col
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
  cat("⚠️  Input file has 0 rows. Writing empty output and filter log.\n")
  file.create(output_file)
  write_csv(tibble(SNP=character(), removed=logical(), reason=character()), filter_log_file)
  quit(status = 0)
}

snps <- unique(gwas$SNP)
if (length(snps) > 10000) {
  cat(sprintf("⚠️  WARNING: You are querying %d SNPs. PhenoScanner batch may fail above 10,000. Consider splitting.\n", length(snps)))
}

cat("🔬 Querying", length(snps), "SNPs via phenoscanner...\n")
cat("🎯 Target trait pattern:", target_keyword, "\n\n")

ps_res <- tryCatch(
  phenoscanner::phenoscanner(snpquery = snps, pvalue = 1),
  error = function(e) {
    cat("❌ ERROR: Phenoscanner batch query failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (is.null(ps_res) || is.null(ps_res$results)) {
  cat("⚠️  No PhenoScanner data for any SNP. Keeping all SNPs (check connection/logs).\n")
  write_csv(gwas, output_file)
  write_csv(tibble(SNP=snps, removed=FALSE, reason="PhenoScanner query failed or no hits"), filter_log_file)
  quit(status = 0)
}

ps_tab <- ps_res$results
ps_tab$trait_lc <- tolower(ps_tab$trait)

# Map: for each SNP, does ANY trait fail to match keyword?
snps_flagged <- ps_tab %>%
  group_by(snp) %>%
  summarize(off_target = any(!grepl(target_keyword, trait_lc)),
            reason = paste(unique(trait_lc[!grepl(target_keyword, trait_lc)]), collapse="; ")) %>%
  ungroup()

snps_to_remove <- snps_flagged %>% filter(off_target) %>% pull(snp)
cat("🧹 SNPs flagged for removal (pleiotropy):", length(snps_to_remove), "\n")

# Filtered output
gwas_filtered <- gwas %>% filter(!SNP %in% snps_to_remove)
write_csv(gwas_filtered, output_file)

# Write filter log for *all* input SNPs
filter_log <- tibble(
  SNP = snps,
  removed = snps %in% snps_to_remove,
  reason = ifelse(snps %in% snps_flagged$snp,
                  snps_flagged$reason[match(snps, snps_flagged$snp)],
                  "No PhenoScanner hit or all traits match")
)
write_csv(filter_log, filter_log_file)

cat("✅ SNPs retained:", nrow(gwas_filtered), "\n")
cat("💾 Output written to:", output_file, "\n")
cat("📝 Filter log written to:", filter_log_file, "\n")
cat("========================================================================\n")
