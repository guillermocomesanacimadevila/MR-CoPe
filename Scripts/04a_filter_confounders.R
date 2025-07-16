#!/usr/bin/env Rscript

# ================================================================
# MR-CoPe | Confounder Filtering via PhenoScanner
# ================================================================
# Author: Guillermo Comesaña & Christian Pepler
#
# Usage:
#   Rscript 04a_filter_confounders.R <input_csv> <output_csv> '<keyword1, keyword2, ...>'
#
# Description:
#   Filters out SNPs that are associated with traits outside of the
#   specified confounder-related keyword list. SNPs with any association
#   to a trait not matching one of the keywords are removed.
#
# Inputs:
# - input_csv:       GWAS file with a 'SNP' or 'rsid' column
# - output_csv:      Filtered GWAS file (SNPs retained only if clean)
# - trait keywords:  One or more trait keywords (comma-separated)
#
# Outputs:
# - Filtered GWAS summary statistics
# - Per-SNP filtering log
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(phenoscanner)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("❌ ERROR: Invalid number of arguments.\n")
  cat("Usage: Rscript 04a_filter_confounders.R <input_csv> <output_csv> '<keyword1, keyword2, ...>'\n")
  quit(status = 1)
}

input_file       <- args[1]
output_file      <- args[2]
raw_keywords     <- tolower(args[3])
filter_log_file  <- sub("\\.csv$", "_filter_log.csv", output_file)

# === Handle Skip Mode ===
if (raw_keywords == ".") {
  cat("⚠️  No confounder keywords provided. Skipping PhenoScanner filtering...\n")

  if (!file.exists(input_file)) {
    cat("❌ ERROR: Input file not found:", input_file, "\n")
    quit(status = 1)
  }

  gwas <- read_csv(input_file, show_col_types = FALSE)

  if (!"SNP" %in% colnames(gwas)) {
    if ("rsid" %in% colnames(gwas)) {
      gwas <- gwas %>% rename(SNP = rsid)
    } else {
      cat("❌ ERROR: Input must have a 'SNP' or 'rsid' column.\n")
      quit(status = 1)
    }
  }

  write_csv(gwas, output_file)
  write_csv(tibble(SNP = gwas$SNP, removed = FALSE, reason = "Filtering skipped ('.')"), filter_log_file)
  cat("✅ SNPs retained:", nrow(gwas), "\n")
  cat("💾 Output written to:", output_file, "\n")
  cat("📝 Filter log written to:", filter_log_file, "\n")
  cat("====================================================================\n")
  quit(status = 0)
}

# === Continue with Filtering Mode ===
keywords <- unlist(strsplit(raw_keywords, ",\\s*"))

cat(sprintf("\n=== MR-CoPe PhenoScanner Confounder Filter ===\n"))
cat(sprintf("Input:       %s\n", input_file))
cat(sprintf("Output:      %s\n", output_file))
cat(sprintf("Keywords:    %s\n", paste(keywords, collapse=", ")))
cat(sprintf("\nFilter log:  %s\n\n", filter_log_file))

if (!file.exists(input_file)) {
  cat("❌ ERROR: Input file not found:", input_file, "\n")
  quit(status = 1)
}

gwas <- read_csv(input_file, show_col_types = FALSE)

if (!"SNP" %in% colnames(gwas)) {
  if ("rsid" %in% colnames(gwas)) {
    cat("🛠️  Renaming 'rsid' to 'SNP'...\n")
    gwas <- gwas %>% rename(SNP = rsid)
  } else {
    cat("❌ ERROR: Input must have a 'SNP' or 'rsid' column.\n")
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
  cat(sprintf("⚠️  WARNING: You are querying %d SNPs. PhenoScanner may fail above 10,000.\n", length(snps)))
}

cat("🔬 Querying", length(snps), "SNPs via PhenoScanner...\n\n")

ps_res <- tryCatch(
  phenoscanner::phenoscanner(snpquery = snps, pvalue = 1),
  error = function(e) {
    cat("❌ ERROR: PhenoScanner query failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (is.null(ps_res) || is.null(ps_res$results) || nrow(ps_res$results) == 0) {
  cat("⚠️  No PhenoScanner data returned. Keeping all SNPs.\n")
  write_csv(gwas, output_file)
  write_csv(tibble(SNP=snps, removed=FALSE, reason="No PhenoScanner results"), filter_log_file)
  quit(status = 0)
}

ps_tab <- ps_res$results
ps_tab$trait_lc <- tolower(ps_tab$trait)

# ---- Check if each trait matches ANY keyword ----
ps_tab$matched <- sapply(ps_tab$trait_lc, function(tr) {
  any(sapply(keywords, function(kw) grepl(kw, tr, ignore.case = TRUE)))
})

# ---- Flag SNPs with any off-target trait ----
snps_flagged <- ps_tab %>%
  group_by(snp) %>%
  summarize(
    off_target = any(!matched),
    reason = paste(unique(trait_lc[!matched]), collapse = "; ")
  ) %>%
  ungroup()

snps_to_remove <- snps_flagged %>% filter(off_target) %>% pull(snp)
cat("🧹 SNPs flagged for removal (off-target traits):", length(snps_to_remove), "\n")

# ---- Write filtered output ----
gwas_filtered <- gwas %>% filter(!SNP %in% snps_to_remove)
write_csv(gwas_filtered, output_file)

# ---- Write full filter log ----
filter_log <- tibble(
  SNP = snps,
  removed = snps %in% snps_to_remove,
  reason = ifelse(snps %in% snps_flagged$snp,
                  snps_flagged$reason[match(snps, snps_flagged$snp)],
                  "No PhenoScanner match or all traits matched")
)
write_csv(filter_log, filter_log_file)

cat("✅ SNPs retained:", nrow(gwas_filtered), "\n")
cat("💾 Output written to:", output_file, "\n")
cat("📝 Filter log written to:", filter_log_file, "\n")
cat("====================================================================\n")
