#!/usr/bin/env Rscript

# ================================================================
# MR-CoPe | Publication-Ready MR Scatter Plots
# ================================================================
# Author: Guillermo Comesaña & Christian Pepler
# Date: 2025
#
# Usage:
# Rscript 06_mr_scatter_plots.R <harmonised_data.csv> <output_dir>
#
# Description:
# Generates scatter plots for MR methods:
# - IVW
# - MR Egger
# - Weighted Median
#
# Outputs:
# - MR_Scatter_IVW.png
# - MR_Scatter_Egger.png
# - MR_Scatter_WME.png
# ================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# ---- Handle Arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript 07_mr_scatter_plots.R <harmonised_data.csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]

cat("\n==============================================================\n")
cat("MR-CoPe | Scatter Plots of MR Results\n")
cat("==============================================================\n\n")

if (!file.exists(input_file)) {
  stop(paste("❌ ERROR: Input file not found:", input_file))
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("📥 Loading harmonised data:", input_file, "\n")
harmonised <- read.csv(input_file)
cat("📊 SNPs available for plotting:", nrow(harmonised), "\n\n")


# ---- General Scatter Plot Function ----
create_plot <- function(data, method, intercept, slope, color, linetype, filename) {
  p <- ggplot(data, aes(x = beta.exposure, y = beta.outcome)) +
    geom_point(size = 4, color = "black") +
    geom_text_repel(aes(label = SNP), size = 4, max.overlaps = 100) +
    geom_abline(intercept = intercept, slope = slope, color = color, linetype = linetype, size = 1.2) +
    labs(
      title = paste("MR Scatter Plot -", method),
      x = "SNP Effect on Exposure",
      y = "SNP Effect on Outcome"
    ) +
    theme_minimal(base_size = 16)

  out_path <- file.path(output_dir, filename)
  ggsave(out_path, p, width = 8, height = 6, dpi = 300)
  cat("✅ Saved:", out_path, "\n")
}


# ---- IVW Plot ----
cat("📊 Creating IVW Scatter Plot...\n")
ivw_model <- lm(beta.outcome ~ beta.exposure, data = harmonised)
create_plot(harmonised, "IVW", coef(ivw_model)[1], coef(ivw_model)[2],
            color = "blue", linetype = "solid", filename = "MR_Scatter_IVW.png")


# ---- MR Egger Plot ----
cat("📊 Creating Egger Scatter Plot...\n")
egger_model <- lm(beta.outcome ~ beta.exposure, data = harmonised)
create_plot(harmonised, "Egger", coef(egger_model)[1], coef(egger_model)[2],
            color = "red", linetype = "dashed", filename = "MR_Scatter_Egger.png")


# ---- Weighted Median Plot ----
cat("📊 Creating Weighted Median Scatter Plot...\n")
wme_model <- lm(beta.outcome ~ beta.exposure + 0, data = harmonised)
create_plot(harmonised, "Weighted Median", intercept = 0, coef(wme_model)[1],
            color = "darkgreen", linetype = "dotdash", filename = "MR_Scatter_WME.png")


cat("\n🎉 Scatter plot generation completed!")
cat("\nAll outputs saved in:", output_dir)
cat("\n==============================================================\n\n")