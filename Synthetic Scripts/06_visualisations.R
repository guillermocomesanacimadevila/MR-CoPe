#!/usr/bin/env Rscript

# === Publication-ready MR Scatter Plots === #

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# --- Load harmonised data --- #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 06_mr_scatter_plots.R <harmonised_data.csv>")
}
harmonised <- read.csv(args[1])

# --- General plot function --- #
create_plot <- function(data, method = "IVW", model, intercept = 0, slope, color = "black", linetype = "solid", filename) {
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

  ggsave(filename, p, width = 8, height = 6, dpi = 300)
  cat("✅ Saved:", filename, "\n")
}

# --- IVW Plot --- #
ivw_model <- lm(beta.outcome ~ beta.exposure, data = harmonised)
create_plot(
  harmonised,
  method = "IVW",
  model = ivw_model,
  intercept = coef(ivw_model)[1],
  slope = coef(ivw_model)[2],
  color = "blue",
  linetype = "solid",
  filename = "MR_Scatter_IVW.png"
)

# --- Egger Plot --- #
egger_model <- lm(beta.outcome ~ beta.exposure, data = harmonised)
create_plot(
  harmonised,
  method = "Egger",
  model = egger_model,
  intercept = coef(egger_model)[1],
  slope = coef(egger_model)[2],
  color = "red",
  linetype = "dashed",
  filename = "MR_Scatter_Egger.png"
)

# --- Weighted Median Plot (zero-intercept) --- #
wme_model <- lm(beta.outcome ~ beta.exposure + 0, data = harmonised)
create_plot(
  harmonised,
  method = "Weighted Median",
  model = wme_model,
  intercept = 0,
  slope = coef(wme_model)[1],
  color = "darkgreen",
  linetype = "dotdash",
  filename = "MR_Scatter_WME.png"
)