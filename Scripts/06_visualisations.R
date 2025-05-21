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
# - MR_LeaveOneOut.png
# ================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(TwoSampleMR)
})

# ---- Handle Arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript 06_mr_scatter_plots.R <harmonised_data.csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]

cat("\n==============================================================\n")
cat("MR-CoPe | Combined MR Scatter + Leave-One-Out\n")
cat("==============================================================\n\n")

# ---- Empty File Check ----
if (!file.exists(input_file)) {
  stop(paste("❌ ERROR: Input file not found:", input_file))
}

if (file.info(input_file)$size == 0) {
  cat("⚠️  Skipping scatter plot — harmonised data is empty.\n")
  quit(status = 0)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("📥 Loading harmonised data:", input_file, "\n")
harmonised <- read.csv(input_file)
cat("📊 SNPs available for plotting:", nrow(harmonised), "\n\n")

# ---- Fit MR Models ----
cat("📊 Fitting MR models (IVW, Egger, WME)...\n")
ivw_fit <- lm(beta.outcome ~ beta.exposure, data = harmonised)
egger_fit <- lm(beta.outcome ~ beta.exposure, data = harmonised)
wme_fit <- lm(beta.outcome ~ beta.exposure + 0, data = harmonised)

# ---- Legend Data ----
legend_df <- tibble(
  x = c(NA, NA, NA),
  y = c(NA, NA, NA),
  method = factor(c("IVW", "Egger", "WME"), levels = c("IVW", "Egger", "WME")),
  slope = c(coef(ivw_fit)[2], coef(egger_fit)[2], coef(wme_fit)[1]),
  intercept = c(coef(ivw_fit)[1], coef(egger_fit)[1], 0)
)

# ---- Combined Scatter Plot ----
cat("📊 Creating combined MR scatter plot...\n")

scatter_plot <- ggplot(harmonised, aes(x = beta.exposure, y = beta.outcome)) +
  geom_errorbar(aes(ymin = beta.outcome - se.outcome,
                    ymax = beta.outcome + se.outcome),
                width = 0, color = "gray70", alpha = 0.7) +
  geom_point(size = 1.7, shape = 21, fill = "black", color = "black", stroke = 0.2, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_abline(aes(slope = slope, intercept = intercept, color = method, linetype = method),
              data = legend_df, size = 1.1, show.legend = TRUE) +
  labs(
    title = "Combined MR Scatter Plot",
    x = expression(beta[exposure]),
    y = expression(beta[outcome]),
    color = "Method",
    linetype = "Method"
  ) +
  scale_color_manual(
    values = c("IVW" = "#0072B2", "Egger" = "#D55E00", "WME" = "#009E73")
  ) +
  scale_linetype_manual(
    values = c("IVW" = "solid", "Egger" = "dashed", "WME" = "dotdash")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# ---- Save Combined Plot ----
combined_path <- file.path(output_dir, "MR_Scatter_Combined.png")
ggsave(combined_path, scatter_plot, width = 8.5, height = 6.5, dpi = 300)
cat("✅ Saved:", combined_path, "\n")

# ---- Leave-One-Out Plot ----
cat("📊 Creating Leave-One-Out plot...\n")

harmonised <- harmonised %>%
  mutate(id.exposure = "Exposure", id.outcome = "Outcome")

res_loo <- mr_leaveoneout(harmonised) %>%
  filter(method == "Inverse variance weighted")

p_loo <- ggplot(res_loo, aes(x = b, y = reorder(SNP, b))) +
  geom_point(size = 2.5, shape = 21, fill = "black", color = "black", stroke = 0.3) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                 height = 0.2, color = "gray40") +
  labs(
    title = "Leave-One-Out Analysis (IVW)",
    x = expression("Causal Estimate " * beta),
    y = "SNP"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  )

loo_path <- file.path(output_dir, "MR_LeaveOneOut.png")
ggsave(loo_path, p_loo, width = 9, height = 7, dpi = 300)
cat("✅ Saved:", loo_path, "\n")

cat("\n🎉 All visualisations complete!\n")
cat("All outputs saved in:", output_dir)
cat("==============================================================\n\n"
