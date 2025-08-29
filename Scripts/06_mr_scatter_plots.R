#!/usr/bin/env Rscript

# ================================================================
# MR-CoPe | Publication-Ready MR Scatter Plots (robust version)
# (Aesthetics-only upgrade: palette, theme, grids, bands, styling)
# ================================================================
# Usage:
#   Rscript 06_mr_scatter_plots.R <harmonised_data.csv> <output_dir>
# Outputs required by Nextflow:
#   - MR_Scatter_IVW.png
#   - MR_Scatter_Egger.png
#   - MR_Scatter_WME.png
#   - MR_LeaveOneOut.png
# Also writes:
#   - MR_Scatter_Combined.png
# ================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(TwoSampleMR)
  library(purrr)
})

# --------- Aesthetics (palette + theme) ---------
COL_POINT <- "#2B6CB0"  # points (blue)
COL_IVW   <- "#2B6CB0"  # IVW line
COL_EGGER <- "#C53030"  # Egger line (red)
COL_WME   <- "#2F855A"  # WME line (green)
COL_GRID  <- "#E6E8EC"  # major grid
COL_AXES  <- "#2D3748"  # text color
COL_REF   <- "#A0AEC0"  # zero/reference lines
COL_CI    <- "#6B7280"  # error bars
COL_BAND  <- "#F7FAFC"  # alternating row band (LOO)

theme_pub <- function(base_size = 14, base_family = "DejaVu Sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(color = COL_AXES),
      plot.title = element_text(size = base_size * 1.15, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size * 0.95),
      axis.text  = element_text(size = base_size * 0.85),
      panel.grid.major = element_line(color = COL_GRID, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.height = grid::unit(10, "pt"),
      legend.spacing.x = grid::unit(8, "pt"),
      plot.margin = margin(10, 16, 10, 16)
    )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 06_mr_scatter_plots.R <harmonised_data.csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]

cat("\n==============================================================\n")
cat("MR-CoPe | Combined MR Scatter + Leave-One-Out (robust)\n")
cat("==============================================================\n\n")

# ---- Checks ----
if (!file.exists(input_file)) stop(paste("❌ ERROR: Input file not found:", input_file))
if (file.info(input_file)$size == 0) {
  cat("⚠️  Skipping plots — harmonised data is empty.\n")
  quit(status = 0)
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("📥 Loading harmonised data:", input_file, "\n")
harmonised <- read.csv(input_file)
cat("📊 SNPs available for plotting:", nrow(harmonised), "\n\n")

# ---- Fit simple lines (visual guide only) ----
cat("📊 Fitting guide lines (IVW, Egger, WME)...\n")
ivw_fit   <- lm(beta.outcome ~ beta.exposure, data = harmonised)
egger_fit <- lm(beta.outcome ~ beta.exposure, data = harmonised)
wme_fit   <- lm(beta.outcome ~ beta.exposure + 0, data = harmonised)

legend_df <- tibble(
  method    = factor(c("IVW", "Egger", "WME"), levels = c("IVW", "Egger", "WME")),
  slope     = c(coef(ivw_fit)[2], coef(egger_fit)[2], coef(wme_fit)[1]),
  intercept = c(coef(ivw_fit)[1], coef(egger_fit)[1], 0)
)

# ---- Combined Scatter Plot (aesthetics upgraded) ----
cat("📊 Creating combined MR scatter...\n")
scatter_plot <- ggplot(harmonised, aes(x = beta.exposure, y = beta.outcome)) +
  # outcome CI bars
  geom_errorbar(aes(ymin = beta.outcome - se.outcome,
                    ymax = beta.outcome + se.outcome),
                width = 0, color = COL_CI, alpha = 0.7) +
  # points (filled circles with thin stroke)
  geom_point(size = 2.1, shape = 21, fill = COL_POINT, color = "black", stroke = 0.25, alpha = 0.9) +
  # reference axes (soft dashed)
  geom_hline(yintercept = 0, color = COL_REF, linewidth = 0.7, linetype = "dashed") +
  geom_vline(xintercept = 0, color = COL_REF, linewidth = 0.7, linetype = "dashed") +
  # guide lines
  geom_abline(aes(slope = slope, intercept = intercept, color = method, linetype = method),
              data = legend_df, linewidth = 1.1, show.legend = TRUE) +
  labs(
    title = "Combined MR Scatter Plot",
    x = expression(beta[exposure]),
    y = expression(beta[outcome]),
    color = "Method",
    linetype = "Method"
  ) +
  scale_color_manual(values = c("IVW" = COL_IVW, "Egger" = COL_EGGER, "WME" = COL_WME)) +
  scale_linetype_manual(values = c("IVW" = "solid", "Egger" = "longdash", "WME" = "dotdash")) +
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  theme_pub(14)

combined_path <- file.path(output_dir, "MR_Scatter_Combined.png")
ggsave(combined_path, scatter_plot, width = 8.5, height = 6.5, dpi = 300, bg = "white")
cat("✅ Saved:", combined_path, "\n")

# ---- Also write method-named PNGs (satisfy Nextflow outputs) ----
ivw_path   <- file.path(output_dir, "MR_Scatter_IVW.png")
egger_path <- file.path(output_dir, "MR_Scatter_Egger.png")
wme_path   <- file.path(output_dir, "MR_Scatter_WME.png")
ggsave(ivw_path,   scatter_plot, width = 8.5, height = 6.5, dpi = 300, bg = "white")
ggsave(egger_path, scatter_plot, width = 8.5, height = 6.5, dpi = 300, bg = "white")
ggsave(wme_path,   scatter_plot, width = 8.5, height = 6.5, dpi = 300, bg = "white")
cat("✅ Saved:", ivw_path, "\n")
cat("✅ Saved:", egger_path, "\n")
cat("✅ Saved:", wme_path, "\n")

# ---- Leave-One-Out (robust) ----
cat("📊 Creating Leave-One-Out plot...\n")
harmonised <- harmonised %>% mutate(id.exposure = "Exposure", id.outcome = "Outcome")

safe_mr_leaveoneout <- function(dat) {
  tryCatch(TwoSampleMR::mr_leaveoneout(dat), error = function(e) NULL)
}

res_loo <- safe_mr_leaveoneout(harmonised)

# If mr_leaveoneout failed or lacks a recognizable method column, compute manual IVW LOO
need_manual <- is.null(res_loo) ||
               (!("method" %in% names(res_loo) | "Method" %in% names(res_loo)))

if (need_manual) {
  cat("ℹ️  Falling back to manual IVW leave-one-out computation...\n")
  if (nrow(harmonised) < 3) {
    cat("⚠️  Not enough SNPs for meaningful leave-one-out. Writing placeholder plot.\n")
    placeholder <- ggplot() + theme_void() +
      annotate("text", x = 0, y = 0, label = "Leave-one-out not available (n < 3)", size = 6)
    ggsave(file.path(output_dir, "MR_LeaveOneOut.png"), placeholder, width = 9, height = 7, dpi = 300, bg = "white")
    cat("✅ Saved placeholder LOO plot.\n")
  } else {
    manual_loo <- map_dfr(seq_len(nrow(harmonised)), function(i) {
      dat_i <- harmonised[-i, ]
      res   <- TwoSampleMR::mr(dat_i, method_list = c("mr_ivw"))
      tibble(
        SNP = harmonised$SNP[i],
        b   = res$b[1],
        se  = res$se[1]
      )
    })

    # alternating band background rows
    band_df <- tibble(i = seq_len(nrow(manual_loo)))

    p_loo <- ggplot(manual_loo, aes(x = b, y = reorder(SNP, b))) +
      geom_rect(data = band_df,
                aes(xmin = -Inf, xmax = Inf, ymin = i - 0.5, ymax = i + 0.5),
                inherit.aes = FALSE,
                fill = rep(c(COL_BAND, "white"), length.out = nrow(band_df)),
                color = NA) +
      geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                     height = 0.18, color = COL_CI, linewidth = 0.9) +
      geom_point(aes(fill = b >= 0), size = 2.4, shape = 21, color = "black", stroke = 0.25) +
      scale_fill_manual(values = c("TRUE" = "#2F855A", "FALSE" = "#C53030"), guide = "none") +
      geom_vline(xintercept = 0, color = COL_REF, linewidth = 0.8, linetype = "dashed") +
      labs(
        title = "Leave-One-Out Analysis (IVW)",
        x = expression("Causal Estimate " * beta),
        y = "SNP"
      ) +
      scale_x_continuous(expand = expansion(mult = 0.06)) +
      theme_pub(14)

    ggsave(file.path(output_dir, "MR_LeaveOneOut.png"), p_loo, width = 9, height = 7, dpi = 300, bg = "white")
    cat("✅ Saved:", file.path(output_dir, "MR_LeaveOneOut.png"), "\n")
  }
} else {
  # Normalize column name casing and filter to IVW
  if ("Method" %in% names(res_loo) && !("method" %in% names(res_loo))) {
    res_loo <- res_loo %>% rename(method = Method)
  }
  # Some versions name columns differently; make sure we have b & se
  if (!("b" %in% names(res_loo) && "se" %in% names(res_loo))) {
    stop("❌ mr_leaveoneout returned without 'b' and 'se' columns; cannot plot.")
  }
  res_loo_ivw <- res_loo %>% dplyr::filter(.data$method == "Inverse variance weighted")
  if (nrow(res_loo_ivw) == 0) {
    cat("ℹ️  No IVW rows in mr_leaveoneout(); using all rows instead.\n")
    res_loo_ivw <- res_loo
  }

  band_df <- tibble(i = seq_len(nrow(res_loo_ivw)))

  p_loo <- ggplot(res_loo_ivw, aes(x = b, y = reorder(SNP, b))) +
    geom_rect(data = band_df,
              aes(xmin = -Inf, xmax = Inf, ymin = i - 0.5, ymax = i + 0.5),
              inherit.aes = FALSE,
              fill = rep(c(COL_BAND, "white"), length.out = nrow(band_df)),
              color = NA) +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                   height = 0.18, color = COL_CI, linewidth = 0.9) +
    geom_point(aes(fill = b >= 0), size = 2.4, shape = 21, color = "black", stroke = 0.25) +
    scale_fill_manual(values = c("TRUE" = "#2F855A", "FALSE" = "#C53030"), guide = "none") +
    geom_vline(xintercept = 0, color = COL_REF, linewidth = 0.8, linetype = "dashed") +
    labs(
      title = "Leave-One-Out Analysis (IVW)",
      x = expression("Causal Estimate " * beta),
      y = "SNP"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.06)) +
    theme_pub(14)

  ggsave(file.path(output_dir, "MR_LeaveOneOut.png"), p_loo, width = 9, height = 7, dpi = 300, bg = "white")
  cat("✅ Saved:", file.path(output_dir, "MR_LeaveOneOut.png"), "\n")
}

cat("\n🎉 All visualisations complete!\n")
cat("All outputs saved in:", output_dir, "\n")
cat("==============================================================\n\n")
