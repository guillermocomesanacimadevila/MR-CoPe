#!/usr/bin/env Rscript

# 📦 Auto-install required R packages for MR-CoPe
required_pkgs <- c(
  "devtools", "remotes", "data.table", "dplyr", "ggplot2", "ggrepel",
  "tibble", "readr", "stringr", "magrittr", "purrr", "gridExtra",
  "plyr", "tidyr", "grid", "tools", "optparse", "lattice", "patchwork",
  "RColorBrewer", "qqman"
)

missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "http://cran.us.r-project.org")
}

# MR-specific tools
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  devtools::install_github("MRCIEU/TwoSampleMR")
}

if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  devtools::install_github("rondolab/MR-PRESSO")
}

