#!/usr/bin/env Rscript

# 📦 Auto-install required R packages for MR-CoPe

# --- Define Required CRAN Packages ---
required_pkgs <- c(
  "devtools", "remotes", "data.table", "dplyr", "ggplot2", "ggrepel",
  "tibble", "readr", "stringr", "magrittr", "purrr", "gridExtra",
  "plyr", "tidyr", "grid", "tools", "optparse", "lattice", "patchwork",
  "RColorBrewer", "qqman"
)

# --- Identify Missing Packages ---
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]

# --- Install Missing CRAN Packages Safely ---
if (length(missing_pkgs) > 0) {
  message("📦 Installing missing CRAN packages: ", paste(missing_pkgs, collapse = ", "))
  tryCatch({
    install.packages(missing_pkgs, repos = "http://cran.us.r-project.org")
  }, error = function(e) {
    message("❌ Error installing CRAN packages: ", e$message)
    quit(status = 1)
  })
} else {
  message("✅ All required CRAN packages are already installed.")
}

# --- Load devtools (install if needed) ---
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "http://cran.us.r-project.org")
}
library(devtools)

# --- Install MR-specific GitHub Packages (if needed) ---
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  message("🌐 Installing TwoSampleMR from GitHub...")
  install_github("MRCIEU/TwoSampleMR")
} else {
  message("✅ TwoSampleMR is already installed.")
}

if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  message("🌐 Installing MRPRESSO from GitHub...")
  install_github("rondolab/MR-PRESSO")
} else {
  message("✅ MRPRESSO is already installed.")
}
