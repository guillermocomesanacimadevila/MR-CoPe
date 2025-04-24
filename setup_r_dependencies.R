#!/usr/bin/env Rscript

# 📦 Auto-install required R packages for MR-CoPe

safe_cran_install <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) == 0) {
    message("✅ All required packages already installed.")
    return()
  }

  message("📦 Installing missing packages: ", paste(missing, collapse = ", "))
  tryCatch({
    install.packages(missing, repos = "http://cran.us.r-project.org", type = "source")
  }, error = function(e) {
    message("⚠️  Warning installing some packages: ", e$message)
    # Continue, but log what failed
  })
}

# --- CRAN dependencies ---
cran_packages <- c(
  "devtools", "remotes", "data.table", "dplyr", "ggplot2", "ggrepel",
  "tibble", "readr", "stringr", "magrittr", "purrr", "gridExtra",
  "plyr", "tidyr", "grid", "tools", "optparse", "lattice", "patchwork",
  "RColorBrewer", "qqman", "MASS", "Matrix", "rlang"
)

safe_cran_install(cran_packages)

# --- GitHub packages ---
safe_github_install <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("🌐 Installing ", pkg, " from GitHub...")
    tryCatch({
      devtools::install_github(repo)
    }, error = function(e) {
      message("❌ Failed to install ", pkg, ": ", e$message)
    })
  } else {
    message("✅ ", pkg, " is already installed.")
  }
}

safe_github_install("TwoSampleMR", "MRCIEU/TwoSampleMR")
safe_github_install("MRPRESSO", "rondolab/MR-PRESSO")