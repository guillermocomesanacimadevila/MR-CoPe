# ===========================================
# MR-CoPe: Dockerfile for Mendelian Randomisation Pipeline
# Base image: R 4.2.3 with Linux tools
# ===========================================

FROM rocker/r-ver:4.2.3

# Prevent interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------
# Install system dependencies
# -----------------------------
RUN apt-get update && apt-get install -y \
    python3.10 python3-pip python-is-python3 \
    build-essential cmake git curl wget \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libpng-dev libtiff5-dev libjpeg-dev \
    && apt-get clean

# -----------------------------
# Make R library path writable
# -----------------------------
RUN mkdir -p /usr/local/lib/R/site-library && chmod -R 777 /usr/local/lib/R/site-library

# -----------------------------
# Install Python dependencies
# -----------------------------
RUN pip install --upgrade pip && pip install \
    pandas numpy matplotlib seaborn scipy

# -----------------------------
# Install R CRAN packages
# -----------------------------
RUN Rscript -e "install.packages(c( \
  'devtools', 'optparse', 'tidyverse', 'qqman', 'ggrepel', 'data.table', \
  'patchwork', 'RColorBrewer', 'MASS', 'Matrix', 'gridExtra', 'lattice', \
  'stringr', 'purrr', 'plyr', 'tidyr' \
), repos='https://cloud.r-project.org')"

# -----------------------------
# ✅ Officially recommended way to install TwoSampleMR
# -----------------------------
RUN Rscript -e "install.packages('TwoSampleMR', repos=c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))" && \
    Rscript -e "stopifnot('TwoSampleMR' %in% rownames(installed.packages()))"

# -----------------------------
# Set default working directory
# -----------------------------
WORKDIR /app
