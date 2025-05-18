# ===========================================
# MR-CoPe: Dockerfile for Mendelian Randomisation Pipeline
# Base image: R 4.2.3 with Linux tools
# ===========================================

FROM rocker/r-ver:4.2.3

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=/usr/local/lib/R/site-library

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
# Ensure R library path is writable
# -----------------------------
RUN mkdir -p ${R_LIBS_USER} && chmod -R 777 ${R_LIBS_USER}

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
# Install TwoSampleMR + ieugwasr (required for ld_clump)
# -----------------------------
RUN Rscript -e "install.packages(c('TwoSampleMR', 'ieugwasr'), \
  repos=c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))" && \
  Rscript -e "stopifnot(all(c('TwoSampleMR', 'ieugwasr') %in% rownames(installed.packages())))"

# -----------------------------
# Set default working directory
# -----------------------------
WORKDIR /app
