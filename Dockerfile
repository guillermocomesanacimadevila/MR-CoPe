# ===========================================
# MR-CoPe: Dockerfile for Mendelian Randomisation Pipeline
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
    libglpk-dev libzmq3-dev libfontconfig1-dev \
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
# Install base R packages from CRAN
# -----------------------------
RUN Rscript -e "install.packages(c( \
    'devtools', 'optparse', 'qqman', 'ggrepel', 'data.table', \
    'patchwork', 'RColorBrewer', 'MASS', 'Matrix', 'gridExtra', 'lattice', \
    'stringr', 'purrr', 'plyr', 'tidyr', 'remotes' \
), repos='https://cloud.r-project.org')"

# -----------------------------
# Install tidyverse with verification
# -----------------------------
RUN Rscript -e "install.packages('tidyverse', repos='https://cloud.r-project.org')" && \
    Rscript -e "stopifnot('tidyverse' %in% rownames(installed.packages()))"

# -----------------------------
# Install TwoSampleMR + ieugwasr from r-universe with verification
# -----------------------------
RUN Rscript -e "install.packages(c('TwoSampleMR', 'ieugwasr'), \
    repos=c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))" && \
    Rscript -e "stopifnot(all(c('TwoSampleMR', 'ieugwasr') %in% rownames(installed.packages())))"

# -----------------------------
# Set working directory
# -----------------------------
WORKDIR /app
