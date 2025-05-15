FROM rocker/r-ver:4.2.3

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------
# System dependencies
# -----------------------------
RUN apt-get update && apt-get install -y \
    python3.10 python3-pip python-is-python3 \
    build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
    cmake git wget curl libgit2-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev \
    && apt-get clean

# -----------------------------
# Python dependencies
# -----------------------------
RUN pip install --upgrade pip && pip install \
    pandas numpy matplotlib seaborn scipy

# -----------------------------
# Create and set permissions for R library path
# -----------------------------
ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN mkdir -p ${R_LIBS_USER} && chmod -R 777 ${R_LIBS_USER}

# -----------------------------
# Install core R packages
# -----------------------------
RUN Rscript -e "install.packages(c( \
  'devtools', 'optparse', 'tidyverse', 'qqman', 'ggrepel', 'data.table', \
  'patchwork', 'RColorBrewer', 'MASS', 'Matrix', 'gridExtra', 'lattice', \
  'stringr', 'purrr', 'plyr', 'tidyr' \
), repos = 'https://cloud.r-project.org', lib=Sys.getenv('R_LIBS_USER'))"

# -----------------------------
# Install remotes and TwoSampleMR
# -----------------------------
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org', lib=Sys.getenv('R_LIBS_USER'))" && \
    Rscript -e "remotes::install_github('MRCIEU/TwoSampleMR', upgrade = 'never', lib=Sys.getenv('R_LIBS_USER'))"

# -----------------------------
# Set default working directory
# -----------------------------
WORKDIR /app
