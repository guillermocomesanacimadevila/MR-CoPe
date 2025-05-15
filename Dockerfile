FROM rocker/r-ver:4.2.3

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------
# Install system dependencies
# -----------------------------
RUN apt-get update && apt-get install -y \
    python3.10 python3-pip python-is-python3 \
    build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
    cmake git wget curl libgit2-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev \
    && apt-get clean

# -----------------------------
# Install Python packages
# -----------------------------
RUN pip install --upgrade pip && pip install \
    pandas numpy matplotlib seaborn scipy

# -----------------------------
# Set R library path and permissions
# -----------------------------
ENV R_LIBS_SITE=/usr/local/lib/R/site-library
RUN mkdir -p ${R_LIBS_SITE} && chmod -R 777 ${R_LIBS_SITE}

# -----------------------------
# Install R CRAN packages
# -----------------------------
RUN Rscript -e "install.packages(c( \
  'devtools', 'optparse', 'tidyverse', 'qqman', 'ggrepel', 'data.table', \
  'patchwork', 'RColorBrewer', 'MASS', 'Matrix', 'gridExtra', 'lattice', \
  'stringr', 'purrr', 'plyr', 'tidyr' \
), repos = 'https://cloud.r-project.org')"

# -----------------------------
# Install TwoSampleMR from GitHub
# -----------------------------
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')" && \
    Rscript -e "remotes::install_github('MRCIEU/TwoSampleMR', upgrade = 'never')"

# -----------------------------
# Set working directory
# -----------------------------
WORKDIR /app
