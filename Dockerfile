FROM rocker/r-ver:4.2.3

ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=/usr/local/lib/R/site-library
ENV PATH="${PATH}:/usr/local/bin"

# -----------------------------
# 🧰 System Dependencies
# -----------------------------
RUN apt-get update && apt-get install -y \
    python3.10 python3-pip python-is-python3 \
    build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
    cmake git wget curl libgit2-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# -----------------------------
# 🐍 Python Dependencies
# -----------------------------
RUN pip install --upgrade pip && pip install \
    pandas numpy matplotlib seaborn scipy

# -----------------------------
# 📁 Ensure Writable R Lib Directory
# -----------------------------
RUN mkdir -p /usr/local/lib/R/site-library && \
    chmod -R 777 /usr/local/lib/R/site-library

# -----------------------------
# 📦 R Dependencies (CRAN)
# -----------------------------
RUN Rscript -e "install.packages(c( \
  'devtools', 'optparse', 'tidyverse', 'qqman', 'ggrepel', 'data.table', \
  'patchwork', 'RColorBrewer', 'MASS', 'Matrix', 'gridExtra', 'lattice', \
  'stringr', 'purrr', 'plyr', 'tidyr' \
), repos = 'https://cloud.r-project.org', lib=Sys.getenv('R_LIBS_USER'))"

# -----------------------------
# 🔗 GitHub Dependencies (TwoSampleMR)
# -----------------------------
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org', lib=Sys.getenv('R_LIBS_USER'))" && \
    Rscript -e "remotes::install_github('MRCIEU/TwoSampleMR', upgrade = 'never', lib=Sys.getenv('R_LIBS_USER'))"

# -----------------------------
# 📂 Working Directory
# -----------------------------
WORKDIR /app
