FROM rocker/r-ver:4.2.3

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install Python, R build tools, and system dependencies
RUN apt-get update && apt-get install -y \
    python3.10 python3-pip python-is-python3 \
    build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
    cmake git wget curl libgit2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    && apt-get clean

# Install Python dependencies
RUN pip install --upgrade pip && pip install \
    pandas numpy matplotlib seaborn scipy

# Install R packages
RUN Rscript -e "install.packages(c(\
    'devtools', 'optparse', 'tidyverse', 'qqman', 'ggrepel', 'data.table', \
    'patchwork', 'RColorBrewer', 'MASS', 'Matrix', 'gridExtra', 'lattice', \
    'stringr', 'purrr', 'plyr', 'tidyr'\
), repos='https://cloud.r-project.org')"

# Install TwoSampleMR from GitHub
RUN Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); \
                remotes::install_github('MRCIEU/TwoSampleMR', upgrade = 'never')"

# Set default workdir
WORKDIR /app

# Default command (override in Nextflow)
CMD ["/bin/bash"]
