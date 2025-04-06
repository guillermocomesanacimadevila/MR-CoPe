# ================================================
# Dockerfile for Real GWAS MR Pipeline
# Includes R, Python, and all required packages
# ================================================

# --- Base image with R preinstalled --- #
FROM rocker/r-ver:4.3.1

# --- Install system dependencies --- #
RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxt-dev \
    git \
    curl \
    unzip && \
    apt-get clean

# --- Install required Python packages --- #
RUN pip3 install --no-cache-dir pandas matplotlib seaborn

# --- Install required R packages --- #
RUN R -e "install.packages(c('TwoSampleMR', 'tidyverse', 'optparse'), repos='https://cloud.r-project.org/')"

# --- Set working directory inside container --- #
WORKDIR /pipeline

# --- Copy all pipeline files into the image (optional for build context) --- #
# You can use this if building from a Git repo or local folder
# COPY . /pipeline

# --- Default to interactive shell --- #
CMD [ "bash" ]
