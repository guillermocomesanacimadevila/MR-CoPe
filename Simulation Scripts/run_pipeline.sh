#!/bin/bash

# ----- RUN PIPELINE ----- #
# ----- -> UNIX
# chmod +x run_pipeline.sh && ./run_pipeline.sh

# Stop pipeline if anything fails
set -e

# --- Create Conda environment ---
echo "🔧 Creating Conda environment (if not exists)..."
conda create -y -n mrcope_env python=3.10 r-base=4.2 || echo "⚠️ Conda env 'mrcope_env' already exists. Continuing..."

# --- Activate Conda environment ---
echo "✅ Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mrcope_env

# --- Install Python dependencies ---
echo "📦 Installing Python packages..."
pip install --quiet pandas numpy matplotlib seaborn

# --- Install essential R dependencies ---
echo "📦 Installing core R packages..."
Rscript -e 'install.packages(c("data.table", "optparse", "ggplot2", "tidyverse", "devtools"), repos="https://cloud.r-project.org")'

# Your main R script (04_mr_analyses.R) handles the rest!

# --- Run the pipeline ---
echo "🚀 Running pipeline..."
nextflow run simulation_mr_pipeline.nf -resume

echo "🎉 Pipeline finished successfully!"
