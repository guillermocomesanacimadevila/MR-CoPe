#!/bin/bash

# ----- RUN PIPELINE ----- #
# ----- -> UNIX
# chmod +x run_pipeline.sh && ./run_pipeline.sh

# Stop pipeline if anything fails
set -e

echo "🔧 Creating Conda environment..."
conda create -y -n mrcope_env python=3.10 r-base=4.2 || echo "⚠️ Environment already exists. Skipping..."

echo "✅ Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mrcope_env

echo "📦 Installing Python packages..."
pip install pandas numpy matplotlib seaborn

echo "📦 Installing R packages..."
Rscript -e 'install.packages(c("tidyverse", "ggplot2", "data.table", "optparse", "devtools"), repos="https://cloud.r-project.org")'

echo "🚀 Running pipeline..."
nextflow run simulation_mr_pipeline.nf -resume

echo "🎉 Pipeline finished!"