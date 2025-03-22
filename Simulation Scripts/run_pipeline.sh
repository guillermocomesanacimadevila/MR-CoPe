#!/bin/bash

# ----- RUN PIPELINE ----- #
# ----- -> UNIX
# chmod +x run_pipeline.sh && ./run_pipeline.sh

# Stop pipeline if anything fails
set -e

# --- Create Conda environment --- #
echo "🔧 Creating Conda environment (if not exists)..."
conda create -y -n mrcope_env python=3.10 r-base=4.2 || echo "⚠️ Conda env 'mrcope_env' already exists. Continuing..."

# --- Activate Conda environment --- #
echo "✅ Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mrcope_env

# --- Install Python dependencies --- #
echo "📦 Installing Python libraries..."
pip install --quiet pandas numpy matplotlib seaborn

# --- Check installed Python libraries --- #
echo "🔧 Installed Python libraries:"
pip list

# --- Run Nextflow pipeline --- #
echo "🚀 Running MR pipeline..."
nextflow run 07_simulation_mr_pipeline.nf -resume

# --- Nice! --- #
echo "🎉 Pipeline finished successfully!"
