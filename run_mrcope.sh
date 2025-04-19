#!/bin/bash

###############################################################################
#                        MR-CoPe Pipeline Runner                              #
###############################################################################
# Usage:
#   chmod +x run_mrcope.sh && ./run_mrcope.sh
#
# Description:
#   This script sets up the environment and runs the MR-CoPe pipeline using 
#   user-supplied Exposure and Outcome GWAS summary statistics.
###############################################################################

set -e  # Exit immediately if a command exits with a non-zero status.

echo ""
echo "🧬 Welcome to MR-CoPe: Mendelian Randomisation Pipeline"
echo "--------------------------------------------------------"

# --- Prompt for Exposure file --- #
read -rp "📥 Enter path to Exposure GWAS summary statistics (.csv or .tsv): " EXPOSURE_PATH
if [[ ! -f "$EXPOSURE_PATH" ]]; then
  echo "❌ ERROR: Exposure file not found at: $EXPOSURE_PATH"
  exit 1
fi

# --- Prompt for Outcome file --- #
read -rp "📥 Enter path to Outcome GWAS summary statistics (.csv or .tsv): " OUTCOME_PATH
if [[ ! -f "$OUTCOME_PATH" ]]; then
  echo "❌ ERROR: Outcome file not found at: $OUTCOME_PATH"
  exit 1
fi

echo ""
echo "--------------------------------------"
echo "🔗 Exposure file: $EXPOSURE_PATH"
echo "🔗 Outcome file : $OUTCOME_PATH"
echo "--------------------------------------"

# --- Conda Environment Setup --- #
echo ""
echo "🔧 Checking Conda environment..."

if conda info --envs | grep -q "mrcope_env"; then
  echo "⚠️  Conda env 'mrcope_env' already exists. Skipping creation."
else
  echo "🛠️  Creating Conda env 'mrcope_env'..."
  conda create -y -n mrcope_env python=3.10 r-base=4.2
fi

echo ""
echo "✅ Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mrcope_env

# --- Install Python dependencies --- #
echo ""
echo "📦 Installing Python packages (if needed)..."
pip install --quiet pandas numpy matplotlib seaborn scipy

# --- Execute Nextflow Pipeline --- #
echo ""
echo "🚀 Launching MR-CoPe Pipeline..."
nextflow run main.nf --exposure "$EXPOSURE_PATH" --outcome "$OUTCOME_PATH" -resume

echo ""
echo "🎉 MR-CoPe Pipeline completed successfully!"
echo "✨ All outputs are located in: ./results/"
echo "--------------------------------------------------------"
