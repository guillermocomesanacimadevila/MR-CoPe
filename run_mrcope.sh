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
read -rp "📥 Enter path to Exposure GWAS summary statistics (.csv): " EXPOSURE_PATH
if [[ ! -f "$EXPOSURE_PATH" ]]; then
  echo "❌ ERROR: Exposure file not found at: $EXPOSURE_PATH"
  exit 1
fi

# --- Prompt for Outcome file --- #
read -rp "📥 Enter path to Outcome GWAS summary statistics (.csv): " OUTCOME_PATH
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
echo "🔧 Creating Conda environment (if needed)..."
conda create -y -n mrcope_env python=3.10 r-base=4.2 >/dev/null 2>&1 || echo "⚠️ Conda env 'mrcope_env' already exists."

echo "✅ Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mrcope_env

# --- Install Python dependencies --- #
echo "📦 Installing Python packages..."
pip install --quiet pandas numpy matplotlib seaborn

# --- Track Start Time --- #
START_TIME=$(date +%s)

# --- Execute Nextflow Pipeline --- #
echo ""
echo "🚀 Launching MR-CoPe Pipeline..."
nextflow run main.nf --exposure "$EXPOSURE_PATH" --outcome "$OUTCOME_PATH" -resume

# --- Track End Time --- #
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo ""
echo "🎉 MR-CoPe Pipeline completed successfully!"
echo "⏱️ Total runtime: $((RUNTIME / 60)) minutes and $((RUNTIME % 60)) seconds"
echo "✨ Outputs are located in: ./results/"
echo "--------------------------------------------------------"
