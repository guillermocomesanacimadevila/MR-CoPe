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

set -e  # Exit on error

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

# --- Ensure CMake is installed --- #
echo ""
echo "🔍 Checking for CMake..."
if ! command -v cmake &> /dev/null; then
  echo "🛠️  CMake not found. Installing via Homebrew..."
  if ! command -v brew &> /dev/null; then
    echo "❌ Homebrew is not installed. Please install it from https://brew.sh and rerun this script."
    exit 1
  fi
  brew install cmake || { echo "❌ Failed to install CMake. Exiting."; exit 1; }
else
  echo "✅ CMake is already installed."
fi

# --- Conda Environment Setup --- #
echo ""
echo "🔧 Checking Conda environment..."

source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | grep -q "mrcope_env"; then
  echo "⚠️  Conda env 'mrcope_env' already exists. Skipping creation."
else
  echo "🛠️  Creating Conda env from environment.yml..."
  conda env create -f environment.yml
fi

echo "📦 Ensuring R packages (e.g. TwoSampleMR) are available..."
conda run -n mrcope_env Rscript post_setup.R

echo ""
echo "✅ Activating environment..."
conda activate mrcope_env

# --- Execute Nextflow Pipeline --- #
echo ""
echo "🚀 Launching MR-CoPe Pipeline..."
nextflow run main.nf --exposure "$EXPOSURE_PATH" --outcome "$OUTCOME_PATH" -resume

echo ""
echo "🎉 MR-CoPe Pipeline completed successfully!"
echo "✨ All outputs are located in: ./results/"
echo "--------------------------------------------------------"
