#!/bin/bash

###############################################################################
#                        MR-CoPe Pipeline Runner (Docker Version)             #
###############################################################################
# Usage:
#   chmod +x run_mrcope.sh && ./run_mrcope.sh
#
# Description:
#   This script launches the MR-CoPe pipeline using Docker and Nextflow.
#   It prompts the user for GWAS input files and executes the full workflow.
###############################################################################

set -e  # Exit on any error

echo ""
echo "🧬 Welcome to MR-CoPe: Mendelian Randomisation Pipeline"
echo "--------------------------------------------------------"

# ------------------------ Input Collection ------------------------

read -rp "📥 Enter path to Exposure GWAS summary statistics (.csv or .tsv): " EXPOSURE_PATH
if [[ ! -f "$EXPOSURE_PATH" ]]; then
  echo "❌ ERROR: Exposure file not found at: $EXPOSURE_PATH"
  exit 1
fi

read -rp "📥 Enter path to Outcome GWAS summary statistics (.csv or .tsv): " OUTCOME_PATH
if [[ ! -f "$OUTCOME_PATH" ]]; then
  echo "❌ ERROR: Outcome file not found at: $OUTCOME_PATH"
  exit 1
fi

echo ""
echo "📁 Exposure: $EXPOSURE_PATH"
echo "📁 Outcome : $OUTCOME_PATH"
echo "--------------------------------------------------------"

# ------------------------ Dependency Checks ------------------------

echo "🔍 Checking for Docker..."
if ! command -v docker &> /dev/null; then
  echo "❌ Docker is not installed. Please install it from https://www.docker.com and rerun this script."
  exit 1
fi
echo "✅ Docker is installed."

echo "🔍 Checking for Nextflow..."
if ! command -v nextflow &> /dev/null; then
  echo "❌ Nextflow is not installed. Please install it from https://www.nextflow.io and rerun this script."
  exit 1
fi
echo "✅ Nextflow is installed."

# ------------------------ Docker Image Check ------------------------

IMAGE_NAME="mrcope:latest"
if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
  echo ""
  echo "🐳 Docker image '$IMAGE_NAME' not found locally."
  read -rp "🔧 Would you like to build it now from Dockerfile? [y/N]: " BUILD_CHOICE
  if [[ "$BUILD_CHOICE" =~ ^[Yy]$ ]]; then
    docker build -t "$IMAGE_NAME" .
    echo "✅ Docker image '$IMAGE_NAME' built successfully."
  else
    echo "❌ Cannot continue without Docker image '$IMAGE_NAME'. Exiting."
    exit 1
  fi
fi

# ------------------------ Run Pipeline ------------------------

echo ""
echo "🚀 Launching MR-CoPe Pipeline with Docker..."
nextflow run main.nf -with-docker "$IMAGE_NAME" \
    --exposure "$EXPOSURE_PATH" \
    --outcome "$OUTCOME_PATH" \
    --output_dir "./results" \
    -resume

# ------------------------ Completion Message ------------------------

echo ""
echo "🎉 MR-CoPe Pipeline completed successfully!"
echo "📦 Results are available in: ./results/"
echo "--------------------------------------------------------"
