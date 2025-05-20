#!/bin/bash

###############################################################################
#                        MR-CoPe Pipeline Runner (Docker Version)             #
###############################################################################
# Usage:
#   chmod +x run_mrcope.sh && ./run_mrcope.sh
#
# Description:
#   This script launches the MR-CoPe pipeline using Docker and Nextflow.
#   It supports CSV/TSV or VCF (gzipped or not) input formats, with p-value logic.
###############################################################################

set -e  # Exit on any error

echo ""
echo "🧬 Welcome to MR-CoPe: Mendelian Randomisation Pipeline"
echo "--------------------------------------------------------"

# ------------------------ Input Format ------------------------

read -rp "🧬 Do you have VCF files or CSV/TSV summary stats? [vcf/csv]: " INPUT_TYPE
INPUT_TYPE=$(echo "$INPUT_TYPE" | tr '[:upper:]' '[:lower:]')

if [[ "$INPUT_TYPE" != "vcf" && "$INPUT_TYPE" != "csv" ]]; then
  echo "❌ ERROR: Unknown input format '$INPUT_TYPE'. Please enter 'vcf' or 'csv'."
  exit 1
fi

# ------------------------ Handle VCF ------------------------

if [[ "$INPUT_TYPE" == "vcf" ]]; then
  read -rp "📥 Enter path to Exposure VCF (.vcf or .vcf.gz): " EXPOSURE_PATH
  read -rp "📥 Enter path to Outcome VCF (.vcf or .vcf.gz): " OUTCOME_PATH

  if [[ ! -f "$EXPOSURE_PATH" || ! -f "$OUTCOME_PATH" ]]; then
    echo "❌ ERROR: One or both VCF files not found."
    exit 1
  fi

  read -rp "🧪 Are the p-values in the VCF encoded as -log10(p)? [y/n]: " LOG10_INPUT
  IS_LOG10="false"
  if [[ "$LOG10_INPUT" == [Yy] ]]; then
    IS_LOG10="true"
  fi

  echo "🔄 Parsing VCF: $EXPOSURE_PATH → ./tmp_exposure.csv"
  python3 Scripts/parse_vcf.py "$EXPOSURE_PATH" "./tmp_exposure.csv" "$IS_LOG10"

  echo "🔄 Parsing VCF: $OUTCOME_PATH → ./tmp_outcome.csv"
  python3 Scripts/parse_vcf.py "$OUTCOME_PATH" "./tmp_outcome.csv" "$IS_LOG10"

  EXPOSURE_PATH="./tmp_exposure.csv"
  OUTCOME_PATH="./tmp_outcome.csv"

# ------------------------ Handle CSV/TSV ------------------------

else
  read -rp "📥 Enter path to Exposure GWAS summary statistics (.csv or .tsv): " EXPOSURE_PATH
  read -rp "📥 Enter path to Outcome GWAS summary statistics (.csv or .tsv): " OUTCOME_PATH

  if [[ ! -f "$EXPOSURE_PATH" || ! -f "$OUTCOME_PATH" ]]; then
    echo "❌ ERROR: One or both summary stat files not found."
    exit 1
  fi
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
