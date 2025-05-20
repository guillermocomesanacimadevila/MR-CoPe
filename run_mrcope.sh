#!/bin/bash

###############################################################################
#                        MR-CoPe Pipeline Runner (Docker Version)             #
###############################################################################
# Usage:
#   chmod +x run_mrcope.sh && ./run_mrcope.sh
#
# Description:
#   This script launches the MR-CoPe pipeline using Docker and Nextflow.
#   It supports both CSV/TSV and VCF (gzipped or not) as inputs.
###############################################################################

set -e  # Exit on any error

echo ""
echo "🧬 Welcome to MR-CoPe: Mendelian Randomisation Pipeline"
echo "--------------------------------------------------------"

# ------------------------ Input Format Prompt ------------------------

read -rp "🧬 Do you have a VCF file or CSV/TSV summary stats? [vcf/csv]: " FORMAT
FORMAT=$(echo "$FORMAT" | tr '[:upper:]' '[:lower:]')

convert_vcf_to_csv() {
  local input_vcf="$1"
  local output_csv="$2"
  echo "🔄 Converting VCF to CSV: $input_vcf → $output_csv"

  python3 - <<EOF
import gzip
import csv

vcf_path = "$input_vcf"
csv_path = "$output_csv"

with gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r') as vcf_in, \
     open(csv_path, 'w', newline='') as csv_out:

    writer = None
    for line in vcf_in:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            headers = line.strip().lstrip("#").split('\t')
            writer = csv.DictWriter(csv_out, fieldnames=headers)
            writer.writeheader()
            continue
        fields = line.strip().split('\t')
        record = dict(zip(headers, fields))
        writer.writerow(record)

print("✅ VCF successfully converted to CSV:", csv_path)
EOF
}

# ------------------------ Input Collection ------------------------

if [[ "$FORMAT" == "vcf" ]]; then
  read -rp "📥 Enter path to Exposure VCF (.vcf or .vcf.gz): " VCF_EXPOSURE
  read -rp "📥 Enter path to Outcome VCF (.vcf or .vcf.gz): " VCF_OUTCOME

  [[ ! -f "$VCF_EXPOSURE" ]] && echo "❌ ERROR: Exposure VCF not found: $VCF_EXPOSURE" && exit 1
  [[ ! -f "$VCF_OUTCOME" ]] && echo "❌ ERROR: Outcome VCF not found: $VCF_OUTCOME" && exit 1

  EXPOSURE_PATH="./tmp_exposure.csv"
  OUTCOME_PATH="./tmp_outcome.csv"

  convert_vcf_to_csv "$VCF_EXPOSURE" "$EXPOSURE_PATH"
  convert_vcf_to_csv "$VCF_OUTCOME" "$OUTCOME_PATH"

elif [[ "$FORMAT" == "csv" ]]; then
  read -rp "📥 Enter path to Exposure GWAS summary statistics (.csv or .tsv): " EXPOSURE_PATH
  read -rp "📥 Enter path to Outcome GWAS summary statistics (.csv or .tsv): " OUTCOME_PATH

  [[ ! -f "$EXPOSURE_PATH" ]] && echo "❌ ERROR: Exposure file not found at: $EXPOSURE_PATH" && exit 1
  [[ ! -f "$OUTCOME_PATH" ]] && echo "❌ ERROR: Outcome file not found at: $OUTCOME_PATH" && exit 1

else
  echo "❌ ERROR: Unknown input format '$FORMAT'. Please enter 'vcf' or 'csv'."
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
