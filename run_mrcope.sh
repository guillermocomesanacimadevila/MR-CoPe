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

# ------------------------ Input Type Selection ------------------------

read -rp "🧬 Do you have VCF files or CSV/TSV summary stats? [vcf/csv]: " INPUT_TYPE
INPUT_TYPE=$(echo "$INPUT_TYPE" | tr '[:upper:]' '[:lower:]')

if [[ "$INPUT_TYPE" != "vcf" && "$INPUT_TYPE" != "csv" ]]; then
  echo "❌ ERROR: Unknown input format '$INPUT_TYPE'. Please enter 'vcf' or 'csv'."
  exit 1
fi

# ------------------------ File Input Collection ------------------------

if [[ "$INPUT_TYPE" == "vcf" ]]; then
  read -rp "📥 Enter path to Exposure VCF (.vcf or .vcf.gz): " EXPOSURE_VCF
  read -rp "📥 Enter path to Outcome VCF (.vcf or .vcf.gz): " OUTCOME_VCF

  for file in "$EXPOSURE_VCF" "$OUTCOME_VCF"; do
    if [[ ! -f "$file" ]]; then
      echo "❌ ERROR: File not found: $file"
      exit 1
    fi
  done

  read -rp "🧪 Are the p-values in the VCF encoded as -log10(p)? [y/n]: " LOG10_INPUT
  LOG10_INPUT=$(echo "$LOG10_INPUT" | tr '[:upper:]' '[:lower:]')

  echo "🔄 Parsing VCF: $EXPOSURE_VCF → ./tmp_exposure.csv"
  python3 - <<EOF
import gzip, csv, sys

vcf_path = "$EXPOSURE_VCF"
csv_path = "./tmp_exposure.csv"
is_log10 = "$LOG10_INPUT" == "y"

def parse_sample_field(field):
    es, se, lp, af, rsid = field.split(":")
    try:
        pval = 10 ** (-float(lp)) if is_log10 else float(lp)
    except:
        pval = None
    return {
        "BETA": float(es),
        "SE": float(se),
        "PVALUE": pval,
        "EAF": float(af),
        "SNP": rsid
    }

with (gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r')) as vcf_in, \
     open(csv_path, 'w', newline='') as csv_out:

    writer = None
    for line in vcf_in:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.strip().lstrip("#").split("\t")
            format_col = header.index("FORMAT")
            sample_col = format_col + 1
            writer = csv.DictWriter(csv_out, fieldnames=[
                "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "PVALUE", "EAF"
            ])
            writer.writeheader()
            continue

        fields = line.strip().split("\t")
        if len(fields) <= sample_col:
            continue
        chrom, pos, snp_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        try:
            stats = parse_sample_field(fields[sample_col])
            writer.writerow({
                "SNP": stats["SNP"],
                "CHR": chrom,
                "BP": pos,
                "A1": alt,
                "A2": ref,
                "BETA": stats["BETA"],
                "SE": stats["SE"],
                "PVALUE": stats["PVALUE"],
                "EAF": stats["EAF"]
            })
        except Exception as e:
            print(f"⚠️ Skipping line at position {pos}: {e}")
EOF

  echo "🔄 Parsing VCF: $OUTCOME_VCF → ./tmp_outcome.csv"
  python3 - <<EOF
import gzip, csv, sys

vcf_path = "$OUTCOME_VCF"
csv_path = "./tmp_outcome.csv"
is_log10 = "$LOG10_INPUT" == "y"

def parse_sample_field(field):
    es, se, lp, af, rsid = field.split(":")
    try:
        pval = 10 ** (-float(lp)) if is_log10 else float(lp)
    except:
        pval = None
    return {
        "BETA": float(es),
        "SE": float(se),
        "PVALUE": pval,
        "EAF": float(af),
        "SNP": rsid
    }

with (gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r')) as vcf_in, \
     open(csv_path, 'w', newline='') as csv_out:

    writer = None
    for line in vcf_in:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.strip().lstrip("#").split("\t")
            format_col = header.index("FORMAT")
            sample_col = format_col + 1
            writer = csv.DictWriter(csv_out, fieldnames=[
                "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "PVALUE", "EAF"
            ])
            writer.writeheader()
            continue

        fields = line.strip().split("\t")
        if len(fields) <= sample_col:
            continue
        chrom, pos, snp_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        try:
            stats = parse_sample_field(fields[sample_col])
            writer.writerow({
                "SNP": stats["SNP"],
                "CHR": chrom,
                "BP": pos,
                "A1": alt,
                "A2": ref,
                "BETA": stats["BETA"],
                "SE": stats["SE"],
                "PVALUE": stats["PVALUE"],
                "EAF": stats["EAF"]
            })
        except Exception as e:
            print(f"⚠️ Skipping line at position {pos}: {e}")
EOF

  EXPOSURE_PATH="./tmp_exposure.csv"
  OUTCOME_PATH="./tmp_outcome.csv"
  LOG10_FLAG="$LOG10_INPUT"

else
  read -rp "📥 Enter path to Exposure GWAS summary statistics (.csv or .tsv): " EXPOSURE_PATH
  read -rp "📥 Enter path to Outcome GWAS summary statistics (.csv or .tsv): " OUTCOME_PATH

  for file in "$EXPOSURE_PATH" "$OUTCOME_PATH"; do
    if [[ ! -f "$file" ]]; then
      echo "❌ ERROR: File not found: $file"
      exit 1
    fi
  done

  read -rp "🧪 Are the p-values in your input files already -log10(p)? [y/n]: " LOG10_FLAG
  LOG10_FLAG=$(echo "$LOG10_FLAG" | tr '[:upper:]' '[:lower:]')
fi

# ------------------------ LD Clumping Settings ------------------------

echo ""
echo "📊 LD Clumping Settings:"
echo "   - These affect SNP pruning based on correlation (LD)."
echo "   - More relaxed values will keep more SNPs."
read -rp "📏 LD clumping window size in kb (default: 10000): " CLUMP_KB
read -rp "🔗 LD clumping r² threshold (default: 0.001): " CLUMP_R2

CLUMP_KB="${CLUMP_KB:-10000}"
CLUMP_R2="${CLUMP_R2:-0.001}"

if ! [[ "$CLUMP_KB" =~ ^[0-9]+$ ]]; then
  echo "❌ Invalid clump_kb value. Must be an integer."
  exit 1
fi
if ! awk "BEGIN {exit !($CLUMP_R2 > 0 && $CLUMP_R2 <= 1)}"; then
  echo "❌ Invalid clump_r2 value. Must be a number > 0 and ≤ 1 (e.g., 0.01, 1.0)."
  exit 1
fi

echo ""
echo "📁 Exposure file : $EXPOSURE_PATH"
echo "📁 Outcome file  : $OUTCOME_PATH"
echo "📏 LD window (kb): $CLUMP_KB"
echo "🔗 LD r² cutoff  : $CLUMP_R2"
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
  read -rp "🔧 Build image from Dockerfile now? [y/N]: " BUILD_CHOICE
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
    --log10_flag "$LOG10_FLAG" \
    --clump_kb "$CLUMP_KB" \
    --clump_r2 "$CLUMP_R2" \
    --output_dir "./results" \
    -resume

echo ""
echo "🎉 MR-CoPe Pipeline completed successfully!"
echo "📦 Results are available in: ./results/"
echo "--------------------------------------------------------"
