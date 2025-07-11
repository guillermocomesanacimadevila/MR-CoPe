#!/usr/bin/env bash

###############################################################################
#                        MR-CoPe Pipeline Runner (Conda Version)              #
###############################################################################
# Usage:
#   chmod +x run_mrcope_conda.sh && ./run_mrcope_conda.sh
###############################################################################

set -e

# --- Color Output (Professional Touch) ---
RED='\033[0;31m'
GRN='\033[0;32m'
YEL='\033[1;33m'
NC='\033[0m'

cleanup() {
  echo -e "\n${RED}🛑 Pipeline interrupted or failed. Exiting...${NC}"
  exit 1
}
trap cleanup INT TERM

echo ""
echo -e "${GRN}🧬 Welcome to MR-CoPe: Mendelian Randomisation Pipeline${NC}"
echo "--------------------------------------------------------"

# --- Detect and enforce running from script dir (avoids path confusion) ---
MYDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$PWD" != "$MYDIR" ]]; then
  echo -e "${YEL}⚠️  Please run this script from its own directory: $MYDIR${NC}"
  exit 1
fi

# --- Version info (from VERSION file or git) ---
if [[ -f VERSION ]]; then
  VERSION=$(cat VERSION)
elif git rev-parse --short HEAD 2>/dev/null; then
  VERSION=$(git rev-parse --short HEAD)
else
  VERSION="(unknown)"
fi
echo -e "${GRN}MR-CoPe Version: $VERSION${NC}"

# --- Conda Check/Env Setup ---
ENV_NAME="mrcope_env"
PYTHON_VERSION="3.10"
R_VERSION="4.2.3"

if ! command -v conda &> /dev/null; then
  echo -e "${RED}❌ Conda is not installed. Please install Miniconda/Anaconda and rerun.${NC}"
  exit 1
fi

ENV_EXISTS=$(conda info --envs | awk '{print $1}' | grep -Fx "$ENV_NAME" || true)
if [[ -n "$ENV_EXISTS" ]]; then
  echo -e "${YEL}⚠️  Conda environment '$ENV_NAME' already exists.${NC}"
  read -rp "❓ Do you want to DELETE and REBUILD the environment? [y/N]: " REBUILD_ENV
  REBUILD_ENV=$(echo "$REBUILD_ENV" | tr '[:upper:]' '[:lower:]')
  if [[ "$REBUILD_ENV" == "y" ]]; then
    echo -e "${YEL}🧹 Removing existing Conda environment '$ENV_NAME'...${NC}"
    conda deactivate || true
    conda env remove -n $ENV_NAME -y
    echo -e "${GRN}📦 Creating new Conda environment '$ENV_NAME'...${NC}"
    conda create -y -n $ENV_NAME python=$PYTHON_VERSION r-base=$R_VERSION r-essentials r-devtools
  else
    echo -e "${GRN}✅ Using existing Conda environment '$ENV_NAME'.${NC}"
  fi
else
  echo -e "${GRN}📦 Creating new Conda environment '$ENV_NAME'...${NC}"
  conda create -y -n $ENV_NAME python=$PYTHON_VERSION r-base=$R_VERSION r-essentials r-devtools
fi

eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

echo -e "${GRN}🐍 Ensuring Python packages are installed...${NC}"
pip install --upgrade pip
pip install pandas numpy matplotlib seaborn scipy plotly jinja2

echo -e "${GRN}📦 Installing/Checking R dependencies...${NC}"
cat > install_mrcope_rpkgs.R <<'EOF'
options(repos = c(
    MRCIEU = "https://mrcieu.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
))
pkgs <- c(
    "devtools", "optparse", "qqman", "ggrepel", "data.table", "patchwork",
    "RColorBrewer", "MASS", "Matrix", "gridExtra", "lattice",
    "stringr", "purrr", "plyr", "tidyr", "remotes", "httr",
    "tidyverse", "TwoSampleMR", "ieugwasr"
)
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies=TRUE)
cat("✅ R dependencies ready!\n")
EOF

Rscript install_mrcope_rpkgs.R

# --- Nextflow Check/Install ---
echo "🔍 Checking for Nextflow..."
NF_CMD=nextflow
if ! command -v nextflow &> /dev/null; then
  if [[ -f "$MYDIR/nextflow" ]]; then
    chmod +x "$MYDIR/nextflow"
    export PATH="$MYDIR:$PATH"
    NF_CMD="$MYDIR/nextflow"
  else
    echo -e "${YEL}⚠️  Nextflow not found, installing locally...${NC}"
    curl -fsSL https://get.nextflow.io | bash
    chmod +x nextflow
    mv nextflow "$MYDIR/"
    export PATH="$MYDIR:$PATH"
    NF_CMD="$MYDIR/nextflow"
  fi
fi
echo -e "${GRN}✅ Nextflow is available: $NF_CMD${NC}"

# --- Create Results Directory If Needed ---
mkdir -p ./results

# ------------------------ Input Type Selection ------------------------
read -rp "🧬 Do you have VCF files or CSV/TSV summary stats? [vcf/csv]: " INPUT_TYPE
INPUT_TYPE=$(echo "$INPUT_TYPE" | tr '[:upper:]' '[:lower:]')

if [[ "$INPUT_TYPE" != "vcf" && "$INPUT_TYPE" != "csv" ]]; then
  echo -e "${RED}❌ ERROR: Unknown input format '$INPUT_TYPE'. Please enter 'vcf' or 'csv'.${NC}"
  exit 1
fi

# ------------------------ File Input Collection ------------------------

if [[ "$INPUT_TYPE" == "vcf" ]]; then
  read -rp "📥 Enter path to Exposure VCF (.vcf or .vcf.gz): " EXPOSURE_VCF
  read -rp "📥 Enter path to Outcome VCF (.vcf or .vcf.gz): " OUTCOME_VCF

  for file in "$EXPOSURE_VCF" "$OUTCOME_VCF"; do
    if [[ ! -f "$file" ]]; then
      echo -e "${RED}❌ ERROR: File not found: $file${NC}"
      exit 1
    fi
  done

  read -rp "🧪 Are the p-values in the VCF encoded as -log10(p)? [y/n]: " LOG10_INPUT
  LOG10_INPUT=$(echo "$LOG10_INPUT" | tr '[:upper:]' '[:lower:]')

  echo "🔄 Parsing VCF: $EXPOSURE_VCF → ./tmp_exposure.csv"
  python3 - <<EOF
import gzip, csv

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
import gzip, csv

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
      echo -e "${RED}❌ ERROR: File not found: $file${NC}"
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
  echo -e "${RED}❌ Invalid clump_kb value. Must be an integer.${NC}"
  exit 1
fi
if ! awk "BEGIN {exit !($CLUMP_R2 > 0 && $CLUMP_R2 <= 1)}"; then
  echo -e "${RED}❌ Invalid clump_r2 value. Must be a number > 0 and ≤ 1 (e.g., 0.01, 1.0).${NC}"
  exit 1
fi

# -------------------- Print Summary Table Before Run -----------------------

echo ""
echo -e "${GRN}🔎 Pipeline launch summary:${NC}"
printf "  %-22s %s\n" "Exposure file:" "$EXPOSURE_PATH"
printf "  %-22s %s\n" "Outcome file:" "$OUTCOME_PATH"
printf "  %-22s %s\n" "LD window (kb):" "$CLUMP_KB"
printf "  %-22s %s\n" "LD r² cutoff:" "$CLUMP_R2"
printf "  %-22s %s\n" "Output dir:" "./results"
printf "  %-22s %s\n" "MR-CoPe version:" "$VERSION"
echo "--------------------------------------------------------"

# --------- Save command to logfile for reproducibility --------------
CMD="$NF_CMD run main.nf \
    --exposure \"$EXPOSURE_PATH\" \
    --outcome \"$OUTCOME_PATH\" \
    --log10_flag \"$LOG10_FLAG\" \
    --clump_kb \"$CLUMP_KB\" \
    --clump_r2 \"$CLUMP_R2\" \
    --output_dir \"./results\" \
    -resume"
echo "$CMD" > ./results/mrcope_command.log

# ------------------------ Run Pipeline ------------------------

echo ""
echo -e "${GRN}🚀 Launching MR-CoPe Pipeline with Conda...${NC}"
eval $CMD

echo ""
echo -e "${GRN}🎉 MR-CoPe Pipeline completed successfully!${NC}"
echo -e "📦 Results are available in: ${YEL}./results/${NC}"
echo -e "🔍 Review the HTML report or output files for your MR results."
echo "--------------------------------------------------------"
