#!/usr/bin/env bash

###############################################################################
#                        MR-CoPe Pipeline Runner (AWS Auto)                   #
###############################################################################
# Usage:
#   chmod +x run_mrcope.sh && ./run_mrcope.sh
#
# Description:
#   This script launches the MR-CoPe pipeline using Docker and Nextflow.
#   On AWS/Ubuntu: auto-installs Docker, Java, Nextflow, Python3, pip, R.
#   On Mac: as before, user action may be needed for Docker Desktop.
###############################################################################

set -e  # Exit on any error

# --- Color Output (Professional Touch) ---
RED='\033[0;31m'
GRN='\033[0;32m'
YEL='\033[1;33m'
NC='\033[0m' # No Color

cleanup() {
  echo -e "\n${RED}🛑 Pipeline interrupted or failed. Exiting...${NC}"
  echo -e "${YEL}💡 You can resume where you left off by running this script again (Nextflow will pick up from the last successful step).${NC}"
  exit 1
}
trap cleanup INT TERM

echo ""
echo -e "${GRN}🧬 Welcome to MR-CoPe: Mendelian Randomisation Pipeline${NC}"
echo "--------------------------------------------------------"

# --- Enforce running from script dir ---
MYDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$PWD" != "$MYDIR" ]]; then
  echo -e "${YEL}⚠️  Please run this script from its own directory: $MYDIR${NC}"
  exit 1
fi

# --- Version info ---
if [[ -f VERSION ]]; then
  VERSION=$(cat VERSION)
elif git rev-parse --short HEAD 2>/dev/null; then
  VERSION=$(git rev-parse --short HEAD)
else
  VERSION="(unknown)"
fi
echo -e "${GRN}MR-CoPe Version: $VERSION${NC}"

###############################################################################
#                  TOOLCHAIN CHECK AND AUTO-INSTALL SECTION                   #
###############################################################################

install_msg() { echo -e "${YEL}🔧 $1${NC}"; }

if ! command -v docker &> /dev/null; then
  install_msg "Docker not found. Installing Docker (this may take a minute)..."
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    curl -fsSL https://get.docker.com | sudo bash || { echo -e "${RED}❌ Failed to install Docker. Install manually.${NC}"; exit 1; }
    sudo usermod -aG docker $USER
    echo -e "${YEL}⚠️  Please LOG OUT and log back in for Docker permissions to update."
    echo -e "   Then rerun this script.${NC}"
    exit 0
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo -e "${RED}❌ Docker Desktop must be installed manually on macOS (https://www.docker.com/products/docker-desktop/).${NC}"
    exit 1
  fi
else
  echo -e "${GRN}✅ Docker is installed.${NC}"
fi

# -------- Python3 and pip check/install ---------
if ! command -v python3 &> /dev/null; then
  install_msg "Python3 not found. Installing Python3..."
  sudo apt-get update && sudo apt-get install -y python3 || { echo -e "${RED}❌ Failed to install Python3. Install manually.${NC}"; exit 1; }
fi
if ! command -v pip3 &> /dev/null; then
  install_msg "pip3 not found. Installing pip3..."
  sudo apt-get install -y python3-pip || { echo -e "${RED}❌ Failed to install pip3. Install manually.${NC}"; exit 1; }
fi

# -------- Java check/install (needed for Nextflow) ---------
if ! command -v java &> /dev/null; then
  install_msg "Java not found. Installing OpenJDK 11..."
  sudo apt-get update && sudo apt-get install -y openjdk-11-jre-headless || { echo -e "${RED}❌ Failed to install Java. Install manually.${NC}"; exit 1; }
fi

# -------- Nextflow check/install ---------
if ! command -v nextflow &> /dev/null; then
  install_msg "Nextflow not found. Installing Nextflow..."
  curl -fsSL https://get.nextflow.io | bash || { echo -e "${RED}❌ Failed to install Nextflow. Install manually.${NC}"; exit 1; }
  sudo mv nextflow /usr/local/bin/
  sudo chmod +x /usr/local/bin/nextflow
fi

# -------- R check/install (needed for MR-CoPe) ---------
if ! command -v Rscript &> /dev/null; then
  install_msg "R not found. Installing R (this may take a minute)..."
  sudo apt-get update && sudo apt-get install -y r-base || { echo -e "${RED}❌ Failed to install R. Install manually.${NC}"; exit 1; }
fi

# --- Essential Python packages (for data prep) ---
pip3 install --user numpy pandas

# --- R packages (for MR-CoPe, e.g. data.table) ---
echo -e "${YEL}🔧 Checking/Installing required R packages...${NC}"
Rscript -e 'pkgs<-c("data.table","dplyr","ggplot2","optparse"); new<-pkgs[!pkgs %in% rownames(installed.packages())]; if(length(new)) install.packages(new, repos="https://cloud.r-project.org")'

###############################################################################

# --- Results Dir, Logging & Debugging ---
mkdir -p ./results

FILTER_SUMMARY="./results/filter_summary.csv"
PARAM_LOG="./results/parameter_log.txt"
echo "step,file,kept,removed,message" > "$FILTER_SUMMARY"
echo "MR-CoPe parameter log" > "$PARAM_LOG"
echo "date: $(date)" >> "$PARAM_LOG"
echo "MR-CoPe version: $VERSION" >> "$PARAM_LOG"

append_filter_summary() {
  echo "$1,$2,$3,$4,$5" >> "$FILTER_SUMMARY"
}
append_param_log() {
  echo "$1: $2" >> "$PARAM_LOG"
}

# ------------------------ Input Type Selection ------------------------
read -rp "🧬 Do you have VCF files or CSV/TSV summary stats? [vcf/csv]: " INPUT_TYPE
INPUT_TYPE=$(echo "$INPUT_TYPE" | tr '[:upper:]' '[:lower:]')
append_param_log "input_type" "$INPUT_TYPE"

if [[ "$INPUT_TYPE" != "vcf" && "$INPUT_TYPE" != "csv" ]]; then
  echo -e "${RED}❌ ERROR: Unknown input format '$INPUT_TYPE'. Please enter 'vcf' or 'csv'.${NC}"
  exit 1
fi

# ------------------------ File Input Collection + Validation ------------------------
if [[ "$INPUT_TYPE" == "vcf" ]]; then
  read -rp "📥 Enter path to Exposure VCF (.vcf or .vcf.gz): " EXPOSURE_VCF
  read -rp "📥 Enter path to Outcome VCF (.vcf or .vcf.gz): " OUTCOME_VCF

  for file in "$EXPOSURE_VCF" "$OUTCOME_VCF"; do
    if [[ ! -f "$file" ]]; then
      echo -e "${RED}❌ ERROR: File not found: $file${NC}"
      exit 1
    elif [[ ! -s "$file" ]]; then
      echo -e "${RED}❌ ERROR: File is empty: $file${NC}"
      exit 1
    fi
  done
  append_param_log "exposure_vcf" "$EXPOSURE_VCF"
  append_param_log "outcome_vcf" "$OUTCOME_VCF"

  read -rp "🧪 Are the p-values in the VCF encoded as -log10(p)? [y/n]: " LOG10_INPUT
  LOG10_INPUT=$(echo "$LOG10_INPUT" | tr '[:upper:]' '[:lower:]')
  append_param_log "log10_flag" "$LOG10_INPUT"

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
        "BETA": float(es), "SE": float(se), "PVALUE": pval,
        "EAF": float(af), "SNP": rsid
    }
with (gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r')) as vcf_in, open(csv_path, 'w', newline='') as csv_out:
    writer = None; rows = 0
    for line in vcf_in:
        if line.startswith("##"): continue
        if line.startswith("#CHROM"):
            header = line.strip().lstrip("#").split("\t")
            format_col = header.index("FORMAT")
            sample_col = format_col + 1
            writer = csv.DictWriter(csv_out, fieldnames=["SNP","CHR","BP","A1","A2","BETA","SE","PVALUE","EAF"])
            writer.writeheader(); continue
        fields = line.strip().split("\t")
        if len(fields) <= sample_col: continue
        chrom, pos, snp_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        try:
            stats = parse_sample_field(fields[sample_col])
            writer.writerow({
                "SNP": stats["SNP"], "CHR": chrom, "BP": pos, "A1": alt, "A2": ref,
                "BETA": stats["BETA"], "SE": stats["SE"], "PVALUE": stats["PVALUE"], "EAF": stats["EAF"]
            }); rows += 1
        except Exception as e:
            print(f"⚠️ Skipping line at position {pos}: {e}", file=sys.stderr)
    print(rows, file=sys.stderr)
EOF
  NROWS=$(awk 'END{print NR-1}' ./tmp_exposure.csv)
  append_filter_summary "vcf_parse_exposure" "./tmp_exposure.csv" "$NROWS" "NA" "VCF→CSV conversion"

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
        "BETA": float(es), "SE": float(se), "PVALUE": pval,
        "EAF": float(af), "SNP": rsid
    }
with (gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r')) as vcf_in, open(csv_path, 'w', newline='') as csv_out:
    writer = None; rows = 0
    for line in vcf_in:
        if line.startswith("##"): continue
        if line.startswith("#CHROM"):
            header = line.strip().lstrip("#").split("\t")
            format_col = header.index("FORMAT")
            sample_col = format_col + 1
            writer = csv.DictWriter(csv_out, fieldnames=["SNP","CHR","BP","A1","A2","BETA","SE","PVALUE","EAF"])
            writer.writeheader(); continue
        fields = line.strip().split("\t")
        if len(fields) <= sample_col: continue
        chrom, pos, snp_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        try:
            stats = parse_sample_field(fields[sample_col])
            writer.writerow({
                "SNP": stats["SNP"], "CHR": chrom, "BP": pos, "A1": alt, "A2": ref,
                "BETA": stats["BETA"], "SE": stats["SE"], "PVALUE": stats["PVALUE"], "EAF": stats["EAF"]
            }); rows += 1
        except Exception as e:
            print(f"⚠️ Skipping line at position {pos}: {e}", file=sys.stderr)
    print(rows, file=sys.stderr)
EOF
  NROWS=$(awk 'END{print NR-1}' ./tmp_outcome.csv)
  append_filter_summary "vcf_parse_outcome" "./tmp_outcome.csv" "$NROWS" "NA" "VCF→CSV conversion"

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
    elif [[ ! -s "$file" ]]; then
      echo -e "${RED}❌ ERROR: File is empty: $file${NC}"
      exit 1
    fi
  done
  append_param_log "exposure_csv" "$EXPOSURE_PATH"
  append_param_log "outcome_csv" "$OUTCOME_PATH"

  read -rp "🧪 Are the p-values in your input files already -log10(p)? [y/n]: " LOG10_FLAG
  LOG10_FLAG=$(echo "$LOG10_FLAG" | tr '[:upper:]' '[:lower:]')
  append_param_log "log10_flag" "$LOG10_FLAG"
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

# ------------------------ Edge Case Warnings ------------------------
if [[ "$CLUMP_KB" -gt 100000 ]]; then
  echo -e "${YEL}⚠️ WARNING: You entered a very large window size (${CLUMP_KB} kb). This may remove too many SNPs.${NC}"
fi
if [[ "$CLUMP_KB" -lt 50 ]]; then
  echo -e "${YEL}⚠️ WARNING: You entered a very small window size (${CLUMP_KB} kb). This may retain highly correlated SNPs.${NC}"
fi
if (( $(echo "$CLUMP_R2 == 1" | bc -l) )); then
  echo -e "${YEL}⚠️ WARNING: r² = 1 means no LD pruning will be performed.${NC}"
fi
if (( $(echo "$CLUMP_R2 < 0.0001" | bc -l) )); then
  echo -e "${YEL}⚠️ WARNING: Very low r² will keep almost no SNPs.${NC}"
fi

if ! [[ "$CLUMP_KB" =~ ^[0-9]+$ ]]; then
  echo -e "${RED}❌ Invalid clump_kb value. Must be an integer.${NC}"
  exit 1
fi
if ! awk "BEGIN {exit !($CLUMP_R2 > 0 && $CLUMP_R2 <= 1)}"; then
  echo -e "${RED}❌ Invalid clump_r2 value. Must be a number > 0 and ≤ 1 (e.g., 0.01, 1.0).${NC}"
  exit 1
fi

# ------------------------ LD Population Selection ------------------------
echo ""
echo "🌍 LD Reference Population:"
echo "   - EUR (European, default)"
echo "   - AFR (African)"
echo "   - EAS (East Asian)"
echo "   - SAS (South Asian)"
echo "   - AMR (Admixed American/Latino)"
read -rp "LD reference population for clumping? [EUR/AFR/EAS/SAS/AMR] (default: EUR): " LD_POP
LD_POP="${LD_POP:-EUR}"
append_param_log "ld_pop" "$LD_POP"

# ------------------------ Clumping Sensitivity Panel ------------------------
echo ""
read -rp "🧪 Run clumping sensitivity panel (try several kb/r² combinations in parallel)? [y/N]: " SENS_PANEL
SENS_PANEL=$(echo "$SENS_PANEL" | tr '[:upper:]' '[:lower:]')

if [[ "$SENS_PANEL" == "y" ]]; then
  echo -e "${YEL}Running clumping with multiple parameter combinations...${NC}"
  PANEL_KB=(500 1000 5000 10000)
  PANEL_R2=(0.01 0.05 0.1 0.5)
  PANEL_LOG="./results/clumping_sensitivity_panel.csv"
  echo "window_kb,r2,num_SNPs_retained" > "$PANEL_LOG"
  for KB in "${PANEL_KB[@]}"; do
    for R2 in "${PANEL_R2[@]}"; do
      OUT="./results/tmp_clump_${KB}_${R2}.csv"
      docker run --rm -v "$PWD":/data mrcope:latest \
        Rscript /app/03_linkage_disequillibrium.R "$EXPOSURE_PATH" "$OUT" "$KB" "$R2" "$LD_POP"
      COUNT=$(awk 'END{print NR-1}' "$OUT")
      echo "$KB,$R2,$COUNT" >> "$PANEL_LOG"
      rm -f "$OUT"
    done
  done
  echo -e "${GRN}Panel complete. Results in $PANEL_LOG${NC}"
fi

# ------------------------ Confounder Keyword Input ------------------------

echo ""
echo -e "${GRN}🧠 Confounder Filtering with PhenoScanner${NC}"
echo "--------------------------------------------------------"
echo "👉 You can enter one or more keywords to RETAIN SNPs associated only with your exposure of interest."
echo "   - Examples: smoking, alcohol, education, socioeconomic, physical activity"
echo "   - Format:  Type keywords separated by commas (NO quotes or spaces)"
echo "   - Case-insensitive; partial matches allowed"
echo -e "   - ${YEL}⚠️  Leaving this blank disables confounder filtering entirely!${NC}"
echo ""

read -rp "📌 Enter confounder keyword(s) to retain (comma-separated): " TRAIT_KEYWORDS

if [[ -z "$TRAIT_KEYWORDS" ]]; then
  echo -e "${YEL}⚠️  No keywords entered. PhenoScanner confounder filtering will be skipped.${NC}"
  TRAIT_KEYWORDS="."  # Placeholder to indicate 'skip'
else
  echo -e "${GRN}🎯 You entered keyword(s):${NC} ${TRAIT_KEYWORDS}"
fi

# -------------------- Print Summary Table Before Run -----------------------
echo ""
echo -e "${GRN}🔎 Pipeline launch summary:${NC}"
printf "  %-22s %s\n" "Exposure file:" "$EXPOSURE_PATH"
printf "  %-22s %s\n" "Outcome file:" "$OUTCOME_PATH"
printf "  %-22s %s\n" "LD window (kb):" "$CLUMP_KB"
printf "  %-22s %s\n" "LD r² cutoff:" "$CLUMP_R2"
printf "  %-22s %s\n" "LD population:" "$LD_POP"
printf "  %-22s %s\n" "Output dir:" "./results"
printf "  %-22s %s\n" "MR-CoPe version:" "$VERSION"
echo "--------------------------------------------------------"

# -------------------- Docker Image Check --------------------------
IMAGE_NAME="mrcope:latest"
if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
  echo ""
  echo -e "${YEL}🐳 Docker image '$IMAGE_NAME' not found locally.${NC}"
  read -rp "🔧 Build image from Dockerfile now? [y/N]: " BUILD_CHOICE
  if [[ "$BUILD_CHOICE" =~ ^[Yy]$ ]]; then
    docker build -t "$IMAGE_NAME" .
    echo -e "${GRN}✅ Docker image '$IMAGE_NAME' built successfully.${NC}"
  else
    echo -e "${RED}❌ Cannot continue without Docker image '$IMAGE_NAME'. Exiting.${NC}"
    exit 1
  fi
fi

# --------- Save command to logfile for reproducibility --------------
NF_CMD=nextflow
CMD="$NF_CMD run main.nf -with-docker \"$IMAGE_NAME\" \
    --exposure \"$EXPOSURE_PATH\" \
    --outcome \"$OUTCOME_PATH\" \
    --log10_flag \"$LOG10_FLAG\" \
    --clump_kb \"$CLUMP_KB\" \
    --clump_r2 \"$CLUMP_R2\" \
    --ld_pop \"$LD_POP\" \
    --trait_keyword \"$TRAIT_KEYWORDS\" \
    --output_dir \"./results\" \
    -resume"
echo "$CMD" > ./results/mrcope_command.log
append_param_log "run_command" "$CMD"

echo ""
echo -e "${YEL}💡 All key pipeline parameters saved in ${PARAM_LOG}."
echo -e "${YEL}💡 Filter summaries will be written to ${FILTER_SUMMARY} after each major step."
echo -e "${YEL}💡 If the pipeline is interrupted, just re-run with '-resume'.${NC}"

# ------------------------ Run Pipeline ------------------------

echo ""
echo -e "${GRN}🚀 Launching MR-CoPe Pipeline with Docker...${NC}"
if ! eval $CMD; then
  echo -e "${RED}❌ Pipeline failed!${NC}"
  echo -e "${YEL}💡 You can resume from where you left off using:${NC}"
  echo -e "${YEL}   $CMD${NC}"
  exit 1
fi

echo ""
echo -e "${GRN}🎉 MR-CoPe Pipeline completed successfully!${NC}"
echo -e "📦 Results are available in: ${YEL}./results/${NC}"
echo -e "🔍 Review the HTML report or output files for your MR results."
echo -e "${YEL}💡 See filter_summary.csv and parameter_log.txt for reproducibility & filtering stats.${NC}"
echo "--------------------------------------------------------"
