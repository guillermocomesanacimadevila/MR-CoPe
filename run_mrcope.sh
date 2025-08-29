#!/usr/bin/env bash

###############################################################################
#                        MR-CoPe Pipeline Runner (Docker Version)             #
###############################################################################
# Usage:
#   chmod +x run_mrcope.sh && ./run_mrcope.sh
#
# Description:
#   This script launches the MR-CoPe pipeline using Docker and Nextflow.
#   It checks if Docker is running and starts it if needed.
#   It also checks/loads or installs Nextflow before running the main workflow.
#   + fzf integration for interactive file selection (optional).
###############################################################################

set -e  # Exit on any error

# --- Color Output (Professional Touch) ---
RED='\033[0;31m'
GRN='\033[0;32m'
YEL='\033[1;33m'
NC='\033[0m' # No Color

# Debug toggle (set MRCOPE_DEBUG=1 in your env to see extra info)
: "${MRCOPE_DEBUG:=0}"

cleanup() {
  echo -e "\n${RED}🛑 Pipeline interrupted or failed. Exiting...${NC}"
  echo -e "${YEL}💡 You can resume where you left off by running this script again (Nextflow will pick up from the last successful step).${NC}"
  exit 1
}
trap cleanup INT TERM

echo ""
echo -e "${GRN}🧬 Welcome to MR-CoPe: Mendelian Randomisation Pipeline${NC}"
echo "--------------------------------------------------------"

# --- Docker runtime helpers (daemon up / start if needed) ---
docker_running() {
  docker info >/dev/null 2>&1
}

start_docker_if_needed() {
  if docker_running; then
    echo -e "${GRN}✅ Docker daemon is running.${NC}"
    return 0
  fi
  echo -e "${YEL}🟡 Docker is installed but not running. Attempting to start it...${NC}"
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    if command -v systemctl >/dev/null 2>&1; then
      sudo systemctl start docker || true
    fi
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    # Start Docker Desktop (no error if already running)
    open -ga Docker || open -a Docker || true
  fi
  # Wait up to ~2 minutes for the daemon to come up
  for _ in {1..60}; do
    if docker_running; then
      echo -e "${GRN}✅ Docker daemon started.${NC}"
      return 0
    fi
    sleep 2
  done
  echo -e "${RED}❌ Docker is not running and could not be started automatically. Please start Docker and re-run.${NC}"
  exit 1
}

# --- macOS: ensure Docker CLI on PATH if Docker Desktop is installed ---
ensure_docker_cli_macos() {
  local candidates=(
    "/Applications/Docker.app/Contents/Resources/bin"
    "$HOME/Applications/Docker.app/Contents/Resources/bin"
    "/opt/homebrew/bin"
    "/usr/local/bin"
  )
  for d in "${candidates[@]}"; do
    if [[ -x "$d/docker" ]]; then
      case ":$PATH:" in
        *":$d:"*) : ;;
        *) export PATH="$d:$PATH" ;;
      esac
      return 0
    fi
  done
  return 1
}

# --- Enforce running from script dir ---
MYDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$PWD" != "$MYDIR" ]]; then
  echo -e "${YEL}⚠️  Please run this script from its own directory: $MYDIR${NC}"
  exit 1
fi

# --- Version info ---
if [[ -f VERSION ]]; then
  VERSION=$(cat VERSION)
elif git rev-parse --short HEAD >/dev/null 2>&1; then
  VERSION=$(git rev-parse --short HEAD)
else
  VERSION="(unknown)"
fi
echo -e "${GRN}MR-CoPe Version: $VERSION${NC}"

# --- Data root detection: prefer ./Data in repo, then /Data, else current dir ---
if [ -n "${DATA_DIR:-}" ]; then
  : # respect user override
elif [ -d "./Data" ]; then
  DATA_DIR="./Data"
elif [ -d "/Data" ]; then
  DATA_DIR="/Data"
else
  DATA_DIR="."
fi
echo -e "${YEL}🔎 Searching for inputs under: ${DATA_DIR}${NC}"

# Force vertical list + border consistently (single definition)
export FZF_DEFAULT_OPTS="--layout=default --border ${FZF_DEFAULT_OPTS}"

###############################################################################
#                  TOOLCHAIN CHECK AND AUTO-INSTALL SECTION                   #
###############################################################################

install_msg() { echo -e "${YEL}🔧 $1${NC}"; }

# -------- Docker check/install ---------
# On macOS, try to discover the Docker CLI even if it's not on PATH
if [[ "$OSTYPE" == "darwin"* ]] && ! command -v docker >/dev/null 2>&1; then
  if ensure_docker_cli_macos && command -v docker >/dev/null 2>&1; then
    echo -e "${GRN}✅ Found Docker CLI via Docker.app.${NC}"
  fi
fi

if ! command -v docker &> /dev/null; then
    install_msg "Docker not found. Attempting to install Docker..."
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        curl -fsSL https://get.docker.com | sudo bash || { echo -e "${RED}❌ Failed to install Docker. Install manually.${NC}"; exit 1; }
        sudo usermod -aG docker "$USER" || true
        echo -e "${YEL}⚠️  Please log out and log back in for Docker permissions to update.${NC}"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        if command -v brew &> /dev/null; then
            brew install --cask docker || { echo -e "${RED}❌ Install Docker Desktop from https://www.docker.com/get-started/ (Mac).${NC}"; exit 1; }
            echo "📝 Please start Docker Desktop manually from Applications."
        else
            echo -e "${RED}❌ Homebrew not found. Install Homebrew or Docker Desktop manually.${NC}"
            exit 1
        fi
    else
        echo -e "${RED}❌ Docker install: Unrecognized OS, please install manually.${NC}"
        exit 1
    fi
else
    echo -e "${GRN}✅ Docker is installed.${NC}"
fi

# Ensure Docker daemon is up (Desktop/Service running) before proceeding
start_docker_if_needed

# -------- Java check/install (needed for Nextflow) ---------
if ! command -v java &> /dev/null; then
    install_msg "Java not found. Attempting to install OpenJDK 11..."
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        sudo apt-get update && sudo apt-get install -y openjdk-11-jre-headless || { echo -e "${RED}❌ Failed to install Java. Install manually.${NC}"; exit 1; }
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        if command -v brew &> /dev/null; then
            brew install openjdk@11 || { echo -e "${RED}❌ Failed to install Java (brew). Install manually.${NC}"; exit 1; }
            sudo ln -sfn /usr/local/opt/openjdk@11/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk-11.jdk || true
            export PATH="/usr/local/opt/openjdk@11/bin:$PATH"
        else
            echo -e "${RED}❌ Homebrew not found. Install Homebrew or Java manually.${NC}"
            exit 1
        fi
    else
        echo -e "${RED}❌ Java install: Unrecognized OS, please install manually.${NC}"
        exit 1
    fi
else
    echo -e "${GRN}✅ Java is installed.${NC}"
fi

# -------- Nextflow check/install ---------
if ! command -v nextflow &> /dev/null; then
    install_msg "Nextflow not found. Installing Nextflow..."
    curl -fsSL https://get.nextflow.io | bash || { echo -e "${RED}❌ Failed to install Nextflow. Install manually.${NC}"; exit 1; }
    sudo mv nextflow /usr/local/bin/
    chmod +x /usr/local/bin/nextflow
    export PATH="/usr/local/bin:$PATH"
else
    echo -e "${GRN}✅ Nextflow is installed.${NC}"
fi

# -------- fzf check / optional install ---------
FZF_AVAILABLE="false"
if command -v fzf >/dev/null 2>&1; then
  FZF_AVAILABLE="true"
  echo -e "${GRN}✅ fzf is installed. Interactive file picker enabled.${NC}"
else
  echo -e "${YEL}ℹ️  fzf (interactive picker) not found.${NC}"
  read -rp "🧭 Install fzf for interactive file selection? [y/N]: " INSTALL_FZF
  INSTALL_FZF=$(echo "$INSTALL_FZF" | tr '[:upper:]' '[:lower:]')
  if [[ "$INSTALL_FZF" == "y" ]]; then
    install_msg "Installing fzf..."
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
      if command -v apt-get >/dev/null 2>&1; then
        sudo apt-get update && sudo apt-get install -y fzf || { echo -e "${RED}❌ fzf install failed. Proceeding without it.${NC}"; }
      elif command -v dnf >/dev/null 2>&1; then
        sudo dnf install -y fzf || { echo -e "${RED}❌ fzf install failed. Proceeding without it.${NC}"; }
      elif command -v yum >/dev/null 2>&1; then
        sudo yum install -y fzf || { echo -e "${RED}❌ fzf install failed. Proceeding without it.${NC}"; }
      else
        echo -e "${YEL}⚠️  Unknown Linux package manager. Consider installing fzf manually.${NC}"
      fi
    elif [[ "$OSTYPE" == "darwin"* ]]; then
      if command -v brew >/dev/null 2>&1; then
        brew install fzf || true
        # Install key-bindings/completion (safe if it fails)
        "$(brew --prefix)"/opt/fzf/install --key-bindings --completion --no-bashrc --no-fish --no-zsh || true
      else
        echo -e "${YEL}⚠️  Homebrew not found. Install fzf manually (https://github.com/junegunn/fzf).${NC}"
      fi
    else
      echo -e "${YEL}⚠️  Unrecognized OS. Install fzf manually if desired.${NC}"
    fi

    if command -v fzf >/dev/null 2>&1; then
      FZF_AVAILABLE="true"
      echo -e "${GRN}✅ fzf installed. Interactive picker enabled.${NC}"
    else
      echo -e "${YEL}➡️  Continuing without fzf. We'll use standard prompts.${NC}"
    fi
  else
    echo -e "${YEL}➡️  Skipping fzf install. We'll use standard prompts.${NC}"
  fi
fi

# -------- findutils (gfind) check/install for macOS ---------
if [[ "$OSTYPE" == "darwin"* ]]; then
  if ! command -v gfind >/dev/null 2>&1; then
    echo -e "${YEL}🔧 GNU findutils (gfind) not found. Installing via Homebrew...${NC}"
    if command -v brew &> /dev/null; then
      brew install findutils || {
        echo -e "${RED}❌ Failed to install findutils. Please install manually: brew install findutils${NC}"
        exit 1
      }
      echo -e "${GRN}✅ gfind installed successfully.${NC}"
    else
      echo -e "${RED}❌ Homebrew not found. Please install Homebrew (https://brew.sh) and rerun this script.${NC}"
      exit 1
    fi
  else
    echo -e "${GRN}✅ gfind is installed.${NC}"
  fi
fi

###############################################################################
#                         Helpers (fzf-based pickers)                         #
###############################################################################

# choose_file <prompt> <extensions_csv> <search_dir>
#   examples: "vcf,vcf.gz" or "csv,tsv"
choose_file() {
  local prompt="$1"
  local exts_csv="$2"
  local search_dir="$3"

  # If fzf isn't available, do plain prompt
  if ! command -v fzf >/dev/null 2>&1; then
    read -rp "$prompt " path
    printf '%s\n' "$path"
    return
  fi

  # Use gfind on macOS if available; otherwise find
  local FIND_BIN="find"
  if command -v gfind >/dev/null 2>&1; then
    FIND_BIN="gfind"
  fi

  # Build case-insensitive -iname patterns safely
  local IFS=,
  local exts=()
  read -r -a exts <<< "$exts_csv"

  local args=( "$search_dir" -type f \( )
  local first=1
  local e
  for e in "${exts[@]}"; do
    if [ $first -eq 0 ]; then args+=( -o ); fi
    args+=( -iname "*.${e}" )
    first=0
  done
  args+=( \) )

  # Get the candidate list FIRST; if empty, skip fzf
  local files
  files="$("$FIND_BIN" "${args[@]}" 2>/dev/null | LC_ALL=C sort)"

  if [ -z "$files" ]; then
    echo "No matches under '$search_dir' for extensions: ${exts[*]}"
    read -rp "$prompt " path
    printf '%s\n' "$path"
    return
  fi

  # Always try fzf when available
  local selection
  selection="$(printf '%s\n' "$files" | fzf --prompt="$prompt " --height=20 --border --layout=default || true)"

  if [ -n "$selection" ]; then
    printf '%s\n' "$selection"
  else
    # user pressed ESC/Enter with no selection — fall back to manual input
    read -rp "$prompt " path
    printf '%s\n' "$path"
  fi
}

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
  # Use fzf picker (or fallback) for VCF paths
  EXPOSURE_VCF="$(choose_file "📥 Select Exposure VCF (.vcf or .vcf.gz):" "vcf,vcf.gz" "$DATA_DIR")"
  OUTCOME_VCF="$(choose_file  "📥 Select Outcome VCF (.vcf or .vcf.gz):"  "vcf,vcf.gz" "$DATA_DIR")"

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

  # --- Export BEFORE calling Python so subprocesses can read them ---
  export EXPOSURE_VCF OUTCOME_VCF LOG10_INPUT

  echo "🔄 Parsing VCF: $EXPOSURE_VCF → ./tmp_exposure.csv"
  python3 - <<'EOF'
import gzip, csv, sys, os
vcf_path = os.environ.get("EXPOSURE_VCF")
csv_path = "./tmp_exposure.csv"
is_log10 = (os.environ.get("LOG10_INPUT","n").lower() == "y")
def parse_sample_field(field):
    es, se, lp, af, rsid = field.split(":")
    try:
        pval = 10 ** (-float(lp)) if is_log10 else float(lp)
    except Exception:
        pval = None
    return {"BETA": float(es), "SE": float(se), "PVALUE": pval, "EAF": float(af), "SNP": rsid}
open_vcf = gzip.open if vcf_path.endswith('.gz') else open
with (open_vcf(vcf_path, 'rt')) as vcf_in, open(csv_path, 'w', newline='') as csv_out:
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
            writer.writerow({"SNP": stats["SNP"], "CHR": chrom, "BP": pos, "A1": alt, "A2": ref,
                             "BETA": stats["BETA"], "SE": stats["SE"], "PVALUE": stats["PVALUE"], "EAF": stats["EAF"]})
            rows += 1
        except Exception as e:
            print(f"⚠️ Skipping line at position {pos}: {e}", file=sys.stderr)
    print(rows, file=sys.stderr)
EOF
  NROWS=$(awk 'END{print NR-1}' ./tmp_exposure.csv)
  append_filter_summary "vcf_parse_exposure" "./tmp_exposure.csv" "$NROWS" "NA" "VCF→CSV conversion"

  echo "🔄 Parsing VCF: $OUTCOME_VCF → ./tmp_outcome.csv"
  python3 - <<'EOF'
import gzip, csv, sys, os
vcf_path = os.environ.get("OUTCOME_VCF")
csv_path = "./tmp_outcome.csv"
is_log10 = (os.environ.get("LOG10_INPUT","n").lower() == "y")
def parse_sample_field(field):
    es, se, lp, af, rsid = field.split(":")
    try:
        pval = 10 ** (-float(lp)) if is_log10 else float(lp)
    except Exception:
        pval = None
    return {"BETA": float(es), "SE": float(se), "PVALUE": pval, "EAF": float(af), "SNP": rsid}
open_vcf = gzip.open if vcf_path.endswith('.gz') else open
with (open_vcf(vcf_path, 'rt')) as vcf_in, open(csv_path, 'w', newline='') as csv_out:
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
            writer.writerow({"SNP": stats["SNP"], "CHR": chrom, "BP": pos, "A1": alt, "A2": ref,
                             "BETA": stats["BETA"], "SE": stats["SE"], "PVALUE": stats["PVALUE"], "EAF": stats["EAF"]})
            rows += 1
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
  # Use fzf picker (or fallback) for CSV/TSV paths
  EXPOSURE_PATH="$(choose_file "📥 Select Exposure GWAS (.csv or .tsv):" "csv,tsv" "$DATA_DIR")"
  OUTCOME_PATH="$(choose_file  "📥 Select Outcome GWAS (.csv or .tsv):"  "csv,tsv" "$DATA_DIR")"

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
echo -e "${GRN}🧩 LD Pruning Step${NC}"
echo "--------------------------------------------------------"
echo "You can SKIP LD clumping if you wish (not generally recommended for most MR)."
read -rp "⏭️  Skip LD pruning? [y/N]: " SKIP_LD
SKIP_LD=$(echo "$SKIP_LD" | tr '[:upper:]' '[:lower:]')

if [[ "$SKIP_LD" == "y" ]]; then
  DO_LD="false"
  echo -e "${YEL}⚠️  LD pruning will be skipped. All SNPs will be retained for MR.${NC}"
else
  DO_LD="true"
fi
append_param_log "skip_ld" "$SKIP_LD"

if [[ "$SKIP_LD" != "y" ]]; then
  echo ""
  echo "📊 LD Clumping Settings:"
  echo "   - These affect SNP pruning based on correlation (LD)."
  echo "   - More relaxed values will keep more SNPs."
  read -rp "📏 LD clumping window size in kb (default: 10000): " CLUMP_KB
  read -rp "🔗 LD clumping r² threshold (default: 0.001): " CLUMP_R2

  CLUMP_KB="${CLUMP_KB:-10000}"
  CLUMP_R2="${CLUMP_R2:-0.001}"

  # Edge case warnings (portable, avoid bc)
  awk -v kb="$CLUMP_KB" 'BEGIN{ if (kb+0>100000) print "WARN_LARGE"; else if (kb+0<50) print "WARN_SMALL"; }' | while read -r w; do
    case "$w" in
      WARN_LARGE) echo -e "${YEL}⚠️ WARNING: You entered a very large window size (${CLUMP_KB} kb). This may remove too many SNPs.${NC}";;
      WARN_SMALL) echo -e "${YEL}⚠️ WARNING: You entered a very small window size (${CLUMP_KB} kb). This may retain highly correlated SNPs.${NC}";;
    esac
  done

  awk -v r2="$CLUMP_R2" 'BEGIN{
    if (r2 == 1) print "WARN_EQ1";
    if (r2+0 < 0.0001) print "WARN_VLOW";
    if (!(r2 > 0 && r2 <= 1)) print "ERR_BAD";
  }' | while read -r w; do
    case "$w" in
      WARN_EQ1)  echo -e "${YEL}⚠️ WARNING: r² = 1 means no LD pruning will be performed.${NC}";;
      WARN_VLOW) echo -e "${YEL}⚠️ WARNING: Very low r² will keep almost no SNPs.${NC}";;
      ERR_BAD)   echo -e "${RED}❌ Invalid clump_r2 value. Must be a number > 0 and ≤ 1 (e.g., 0.01, 1.0).${NC}"; exit 1;;
    esac
  done

  if ! [[ "$CLUMP_KB" =~ ^[0-9]+$ ]]; then
    echo -e "${RED}❌ Invalid clump_kb value. Must be an integer.${NC}"
    exit 1
  fi

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
else
  CLUMP_KB=""
  CLUMP_R2=""
  LD_POP=""
fi

# ------------------------ Clumping Sensitivity Panel ------------------------
if [[ "$SKIP_LD" != "y" ]]; then
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
if [[ "$SKIP_LD" == "y" ]]; then
  printf "  %-22s %s\n" "LD pruning:" "SKIPPED"
else
  printf "  %-22s %s\n" "LD window (kb):" "$CLUMP_KB"
  printf "  %-22s %s\n" "LD r² cutoff:" "$CLUMP_R2"
  printf "  %-22s %s\n" "LD population:" "$LD_POP"
fi
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
    --skip_ld \"$SKIP_LD\" \
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
# shellcheck disable=SC2086
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

# --- Auto-open HTML report (best effort) ---
FIND_BIN="find"; command -v gfind >/dev/null 2>&1 && FIND_BIN="gfind"
# Try to locate the most recent HTML report in ./results (including subfolders)
REPORT_FILE="$(
  if [[ "$FIND_BIN" = "gfind" ]]; then
    $FIND_BIN ./results -type f -iname '*.html' -printf '%T@ %p\n' 2>/dev/null | sort -nr | awk 'NR==1{$1=""; sub(/^ /,""); print}'
  else
    ls -t ./results/*.html 2>/dev/null | head -n1
  fi
)"
if [[ -n "$REPORT_FILE" && -f "$REPORT_FILE" ]]; then
  echo -e "${GRN}📖 Opening report:${NC} $REPORT_FILE"
  if [[ "$OSTYPE" == "darwin"* ]]; then
    open "$REPORT_FILE" >/dev/null 2>&1 || true
  elif command -v xdg-open >/dev/null 2>&1; then
    xdg-open "$REPORT_FILE" >/dev/null 2>&1 || true
  fi
else
  echo -e "${YEL}ℹ️  No HTML report found to open automatically.${NC}"
fi
