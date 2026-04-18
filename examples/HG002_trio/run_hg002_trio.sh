#!/usr/bin/env bash
# =============================================================================
# run_hg002_trio.sh
#
# End-to-end de novo variant filtering for the GIAB HG002 Ashkenazi trio.
# Targets HPC systems using Apptainer (Singularity) via SLURM.
#
# Trio
# ----
#   HG002 / NA24385  –  Son   (child / proband)
#   HG003 / NA24149  –  Father
#   HG004 / NA24143  –  Mother
#
# Pipeline
# --------
#   1. Download GIAB trio alignment files and benchmark VCFs via Aspera
#      (falls back to HTTPS/wget when Aspera is unavailable).
#   2. Identify putative de novo variants by finding child-private sites
#      in the trio VCFs using bcftools.
#   3. Subset the child VCF to putative de novo variants.
#   4. Run kmer-denovo via Apptainer to annotate each candidate with
#      k-mer–based evidence of de novo origin.
#   5. Extract "mini" alignment files (CRAM or BAM) around each candidate
#      for IGV review (±1 kb by default).
#
# Usage
# -----
#   # SLURM submission (recommended):
#   sbatch [--partition=<name>] [--account=<acct>] \
#       examples/HG002_trio/run_hg002_trio.sh \
#       --data-dir /scratch/$USER/hg002_data \
#       --results-dir /scratch/$USER/hg002_results
#
#   # Interactive (e.g. on a login node for testing):
#   bash examples/HG002_trio/run_hg002_trio.sh \
#       --data-dir /scratch/$USER/hg002_data \
#       --results-dir /scratch/$USER/hg002_results
#
# Disk & Time Estimates
# ---------------------
#   • Downloads : ~500 GB (three ~160 GB BAMs + VCFs + indices)
#   • Working   : ~200 GB (jellyfish hashes, intermediate files)
#   • Wall time : 6–24 h depending on network speed and cluster load
#
# Prerequisites
# -------------
#   • Apptainer ≥ 1.1 (or Singularity ≥ 3.8) on PATH
#   • bcftools ≥ 1.10 on PATH (for de novo candidate identification)
#   • Aspera CLI (ascp) on PATH – or wget for HTTPS fallback
#   • tabix and bgzip on PATH (usually bundled with htslib/bcftools)
#
# =============================================================================

# ── SLURM directives (override with sbatch flags) ───────────────────────────
#SBATCH --job-name=hg002-kmer-denovo
#SBATCH --output=hg002_kmer_denovo_%j.log
#SBATCH --error=hg002_kmer_denovo_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

module load bcftools
module load aspera
set -euo pipefail

# ── Constant: in-container path to helper scripts ───────────────────────────
# All companion scripts are packaged in the Docker/Apptainer image at this
# well-known path, guaranteeing they are available regardless of where this
# script is staged (e.g. SLURM spool directories).
CONTAINER_SCRIPT_DIR="/app/examples/HG002_trio"

# ── Configurable defaults ───────────────────────────────────────────────────
# All defaults can be overridden via command-line arguments or environment
# variables (the CLI argument takes precedence).
DATA_DIR="${DATA_DIR:-./hg002_data}"
RESULTS_DIR="${RESULTS_DIR:-./hg002_results}"
TMP_DIR="${TMP_DIR:-}"  # empty = auto (RESULTS_DIR/tmp)
THREADS="${THREADS:-16}"
MEMORY_GB="${MEMORY_GB:-64}"
KMER_SIZE="${KMER_SIZE:-31}"
CONTAINER_URI="${CONTAINER_URI:-docker://ghcr.io/jlanej/kmer_denovo_filter:latest}"
SIF_CACHE="${SIF_CACHE:-./containers}"
ASPERA_KEY="${ASPERA_KEY:-}"        # auto-discovered if empty
ASPERA_MAX_RATE="${ASPERA_MAX_RATE:-500m}"
SKIP_DOWNLOAD="${SKIP_DOWNLOAD:-0}"
FORCE_DOWNLOAD="${FORCE_DOWNLOAD:-0}"
REPORT_ONLY="${REPORT_ONLY:-0}"
REF_FASTA="${REF_FASTA:-}"         # optional; required only for CRAM input
VARIANT_TYPES="${VARIANT_TYPES:-}" # e.g. "snps" or "snps,indels"; empty = all
PROBAND_ID="${PROBAND_ID:-HG002}"
KRAKEN2_DB="${KRAKEN2_DB:-}"       # optional Kraken2 database path
EXTRA_ARGS="${EXTRA_ARGS:-}"       # additional kmer-denovo arguments
MINI_CRAM_PADDING="${MINI_CRAM_PADDING:-1000}"  # ±bp for mini CRAM extraction

# ── GIAB data paths (NCBI FTP) ─────────────────────────────────────────────
NCBI_FTP_HOST="anonftp@ftp.ncbi.nlm.nih.gov"
NCBI_FTP_PORT=33001
GIAB_FTP_BASE="/ReferenceSamples/giab"
GIAB_HTTPS_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"

# Illumina 2×250 bp WGS BAMs (GRCh38, novoalign)
BAM_FTP_BASE="${GIAB_FTP_BASE}/data/AshkenazimTrio"
HG002_BAM_PATH="${BAM_FTP_BASE}/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
HG003_BAM_PATH="${BAM_FTP_BASE}/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.GRCh38.2x250.bam"
HG004_BAM_PATH="${BAM_FTP_BASE}/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.GRCh38.2x250.bam"

# GIAB v4.2.1 benchmark VCFs (GRCh38)
VCF_FTP_BASE="${GIAB_FTP_BASE}/release/AshkenazimTrio"
HG002_VCF_PATH="${VCF_FTP_BASE}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG003_VCF_PATH="${VCF_FTP_BASE}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG004_VCF_PATH="${VCF_FTP_BASE}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

# ── Helper functions ────────────────────────────────────────────────────────

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
die()  { log "ERROR: $*"; exit 1; }

usage() {
    cat <<'EOF'
Usage: run_hg002_trio.sh [OPTIONS]

End-to-end de novo variant filtering for the GIAB HG002 Ashkenazi trio.

Data & Output:
  --data-dir DIR          Download directory for BAMs/VCFs (default: ./hg002_data)
  --results-dir DIR       Results directory (default: ./hg002_results)
  --tmp-dir DIR           Temp directory for jellyfish hashes; avoid RAM-backed
                          filesystems (default: RESULTS_DIR/tmp)

Compute:
  --threads N             Thread count for jellyfish/kmer-denovo (default: 16)
  --memory N              Available memory in GB (default: 64)
  --kmer-size N           K-mer size (default: 31)

Container:
  --container URI         Apptainer image URI or local .sif path
                          (default: docker://ghcr.io/jlanej/kmer_denovo_filter:latest)
  --sif-cache DIR         Directory to cache pulled .sif files (default: ./containers)

Download:
  --aspera-key PATH       Path to Aspera SSH key (auto-discovered if omitted)
  --aspera-max-rate RATE  Aspera max transfer rate (default: 500m)
  --skip-download         Skip all downloads (use pre-existing files)
  --force-download        Force re-download even if files exist

Report:
  --report-only           Skip all pipeline steps and only regenerate the HTML
                          report from existing output files in --results-dir.
                          Requires a prior successful pipeline run.
                          Equivalent to running kmer-report directly.

Analysis:
  --ref-fasta PATH        Reference FASTA (required for CRAM; optional for BAM)
  --variant-types TYPES   Comma-separated variant types to include in de novo
                          scan (e.g. "snps", "snps,indels"; default: all)
  --proband-id ID         Proband sample ID in VCF (default: HG002)
  --kraken2-db PATH       Optional Kraken2 database path for non-human
                          content annotations
  --extra-args "ARGS"     Additional arguments passed to kmer-denovo
  --mini-cram-padding N   Padding in bp for mini CRAM extraction (default: 1000)

General:
  -h, --help              Show this help
EOF
    exit 0
}

check_tool() {
    local tool="$1" required="${2:-1}"
    if command -v "$tool" >/dev/null 2>&1; then
        return 0
    elif [[ "$required" -eq 1 ]]; then
        die "$tool is required but not found on PATH"
    else
        return 1
    fi
}

# find_apptainer – locate apptainer or singularity
find_apptainer() {
    if command -v apptainer >/dev/null 2>&1; then
        echo "apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        echo "singularity"
    else
        die "Neither apptainer nor singularity found on PATH. "\
            "Load the appropriate module (e.g. 'module load apptainer')."
    fi
}

# find_aspera_key – search common locations for the Aspera SSH key
find_aspera_key() {
    local search_paths=(
        "${ASPERA_KEY:-}"
        "${CONDA_PREFIX:-}/etc/asperaweb_id_dsa.openssh"
        "$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh"
        "/opt/aspera/connect/etc/asperaweb_id_dsa.openssh"
        "/usr/local/etc/asperaweb_id_dsa.openssh"
    )
    for p in "${search_paths[@]}"; do
        [[ -n "$p" && -f "$p" ]] && { echo "$p"; return 0; }
    done
    return 1
}

# has_required_kraken2_files – validate core Kraken2 DB artifacts
has_required_kraken2_files() {
    local db_dir="$1"
    local req
    for req in hash.k2d opts.k2d taxo.k2d; do
        [[ -f "$db_dir/$req" ]] || return 1
    done
    return 0
}

# resolve_kraken2_db_dir – handle DB root or nested extracted subdirectory
resolve_kraken2_db_dir() {
    local db_path="$1"
    local resolved=""
    local candidate=""
    local -a _db_candidates=()
    local -a _matching_candidates=()

    [[ -d "$db_path" ]] || die "Kraken2 database path is not a directory: $db_path"

    if has_required_kraken2_files "$db_path"; then
        (cd "$db_path" && pwd)
        return 0
    fi

    mapfile -t _db_candidates < <(
        find "$db_path" -mindepth 1 -maxdepth 3 -type f -name "hash.k2d" \
            -exec dirname {} + | sort -u
    )

    for candidate in "${_db_candidates[@]}"; do
        if has_required_kraken2_files "$candidate"; then
            _matching_candidates+=("$candidate")
        fi
    done

    if [[ ${#_matching_candidates[@]} -eq 1 ]]; then
        resolved="${_matching_candidates[0]}"
        (cd "$resolved" && pwd)
        return 0
    fi

    if [[ ${#_matching_candidates[@]} -gt 1 ]]; then
        log "Multiple Kraken2 database directories detected under: $db_path"
        for candidate in "${_matching_candidates[@]}"; do
            log "  - $candidate"
        done
        die "Please set --kraken2-db to the specific directory containing hash.k2d/opts.k2d/taxo.k2d"
    fi

    die "Could not find Kraken2 DB files (hash.k2d, opts.k2d, taxo.k2d) under: $db_path"
}

# download_file – idempotent download via Aspera with HTTPS/wget fallback
#   $1 = FTP path (relative to NCBI FTP root)
#   $2 = destination file path
download_file() {
    local ftp_path="$1" dest="$2"
    local filename
    filename="$(basename "$dest")"

    # Idempotency: skip if file exists and is non-empty
    if [[ "$FORCE_DOWNLOAD" -ne 1 && -f "$dest" && -s "$dest" ]]; then
        log "  [skip] $filename (already exists)"
        return 0
    fi

    mkdir -p "$(dirname "$dest")"

    # Try Aspera first
    if [[ "$USE_ASPERA" -eq 1 ]]; then
        log "  [aspera] Downloading $filename ..."
        if ascp -v -i "$RESOLVED_ASPERA_KEY" \
                -k 1 -T -l "$ASPERA_MAX_RATE" -P "$NCBI_FTP_PORT" \
                "${NCBI_FTP_HOST}:${ftp_path}" \
                "$dest" 2>&1 | tail -5; then
            log "  [aspera] $filename complete"
            return 0
        else
            log "  [aspera] Failed – falling back to HTTPS for $filename"
        fi
    fi

    # HTTPS fallback via wget
    local url="${GIAB_HTTPS_BASE}${ftp_path#"$GIAB_FTP_BASE"}"
    log "  [https] Downloading $filename from $url ..."
    wget -c -q --show-progress -O "$dest" "$url" \
        || die "Failed to download $filename via HTTPS"
    log "  [https] $filename complete"
}

# ── Parse command-line arguments ────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)         DATA_DIR="${2:-}";         shift 2 ;;
        --results-dir)      RESULTS_DIR="${2:-}";      shift 2 ;;
        --tmp-dir)          TMP_DIR="${2:-}";           shift 2 ;;
        --threads)          THREADS="${2:-}";           shift 2 ;;
        --memory)           MEMORY_GB="${2:-}";         shift 2 ;;
        --kmer-size)        KMER_SIZE="${2:-}";         shift 2 ;;
        --container)        CONTAINER_URI="${2:-}";     shift 2 ;;
        --sif-cache)        SIF_CACHE="${2:-}";         shift 2 ;;
        --aspera-key)       ASPERA_KEY="${2:-}";        shift 2 ;;
        --aspera-max-rate)  ASPERA_MAX_RATE="${2:-}";   shift 2 ;;
        --skip-download)    SKIP_DOWNLOAD=1;           shift ;;
        --force-download)   FORCE_DOWNLOAD=1;          shift ;;
        --report-only)      REPORT_ONLY=1;             shift ;;
        --ref-fasta)        REF_FASTA="${2:-}";         shift 2 ;;
        --variant-types)    VARIANT_TYPES="${2:-}";     shift 2 ;;
        --proband-id)       PROBAND_ID="${2:-}";        shift 2 ;;
        --kraken2-db)       KRAKEN2_DB="${2:-}";        shift 2 ;;
        --extra-args)       EXTRA_ARGS="${2:-}";        shift 2 ;;
        --mini-cram-padding) MINI_CRAM_PADDING="${2:-}"; shift 2 ;;
        -h|--help)          usage ;;
        *)                  die "Unknown argument: $1" ;;
    esac
done

# Resolve TMP_DIR default
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${RESULTS_DIR}/tmp"
fi

# ── Banner ──────────────────────────────────────────────────────────────────
log "========================================================================"
log "  HG002 Trio – End-to-End De Novo Variant Filtering"
log "========================================================================"
log "  Data dir     : $DATA_DIR"
log "  Results dir  : $RESULTS_DIR"
log "  Tmp dir      : $TMP_DIR"
log "  Threads      : $THREADS"
log "  Memory       : ${MEMORY_GB} GB"
log "  K-mer size   : $KMER_SIZE"
log "  Container    : $CONTAINER_URI"
log "  Proband ID   : $PROBAND_ID"
log "  Kraken2 DB   : ${KRAKEN2_DB:-"(disabled)"}"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    log "  SLURM job    : $SLURM_JOB_ID"
    log "  SLURM node   : ${SLURM_NODELIST:-unknown}"
fi
log "========================================================================"

# ── Preflight checks ───────────────────────────────────────────────────────
log "Checking prerequisites ..."

APPTAINER_CMD="$(find_apptainer)"
log "  Container runtime : $APPTAINER_CMD ($(command -v "$APPTAINER_CMD"))"

check_tool bcftools
log "  bcftools          : $(bcftools --version | head -1)"

check_tool tabix
check_tool bgzip

# Aspera: optional – fall back to wget
USE_ASPERA=0
RESOLVED_ASPERA_KEY=""
if check_tool ascp 0; then
    if RESOLVED_ASPERA_KEY="$(find_aspera_key)"; then
        USE_ASPERA=1
        log "  Aspera            : ascp (key: $RESOLVED_ASPERA_KEY)"
    else
        log "  Aspera            : ascp found but no SSH key; using HTTPS fallback"
    fi
else
    log "  Aspera            : not found; using HTTPS/wget fallback"
    check_tool wget
fi

RESOLVED_KRAKEN2_DB=""
if [[ -n "$KRAKEN2_DB" ]]; then
    RESOLVED_KRAKEN2_DB="$(resolve_kraken2_db_dir "$KRAKEN2_DB")"
    log "  Kraken2 DB        : $RESOLVED_KRAKEN2_DB"
fi

# ── Create directories ─────────────────────────────────────────────────────
mkdir -p "$DATA_DIR/bams" "$DATA_DIR/vcfs" "$RESULTS_DIR" "$TMP_DIR" "$SIF_CACHE"

# ============================================================================
# REPORT-ONLY MODE – regenerate the HTML report from existing output files
# ============================================================================
if [[ "$REPORT_ONLY" -eq 1 ]]; then
    log ""
    log "Report-only mode: regenerating HTML report from existing outputs ..."
    METRICS_JSON="$RESULTS_DIR/HG002_metrics.json"
    SUMMARY_TXT="$RESULTS_DIR/HG002_summary.txt"
    OUTPUT_VCF="$RESULTS_DIR/HG002_denovo_annotated.vcf"
    REPORT_HTML="$RESULTS_DIR/HG002_report.html"

    # Verify required files exist
    [[ -f "$METRICS_JSON" ]] || die "Metrics file not found: $METRICS_JSON (run the full pipeline first)"
    [[ -f "$SUMMARY_TXT"  ]] || die "Summary file not found: $SUMMARY_TXT (run the full pipeline first)"

    # Prepare Apptainer container (pull if needed) so kmer-report is available
    if [[ "$CONTAINER_URI" == *.sif ]]; then
        SIF_FILE="$CONTAINER_URI"
        [[ -f "$SIF_FILE" ]] || die "Container .sif not found: $SIF_FILE"
    else
        APPTAINER_CMD="$(find_apptainer)"
        SIF_NAME="kmer_denovo_filter_${CONTAINER_URI//[:\/@]/_}.sif"
        SIF_FILE="${SIF_CACHE}/${SIF_NAME}"
        if [[ ! -f "$SIF_FILE" ]]; then
            log "  Pulling container: $CONTAINER_URI"
            "$APPTAINER_CMD" pull "$SIF_FILE" "$CONTAINER_URI" \
                || die "Failed to pull container: $CONTAINER_URI"
        fi
    fi

    APPTAINER_CMD="${APPTAINER_CMD:-$(find_apptainer)}"
    declare -A BIND_DIRS_RO
    BIND_DIRS_RO["$(cd "$RESULTS_DIR" && pwd)"]=1
    BIND_ARGS_RO=""
    for d in "${!BIND_DIRS_RO[@]}"; do
        BIND_ARGS_RO="${BIND_ARGS_RO:+${BIND_ARGS_RO},}${d}"
    done

    REPORT_CMD=(
        kmer-report
        --output        "$REPORT_HTML"
        --vcf-metrics   "$METRICS_JSON"
        --vcf-summary   "$SUMMARY_TXT"
    )
    if [[ -f "$OUTPUT_VCF" || -f "${OUTPUT_VCF}.gz" ]]; then
        report_vcf="$OUTPUT_VCF"
        [[ -f "${OUTPUT_VCF}.gz" ]] && report_vcf="${OUTPUT_VCF}.gz"
        REPORT_CMD+=(--vcf "$report_vcf")
    fi

    log "  Command: ${REPORT_CMD[*]}"
    "$APPTAINER_CMD" exec \
        --bind "$BIND_ARGS_RO" \
        "$SIF_FILE" \
        "${REPORT_CMD[@]}"

    log ""
    log "Report regenerated: $REPORT_HTML"
    log "Done."
    exit 0
fi

# ============================================================================
# STEP 1 – Download GIAB trio data
# ============================================================================
if [[ "$SKIP_DOWNLOAD" -eq 1 ]]; then
    log ""
    log "Step 1: SKIPPED (--skip-download)"
else
    log ""
    log "Step 1: Downloading GIAB HG002 trio data ..."
    log "  BAMs: NIST Illumina 2×250 bp WGS (GRCh38, novoalign)"
    log "  VCFs: GIAB v4.2.1 benchmark (GRCh38, chr1-22)"
    log ""

    # ── BAMs + indices ──────────────────────────────────────────────────────
    log "  --- Alignment files (BAMs) ---"
    download_file "$HG002_BAM_PATH"         "$DATA_DIR/bams/HG002.GRCh38.2x250.bam"
    download_file "${HG002_BAM_PATH}.bai"   "$DATA_DIR/bams/HG002.GRCh38.2x250.bam.bai"

    download_file "$HG003_BAM_PATH"         "$DATA_DIR/bams/HG003.GRCh38.2x250.bam"
    download_file "${HG003_BAM_PATH}.bai"   "$DATA_DIR/bams/HG003.GRCh38.2x250.bam.bai"

    download_file "$HG004_BAM_PATH"         "$DATA_DIR/bams/HG004.GRCh38.2x250.bam"
    download_file "${HG004_BAM_PATH}.bai"   "$DATA_DIR/bams/HG004.GRCh38.2x250.bam.bai"

    # ── VCFs + indices ──────────────────────────────────────────────────────
    log ""
    log "  --- Benchmark VCFs ---"
    download_file "$HG002_VCF_PATH"         "$DATA_DIR/vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    download_file "${HG002_VCF_PATH}.tbi"   "$DATA_DIR/vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"

    download_file "$HG003_VCF_PATH"         "$DATA_DIR/vcfs/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    download_file "${HG003_VCF_PATH}.tbi"   "$DATA_DIR/vcfs/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"

    download_file "$HG004_VCF_PATH"         "$DATA_DIR/vcfs/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    download_file "${HG004_VCF_PATH}.tbi"   "$DATA_DIR/vcfs/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"

    log ""
    log "Step 1: Downloads complete."
fi

# Resolve local file paths
CHILD_BAM="$DATA_DIR/bams/HG002.GRCh38.2x250.bam"
FATHER_BAM="$DATA_DIR/bams/HG003.GRCh38.2x250.bam"
MOTHER_BAM="$DATA_DIR/bams/HG004.GRCh38.2x250.bam"
CHILD_VCF="$DATA_DIR/vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
FATHER_VCF="$DATA_DIR/vcfs/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
MOTHER_VCF="$DATA_DIR/vcfs/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

# Verify files exist
for f in "$CHILD_BAM" "$FATHER_BAM" "$MOTHER_BAM" \
         "$CHILD_VCF" "$FATHER_VCF" "$MOTHER_VCF"; do
    [[ -f "$f" ]] || die "Required file not found: $f"
done

# ============================================================================
# STEP 2 – Pull / cache Apptainer container
# ============================================================================
log ""
log "Step 2: Preparing Apptainer container ..."

if [[ "$CONTAINER_URI" == *.sif ]]; then
    # User provided a local .sif file
    SIF_FILE="$CONTAINER_URI"
    [[ -f "$SIF_FILE" ]] || die "Container .sif not found: $SIF_FILE"
else
    # Pull from registry and cache locally
    SIF_NAME="kmer_denovo_filter_${CONTAINER_URI//[:\/@]/_}.sif"
    SIF_FILE="${SIF_CACHE}/${SIF_NAME}"
    if [[ -f "$SIF_FILE" ]]; then
        log "  Container already cached: $SIF_FILE"
    else
        log "  Pulling container: $CONTAINER_URI"
        "$APPTAINER_CMD" pull "$SIF_FILE" "$CONTAINER_URI" \
            || die "Failed to pull container: $CONTAINER_URI"
        log "  Cached at: $SIF_FILE"
    fi
fi

# ── Build bind mount list ────────────────────────────────────────────────────
# Compute once here so it is available for all apptainer exec calls below.
# Apptainer needs --bind for any path not in its default bind list.
declare -A BIND_DIRS
for d in "$(cd "$DATA_DIR" && pwd)" \
         "$(cd "$RESULTS_DIR" && pwd)" \
         "$(cd "$TMP_DIR" && pwd)"; do
    BIND_DIRS["$d"]=1
done
if [[ -n "$REF_FASTA" && -f "$REF_FASTA" ]]; then
    BIND_DIRS["$(cd "$(dirname "$REF_FASTA")" && pwd)"]=1
fi
if [[ -n "$RESOLVED_KRAKEN2_DB" ]]; then
    BIND_DIRS["$RESOLVED_KRAKEN2_DB"]=1
fi

BIND_ARGS=""
for d in "${!BIND_DIRS[@]}"; do
    BIND_ARGS="${BIND_ARGS:+${BIND_ARGS},}${d}"
done

# ============================================================================
# STEP 3 – Identify putative de novo variants
# ============================================================================
log ""
log "Step 3: Identifying putative de novo variants ..."

DENOVO_VCF="$RESULTS_DIR/putative_denovos.vcf.gz"

HELPER_ARGS=(
    --child-vcf  "$CHILD_VCF"
    --father-vcf "$FATHER_VCF"
    --mother-vcf "$MOTHER_VCF"
    --output     "$DENOVO_VCF"
    --tmp-dir    "$TMP_DIR/identify_denovos"
)
if [[ -n "$VARIANT_TYPES" ]]; then
    HELPER_ARGS+=(--variant-types "$VARIANT_TYPES")
fi

"$APPTAINER_CMD" exec \
    --bind "$BIND_ARGS" \
    "$SIF_FILE" \
    bash "${CONTAINER_SCRIPT_DIR}/identify_putative_denovos.sh" "${HELPER_ARGS[@]}"

DENOVO_COUNT=$(bcftools view -H "$DENOVO_VCF" | wc -l)
log "  Putative de novo variants: $DENOVO_COUNT"

if [[ "$DENOVO_COUNT" -eq 0 ]]; then
    die "No putative de novo variants found. Check input VCFs."
fi

# ============================================================================
# STEP 4 – Run kmer-denovo via Apptainer
# ============================================================================
log ""
log "Step 4: Running kmer-denovo via Apptainer ..."

OUTPUT_VCF="$RESULTS_DIR/HG002_denovo_annotated.vcf"
METRICS_JSON="$RESULTS_DIR/HG002_metrics.json"
SUMMARY_TXT="$RESULTS_DIR/HG002_summary.txt"
INFO_READS_BAM="$RESULTS_DIR/HG002_informative_reads.bam"
REPORT_HTML="$RESULTS_DIR/HG002_report.html"

# Build kmer-denovo command
KMER_CMD=(
    kmer-denovo
    --child   "$CHILD_BAM"
    --mother  "$MOTHER_BAM"
    --father  "$FATHER_BAM"
    --vcf     "$DENOVO_VCF"
    --output  "$OUTPUT_VCF"
    --proband-id "$PROBAND_ID"
    --threads "$THREADS"
    --memory  "$MEMORY_GB"
    --kmer-size "$KMER_SIZE"
    --metrics "$METRICS_JSON"
    --summary "$SUMMARY_TXT"
    --informative-reads "$INFO_READS_BAM"
    --report  "$REPORT_HTML"
    --tmp-dir "$TMP_DIR/kmer_denovo"
)

if [[ -n "$REF_FASTA" ]]; then
    KMER_CMD+=(--ref-fasta "$REF_FASTA")
fi
if [[ -n "$RESOLVED_KRAKEN2_DB" ]]; then
    KMER_CMD+=(--kraken2-db "$RESOLVED_KRAKEN2_DB")
fi

# Append any extra user-supplied arguments
if [[ -n "$EXTRA_ARGS" ]]; then
    # shellcheck disable=SC2206
    KMER_CMD+=($EXTRA_ARGS)
fi

log "  Container : $SIF_FILE"
log "  Bind dirs : $BIND_ARGS"
log "  Command   : ${KMER_CMD[*]}"
log ""

mkdir -p "$TMP_DIR/kmer_denovo"

"$APPTAINER_CMD" exec \
    --bind "$BIND_ARGS" \
    "$SIF_FILE" \
    "${KMER_CMD[@]}"

# ============================================================================
# STEP 5 – Extract mini alignment files for IGV review
# ============================================================================
log ""
log "Step 5: Extracting mini alignment files (±${MINI_CRAM_PADDING} bp) ..."

MINI_DIR="$RESULTS_DIR/mini_crams"

EXTRACT_ARGS=(
    --vcf        "$DENOVO_VCF"
    --child-bam  "$CHILD_BAM"
    --father-bam "$FATHER_BAM"
    --mother-bam "$MOTHER_BAM"
    --output-dir "$MINI_DIR"
    --padding    "$MINI_CRAM_PADDING"
    --prefix     "HG002_trio"
)
if [[ -n "$REF_FASTA" ]]; then
    EXTRACT_ARGS+=(--ref-fasta "$REF_FASTA")
fi

"$APPTAINER_CMD" exec \
    --bind "$BIND_ARGS" \
    "$SIF_FILE" \
    bash "${CONTAINER_SCRIPT_DIR}/extract_mini_crams.sh" "${EXTRACT_ARGS[@]}"

# ============================================================================
# STEP 6 – Create IGV variant review TSV
# ============================================================================
log ""
log "Step 6: Creating IGV variant review TSV ..."

IGV_TSV="$RESULTS_DIR/HG002_igv_review.tsv"

# Use the annotated VCF as the source; the script handles bgzip/tabix if needed
"$APPTAINER_CMD" exec \
    --bind "$BIND_ARGS" \
    "$SIF_FILE" \
    bash "${CONTAINER_SCRIPT_DIR}/create_igv_review_tsv.sh" \
        --vcf         "$OUTPUT_VCF" \
        --mini-dir    "$MINI_DIR"   \
        --prefix      "HG002_trio"  \
        --output      "$IGV_TSV"    \
        --proband-id  "$PROBAND_ID"

# ============================================================================
# STEP 7 – Report results
# ============================================================================
log ""
log "========================================================================"
log "  Pipeline complete!"
log "========================================================================"
log ""
log "  Results directory: $RESULTS_DIR"
log ""
log "  Output files:"
log "    Annotated VCF        : $OUTPUT_VCF"
log "    Metrics (JSON)       : $METRICS_JSON"
log "    Summary              : $SUMMARY_TXT"
log "    Interactive report   : $REPORT_HTML"
log "    Informative reads    : $INFO_READS_BAM"
log "    Putative de novos    : $DENOVO_VCF"
log "    Mini alignments dir  : $MINI_DIR"
log "    IGV review TSV       : $IGV_TSV"
log ""

if [[ -f "$SUMMARY_TXT" ]]; then
    log "  --- Summary excerpt ---"
    head -30 "$SUMMARY_TXT" | while IFS= read -r line; do
        log "  $line"
    done
    log "  --- (see $SUMMARY_TXT for full details) ---"
fi

log ""
log "Done."
