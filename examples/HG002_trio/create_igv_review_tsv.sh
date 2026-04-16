#!/usr/bin/env bash
# =============================================================================
# create_igv_review_tsv.sh
#
# Generates a variant TSV file for the IGV de novo variant review server
# (https://github.com/jlanej/igv.js/tree/master/server).
#
# The output file contains:
#   • Required columns  – chrom, pos, ref, alt
#   • Quality columns   – quality (QUAL), filter (FILTER), child_gt (GT)
#   • kmer-denovo annotations – DKU, DKT, DKA, DKU_DKT, DKA_DKT, and any
#     other DK* FORMAT fields present in the annotated VCF (e.g. kraken2
#     contamination fractions DKU_NHF, DKA_NHF, …)
#   • Inheritance       – hardcoded "de_novo" (all candidates are putative
#     de novos by construction)
#   • Alignment tracks  – child_file/child_index, father_file/father_index,
#     mother_file/mother_index — pointing to the mini CRAM/BAM files produced
#     by extract_mini_crams.sh
#   • VCF track        – child_vcf/child_vcf_index/child_vcf_id pointing to
#     the bgzipped+indexed annotated VCF from kmer-denovo
#
# Population frequency, variant impact, and gene annotations are NOT available
# from the kmer-denovo pipeline and are therefore omitted.  They can be added
# manually or by a downstream annotation step (e.g. Ensembl VEP) by joining on
# the chrom/pos/ref/alt key.
#
# Usage
# -----
#   create_igv_review_tsv.sh \
#       --vcf         /results/HG002_denovo_annotated.vcf.gz \
#       --mini-dir    /results/mini_crams                    \
#       --output      /results/HG002_igv_review.tsv          \
#       [--prefix     HG002_trio]                            \
#       [--proband-id HG002]
#
# Prerequisites
# -------------
#   • bcftools ≥ 1.10  (for VCF querying and header inspection)
#   • bgzip + tabix    (bundled with htslib/bcftools)
#
# IGV Review Server Quick Start
# ------------------------------
#   # Start the server (requires Node.js):
#   node /path/to/igv.js/server/server.js \
#       --variants   /results/HG002_igv_review.tsv \
#       --data-dir   /results \
#       --genome     hg38 \
#       --port       3000
#
#   # Then open in browser: http://127.0.0.1:3000
#
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
ANNOTATED_VCF=""
MINI_DIR=""
OUTPUT_TSV=""
PREFIX="mini"           # must match the --prefix used with extract_mini_crams.sh
PROBAND_ID="HG002"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [igv_tsv] $*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
    cat <<'EOF'
Usage: create_igv_review_tsv.sh [OPTIONS]

Generate a variant TSV for the IGV de novo variant review server.

Required:
  --vcf          FILE    Annotated VCF from kmer-denovo (.vcf or .vcf.gz)
  --mini-dir     DIR     Directory with mini CRAM/BAM files from
                         extract_mini_crams.sh
  --output       FILE    Output TSV path

Optional:
  --prefix       STR     Filename prefix used by extract_mini_crams.sh
                         (default: "mini"); e.g. "HG002_trio" produces
                         HG002_trio_child.cram / HG002_trio_father.bam …
  --proband-id   ID      Proband sample ID in the VCF (default: HG002)
  -h, --help             Show this help

Notes:
  • If --vcf points to a plain .vcf file it will be bgzipped and tabix-indexed
    as <vcf>.gz (adjacent to the original) so the IGV VCF track can read it.
  • Alignment file paths in the TSV are absolute.
  • Population frequency, variant impact, and gene annotations are not
    populated (not available from kmer-denovo output).
EOF
    exit "${1:-0}"
}

check_tool() {
    command -v "$1" >/dev/null 2>&1 \
        || die "$1 is required but not found on PATH"
}

# abs_path – portable absolute path resolution
abs_path() {
    local p="$1"
    if command -v realpath >/dev/null 2>&1; then
        realpath "$p"
    else
        echo "$(cd "$(dirname "$p")" && pwd)/$(basename "$p")"
    fi
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)        ANNOTATED_VCF="${2:-}"; shift 2 ;;
        --mini-dir)   MINI_DIR="${2:-}";      shift 2 ;;
        --output)     OUTPUT_TSV="${2:-}";    shift 2 ;;
        --prefix)     PREFIX="${2:-}";        shift 2 ;;
        --proband-id) PROBAND_ID="${2:-}";    shift 2 ;;
        -h|--help)    usage 0 ;;
        *)            die "Unknown argument: $1" ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
[[ -n "$ANNOTATED_VCF" ]] || die "--vcf is required"
[[ -n "$MINI_DIR" ]]      || die "--mini-dir is required"
[[ -n "$OUTPUT_TSV" ]]    || die "--output is required"

check_tool bcftools
check_tool bgzip
check_tool tabix

[[ -f "$ANNOTATED_VCF" ]] || die "Annotated VCF not found: $ANNOTATED_VCF"
[[ -d "$MINI_DIR" ]]      || die "Mini alignment directory not found: $MINI_DIR"

# ---------------------------------------------------------------------------
# Step 1 – Ensure VCF is bgzipped + tabix-indexed (required for IGV tracks)
# ---------------------------------------------------------------------------
log "=== Creating IGV review TSV ==="
log "  Annotated VCF : $ANNOTATED_VCF"
log "  Mini dir      : $MINI_DIR"
log "  Prefix        : $PREFIX"
log "  Proband ID    : $PROBAND_ID"
log "  Output        : $OUTPUT_TSV"

if [[ "$ANNOTATED_VCF" == *.vcf.gz ]]; then
    VCF_GZ="$ANNOTATED_VCF"
    log "Step 1: VCF is already bgzipped: $(basename "$VCF_GZ")"
elif [[ "$ANNOTATED_VCF" == *.vcf ]]; then
    VCF_GZ="${ANNOTATED_VCF}.gz"
    if [[ ! -f "$VCF_GZ" ]]; then
        log "Step 1: Compressing VCF with bgzip → $(basename "$VCF_GZ") ..."
        bgzip -@ 2 -c "$ANNOTATED_VCF" > "$VCF_GZ"
    else
        log "Step 1: Found existing bgzipped VCF: $(basename "$VCF_GZ")"
    fi
else
    die "Unrecognised VCF extension (expected .vcf or .vcf.gz): $ANNOTATED_VCF"
fi

if [[ ! -f "${VCF_GZ}.tbi" ]]; then
    log "  Indexing with tabix ..."
    tabix -p vcf "$VCF_GZ"
fi

VCF_ABS="$(abs_path "$VCF_GZ")"
VCF_TBI="${VCF_ABS}.tbi"

# ---------------------------------------------------------------------------
# Step 2 – Detect mini alignment file format (CRAM or BAM)
# ---------------------------------------------------------------------------
log "Step 2: Detecting mini alignment file format ..."

CHILD_CRAM="$MINI_DIR/${PREFIX}_child.cram"
CHILD_BAM="$MINI_DIR/${PREFIX}_child.bam"

if [[ -f "$CHILD_CRAM" ]]; then
    EXT="cram"
    IDX_SUFFIX=".crai"
    log "  Format: CRAM"
elif [[ -f "$CHILD_BAM" ]]; then
    EXT="bam"
    IDX_SUFFIX=".bai"
    log "  Format: BAM"
else
    die "No mini alignment files found in $MINI_DIR with prefix '${PREFIX}'" \
        "(expected ${PREFIX}_child.cram or ${PREFIX}_child.bam)"
fi

# Build absolute paths for all trio members
CHILD_FILE="$(abs_path "$MINI_DIR/${PREFIX}_child.${EXT}")"
CHILD_INDEX="${CHILD_FILE}${IDX_SUFFIX}"
FATHER_FILE="$(abs_path "$MINI_DIR/${PREFIX}_father.${EXT}")"
FATHER_INDEX="${FATHER_FILE}${IDX_SUFFIX}"
MOTHER_FILE="$(abs_path "$MINI_DIR/${PREFIX}_mother.${EXT}")"
MOTHER_INDEX="${MOTHER_FILE}${IDX_SUFFIX}"

# Verify files exist
for f in "$CHILD_FILE" "$CHILD_INDEX" \
         "$FATHER_FILE" "$FATHER_INDEX" \
         "$MOTHER_FILE" "$MOTHER_INDEX"; do
    [[ -f "$f" ]] || die "Required mini alignment file not found: $f"
done

log "  Child  : $CHILD_FILE"
log "  Father : $FATHER_FILE"
log "  Mother : $MOTHER_FILE"

# ---------------------------------------------------------------------------
# Step 3 – Discover kmer-denovo FORMAT fields in VCF header
#
# Include all FORMAT fields with the "DK" prefix:
#   Always present:  DKU  DKT  DKA  DKU_DKT  DKA_DKT
#   With Kraken2:    DKU_BF  DKA_BF  DKU_AF  DKA_AF  DKU_FF  DKA_FF
#                    DKU_PF  DKA_PF  DKU_VF  DKA_VF  DKU_UCF DKA_UCF
#                    DKU_NHF DKA_NHF DKU_UF  DKA_UF  DKU_HLF DKA_HLF
# ---------------------------------------------------------------------------
log "Step 3: Discovering kmer-denovo FORMAT fields ..."

DENOVO_TAGS=()
while IFS= read -r fmt_line; do
    # Extract the ID value from lines like:  ##FORMAT=<ID=DKU,Number=...>
    tag="${fmt_line#*<ID=}"
    tag="${tag%%,*}"
    [[ "$tag" == DK* ]] && DENOVO_TAGS+=("$tag")
done < <(bcftools view -h "$VCF_GZ" 2>/dev/null | grep "^##FORMAT=")

if [[ ${#DENOVO_TAGS[@]} -eq 0 ]]; then
    log "  WARNING: No DK* FORMAT fields found in VCF header."
    log "           TSV will be written without kmer-denovo annotation columns."
else
    log "  Found ${#DENOVO_TAGS[@]} DK* FORMAT fields: ${DENOVO_TAGS[*]}"
fi

# ---------------------------------------------------------------------------
# Step 4 – Build bcftools query format string
#
# Columns extracted per variant:
#   chrom  pos  ref  alt  quality  filter  child_gt  [DK*…]
# ---------------------------------------------------------------------------
# Always-present basic fields (no FORMAT dependency)
QUERY_FMT='%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]'

# Append each DK* FORMAT field using sample-iteration syntax [%TAG]
# For a single-sample VCF this produces a scalar value per field.
for tag in "${DENOVO_TAGS[@]}"; do
    QUERY_FMT+="\\t[%${tag}]"
done
QUERY_FMT+='\n'

# ---------------------------------------------------------------------------
# Step 5 – Build TSV header row
# ---------------------------------------------------------------------------
# Required IGV review columns: chrom, pos, ref, alt
# Recommended: quality, filter, child_gt, inheritance, file/vcf columns
# Extra kmer-denovo annotation columns: lowercased DK* tag names
HEADER="chrom\tpos\tref\talt\tquality\tfilter\tchild_gt"
for tag in "${DENOVO_TAGS[@]}"; do
    HEADER+="\t$(printf '%s' "$tag" | tr '[:upper:]' '[:lower:]')"
done
HEADER+="\tinheritance"
HEADER+="\tchild_file\tchild_index"
HEADER+="\tfather_file\tfather_index"
HEADER+="\tmother_file\tmother_index"
HEADER+="\tchild_vcf\tchild_vcf_index\tchild_vcf_id"

# ---------------------------------------------------------------------------
# Step 4 – Write TSV
# ---------------------------------------------------------------------------
log "Step 4: Writing TSV to $(basename "$OUTPUT_TSV") ..."

mkdir -p "$(dirname "$OUTPUT_TSV")"

# Write header
printf '%b\n' "$HEADER" > "$OUTPUT_TSV"

# Extract per-variant rows from VCF and append fixed path columns.
# awk receives:
#   bcftools query output  (tab-separated variant columns, one line per variant)
#   appended via print:    inheritance + all 9 file/vcf path columns
bcftools query -f "$QUERY_FMT" "$VCF_GZ" \
    | awk \
        -v inherit="de_novo" \
        -v cf="$CHILD_FILE"  -v ci="$CHILD_INDEX"  \
        -v ff="$FATHER_FILE" -v fi="$FATHER_INDEX" \
        -v mf="$MOTHER_FILE" -v mi="$MOTHER_INDEX"  \
        -v vf="$VCF_ABS"     -v vt="$VCF_TBI"       \
        -v pi="$PROBAND_ID"  \
        'BEGIN{OFS="\t"} {
            print $0, inherit, cf, ci, ff, fi, mf, mi, vf, vt, pi
        }' \
    >> "$OUTPUT_TSV"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
N_VARS=$(bcftools view -H "$VCF_GZ" | wc -l)
N_COLS=$(head -1 "$OUTPUT_TSV" | awk -F'\t' '{print NF}')

log ""
log "=== IGV review TSV complete ==="
log ""
log "  Output file  : $OUTPUT_TSV"
log "  Variants     : $N_VARS"
log "  Columns      : $N_COLS"
log ""
log "  Alignment tracks:"
log "    child  : $(basename "$CHILD_FILE")"
log "    father : $(basename "$FATHER_FILE")"
log "    mother : $(basename "$MOTHER_FILE")"
log ""
log "  VCF track    : $(basename "$VCF_ABS") (sample: $PROBAND_ID)"
log ""
log "  Quick start (requires Node.js + igv.js/server):"
log "    node server.js \\"
log "      --variants   $OUTPUT_TSV \\"
log "      --data-dir   $(dirname "$OUTPUT_TSV") \\"
log "      --genome     hg38 \\"
log "      --port       3000"
log "    # Then open: http://127.0.0.1:3000"
log ""
log "=== Done ==="
