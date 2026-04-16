#!/usr/bin/env bash
# =============================================================================
# identify_putative_denovos.sh
#
# Identifies putative de novo variants from a trio of VCFs by finding
# child-private sites: variants present in the child but absent from
# both parents.
#
# This script uses bcftools isec to perform the comparison in two passes:
#   1. Remove father variants from the child VCF
#   2. Remove mother variants from the result
#
# The output is a bgzipped, tabix-indexed VCF of putative de novo variants
# suitable for input to kmer-denovo.
#
# Usage
# -----
#   identify_putative_denovos.sh \
#       --child-vcf   child.vcf.gz  \
#       --father-vcf  father.vcf.gz \
#       --mother-vcf  mother.vcf.gz \
#       --output      putative_denovos.vcf.gz \
#       [--tmp-dir /scratch/tmp] \
#       [--variant-types snps,indels]
#
# Prerequisites
# -------------
#   • bcftools ≥ 1.10
#   • htslib (bgzip, tabix) – usually bundled with bcftools
#
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
CHILD_VCF=""
FATHER_VCF=""
MOTHER_VCF=""
OUTPUT_VCF=""
TMP_DIR=""
VARIANT_TYPES=""   # e.g. "snps" or "snps,indels"; empty = all types

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [identify_denovos] $*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
    cat <<'EOF'
Usage: identify_putative_denovos.sh [OPTIONS]

Identify child-private variants (putative de novos) from trio VCFs.

Required:
  --child-vcf   FILE    Child VCF (bgzipped + tabix-indexed)
  --father-vcf  FILE    Father VCF (bgzipped + tabix-indexed)
  --mother-vcf  FILE    Mother VCF (bgzipped + tabix-indexed)
  --output      FILE    Output VCF path (.vcf.gz)

Optional:
  --tmp-dir     DIR     Directory for intermediate files (default: auto)
  --variant-types TYPE  Comma-separated variant types to include
                        (e.g. "snps", "snps,indels"; default: all)
  -h, --help            Show this help
EOF
    exit "${1:-0}"
}

check_tool() {
    command -v "$1" >/dev/null 2>&1 \
        || die "$1 is required but not found on PATH"
}

check_indexed_vcf() {
    local vcf="$1" label="$2"
    [[ -f "$vcf" ]] || die "${label} not found: ${vcf}"
    # Check for tabix index (.tbi) or CSI index (.csi)
    if [[ ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]]; then
        log "WARNING: No index found for ${label} (${vcf}). Attempting to create one ..."
        tabix -p vcf "$vcf" \
            || die "Failed to index ${label}: ${vcf}. Ensure it is bgzipped."
    fi
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --child-vcf)    CHILD_VCF="${2:-}";    shift 2 ;;
        --father-vcf)   FATHER_VCF="${2:-}";   shift 2 ;;
        --mother-vcf)   MOTHER_VCF="${2:-}";   shift 2 ;;
        --output)       OUTPUT_VCF="${2:-}";    shift 2 ;;
        --tmp-dir)      TMP_DIR="${2:-}";       shift 2 ;;
        --variant-types) VARIANT_TYPES="${2:-}"; shift 2 ;;
        -h|--help)      usage 0 ;;
        *)              die "Unknown argument: $1" ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
[[ -n "$CHILD_VCF" ]]  || die "--child-vcf is required"
[[ -n "$FATHER_VCF" ]] || die "--father-vcf is required"
[[ -n "$MOTHER_VCF" ]] || die "--mother-vcf is required"
[[ -n "$OUTPUT_VCF" ]] || die "--output is required"

check_tool bcftools
check_tool bgzip
check_tool tabix

check_indexed_vcf "$CHILD_VCF"  "child VCF"
check_indexed_vcf "$FATHER_VCF" "father VCF"
check_indexed_vcf "$MOTHER_VCF" "mother VCF"

# ---------------------------------------------------------------------------
# Set up working directory
# ---------------------------------------------------------------------------
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="$(mktemp -d -t identify_denovos_XXXXXX)"
    CLEANUP_TMP=1
else
    mkdir -p "$TMP_DIR"
    CLEANUP_TMP=0
fi

cleanup() {
    if [[ "${CLEANUP_TMP:-0}" -eq 1 && -d "${TMP_DIR:-}" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT

log "=== Identifying putative de novo variants ==="
log "  Child VCF  : $CHILD_VCF"
log "  Father VCF : $FATHER_VCF"
log "  Mother VCF : $MOTHER_VCF"
log "  Output     : $OUTPUT_VCF"
log "  Tmp dir    : $TMP_DIR"

# Optional: filter to specific variant types before comparison
INPUT_VCF="$CHILD_VCF"
if [[ -n "$VARIANT_TYPES" ]]; then
    log "  Filtering child VCF to variant types: $VARIANT_TYPES"
    FILTERED="$TMP_DIR/child_filtered.vcf.gz"
    bcftools view -v "$VARIANT_TYPES" "$CHILD_VCF" -Oz -o "$FILTERED"
    tabix -p vcf "$FILTERED"
    INPUT_VCF="$FILTERED"
fi

# ---------------------------------------------------------------------------
# Step 1 – Remove variants present in father
# ---------------------------------------------------------------------------
log "Step 1: Removing variants present in father ..."

STEP1_OUT="$TMP_DIR/child_not_father.vcf.gz"

bcftools isec -C \
    "$INPUT_VCF" "$FATHER_VCF" \
    -Oz -o "$STEP1_OUT"

tabix -p vcf "$STEP1_OUT"

CHILD_COUNT=$(bcftools view -H "$INPUT_VCF" | wc -l)
STEP1_COUNT=$(bcftools view -H "$STEP1_OUT" | wc -l)
log "  Child variants         : $CHILD_COUNT"
log "  After removing father  : $STEP1_COUNT"

# ---------------------------------------------------------------------------
# Step 2 – Remove variants present in mother
# ---------------------------------------------------------------------------
log "Step 2: Removing variants present in mother ..."

bcftools isec -C \
    "$STEP1_OUT" "$MOTHER_VCF" \
    -Oz -o "$OUTPUT_VCF"

tabix -p vcf "$OUTPUT_VCF"

FINAL_COUNT=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
log "  After removing mother  : $FINAL_COUNT"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
FATHER_REMOVED=$((CHILD_COUNT - STEP1_COUNT))
MOTHER_REMOVED=$((STEP1_COUNT - FINAL_COUNT))
TOTAL_REMOVED=$((CHILD_COUNT - FINAL_COUNT))

log "=== Summary ==="
log "  Total child variants       : $CHILD_COUNT"
log "  Shared with father         : $FATHER_REMOVED"
log "  Shared with mother (only)  : $MOTHER_REMOVED"
log "  Total inherited (removed)  : $TOTAL_REMOVED"
log "  Putative de novos          : $FINAL_COUNT"
log "  Output: $OUTPUT_VCF"
log "=== Done ==="
