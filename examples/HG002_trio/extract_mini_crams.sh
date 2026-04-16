#!/usr/bin/env bash
# =============================================================================
# extract_mini_crams.sh
#
# Extracts small alignment files (CRAM or BAM) for each trio member
# containing only the reads within ±padding of candidate de novo variant
# sites.  These "mini" files are ideal for IGV review without requiring
# hundreds of gigabytes of full-genome alignment data.
#
# Algorithm
# ---------
#   1. Parse variant positions from the input VCF (CHROM + POS).
#   2. Build a BED file of extraction regions: each variant ± padding bp.
#   3. For each trio member (child, father, mother):
#      a. Extract reads overlapping the regions (samtools view -L).
#      b. Sort the output.
#      c. Index the output.
#   4. Report file sizes and read counts.
#
# Output Format
# -------------
#   When --ref-fasta is provided, output is CRAM (highly compressed).
#   Otherwise output is BAM (no reference required).
#   The --format flag overrides the automatic selection.
#
# Usage
# -----
#   extract_mini_crams.sh \
#       --vcf         candidates.vcf.gz     \
#       --child-bam   child.bam             \
#       --father-bam  father.bam            \
#       --mother-bam  mother.bam            \
#       --output-dir  mini_crams/           \
#       [--ref-fasta  GRCh38.fa]            \
#       [--padding    1000]                 \
#       [--format     cram|bam]             \
#       [--prefix     HG002_trio]
#
# Prerequisites
# -------------
#   • samtools ≥ 1.10
#   • bcftools ≥ 1.10 (or any tool that can read VCF)
#
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
VCF=""
CHILD_BAM=""
FATHER_BAM=""
MOTHER_BAM=""
OUTPUT_DIR=""
REF_FASTA=""
PADDING=1000
FORMAT=""       # auto: cram if ref provided, else bam
PREFIX="mini"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [extract_mini] $*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
    cat <<'EOF'
Usage: extract_mini_crams.sh [OPTIONS]

Extract small alignment files around candidate de novo variant sites.

Required:
  --vcf          FILE    Input VCF with candidate variants (bgzipped or plain)
  --child-bam    FILE    Child BAM/CRAM (indexed)
  --father-bam   FILE    Father BAM/CRAM (indexed)
  --mother-bam   FILE    Mother BAM/CRAM (indexed)
  --output-dir   DIR     Output directory for mini alignment files

Optional:
  --ref-fasta    FILE    Reference FASTA (enables CRAM output; required for
                         CRAM input files)
  --padding      N       Base-pairs of padding around each variant (default: 1000)
  --format       FMT     Output format: "cram" or "bam" (default: cram if
                         --ref-fasta provided, bam otherwise)
  --prefix       STR     Output filename prefix (default: "mini")
  -h, --help             Show this help
EOF
    exit "${1:-0}"
}

check_tool() {
    command -v "$1" >/dev/null 2>&1 \
        || die "$1 is required but not found on PATH"
}

# human_size – format byte count for display
human_size() {
    local bytes="$1"
    if [[ "$bytes" -ge 1073741824 ]]; then
        awk "BEGIN{printf \"%.1f GB\", $bytes/1073741824}"
    elif [[ "$bytes" -ge 1048576 ]]; then
        awk "BEGIN{printf \"%.1f MB\", $bytes/1048576}"
    elif [[ "$bytes" -ge 1024 ]]; then
        awk "BEGIN{printf \"%.1f KB\", $bytes/1024}"
    else
        echo "${bytes} B"
    fi
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)        VCF="${2:-}";        shift 2 ;;
        --child-bam)  CHILD_BAM="${2:-}";  shift 2 ;;
        --father-bam) FATHER_BAM="${2:-}"; shift 2 ;;
        --mother-bam) MOTHER_BAM="${2:-}"; shift 2 ;;
        --output-dir) OUTPUT_DIR="${2:-}"; shift 2 ;;
        --ref-fasta)  REF_FASTA="${2:-}";  shift 2 ;;
        --padding)    PADDING="${2:-}";    shift 2 ;;
        --format)     FORMAT="${2:-}";     shift 2 ;;
        --prefix)     PREFIX="${2:-}";     shift 2 ;;
        -h|--help)    usage 0 ;;
        *)            die "Unknown argument: $1" ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
[[ -n "$VCF" ]]        || die "--vcf is required"
[[ -n "$CHILD_BAM" ]]  || die "--child-bam is required"
[[ -n "$FATHER_BAM" ]] || die "--father-bam is required"
[[ -n "$MOTHER_BAM" ]] || die "--mother-bam is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"

check_tool samtools
check_tool bcftools

[[ -f "$VCF" ]]        || die "VCF not found: $VCF"
[[ -f "$CHILD_BAM" ]]  || die "Child BAM/CRAM not found: $CHILD_BAM"
[[ -f "$FATHER_BAM" ]] || die "Father BAM/CRAM not found: $FATHER_BAM"
[[ -f "$MOTHER_BAM" ]] || die "Mother BAM/CRAM not found: $MOTHER_BAM"

if [[ -n "$REF_FASTA" ]]; then
    [[ -f "$REF_FASTA" ]] || die "Reference FASTA not found: $REF_FASTA"
fi

# Determine output format
if [[ -z "$FORMAT" ]]; then
    if [[ -n "$REF_FASTA" ]]; then
        FORMAT="cram"
    else
        FORMAT="bam"
    fi
fi

case "$FORMAT" in
    cram)
        EXT="cram"
        if [[ -z "$REF_FASTA" ]]; then
            die "CRAM output requires --ref-fasta"
        fi
        ;;
    bam)
        EXT="bam"
        ;;
    *)
        die "Unknown format: $FORMAT (use 'cram' or 'bam')"
        ;;
esac

mkdir -p "$OUTPUT_DIR"

log "=== Extracting mini alignment files ==="
log "  VCF          : $VCF"
log "  Child BAM    : $CHILD_BAM"
log "  Father BAM   : $FATHER_BAM"
log "  Mother BAM   : $MOTHER_BAM"
log "  Output dir   : $OUTPUT_DIR"
log "  Padding      : ±${PADDING} bp"
log "  Format       : $FORMAT"
if [[ -n "$REF_FASTA" ]]; then
    log "  Reference    : $REF_FASTA"
fi

# ---------------------------------------------------------------------------
# Step 1 – Build extraction BED from VCF positions
# ---------------------------------------------------------------------------
log "Step 1: Building extraction regions from VCF ..."

REGIONS_BED="$OUTPUT_DIR/${PREFIX}_regions.bed"

# Extract CHROM and POS from VCF, convert to 0-based BED with padding.
# Clamp start to 0 to avoid negative coordinates.
bcftools query -f '%CHROM\t%POS\n' "$VCF" \
    | awk -v pad="$PADDING" 'BEGIN{OFS="\t"} {
        chrom = $1
        pos   = $2
        start = pos - pad - 1
        if (start < 0) start = 0
        end   = pos + pad
        print chrom, start, end
    }' \
    | sort -k1,1 -k2,2n \
    > "$REGIONS_BED"

NUM_REGIONS=$(wc -l < "$REGIONS_BED" | tr -d ' ')
log "  Extraction regions: $NUM_REGIONS"

if [[ "$NUM_REGIONS" -eq 0 ]]; then
    die "No regions found in VCF: $VCF"
fi

# Merge overlapping regions to minimize redundant extraction
MERGED_BED="$OUTPUT_DIR/${PREFIX}_regions_merged.bed"
bedtools_available=0
if command -v bedtools >/dev/null 2>&1; then
    bedtools_available=1
fi

if [[ "$bedtools_available" -eq 1 ]]; then
    bedtools merge -i "$REGIONS_BED" > "$MERGED_BED"
    MERGED_COUNT=$(wc -l < "$MERGED_BED" | tr -d ' ')
    log "  Merged regions (bedtools): $MERGED_COUNT"
else
    # Fallback: use awk to merge overlapping regions (assumes sorted input)
    awk 'BEGIN{OFS="\t"}
    NR==1 { chr=$1; s=$2; e=$3; next }
    $1==chr && $2<=e { if($3>e) e=$3; next }
    { print chr, s, e; chr=$1; s=$2; e=$3 }
    END { if(NR>0) print chr, s, e }' "$REGIONS_BED" > "$MERGED_BED"
    MERGED_COUNT=$(wc -l < "$MERGED_BED" | tr -d ' ')
    log "  Merged regions (awk): $MERGED_COUNT"
fi

# Use merged BED for extraction
EXTRACT_BED="$MERGED_BED"

# Log total span
TOTAL_SPAN=$(awk '{sum += $3 - $2} END{print sum}' "$EXTRACT_BED")
log "  Total extraction span: ${TOTAL_SPAN} bp"

# ---------------------------------------------------------------------------
# Step 2 – Extract reads for each trio member
# ---------------------------------------------------------------------------
log "Step 2: Extracting reads ..."

# extract_reads – extract, sort, and index reads for one sample
#   $1 = label (e.g. "child")
#   $2 = input BAM/CRAM path
#   $3 = output path (without extension – added automatically)
extract_reads() {
    local label="$1" input="$2" output_base="$3"
    local output="${output_base}.${EXT}"
    local view_args=(-h -L "$EXTRACT_BED" -@ 2)

    if [[ "$FORMAT" == "cram" ]]; then
        view_args+=(--output-fmt CRAM --reference "$REF_FASTA")
    else
        view_args+=(--output-fmt BAM)
    fi

    log "  Extracting ${label} reads ..."
    samtools view "${view_args[@]}" "$input" \
        | samtools sort -@ 2 -o "$output" ${REF_FASTA:+--reference "$REF_FASTA"}

    samtools index "$output"

    local read_count file_size file_size_h
    read_count=$(samtools view -c "$output")
    file_size=$(stat -c %s "$output" 2>/dev/null || stat -f %z "$output")
    file_size_h=$(human_size "$file_size")

    log "    ${label}: ${read_count} reads, ${file_size_h} -> $(basename "$output")"
}

extract_reads "child"  "$CHILD_BAM"  "$OUTPUT_DIR/${PREFIX}_child"
extract_reads "father" "$FATHER_BAM" "$OUTPUT_DIR/${PREFIX}_father"
extract_reads "mother" "$MOTHER_BAM" "$OUTPUT_DIR/${PREFIX}_mother"

# ---------------------------------------------------------------------------
# Step 3 – Summary
# ---------------------------------------------------------------------------
log ""
log "=== Extraction complete ==="
log ""
log "  Output directory: $OUTPUT_DIR"
log ""
log "  Files:"
for f in "$OUTPUT_DIR/${PREFIX}_"*."${EXT}"; do
    [[ -f "$f" ]] || continue
    fsize=$(stat -c %s "$f" 2>/dev/null || stat -f %z "$f")
    fsize_h=$(human_size "$fsize")
    rcount=$(samtools view -c "$f")
    log "    $(basename "$f"): ${rcount} reads, ${fsize_h}"
done
log ""
log "  Regions BED   : $REGIONS_BED"
log "  Merged BED    : $MERGED_BED"
log ""

# Compare to original BAM sizes
log "  Size comparison (mini vs original):"
for pair in "child:$CHILD_BAM:$OUTPUT_DIR/${PREFIX}_child.${EXT}" \
            "father:$FATHER_BAM:$OUTPUT_DIR/${PREFIX}_father.${EXT}" \
            "mother:$MOTHER_BAM:$OUTPUT_DIR/${PREFIX}_mother.${EXT}"; do
    IFS=':' read -r plabel porig pmini <<< "$pair"
    if [[ -f "$porig" && -f "$pmini" ]]; then
        orig_size=$(stat -c %s "$porig" 2>/dev/null || stat -f %z "$porig")
        mini_size=$(stat -c %s "$pmini" 2>/dev/null || stat -f %z "$pmini")
        orig_h=$(human_size "$orig_size")
        mini_h=$(human_size "$mini_size")
        if [[ "$orig_size" -gt 0 ]]; then
            pct=$(awk "BEGIN{printf \"%.2f\", ($mini_size/$orig_size)*100}")
        else
            pct="N/A"
        fi
        log "    ${plabel}: ${orig_h} -> ${mini_h} (${pct}%)"
    fi
done

log ""
log "=== Done ==="
