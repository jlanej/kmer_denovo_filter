#!/usr/bin/env bash
# =============================================================================
# download_giab_dnm_testdata.sh
#
# Discovers child-private SNVs from the GIAB HG002 (NA24385) Ashkenazi trio
# and extracts small BAM regions around them for integration testing.
#
# Overview
# --------
# The Genome in a Bottle (GIAB) consortium provides extensively characterized
# benchmark data for the Ashkenazi Jewish trio:
#
#   HG002 / NA24385  –  Son   (child / proband)
#   HG003 / NA24149  –  Father
#   HG004 / NA24143  –  Mother
#
# This script dynamically discovers SNVs that are present in the HG002
# benchmark VCF (v4.2.1) but absent from both parents' benchmark VCFs.
# These child-private variants serve as realistic candidates for testing
# de novo mutation filtering tools.
#
# NO large file downloads are performed.  All VCF and BAM queries use
# htslib HTTPS random access against the public GIAB FTP, transferring only
# a few MB of data in total.
#
# Algorithm
# ---------
#   1. Stream small windows (~50 kb) from the HG002 benchmark VCF across
#      multiple chromosomes, collecting SNVs via HTTPS random access.
#   2. For each HG002 SNV, verify via bcftools that it is absent from
#      both HG003 and HG004 benchmark VCFs (child-private).
#   3. Select a subset of verified candidates (default 5).
#   4. Extract ±500 bp BAM slices around each candidate from the public GIAB
#      Illumina 2×250 bp WGS BAMs (remote HTTPS via htslib).
#   5. Write a candidate VCF containing the selected child-private SNVs.
#
# Prerequisites
# -------------
#   • samtools  ≥ 1.10  (compiled with htslib HTTP/S support)
#   • bcftools  ≥ 1.10
#
# Usage
# -----
#   ./scripts/download_giab_dnm_testdata.sh \
#       [-o output_dir]                     \
#       [-n num_variants]                   \
#       [-p padding_bp]
#
# Options
# -------
#   -o OUTPUT_DIR  Output directory  (default: tests/data/giab)
#   -n NUM         Number of DNM candidates to select  (default: 5)
#   -p PADDING     Base-pairs of padding around each variant  (default: 500)
#
# Data sources
# ------------
#   BAMs  :  NIST Illumina 2×250 bp WGS, novoalign, GRCh38
#   VCFs  :  GIAB v4.2.1 benchmark  (GRCh38, chr1-22)  [queried over HTTPS]
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_DIR="${REPO_DIR}/tests/data/giab"
NUM_VARIANTS=20
PADDING=500

# ---------------------------------------------------------------------------
# GIAB data URLs (all accessed via HTTPS random access – no bulk downloads)
# ---------------------------------------------------------------------------
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"

# Illumina 2×250 bp WGS BAMs (GRCh38, novoalign)
BAM_BASE="${GIAB_BASE}/data/AshkenazimTrio"
HG002_BAM="${BAM_BASE}/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
HG003_BAM="${BAM_BASE}/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.GRCh38.2x250.bam"
HG004_BAM="${BAM_BASE}/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.GRCh38.2x250.bam"

# GIAB v4.2.1 benchmark VCFs (GRCh38) – queried over HTTPS, never downloaded
# NOTE: Use explicit version directory (NISTv4.2.1), NOT the 'latest' symlink,
# because 'latest' may point to a newer release that uses different filenames.
BENCH_BASE="${GIAB_BASE}/release/AshkenazimTrio"
HG002_VCF="${BENCH_BASE}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG003_VCF="${BENCH_BASE}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG004_VCF="${BENCH_BASE}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

# ---------------------------------------------------------------------------
# Discovery windows – small regions on different chromosomes to stream
# from HG002's benchmark VCF.  Each ~50 kb window typically contains
# 100–200 HG002 benchmark SNVs, of which a few are child-private.
# ---------------------------------------------------------------------------
DISCOVERY_WINDOWS=(
    # Round 1 – one window per chromosome
    "chr1:5000000-5050000"
    "chr2:10000000-10050000"
    "chr3:15000000-15050000"
    "chr4:20000000-20050000"
    "chr5:25000000-25050000"
    "chr6:30000000-30050000"
    "chr7:35000000-35050000"
    "chr8:40000000-40050000"
    "chr9:45000000-45050000"
    "chr10:50000000-50050000"
    "chr11:55000000-55050000"
    "chr12:60000000-60050000"
    "chr13:40000000-40050000"
    "chr14:50000000-50050000"
    "chr15:35000000-35050000"
    "chr16:20000000-20050000"
    "chr17:25000000-25050000"
    "chr18:30000000-30050000"
    "chr19:15000000-15050000"
    "chr20:10000000-10050000"
    "chr21:20000000-20050000"
    "chr22:25000000-25050000"
    # Round 2 – different offsets as fallback
    "chr1:100000000-100050000"
    "chr2:50000000-50050000"
    "chr3:80000000-80050000"
    "chr4:90000000-90050000"
    "chr5:70000000-70050000"
    "chr6:60000000-60050000"
    "chr7:80000000-80050000"
    "chr8:70000000-70050000"
    "chr9:30000000-30050000"
    "chr10:80000000-80050000"
    "chr11:30000000-30050000"
    "chr12:40000000-40050000"
    "chr13:60000000-60050000"
    "chr14:30000000-30050000"
    "chr15:50000000-50050000"
    "chr16:40000000-40050000"
    "chr17:45000000-45050000"
    "chr18:50000000-50050000"
    "chr19:30000000-30050000"
    "chr20:30000000-30050000"
    "chr21:25000000-25050000"
    "chr22:30000000-30050000"
    # Round 3 – further fallback at yet another offset
    "chr1:150000000-150050000"
    "chr2:120000000-120050000"
    "chr3:120000000-120050000"
    "chr4:130000000-130050000"
    "chr5:110000000-110050000"
    "chr6:100000000-100050000"
    "chr7:100000000-100050000"
    "chr8:100000000-100050000"
    "chr9:80000000-80050000"
    "chr10:100000000-100050000"
    "chr11:80000000-80050000"
    "chr12:90000000-90050000"
    "chr13:70000000-70050000"
    "chr14:60000000-60050000"
    "chr15:60000000-60050000"
    "chr16:50000000-50050000"
    "chr17:50000000-50050000"
    "chr18:40000000-40050000"
    "chr19:20000000-20050000"
    "chr20:40000000-40050000"
    "chr21:30000000-30050000"
    "chr22:35000000-35050000"
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
die()  { log "ERROR: $*"; exit 1; }

usage() {
    sed -n '/^# Usage/,/^# Data sources/p' "$0" | sed '$d' | sed 's/^# \?//'
    exit 1
}

check_tool() {
    command -v "$1" >/dev/null 2>&1 || die "$1 is required but not found on PATH"
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while getopts ":o:n:p:h" opt; do
    case ${opt} in
        o) OUTPUT_DIR="${OPTARG}" ;;
        n) NUM_VARIANTS="${OPTARG}" ;;
        p) PADDING="${OPTARG}" ;;
        h) usage ;;
        :) die "Option -${OPTARG} requires an argument" ;;
        *) die "Unknown option -${OPTARG}" ;;
    esac
done

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
check_tool samtools
check_tool bcftools

log "=== GIAB HG002 trio – child-private variant test-data extractor ==="
log "Output directory: ${OUTPUT_DIR}"
log "Num variants    : ${NUM_VARIANTS}"
log "Padding (bp)    : ${PADDING}"

# Quick connectivity check: fetch the VCF header to confirm remote access works
log "Checking GIAB FTP connectivity ..."
if ! bcftools view -h "${HG002_VCF}" >/dev/null 2>&1; then
    die "Cannot access HG002 benchmark VCF at: ${HG002_VCF}
         Check network connectivity and that the URL is still valid.
         Browse: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/"
fi
log "  OK – remote VCF accessible"

# ---------------------------------------------------------------------------
# Set up directories
# ---------------------------------------------------------------------------
WORK_DIR="$(mktemp -d -t giab_dnm_XXXXXX)"
trap 'rm -rf "${WORK_DIR}"' EXIT
log "Working directory: ${WORK_DIR}"

mkdir -p "${OUTPUT_DIR}"

# ---------------------------------------------------------------------------
# Step 1 – Discover child-private SNVs from GIAB benchmark VCFs
# ---------------------------------------------------------------------------
log "Step 1: Discovering child-private SNVs from GIAB v4.2.1 benchmark VCFs ..."
log "  (streaming small windows via HTTPS – no bulk downloads)"
log "  Each window = 3 HTTP range requests (child + father + mother)"

VERIFIED=0

for window in "${DISCOVERY_WINDOWS[@]}"; do
    [[ "${VERIFIED}" -ge "${NUM_VARIANTS}" ]] && break

    log "  Scanning window ${window} ..."

    # Fetch all three VCFs for this window in one batch (3 HTTP requests).
    # Extract only CHROM and POS to build position sets for comparison.
    bcftools view -H -v snps -r "${window}" "${HG002_VCF}" 2>/dev/null \
        > "${WORK_DIR}/hg002_window.tsv" || continue
    bcftools view -H -v snps -r "${window}" "${HG003_VCF}" 2>/dev/null \
        | awk -F'\t' '{print $1"\t"$2}' > "${WORK_DIR}/hg003_positions.tsv" || true
    bcftools view -H -v snps -r "${window}" "${HG004_VCF}" 2>/dev/null \
        | awk -F'\t' '{print $1"\t"$2}' > "${WORK_DIR}/hg004_positions.tsv" || true

    num_child=$(wc -l < "${WORK_DIR}/hg002_window.tsv" | tr -d ' ')
    num_father=$(wc -l < "${WORK_DIR}/hg003_positions.tsv" | tr -d ' ')
    num_mother=$(wc -l < "${WORK_DIR}/hg004_positions.tsv" | tr -d ' ')
    log "    HG002: ${num_child} SNVs, HG003: ${num_father}, HG004: ${num_mother}"

    if [[ "${num_child}" -eq 0 ]]; then
        continue
    fi

    # Find child-private SNVs: in HG002 but not in either parent.
    # Compare locally – no additional HTTP requests needed.
    while IFS=$'\t' read -r chrom pos _ ref alt _rest; do
        [[ "${VERIFIED}" -ge "${NUM_VARIANTS}" ]] && break

        # Only consider biallelic SNVs (single base ref and alt)
        [[ ${#ref} -ne 1 || ${#alt} -ne 1 ]] && continue
        # Skip multi-allelic
        [[ "${alt}" == *","* ]] && continue

        # Check if this position exists in either parent (local grep)
        if grep -q "^${chrom}	${pos}$" "${WORK_DIR}/hg003_positions.tsv" 2>/dev/null; then
            continue
        fi
        if grep -q "^${chrom}	${pos}$" "${WORK_DIR}/hg004_positions.tsv" 2>/dev/null; then
            continue
        fi

        log "    FOUND child-private SNV: ${chrom}:${pos} ${ref}>${alt}"
        printf "%s\t%s\t%s\t%s\n" "${chrom}" "${pos}" "${ref}" "${alt}" \
            >> "${WORK_DIR}/dnm_verified.tsv"
        # Save the full VCF line for later extraction
        grep "^${chrom}	${pos}	" "${WORK_DIR}/hg002_window.tsv" | head -1 \
            >> "${WORK_DIR}/dnm_vcflines.tsv"
        VERIFIED=$((VERIFIED + 1))
    done < "${WORK_DIR}/hg002_window.tsv"

    log "    ${VERIFIED}/${NUM_VARIANTS} candidates found so far"
done

if [[ "${VERIFIED}" -eq 0 ]]; then
    die "No child-private SNVs found. Check network connectivity to GIAB FTP."
fi

log "  Discovered ${VERIFIED} child-private SNVs"

# ---------------------------------------------------------------------------
# Step 2 – Create extraction regions (±PADDING around each candidate)
# ---------------------------------------------------------------------------
log "Step 2: Creating extraction regions (±${PADDING} bp) ..."

awk -v pad="${PADDING}" 'BEGIN{OFS="\t"} {
    chrom = $1
    pos   = $2         # 1-based VCF POS
    start = pos - pad - 1   # 0-based BED start
    if (start < 0) start = 0
    end   = pos + pad  # BED end (exclusive)
    print chrom, start, end
}' "${WORK_DIR}/dnm_verified.tsv" > "${WORK_DIR}/regions.bed"

log "  Regions:"
while IFS=$'\t' read -r chrom start end; do
    log "    ${chrom}:${start}-${end}"
done < "${WORK_DIR}/regions.bed"

# Build samtools region strings (space-separated)
REGIONS=$(awk 'BEGIN{OFS=""} {printf "%s%s:%d-%d", (NR>1 ? " " : ""), $1, $2+1, $3}' "${WORK_DIR}/regions.bed")

# ---------------------------------------------------------------------------
# Step 3 – Extract BAM regions for each trio member (via HTTPS)
# ---------------------------------------------------------------------------
extract_bam_regions() {
    local url="$1" label="$2" outbam="$3"
    log "  Extracting ${label} reads (HTTPS random access) ..."
    # shellcheck disable=SC2086
    samtools view -b -h -o "${WORK_DIR}/${label}_unsorted.bam" \
        "${url}" ${REGIONS} \
        || die "Failed to extract ${label} BAM regions. Ensure samtools has HTTP support."

    samtools sort -o "${outbam}" "${WORK_DIR}/${label}_unsorted.bam"
    samtools index "${outbam}"
    local count
    count=$(samtools view -c "${outbam}")
    log "    ${count} reads -> ${outbam}"
}

log "Step 3: Extracting BAM regions from GIAB Illumina 2x250bp WGS (HTTPS) ..."
extract_bam_regions "${HG002_BAM}" "HG002_child"  "${OUTPUT_DIR}/HG002_child.bam"
extract_bam_regions "${HG003_BAM}" "HG003_father" "${OUTPUT_DIR}/HG003_father.bam"
extract_bam_regions "${HG004_BAM}" "HG004_mother" "${OUTPUT_DIR}/HG004_mother.bam"

# ---------------------------------------------------------------------------
# Step 4 – Create candidate VCF from GIAB benchmark records
# ---------------------------------------------------------------------------
log "Step 4: Creating candidate VCF ..."

# Grab the full VCF header from HG002's benchmark VCF (includes FORMAT
# definitions, contig lines, FILTER, INFO, and the sample column).
bcftools view -h "${HG002_VCF}" 2>/dev/null > "${WORK_DIR}/giab_header.txt" \
    || die "Failed to fetch HG002 VCF header"

# Sort the saved full-record lines by chromosome + position
sort -k1,1V -k2,2n "${WORK_DIR}/dnm_vcflines.tsv" > "${WORK_DIR}/dnm_vcflines_sorted.tsv"

# Combine GIAB header + sorted records into a valid VCF
{
    cat "${WORK_DIR}/giab_header.txt"
    cat "${WORK_DIR}/dnm_vcflines_sorted.tsv"
} > "${WORK_DIR}/candidates.vcf"

bcftools sort "${WORK_DIR}/candidates.vcf" -Oz -o "${OUTPUT_DIR}/candidates.vcf.gz"
bcftools index -t "${OUTPUT_DIR}/candidates.vcf.gz"
log "  Candidates VCF: ${OUTPUT_DIR}/candidates.vcf.gz"

# Log the sample name from the VCF
SAMPLE_NAME=$(bcftools query -l "${OUTPUT_DIR}/candidates.vcf.gz" | head -1)
log "  Sample: ${SAMPLE_NAME}"

# ---------------------------------------------------------------------------
# Step 5 – Write manifest
# ---------------------------------------------------------------------------
log "Step 5: Writing manifest ..."

cat > "${OUTPUT_DIR}/README.md" << 'MANIFEST_EOF'
# GIAB HG002 Trio – Child-Private Variant Test Data

Small BAM slices around child-private SNVs from the GIAB Ashkenazi trio,
for integration testing of `kmer_denovo_filter`.

Child-private variants are SNVs present in the HG002 (child) GIAB v4.2.1
benchmark VCF but absent from both HG003 (father) and HG004 (mother)
benchmark VCFs.  These serve as realistic candidates for testing de novo
mutation filtering workflows.

## Samples

| File              | Sample | Role   | GIAB ID   |
|-------------------|--------|--------|-----------|
| HG002_child.bam   | HG002  | Child  | NA24385   |
| HG003_father.bam  | HG003  | Father | NA24149   |
| HG004_mother.bam  | HG004  | Mother | NA24143   |

## Data sources

- **BAMs**: NIST Illumina 2×250 bp WGS, aligned with novoalign to GRCh38
  (queried via HTTPS random access – never downloaded in full)
- **Variant discovery**: SNVs streamed from HG002 GIAB v4.2.1 benchmark VCF
  across small chromosomal windows, filtered for absence in both parents
- **Verification**: Each position verified against HG003 and HG004 GIAB
  v4.2.1 benchmark VCFs over HTTPS

## Files

- `HG002_child.bam` / `.bai` – Child reads around variant sites
- `HG003_father.bam` / `.bai` – Father reads around variant sites
- `HG004_mother.bam` / `.bai` – Mother reads around variant sites
- `candidates.vcf.gz` / `.tbi` – VCF with child-private SNVs

## Usage with kmer-denovo

```bash
kmer-denovo \
    --child   tests/data/giab/HG002_child.bam   \
    --father  tests/data/giab/HG003_father.bam  \
    --mother  tests/data/giab/HG004_mother.bam  \
    --vcf     tests/data/giab/candidates.vcf.gz  \
    --output  output.vcf
```
MANIFEST_EOF

log "  Manifest: ${OUTPUT_DIR}/README.md"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log ""
log "=== Done! ==="
log "Output directory: ${OUTPUT_DIR}"
log "Files:"
ls -lh "${OUTPUT_DIR}/" | tail -n +2 | while read -r line; do
    log "  ${line}"
done
log ""
log "Verified DNM candidates:"
awk '{printf "  %s:%s  %s>%s\n", $1, $2, $3, $4}' "${WORK_DIR}/dnm_verified.tsv" >&2
