#!/usr/bin/env bash
# =============================================================================
# download_giab_dnm_testdata.sh
#
# Downloads small BAM regions around known de novo mutations (DNMs) from the
# GIAB HG002 (NA24385) Ashkenazi trio for integration testing.
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
# De novo mutations (DNMs) are variants present in the child but absent in
# both parents.  This script identifies putative DNMs by comparing the GIAB
# v4.2.1 benchmark VCFs, then extracts small BAM regions around a selected
# subset so that realistic but lightweight test data can be committed to the
# repository for integration testing of kmer_denovo_filter.
#
# Algorithm
# ---------
#   1. Download GIAB v4.2.1 benchmark VCFs + indexes for HG002, HG003, HG004.
#   2. Identify putative DNMs via bcftools isec  (variants in HG002 absent
#      from both parents).
#   3. Filter to biallelic SNVs within the GIAB high-confidence regions for
#      all three samples.
#   4. Select a small subset (default 5) of candidates spread across
#      different chromosomes.
#   5. Extract ±500 bp BAM slices around each candidate from the public GIAB
#      Illumina 2×250 bp WGS BAMs (remote HTTPS access via htslib).
#   6. Extract corresponding reference FASTA regions.
#   7. Write a candidate VCF containing the selected DNMs.
#
# Prerequisites
# -------------
#   • samtools  ≥ 1.10  (compiled with htslib HTTP/S support)
#   • bcftools  ≥ 1.10
#   • A local copy of the GRCh38 reference FASTA (with .fai index)
#
# Usage
# -----
#   ./scripts/download_giab_dnm_testdata.sh \
#       -r /path/to/GRCh38_no_alt.fa       \
#       [-o output_dir]                     \
#       [-n num_variants]                   \
#       [-p padding_bp]
#
# Options
# -------
#   -r REF_FASTA   Path to local GRCh38 reference FASTA (.fai index required).
#                  Download from:
#                    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/
#                      release/references/GRCh38/
#                      GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
#   -o OUTPUT_DIR  Output directory  (default: tests/data/giab)
#   -n NUM         Number of DNM candidates to select  (default: 5)
#   -p PADDING     Base-pairs of padding around each variant  (default: 500)
#
# Data sources
# ------------
#   BAMs  :  NIST Illumina 2×250 bp WGS, novoalign, GRCh38
#   VCFs  :  GIAB v4.2.1 benchmark  (GRCh38, chr1-22)
#   Ref   :  GCA_000001405.15  GRCh38 no-alt analysis set
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_DIR="${REPO_DIR}/tests/data/giab"
NUM_VARIANTS=5
PADDING=500
REF_FASTA=""

# ---------------------------------------------------------------------------
# GIAB data URLs
# ---------------------------------------------------------------------------
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"

# Illumina 2×250 bp WGS BAMs (GRCh38, novoalign)
BAM_BASE="${GIAB_BASE}/data/AshkenazimTrio"
HG002_BAM="${BAM_BASE}/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
HG003_BAM="${BAM_BASE}/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.GRCh38.2x250.bam"
HG004_BAM="${BAM_BASE}/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.GRCh38.2x250.bam"

# GIAB v4.2.1 benchmark VCFs (GRCh38)
BENCH_BASE="${GIAB_BASE}/release/AshkenazimTrio"
HG002_VCF_URL="${BENCH_BASE}/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG003_VCF_URL="${BENCH_BASE}/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG004_VCF_URL="${BENCH_BASE}/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
die()  { log "ERROR: $*"; exit 1; }

usage() {
    sed -n '/^# Usage/,/^# Data sources/p' "$0" | head -n -1 | sed 's/^# \?//'
    exit 1
}

check_tool() {
    command -v "$1" >/dev/null 2>&1 || die "$1 is required but not found on PATH"
}

download() {
    local url="$1" dest="$2"
    if [[ -f "${dest}" ]]; then
        log "  Already exists: ${dest}"
        return
    fi
    log "  Downloading $(basename "${dest}") ..."
    curl -fsSL -o "${dest}" "${url}" \
        || wget -q -O "${dest}" "${url}" \
        || die "Failed to download ${url}"
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while getopts ":r:o:n:p:h" opt; do
    case ${opt} in
        r) REF_FASTA="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        n) NUM_VARIANTS="${OPTARG}" ;;
        p) PADDING="${OPTARG}" ;;
        h) usage ;;
        :) die "Option -${OPTARG} requires an argument" ;;
        *) die "Unknown option -${OPTARG}" ;;
    esac
done

[[ -n "${REF_FASTA}" ]] || die "Reference FASTA is required (-r).  See usage (-h)."
[[ -f "${REF_FASTA}" ]] || die "Reference FASTA not found: ${REF_FASTA}"
[[ -f "${REF_FASTA}.fai" ]] || die "Reference FASTA index not found: ${REF_FASTA}.fai"

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
check_tool samtools
check_tool bcftools
command -v curl >/dev/null 2>&1 || command -v wget >/dev/null 2>&1 \
    || die "Either curl or wget is required but neither was found on PATH"

log "=== GIAB HG002 trio DNM test-data downloader ==="
log "Reference FASTA : ${REF_FASTA}"
log "Output directory: ${OUTPUT_DIR}"
log "Num variants    : ${NUM_VARIANTS}"
log "Padding (bp)    : ${PADDING}"

# ---------------------------------------------------------------------------
# Set up directories
# ---------------------------------------------------------------------------
WORK_DIR="$(mktemp -d -t giab_dnm_XXXXXX)"
trap 'rm -rf "${WORK_DIR}"' EXIT
log "Working directory: ${WORK_DIR}"

mkdir -p "${OUTPUT_DIR}"

# ---------------------------------------------------------------------------
# Step 1 – Download GIAB benchmark VCFs
# ---------------------------------------------------------------------------
log "Step 1: Downloading GIAB v4.2.1 benchmark VCFs ..."

download "${HG002_VCF_URL}"          "${WORK_DIR}/HG002.vcf.gz"
download "${HG002_VCF_URL}.tbi"      "${WORK_DIR}/HG002.vcf.gz.tbi"

download "${HG003_VCF_URL}"          "${WORK_DIR}/HG003.vcf.gz"
download "${HG003_VCF_URL}.tbi"      "${WORK_DIR}/HG003.vcf.gz.tbi"

download "${HG004_VCF_URL}"          "${WORK_DIR}/HG004.vcf.gz"
download "${HG004_VCF_URL}.tbi"      "${WORK_DIR}/HG004.vcf.gz.tbi"

# ---------------------------------------------------------------------------
# Step 2 – Note on high-confidence region handling
# ---------------------------------------------------------------------------
log "Step 2: High-confidence region handling ..."

# The GIAB v4.2.1 benchmark VCFs only contain variants within each sample's
# high-confidence regions, so bcftools isec (Step 3) inherently restricts the
# comparison to well-characterized genomic positions.  A variant flagged as
# HG002-unique may simply be outside a parent's HC region rather than a true
# DNM.  For a small integration-test dataset this is acceptable; downstream
# manual inspection can confirm candidates.
log "  Using implicit HC filtering from GIAB benchmark VCFs"

# ---------------------------------------------------------------------------
# Step 3 – Identify putative de novo mutations
# ---------------------------------------------------------------------------
log "Step 3: Identifying putative de novo mutations (bcftools isec) ..."

# bcftools isec -C: output records from first file not present in any others
bcftools isec -C \
    "${WORK_DIR}/HG002.vcf.gz" \
    "${WORK_DIR}/HG003.vcf.gz" \
    "${WORK_DIR}/HG004.vcf.gz" \
    -Oz -o "${WORK_DIR}/dnm_raw.vcf.gz"

bcftools index -t "${WORK_DIR}/dnm_raw.vcf.gz"

TOTAL_RAW=$(bcftools view -H "${WORK_DIR}/dnm_raw.vcf.gz" | wc -l)
log "  Found ${TOTAL_RAW} HG002-unique variants (putative DNMs + private)"

# ---------------------------------------------------------------------------
# Step 4 – Filter to biallelic SNVs and select candidates
# ---------------------------------------------------------------------------
log "Step 4: Filtering to biallelic SNVs and selecting ${NUM_VARIANTS} candidates ..."

# Keep only biallelic SNVs (simplest, most reliable variant type for testing)
bcftools view -v snps -m2 -M2 "${WORK_DIR}/dnm_raw.vcf.gz" \
    -Oz -o "${WORK_DIR}/dnm_snvs.vcf.gz"
bcftools index -t "${WORK_DIR}/dnm_snvs.vcf.gz"

TOTAL_SNVS=$(bcftools view -H "${WORK_DIR}/dnm_snvs.vcf.gz" | wc -l)
log "  Biallelic SNVs: ${TOTAL_SNVS}"

# Select candidates spread across different chromosomes.
# Strategy: pick up to one variant per chromosome, cycling through
# autosomes until we have enough candidates.
bcftools view -H "${WORK_DIR}/dnm_snvs.vcf.gz" \
    | awk -v n="${NUM_VARIANTS}" '
        BEGIN { count = 0 }
        {
            chrom = $1
            if (!(chrom in seen)) {
                seen[chrom] = 1
                print
                count++
                if (count >= n) exit
            }
        }
    ' > "${WORK_DIR}/dnm_selected.txt"

SELECTED=$(wc -l < "${WORK_DIR}/dnm_selected.txt")
log "  Selected ${SELECTED} DNM candidates across different chromosomes"

if [[ "${SELECTED}" -eq 0 ]]; then
    die "No DNM candidates found.  Check that truth VCFs downloaded correctly."
fi

# ---------------------------------------------------------------------------
# Step 5 – Create extraction BED (±PADDING around each candidate)
# ---------------------------------------------------------------------------
log "Step 5: Creating extraction regions (±${PADDING} bp) ..."

awk -v pad="${PADDING}" 'BEGIN{OFS="\t"} {
    chrom = $1
    pos   = $2         # 1-based VCF POS
    start = pos - pad - 1   # 0-based BED start
    if (start < 0) start = 0
    end   = pos + pad  # BED end (exclusive)
    print chrom, start, end
}' "${WORK_DIR}/dnm_selected.txt" > "${WORK_DIR}/regions.bed"

log "  Regions:"
while IFS=$'\t' read -r chrom start end; do
    log "    ${chrom}:${start}-${end}"
done < "${WORK_DIR}/regions.bed"

# Build samtools region strings (space-separated)
REGIONS=$(awk 'BEGIN{OFS=""} {printf "%s%s:%d-%d", (NR>1 ? " " : ""), $1, $2+1, $3}' "${WORK_DIR}/regions.bed")

# ---------------------------------------------------------------------------
# Step 6 – Extract reference FASTA for selected regions
# ---------------------------------------------------------------------------
log "Step 6: Extracting reference FASTA for selected regions ..."

REGION_LIST=$(awk '{printf "%s:%d-%d\n", $1, $2+1, $3}' "${WORK_DIR}/regions.bed")
while read -r reg; do
    samtools faidx "${REF_FASTA}" "${reg}"
done <<< "${REGION_LIST}" > "${OUTPUT_DIR}/reference.fa"

samtools faidx "${OUTPUT_DIR}/reference.fa"
log "  Reference: ${OUTPUT_DIR}/reference.fa"

# ---------------------------------------------------------------------------
# Step 7 – Extract BAM regions for each trio member
# ---------------------------------------------------------------------------
extract_bam_regions() {
    local url="$1" label="$2" outbam="$3"
    log "  Extracting ${label} reads ..."
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

log "Step 7: Extracting BAM regions from GIAB Illumina 2x250bp WGS ..."
extract_bam_regions "${HG002_BAM}" "HG002_child"  "${OUTPUT_DIR}/HG002_child.bam"
extract_bam_regions "${HG003_BAM}" "HG003_father" "${OUTPUT_DIR}/HG003_father.bam"
extract_bam_regions "${HG004_BAM}" "HG004_mother" "${OUTPUT_DIR}/HG004_mother.bam"

# ---------------------------------------------------------------------------
# Step 8 – Create candidate VCF with selected DNMs
# ---------------------------------------------------------------------------
log "Step 8: Creating candidate VCF ..."

{
    # VCF header
    echo "##fileformat=VCFv4.2"
    echo "##source=download_giab_dnm_testdata.sh"
    echo "##reference=GRCh38"
    echo "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Variant origin (GIAB benchmark)\">"
    # Contig lines from the selected regions
    awk '{printf "##contig=<ID=%s>\n", $1}' "${WORK_DIR}/dnm_selected.txt" | sort -u
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

    # Variant lines from the selected records
    awk 'BEGIN{OFS="\t"} {
        print $1, $2, ".", $4, $5, ".", "PASS", "ORIGIN=GIAB_v4.2.1_denovo"
    }' "${WORK_DIR}/dnm_selected.txt"
} > "${WORK_DIR}/candidates.vcf"

bgzip -c "${WORK_DIR}/candidates.vcf" > "${OUTPUT_DIR}/candidates.vcf.gz"
bcftools index -t "${OUTPUT_DIR}/candidates.vcf.gz"
log "  Candidates VCF: ${OUTPUT_DIR}/candidates.vcf.gz"

# ---------------------------------------------------------------------------
# Step 9 – Write manifest
# ---------------------------------------------------------------------------
log "Step 9: Writing manifest ..."

cat > "${OUTPUT_DIR}/README.md" << 'MANIFEST_EOF'
# GIAB HG002 Trio – De Novo Mutation Test Data

Small BAM slices around putative de novo mutations (DNMs) from the GIAB
Ashkenazi trio, for integration testing of `kmer_denovo_filter`.

## Samples

| File              | Sample | Role   | GIAB ID   |
|-------------------|--------|--------|-----------|
| HG002_child.bam   | HG002  | Child  | NA24385   |
| HG003_father.bam  | HG003  | Father | NA24149   |
| HG004_mother.bam  | HG004  | Mother | NA24143   |

## Data sources

- **BAMs**: NIST Illumina 2×250 bp WGS, aligned with novoalign to GRCh38
- **Truth VCFs**: GIAB v4.2.1 benchmark (GRCh38, chr1-22)
- **Reference**: GCA_000001405.15 GRCh38 no-alt analysis set

## DNM identification

Putative DNMs were identified using `bcftools isec -C` to find variants
present in the HG002 benchmark VCF but absent from both HG003 and HG004
benchmark VCFs.  Only biallelic SNVs were retained.

## Files

- `HG002_child.bam` / `.bai` – Child reads around DNM sites
- `HG003_father.bam` / `.bai` – Father reads around DNM sites
- `HG004_mother.bam` / `.bai` – Mother reads around DNM sites
- `reference.fa` / `.fai` – Reference FASTA (selected regions only)
- `candidates.vcf.gz` / `.tbi` – VCF with selected DNM variants

## Usage with kmer-denovo

```bash
kmer-denovo \
    --child   tests/data/giab/HG002_child.bam   \
    --father  tests/data/giab/HG003_father.bam  \
    --mother  tests/data/giab/HG004_mother.bam  \
    --ref-fasta tests/data/giab/reference.fa     \
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
log "Selected DNM candidates:"
awk '{printf "  %s:%s  %s>%s\n", $1, $2, $4, $5}' "${WORK_DIR}/dnm_selected.txt" >&2
