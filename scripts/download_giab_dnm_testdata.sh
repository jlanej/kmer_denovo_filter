#!/usr/bin/env bash
# =============================================================================
# download_giab_dnm_testdata.sh
#
# Extracts small BAM regions around known de novo mutations (DNMs) from the
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
# both parents.  This script uses a curated list of known DNM positions from
# the GIAB HG002 mosaic/de novo benchmark (usnistgov/giab-HG002-mosaic-
# benchmark), verifies them against the GIAB v4.2.1 truth VCFs over HTTPS,
# then extracts small BAM regions via samtools random access.
#
# NO large file downloads are performed.  All VCF and BAM queries use
# htslib HTTPS random access against the public GIAB FTP, transferring only
# a few MB of data in total.
#
# Algorithm
# ---------
#   1. Start with a curated list of published DNM positions (SNVs) sourced
#      from the GIAB HG002 mosaic benchmark (Strelka2 somatic caller,
#      IGV-curated, see usnistgov/giab-HG002-mosaic-benchmark).
#   2. For each candidate, verify via bcftools over HTTPS that the variant
#      is present in the HG002 benchmark VCF and absent from both parents.
#   3. Select a subset of verified candidates (default 5).
#   4. Extract ±500 bp BAM slices around each candidate from the public GIAB
#      Illumina 2×250 bp WGS BAMs (remote HTTPS via htslib).
#   5. Extract corresponding reference FASTA regions.
#   6. Write a candidate VCF containing the selected DNMs.
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
#   DNMs  :  usnistgov/giab-HG002-mosaic-benchmark (Strelka2, IGV-curated)
#   BAMs  :  NIST Illumina 2×250 bp WGS, novoalign, GRCh38
#   VCFs  :  GIAB v4.2.1 benchmark  (GRCh38, chr1-22)  [queried over HTTPS]
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
# GIAB data URLs (all accessed via HTTPS random access – no bulk downloads)
# ---------------------------------------------------------------------------
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"

# Illumina 2×250 bp WGS BAMs (GRCh38, novoalign)
BAM_BASE="${GIAB_BASE}/data/AshkenazimTrio"
HG002_BAM="${BAM_BASE}/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
HG003_BAM="${BAM_BASE}/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.GRCh38.2x250.bam"
HG004_BAM="${BAM_BASE}/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.GRCh38.2x250.bam"

# GIAB v4.2.1 benchmark VCFs (GRCh38) – queried over HTTPS, never downloaded
BENCH_BASE="${GIAB_BASE}/release/AshkenazimTrio"
HG002_VCF="${BENCH_BASE}/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG003_VCF="${BENCH_BASE}/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG004_VCF="${BENCH_BASE}/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

# ---------------------------------------------------------------------------
# Curated DNM candidate positions (GRCh38)
#
# Source: usnistgov/giab-HG002-mosaic-benchmark data/igv_curated_variants_RM.tsv
# These are SNVs identified by Strelka2 somatic calling (HG002 as tumor,
# parents as normal), then IGV-curated.  Selected "candidate" variants with
# high SomaticEVS scores (>17) across different chromosomes.
#
# Format: CHROM  POS  REF  ALT
# ---------------------------------------------------------------------------
CANDIDATE_POSITIONS=(
    "chr1	79144396	T	C"
    "chr2	4482351	G	A"
    "chr3	21601448	G	C"
    "chr4	58670824	G	A"
    "chr5	41022783	A	G"
    "chr6	27836914	A	G"
    "chr7	550099	G	A"
    "chr8	16522616	C	A"
    "chr9	95143456	C	A"
    "chr10	18475692	A	T"
    "chr11	97685827	C	A"
    "chr12	112675498	T	G"
    "chr13	100629230	G	T"
    "chr14	48779433	A	G"
    "chr15	32467857	A	G"
    "chr16	11815294	T	C"
    "chr17	7131851	G	A"
    "chr18	25733623	C	T"
    "chr19	18888816	C	T"
    "chr20	12956555	T	C"
    "chr21	30498752	T	C"
    "chr22	38620848	C	T"
)

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

log "=== GIAB HG002 trio DNM test-data extractor ==="
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
# Step 1 – Verify DNM candidates against GIAB truth VCFs over HTTPS
# ---------------------------------------------------------------------------
log "Step 1: Verifying DNM candidates against GIAB v4.2.1 benchmark VCFs ..."
log "  (querying remote VCFs via HTTPS – no bulk downloads)"

VERIFIED=0

for entry in "${CANDIDATE_POSITIONS[@]}"; do
    [[ "${VERIFIED}" -ge "${NUM_VARIANTS}" ]] && break

    chrom=$(echo "${entry}" | cut -f1)
    pos=$(echo "${entry}" | cut -f2)
    ref=$(echo "${entry}" | cut -f3)
    alt=$(echo "${entry}" | cut -f4)
    region="${chrom}:${pos}-${pos}"

    log "  Checking ${chrom}:${pos} ${ref}>${alt} ..."

    # Verify variant is in HG002 benchmark VCF
    in_child=$(bcftools view -H -r "${region}" "${HG002_VCF}" 2>/dev/null | wc -l)
    if [[ "${in_child}" -eq 0 ]]; then
        log "    SKIP – not in HG002 benchmark VCF"
        continue
    fi

    # Verify variant is absent from both parents
    in_father=$(bcftools view -H -r "${region}" "${HG003_VCF}" 2>/dev/null | wc -l)
    in_mother=$(bcftools view -H -r "${region}" "${HG004_VCF}" 2>/dev/null | wc -l)
    if [[ "${in_father}" -gt 0 || "${in_mother}" -gt 0 ]]; then
        log "    SKIP – present in parent(s) (father=${in_father}, mother=${in_mother})"
        continue
    fi

    log "    VERIFIED – child-only variant"
    printf "%s\t%s\t%s\t%s\n" "${chrom}" "${pos}" "${ref}" "${alt}" \
        >> "${WORK_DIR}/dnm_verified.tsv"
    VERIFIED=$((VERIFIED + 1))
done

if [[ "${VERIFIED}" -eq 0 ]]; then
    die "No DNM candidates could be verified. Check network connectivity."
fi

log "  Verified ${VERIFIED} DNM candidates"

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
# Step 3 – Extract reference FASTA for selected regions
# ---------------------------------------------------------------------------
log "Step 3: Extracting reference FASTA for selected regions ..."

REGION_LIST=$(awk '{printf "%s:%d-%d\n", $1, $2+1, $3}' "${WORK_DIR}/regions.bed")
while read -r reg; do
    samtools faidx "${REF_FASTA}" "${reg}"
done <<< "${REGION_LIST}" > "${OUTPUT_DIR}/reference.fa"

samtools faidx "${OUTPUT_DIR}/reference.fa"
log "  Reference: ${OUTPUT_DIR}/reference.fa"

# ---------------------------------------------------------------------------
# Step 4 – Extract BAM regions for each trio member (via HTTPS)
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

log "Step 4: Extracting BAM regions from GIAB Illumina 2x250bp WGS (HTTPS) ..."
extract_bam_regions "${HG002_BAM}" "HG002_child"  "${OUTPUT_DIR}/HG002_child.bam"
extract_bam_regions "${HG003_BAM}" "HG003_father" "${OUTPUT_DIR}/HG003_father.bam"
extract_bam_regions "${HG004_BAM}" "HG004_mother" "${OUTPUT_DIR}/HG004_mother.bam"

# ---------------------------------------------------------------------------
# Step 5 – Create candidate VCF with verified DNMs
# ---------------------------------------------------------------------------
log "Step 5: Creating candidate VCF ..."

{
    echo "##fileformat=VCFv4.2"
    echo "##source=download_giab_dnm_testdata.sh"
    echo "##reference=GRCh38"
    echo "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Variant origin\">"
    awk '{printf "##contig=<ID=%s>\n", $1}' "${WORK_DIR}/dnm_verified.tsv" | sort -u
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    awk 'BEGIN{OFS="\t"} {
        print $1, $2, ".", $3, $4, ".", "PASS", "ORIGIN=GIAB_mosaic_benchmark_denovo"
    }' "${WORK_DIR}/dnm_verified.tsv"
} > "${WORK_DIR}/candidates.vcf"

bgzip -c "${WORK_DIR}/candidates.vcf" > "${OUTPUT_DIR}/candidates.vcf.gz"
bcftools index -t "${OUTPUT_DIR}/candidates.vcf.gz"
log "  Candidates VCF: ${OUTPUT_DIR}/candidates.vcf.gz"

# ---------------------------------------------------------------------------
# Step 6 – Write manifest
# ---------------------------------------------------------------------------
log "Step 6: Writing manifest ..."

cat > "${OUTPUT_DIR}/README.md" << 'MANIFEST_EOF'
# GIAB HG002 Trio – De Novo Mutation Test Data

Small BAM slices around known de novo mutations (DNMs) from the GIAB
Ashkenazi trio, for integration testing of `kmer_denovo_filter`.

## Samples

| File              | Sample | Role   | GIAB ID   |
|-------------------|--------|--------|-----------|
| HG002_child.bam   | HG002  | Child  | NA24385   |
| HG003_father.bam  | HG003  | Father | NA24149   |
| HG004_mother.bam  | HG004  | Mother | NA24143   |

## Data sources

- **BAMs**: NIST Illumina 2×250 bp WGS, aligned with novoalign to GRCh38
  (queried via HTTPS random access – never downloaded in full)
- **DNM positions**: Curated from usnistgov/giab-HG002-mosaic-benchmark
  (Strelka2 somatic caller, IGV-curated candidate SNVs)
- **Verification**: Each position verified against GIAB v4.2.1 benchmark VCFs
  over HTTPS (present in HG002, absent from HG003 and HG004)
- **Reference**: GCA_000001405.15 GRCh38 no-alt analysis set

## Files

- `HG002_child.bam` / `.bai` – Child reads around DNM sites
- `HG003_father.bam` / `.bai` – Father reads around DNM sites
- `HG004_mother.bam` / `.bai` – Mother reads around DNM sites
- `reference.fa` / `.fai` – Reference FASTA (selected regions only)
- `candidates.vcf.gz` / `.tbi` – VCF with verified DNM variants

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
log "Verified DNM candidates:"
awk '{printf "  %s:%s  %s>%s\n", $1, $2, $3, $4}' "${WORK_DIR}/dnm_verified.tsv" >&2
