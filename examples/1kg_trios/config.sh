#!/usr/bin/env bash
# =============================================================================
# config.sh — shared configuration for the 1000 Genomes trio screen
# =============================================================================
#
# Every script in this directory sources this file. Each setting is a
# variable with a default; set it any of these ways, highest precedence first:
#
#   1. Command-line flags, on any script: --work-dir DIR, --<variable> VALUE
#      for any variable listed in CONFIG_VARIABLES below, spelled in lower
#      case with dashes (--results-dir, --ref-fasta, --n-trios,
#      --sbatch-extra-args "..."), --set NAME=VALUE, and --config FILE.
#      --help prints a script's usage and the effective configuration.
#   2. Environment variables exported before running a script:
#        export WORK_DIR=/scratch/$USER/1kg_trios
#   3. A site config file, sourced before the defaults are applied:
#      --config FILE, else $KMER_TRIOS_CONFIG, else config.local.sh beside
#      this file (gitignored; see config.local.example.sh). Write it with
#      the  VAR="${VAR:-value}"  form so flags and exported variables still
#      win.
#   4. The defaults below.
#
# Jobs the scripts submit get the effective configuration as a snapshot
# (--config $CONFIG_SNAPSHOT), so they see the same paths whether or not the
# scheduler exports the environment.
#
# Nothing here is specific to one cluster: paths default to a scratch tree
# under WORK_DIR, scheduler specifics go in SBATCH_EXTRA_ARGS, environment
# modules in MODULES_TO_LOAD, and every URL points at the public
# 1000 Genomes / ENA / AWS Open Data copies.
# =============================================================================

set -euo pipefail

# Directory holding this file (the scripts address each other through it,
# because a SLURM batch script runs from the spool directory).
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Every setting, in the order --help lists them. A flag --name-with-dashes
# is accepted for each; the snapshot a job receives holds all of them.
CONFIG_VARIABLES="WORK_DIR CRAM_DIR RESULTS_DIR COHORT_DIR LOG_DIR REF_DIR STATE_DIR TMP_ROOT CONFIG_SNAPSHOT
  CONTAINER_URI SIF_IMAGE CONTAINER_EXAMPLE_DIR MODULES_TO_LOAD
  REF_BASE_URL REF_NAME REF_FASTA KMER_SIZE REF_JF
  COLLECTION_URL INDEX_URL_2504 INDEX_URL_698 PED_URL INDEX_FILE_2504 INDEX_FILE_698 PED_FILE EXPECTED_COMPLETE_TRIOS
  TRIOS_TSV TRIO_FILES_TSV N_TRIOS TRIO_SEED TRIO_POPULATIONS TRIO_SUPERPOPULATIONS TRIO_CHILDREN TRIO_ALL_SIBLINGS
  SV_VCF_URL SV_VCF SV_COHORT_VCF USE_SV_CALLSET SV_MAX_AF SV_MIN_GQ SV_PASS_ONLY
  EXPECTED_FTP_PREFIX ENA_HTTPS_BASE S3_HTTPS_BASE S3_PREFIX_2504 S3_PREFIX_698 DOWNLOAD_CONNECTIONS
  TRIO_SLOTS MAX_STAGED_TRIOS TRIO_LIMIT KEEP_CRAMS
  DISCOVERY_CPUS DISCOVERY_MEM DISCOVERY_TIME MIN_CHILD_COUNT DISCOVERY_EXTRA_ARGS DISCOVERY_TMP_DIR KRAKEN2_DB
  REF_INDEX_CPUS REF_INDEX_MEM REF_INDEX_TIME REF_JF_HASH_SIZE
  SBATCH_EXTRA_ARGS JOB_PREFIX RECURRENT_MIN_TRIOS RECURRENT_SLACK_BP"

# ── Command-line flags ───────────────────────────────────────────────────────
# Parsed from the arguments the sourcing script passed (source config.sh "$@").
# Recognised flags become exported variables; everything else is handed back
# in CONFIG_REMAINING_ARGS for the script's own arguments.

config_variable_for_flag() {
  # --results-dir -> RESULTS_DIR, when that is a configuration variable
  local name
  name="$(printf '%s' "${1#--}" | tr 'a-z-' 'A-Z_')"
  case " ${CONFIG_VARIABLES} " in
    *" ${name} "*) echo "${name}" ;;
    *) return 1 ;;
  esac
}

apply_config_assignments() {
  local assignment name value
  for assignment in ${CONFIG_ASSIGNMENTS[@]+"${CONFIG_ASSIGNMENTS[@]}"}; do
    name="${assignment%%=*}"
    value="${assignment#*=}"
    if [[ ! "${name}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then
      echo "ERROR: --set expects NAME=VALUE, got '${assignment}'" >&2
      exit 2
    fi
    export "${name}=${value}"
  done
}

CONFIG_REMAINING_ARGS=()
CONFIG_ASSIGNMENTS=()
CONFIG_SHOW_HELP=0
CONFIG_SITE_FILE="${KMER_TRIOS_CONFIG:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      CONFIG_SHOW_HELP=1; shift ;;
    --set)
      [[ $# -ge 2 ]] || { echo "ERROR: --set needs NAME=VALUE" >&2; exit 2; }
      CONFIG_ASSIGNMENTS+=("$2"); shift 2 ;;
    --set=*)
      CONFIG_ASSIGNMENTS+=("${1#*=}"); shift ;;
    --config)
      [[ $# -ge 2 ]] || { echo "ERROR: --config needs a file" >&2; exit 2; }
      CONFIG_SITE_FILE="$2"; shift 2 ;;
    --config=*)
      CONFIG_SITE_FILE="${1#*=}"; shift ;;
    --*=*)
      if name="$(config_variable_for_flag "${1%%=*}")"; then
        CONFIG_ASSIGNMENTS+=("${name}=${1#*=}"); shift
      else
        CONFIG_REMAINING_ARGS+=("$1"); shift
      fi ;;
    --*)
      if name="$(config_variable_for_flag "$1")"; then
        [[ $# -ge 2 ]] || { echo "ERROR: $1 needs a value" >&2; exit 2; }
        CONFIG_ASSIGNMENTS+=("${name}=$2"); shift 2
      else
        CONFIG_REMAINING_ARGS+=("$1"); shift
      fi ;;
    *)
      CONFIG_REMAINING_ARGS+=("$1"); shift ;;
  esac
done
unset name
apply_config_assignments

# ── Site config file ─────────────────────────────────────────────────────────
if [[ -z "${CONFIG_SITE_FILE}" && -f "${CONFIG_DIR}/config.local.sh" ]]; then
  CONFIG_SITE_FILE="${CONFIG_DIR}/config.local.sh"
fi
if [[ -n "${CONFIG_SITE_FILE}" ]]; then
  if [[ ! -f "${CONFIG_SITE_FILE}" ]]; then
    echo "ERROR: config file not found: ${CONFIG_SITE_FILE}" >&2
    exit 2
  fi
  # shellcheck disable=SC1090
  source "${CONFIG_SITE_FILE}"
  apply_config_assignments   # flags win over the site file
fi

# ── Working tree ─────────────────────────────────────────────────────────────
WORK_DIR="${WORK_DIR:-/scratch/${USER:-$(id -un)}/1kg_trios}"
CRAM_DIR="${CRAM_DIR:-${WORK_DIR}/crams}"           # staged CRAMs (deleted per trio)
RESULTS_DIR="${RESULTS_DIR:-${WORK_DIR}/results}"   # one subdirectory per trio
COHORT_DIR="${COHORT_DIR:-${RESULTS_DIR}/cohort}"   # cross-trio summary
LOG_DIR="${LOG_DIR:-${WORK_DIR}/logs}"
REF_DIR="${REF_DIR:-${WORK_DIR}/reference}"
STATE_DIR="${STATE_DIR:-${WORK_DIR}/state}"         # manager bookkeeping
TMP_ROOT="${TMP_ROOT:-${WORK_DIR}/tmp}"             # kmer-discovery --tmp-dir root
# The effective configuration, written by the submitting script and passed
# to every job it submits as --config.
CONFIG_SNAPSHOT="${CONFIG_SNAPSHOT:-${STATE_DIR}/config.env}"

# ── Container and environment ────────────────────────────────────────────────
# The published image holds kmer-discovery, samtools, bcftools, jellyfish,
# aria2c and this directory's helper scripts (at CONTAINER_EXAMPLE_DIR).
CONTAINER_URI="${CONTAINER_URI:-docker://ghcr.io/jlanej/kmer_denovo_filter:latest}"
SIF_IMAGE="${SIF_IMAGE:-${WORK_DIR}/kmer_denovo_filter.sif}"  # or an existing .sif
CONTAINER_EXAMPLE_DIR="${CONTAINER_EXAMPLE_DIR:-/app/examples/1kg_trios}"
# Environment modules every script loads first (space-separated), e.g.
# "apptainer aria2". Empty = none.
MODULES_TO_LOAD="${MODULES_TO_LOAD:-}"

# ── Reference (GRCh38 as used by the NYGC alignments; needed to decode CRAM) ─
REF_BASE_URL="${REF_BASE_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome}"
REF_NAME="${REF_NAME:-GRCh38_full_analysis_set_plus_decoy_hla.fa}"
REF_FASTA="${REF_FASTA:-${REF_DIR}/${REF_NAME}}"   # an existing copy works too
KMER_SIZE="${KMER_SIZE:-31}"
# Jellyfish index of the reference, built once by 00b_build_ref_index.sh and
# shared by every trio job (kmer-discovery --ref-jf).
REF_JF="${REF_JF:-${REF_FASTA}.k${KMER_SIZE}.jf}"

# ── 1000 Genomes 30x data collection (NYGC, 3,202 samples, 602 trios) ────────
# Byrska-Bishop et al. (2022) Cell 185(18):3426-3440.
# https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
COLLECTION_URL="${COLLECTION_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage}"
INDEX_URL_2504="${INDEX_URL_2504:-${COLLECTION_URL}/1000G_2504_high_coverage.sequence.index}"
INDEX_URL_698="${INDEX_URL_698:-${COLLECTION_URL}/1000G_698_related_high_coverage.sequence.index}"
# Pedigree of all 3,202 samples: FamilyID SampleID FatherID MotherID Sex
# Population Superpopulation. 603 rows name both parents; 602 of those trios
# have CRAMs for all three members.
PED_URL="${PED_URL:-${COLLECTION_URL}/20130606_g1k_3202_samples_ped_population.txt}"
INDEX_FILE_2504="${INDEX_FILE_2504:-${WORK_DIR}/1000G_2504_high_coverage.sequence.index}"
INDEX_FILE_698="${INDEX_FILE_698:-${WORK_DIR}/1000G_698_related_high_coverage.sequence.index}"
PED_FILE="${PED_FILE:-${WORK_DIR}/1kg_3202_samples.ped}"
EXPECTED_COMPLETE_TRIOS="${EXPECTED_COMPLETE_TRIOS:-602}"

# ── Trio selection (00_setup.sh -> select_trios.py) ──────────────────────────
TRIOS_TSV="${TRIOS_TSV:-${WORK_DIR}/trios.tsv}"             # one row per trio
TRIO_FILES_TSV="${TRIO_FILES_TSV:-${WORK_DIR}/trio_files.tsv}"  # three rows per trio
N_TRIOS="${N_TRIOS:-10}"                 # 0 = every complete trio
TRIO_SEED="${TRIO_SEED-1}"               # seeded shuffle before taking N_TRIOS; export TRIO_SEED="" for sorted order
TRIO_POPULATIONS="${TRIO_POPULATIONS:-}"          # e.g. "GBR,CHS" (empty = any)
TRIO_SUPERPOPULATIONS="${TRIO_SUPERPOPULATIONS:-}" # e.g. "EUR,EAS"
TRIO_CHILDREN="${TRIO_CHILDREN:-}"       # explicit child IDs ("HG00405,NA19675") or a file
TRIO_ALL_SIBLINGS="${TRIO_ALL_SIBLINGS:-0}"  # 1 = keep every child of a family

# ── Published SV call set (optional screen of child-only SV genotypes) ───────
# The 3,202-sample Illumina SV call set (GATK-SV + svtools + Absinthe, freeze
# V3 with allele frequencies). Per trio, SVs genotyped non-reference in the
# child and homozygous reference in both parents become a --dnm-regions file,
# so each run reports whether discovery sees them. Most such calls in a
# population call set are parental genotyping misses, which is why the AF and
# GQ filters below are strict; see child_only_svs.py.
SV_VCF_URL="${SV_VCF_URL:-${COLLECTION_URL}/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz}"
SV_VCF="${SV_VCF:-${WORK_DIR}/sv_callset/$(basename "${SV_VCF_URL}")}"
SV_COHORT_VCF="${SV_COHORT_VCF:-${WORK_DIR}/sv_callset/cohort_trios.sv.vcf.gz}"
USE_SV_CALLSET="${USE_SV_CALLSET:-1}"
SV_MAX_AF="${SV_MAX_AF:-0.01}"
SV_MIN_GQ="${SV_MIN_GQ:-20}"
SV_PASS_ONLY="${SV_PASS_ONLY:-1}"

# ── Downloads ────────────────────────────────────────────────────────────────
# CRAM URLs in the sequence indexes point at the ENA FTP. The same files are
# served over HTTPS with range support by ENA and, byte-identical, by the AWS
# Open Data mirror (free egress, S3 throughput) - tried first. Every file
# passes the index MD5 before a job sees it, whichever source served it.
EXPECTED_FTP_PREFIX="${EXPECTED_FTP_PREFIX:-ftp://ftp.sra.ebi.ac.uk/}"
ENA_HTTPS_BASE="${ENA_HTTPS_BASE:-https://ftp.sra.ebi.ac.uk}"
S3_HTTPS_BASE="${S3_HTTPS_BASE:-https://1000genomes.s3.amazonaws.com}"  # empty = ENA only
# Where each release batch lives in the bucket; the ENA run accession and
# file name are appended.
S3_PREFIX_2504="${S3_PREFIX_2504:-1000G_2504_high_coverage/data}"
S3_PREFIX_698="${S3_PREFIX_698:-1000G_2504_high_coverage/additional_698_related/data}"
DOWNLOAD_CONNECTIONS="${DOWNLOAD_CONNECTIONS:-16}"  # HTTPS range streams per file
# The manager stages this many trios at once (one file at a time per trio,
# DOWNLOAD_CONNECTIONS streams each). The WAN is the scarce resource: a few
# fast transfers beat many starved ones.
TRIO_SLOTS="${TRIO_SLOTS:-2}"
# Downloads pause while this many trios sit on disk awaiting or running
# discovery (about 45 GB each), so a stalled queue cannot fill scratch.
MAX_STAGED_TRIOS="${MAX_STAGED_TRIOS:-4}"
# Stop after staging this many trios; 0 = the whole manifest. For smoke tests.
TRIO_LIMIT="${TRIO_LIMIT:-0}"
# 1 = leave the CRAMs in place after discovery (default: a finished trio
# deletes its CRAMs unless another selected trio still needs one of them).
KEEP_CRAMS="${KEEP_CRAMS:-0}"

# ── Per-trio discovery job (01b_discover_trio.sh) ────────────────────────────
# A 30x NYGC trio is ~45 GB of CRAM (about 15 GB each). Module 1 sizes its
# Bloom filter and hash to --memory; 128 GB is comfortable for a 30x child.
# Walltime depends on the node's I/O; the README gives observed figures to
# tune these.
DISCOVERY_CPUS="${DISCOVERY_CPUS:-16}"
DISCOVERY_MEM="${DISCOVERY_MEM:-128G}"
DISCOVERY_TIME="${DISCOVERY_TIME:-36:00:00}"
MIN_CHILD_COUNT="${MIN_CHILD_COUNT:-3}"
DISCOVERY_EXTRA_ARGS="${DISCOVERY_EXTRA_ARGS:-}"   # appended to kmer-discovery
# Scratch for kmer-discovery intermediates (100+ GB per trio). Empty = a
# per-trio directory under TMP_ROOT. Point it at node-local disk only if that
# disk is large enough; never at a RAM-backed /tmp.
DISCOVERY_TMP_DIR="${DISCOVERY_TMP_DIR:-}"
# The non-human fraction screen (kmer-discovery --kraken2-db): a Kraken2
# database directory (hash.k2d, opts.k2d, taxo.k2d and its taxonomy; the
# image has kraken2, and scripts/download_kraken2_db.sh fetches PrackenDB).
# It is bound into the discovery job, which classifies the informative reads
# of every SV and AMBIGUOUS region and writes the BED columns nhf_reads and
# nhf ("." = not screened). Empty = no screen, recorded as disabled in the
# run's metrics; kmer-nhf can add the screen to a finished run afterwards.
# The database is loaded into memory (PrackenDB: tens of GB): size
# DISCOVERY_MEM for it, or add --kraken2-memory-mapping to
# DISCOVERY_EXTRA_ARGS.
KRAKEN2_DB="${KRAKEN2_DB:-}"

# ── Reference index job (00b_build_ref_index.sh) ─────────────────────────────
# jellyfish count -s 4G at k=31 takes about 23 GB of RAM.
REF_INDEX_CPUS="${REF_INDEX_CPUS:-8}"
REF_INDEX_MEM="${REF_INDEX_MEM:-64G}"
REF_INDEX_TIME="${REF_INDEX_TIME:-06:00:00}"
REF_JF_HASH_SIZE="${REF_JF_HASH_SIZE:-4G}"

# ── Scheduler ────────────────────────────────────────────────────────────────
# Site-specific sbatch flags for EVERY submission: partition, account, QoS,
# constraints. Word-split on purpose, so keep it to simple flags.
SBATCH_EXTRA_ARGS="${SBATCH_EXTRA_ARGS:-}"
JOB_PREFIX="${JOB_PREFIX:-1kgtrio}"

# ── Cohort summary (02_summarize.sh) ─────────────────────────────────────────
# A region seen in this many children or more is reported as recurrent and
# dropped from the de novo SV candidate list.
RECURRENT_MIN_TRIOS="${RECURRENT_MIN_TRIOS:-2}"
RECURRENT_SLACK_BP="${RECURRENT_SLACK_BP:-500}"

# ── Helpers shared by the scripts ────────────────────────────────────────────

# The script that sourced this file: the first call frame outside it.
config_caller_script() {
  local frame
  for frame in ${BASH_SOURCE[@]+"${BASH_SOURCE[@]}"}; do
    [[ "${frame}" != "${BASH_SOURCE[0]}" ]] && { echo "${frame}"; return 0; }
  done
  echo "$0"
}

# The sourcing script's header comment, then the effective configuration.
config_usage() {
  local script name
  script="$(config_caller_script)"
  sed -n '2,/^$/p' "${script}" | sed '$d' | sed 's/^# \{0,1\}//'
  echo "Configuration: any flag below, an exported variable, a site file (--config FILE,"
  echo "\$KMER_TRIOS_CONFIG, or config.local.sh), or the defaults in config.sh. Flags also"
  echo "take --set NAME=VALUE. Effective values:"
  for name in ${CONFIG_VARIABLES}; do
    printf '  --%-24s %s\n' "$(printf '%s' "${name}" | tr 'A-Z_' 'a-z-')" "${!name}"
  done
}

# Write the effective configuration as a file --config can read back.
write_config_snapshot() {
  local path="$1" name
  mkdir -p "$(dirname "${path}")"
  {
    echo "# 1000 Genomes trio screen: configuration snapshot, written $(date) by $(config_caller_script)"
    for name in ${CONFIG_VARIABLES}; do
      printf 'export %s=%q\n' "${name}" "${!name}"
    done
  } > "${path}.tmp$$"
  mv "${path}.tmp$$" "${path}"
}

# Load MODULES_TO_LOAD with the environment-modules system, when there is one.
load_modules() {
  [[ -n "${MODULES_TO_LOAD}" ]] || return 0
  local init m
  if ! command -v module >/dev/null 2>&1; then
    for init in "${MODULESHOME:-}/init/bash" /etc/profile.d/modules.sh /etc/profile.d/lmod.sh; do
      if [[ -n "${init}" && -r "${init}" ]]; then
        # shellcheck disable=SC1090
        source "${init}" >/dev/null 2>&1 || true
        break
      fi
    done
  fi
  for m in ${MODULES_TO_LOAD}; do
    module load "${m}" 2>/dev/null || echo "WARNING: could not load module ${m}" >&2
  done
}

# For scripts without positional arguments: anything left over is a typo.
reject_extra_args() {
  (( $# == 0 )) && return 0
  echo "ERROR: unexpected argument(s): $* (see --help for the flags)" >&2
  exit 2
}

# MD5 of a file, with md5sum (Linux) or md5 (macOS) - whichever exists.
file_md5() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

file_size() {
  stat -c %s "$1" 2>/dev/null || stat -f %z "$1" 2>/dev/null
}

# apptainer, or singularity on sites that still ship it under that name.
find_container_runtime() {
  if command -v apptainer >/dev/null 2>&1; then
    echo apptainer
  elif command -v singularity >/dev/null 2>&1; then
    echo singularity
  else
    return 1
  fi
}

# Comma-joined list of distinct directories to bind into the container.
# Usage: bind_list DIR [DIR...]  (missing directories are skipped)
bind_list() {
  local d out=""
  for d in "$@"; do
    [[ -n "${d}" && -d "${d}" ]] || continue
    d="$(cd "${d}" && pwd)"
    case ",${out}," in
      *",${d},"*) ;;
      *) out="${out:+${out},}${d}" ;;
    esac
  done
  echo "${out}"
}

# Submit a job with the site flags, retrying transient sbatch failures.
# Prints the job id (sbatch --parsable). Usage: submit_job [sbatch flags...] SCRIPT [args...]
submit_job() {
  local try jobid
  for try in 1 2 3; do
    # shellcheck disable=SC2086
    if jobid=$(sbatch --parsable ${SBATCH_EXTRA_ARGS} "$@"); then
      echo "${jobid%%;*}"
      return 0
    fi
    echo "  sbatch failed (attempt ${try}/3); retrying in 30 s..." >&2
    sleep 30
  done
  return 1
}

# Memory in whole GB from a SLURM-style size ("128G", "4096M", "2T").
mem_gb() {
  local v="$1" n u
  n="${v%[GgMmTtKk]}"
  u="${v#"${n}"}"
  case "${u}" in
    T|t) echo $(( n * 1024 )) ;;
    M|m) echo $(( n / 1024 )) ;;
    K|k) echo $(( n / 1024 / 1024 )) ;;
    *)   echo "${n}" ;;
  esac
}

# The three (ROLE SAMPLE CRAM_URL CRAI_URL MD5 BATCH) rows of a trio.
trio_rows() {
  awk -F'\t' -v t="$1" 'NR > 1 && $1 == t { print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' "${TRIO_FILES_TSV}"
}

# A trio is done when its job wrote the marker beside the outputs.
trio_done() {
  [[ -s "${RESULTS_DIR}/$1/$1.done" ]]
}

if (( CONFIG_SHOW_HELP )); then
  config_usage
  exit 0
fi
