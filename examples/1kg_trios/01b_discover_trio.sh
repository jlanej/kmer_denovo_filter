#!/usr/bin/env bash
# =============================================================================
# 01b_discover_trio.sh — kmer-discovery on one trio whose CRAMs are staged
# =============================================================================
#
# Spawned by 01_stage_and_discover.sh once a trio's three CRAMs have passed
# their MD5 checks. It
#
#   1. lists the child-only SVs of the published call set for this trio
#      (child_only_svs.py) and passes them as --dnm-regions, so the run
#      reports which of them discovery supports;
#   2. runs kmer-discovery in the container with the shared reference index,
#      sized to the job's allocation, and with --kraken2-db when KRAKEN2_DB
#      is set (the non-human fraction screen of its informative reads);
#   3. records the wall time, writes the done marker, and deletes the CRAMs
#      unless KEEP_CRAMS=1 or another selected trio still needs one.
#
# Outputs land in RESULTS_DIR/<trio>/<trio>.* (BED, BEDPE, informative BAM,
# metrics JSON, summary, HTML report, log, timing).
#
# Usage (normally via the manager; by hand for a trio staged some other way):
#   sbatch 01b_discover_trio.sh [--config FILE | --work-dir DIR ...] TRIO_ID [verify]
#   (verify = check the MD5s first; add --output=$WORK_DIR/logs/... to keep the
#   log with the others)
# =============================================================================

#SBATCH --job-name=1kgtrio_disc
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=36:00:00

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/config.sh" ]]; then
  CONFIG_FILE="${SLURM_SUBMIT_DIR}/config.sh"
elif [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  CONFIG_FILE="${SCRIPT_DIR}/config.sh"
elif [[ -f "$(pwd)/config.sh" ]]; then
  CONFIG_FILE="$(pwd)/config.sh"
else
  echo "ERROR: Could not find config.sh (checked \$SLURM_SUBMIT_DIR, script directory, and current directory)."
  exit 1
fi
source "${CONFIG_FILE}" "$@"
set -- ${CONFIG_REMAINING_ARGS[@]+"${CONFIG_REMAINING_ARGS[@]}"}
load_modules

TRIO="${1:?usage: 01b_discover_trio.sh [--config FILE | --work-dir DIR ...] TRIO_ID [verify]}"
if [[ "${TRIO}" == -* ]]; then
  echo "ERROR: unknown flag ${TRIO} (see --help); usage: 01b_discover_trio.sh TRIO_ID [verify]" >&2
  exit 2
fi
VERIFY="${2:-}"

OUT_DIR="${RESULTS_DIR}/${TRIO}"
PREFIX="${OUT_DIR}/${TRIO}"
mkdir -p "${OUT_DIR}"

# The trio's members, from the manifest.
CHILD=""; FATHER=""; MOTHER=""
CHILD_MD5=""; FATHER_MD5=""; MOTHER_MD5=""
while IFS=$'\t' read -r role sample _cram _crai md5 _batch; do
  case "${role}" in
    child)  CHILD="${sample}";  CHILD_MD5="${md5}" ;;
    father) FATHER="${sample}"; FATHER_MD5="${md5}" ;;
    mother) MOTHER="${sample}"; MOTHER_MD5="${md5}" ;;
  esac
done < <(trio_rows "${TRIO}")
if [[ -z "${CHILD}" || -z "${FATHER}" || -z "${MOTHER}" ]]; then
  echo "ERROR: trio ${TRIO} is not in ${TRIO_FILES_TSV} (or lacks a member)."
  exit 1
fi
CHILD_CRAM="${CRAM_DIR}/${CHILD}.cram"
FATHER_CRAM="${CRAM_DIR}/${FATHER}.cram"
MOTHER_CRAM="${CRAM_DIR}/${MOTHER}.cram"

echo "============================================================"
echo " Discovery: trio ${TRIO} (child ${CHILD}, father ${FATHER}, mother ${MOTHER})"
echo " Started: $(date) on $(hostname)${SLURM_JOB_ID:+ (job ${SLURM_JOB_ID})}"
echo "============================================================"

# Delete the trio's CRAMs, keeping any that another unfinished trio needs.
cleanup_crams() {
  local sample other keep
  [[ "${KEEP_CRAMS}" == "1" ]] && { echo "KEEP_CRAMS=1: CRAMs left in ${CRAM_DIR}."; return 0; }
  for sample in "${CHILD}" "${FATHER}" "${MOTHER}"; do
    keep=""
    while IFS= read -r other; do
      [[ -n "${other}" ]] || continue
      trio_done "${other}" || keep="${other}"
    done < <(awk -F'\t' -v s="${sample}" -v t="${TRIO}" 'NR > 1 && $3 == s && $1 != t { print $1 }' "${TRIO_FILES_TSV}")
    if [[ -n "${keep}" ]]; then
      echo "  Keeping ${sample}.cram: trio ${keep} still needs it."
    else
      rm -f "${CRAM_DIR}/${sample}.cram" "${CRAM_DIR}/${sample}.cram.crai"
      echo "  Removed ${sample}.cram"
    fi
  done
}

# Idempotent under a requeue: finished work only tidies up.
if trio_done "${TRIO}"; then
  echo "SKIP: ${TRIO} is already done ($(cat "${PREFIX}.done"))."
  cleanup_crams
  exit 0
fi

for f in "${CHILD_CRAM}" "${FATHER_CRAM}" "${MOTHER_CRAM}"; do
  if [[ ! -s "${f}" || ! -s "${f}.crai" ]]; then
    echo "ERROR: ${f} (or its .crai) is missing. The manager verifies both before submitting"
    echo "       this job, so something removed them; re-run 01_stage_and_discover.sh."
    exit 1
  fi
done

if [[ "${VERIFY}" == "verify" ]]; then
  echo "Verifying staged CRAMs against the index MD5s..."
  for pair in "${CHILD_CRAM}:${CHILD_MD5}" "${FATHER_CRAM}:${FATHER_MD5}" "${MOTHER_CRAM}:${MOTHER_MD5}"; do
    f="${pair%%:*}"; expected="${pair##*:}"
    actual="$(file_md5 "${f}")"
    if [[ "${actual}" != "${expected}" ]]; then
      echo "ERROR: ${f} failed MD5 (expected ${expected}, got ${actual}); left in place."
      exit 1
    fi
    echo "  $(basename "${f}"): MD5 verified."
  done
fi

RUNTIME="$(find_container_runtime)" || { echo "ERROR: apptainer/singularity not found on PATH."; exit 1; }
[[ -s "${SIF_IMAGE}" ]] || { echo "ERROR: container image missing: ${SIF_IMAGE}"; exit 1; }
[[ -s "${REF_JF}" ]] || { echo "ERROR: reference Jellyfish index missing: ${REF_JF}"; exit 1; }

# ── Resources, from the allocation ───────────────────────────────────────────
THREADS="${SLURM_CPUS_PER_TASK:-${DISCOVERY_CPUS}}"
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
  MEMORY_GB=$(( SLURM_MEM_PER_NODE / 1024 ))
elif [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
  MEMORY_GB=$(( SLURM_MEM_PER_CPU * THREADS / 1024 ))
else
  MEMORY_GB="$(mem_gb "${DISCOVERY_MEM}")"
fi
(( MEMORY_GB >= 1 )) || MEMORY_GB=1
TMP_DIR="${DISCOVERY_TMP_DIR:-${TMP_ROOT}}/${TRIO}"
mkdir -p "${TMP_DIR}"
echo "Threads: ${THREADS}; memory budget: ${MEMORY_GB} GB; tmp: ${TMP_DIR}"

# CONTAINER_EXAMPLE_DIR is bound when it is a host directory (a checkout of
# these helpers in place of the image's copy); the image's own /app path is
# not one and is skipped.
BINDS="$(bind_list "${WORK_DIR}" "${CRAM_DIR}" "${RESULTS_DIR}" "${TMP_DIR}" \
  "$(dirname "${REF_FASTA}")" "$(dirname "${REF_JF}")" "$(dirname "${SV_COHORT_VCF}")" "$(dirname "${SV_VCF}")" \
  "${CONTAINER_EXAMPLE_DIR}" "${KRAKEN2_DB}")"

# ── 1. Child-only SVs from the published call set ────────────────────────────
DNM_ARGS=()
if [[ "${USE_SV_CALLSET}" == "1" ]]; then
  if [[ -s "${SV_COHORT_VCF}" ]]; then
    SV_SOURCE="${SV_COHORT_VCF}"
  elif [[ -s "${SV_VCF}" ]]; then
    SV_SOURCE="${SV_VCF}"
  else
    SV_SOURCE=""
    echo "WARNING: SV call set not found (${SV_COHORT_VCF}); running without --dnm-regions."
  fi
  if [[ -n "${SV_SOURCE}" ]]; then
    echo "[1/3] Child-only SVs for ${CHILD} from $(basename "${SV_SOURCE}")..."
    SV_ARGS=(
      --vcf "${SV_SOURCE}" --child "${CHILD}" --father "${FATHER}" --mother "${MOTHER}"
      --out-regions "${PREFIX}.child_only_svs.tsv"
      --out-table "${PREFIX}.child_only_svs.detail.tsv"
      --out-inherited "${PREFIX}.inherited_svs.tsv"
      --max-af "${SV_MAX_AF}" --min-gq "${SV_MIN_GQ}"
    )
    [[ "${SV_PASS_ONLY}" == "1" ]] || SV_ARGS+=(--no-pass-only)
    "${RUNTIME}" exec --bind "${BINDS}" "${SIF_IMAGE}" \
      python3 "${CONTAINER_EXAMPLE_DIR}/child_only_svs.py" "${SV_ARGS[@]}"
    N_SV=$(grep -cv '^#' "${PREFIX}.child_only_svs.tsv" || true)
    if (( N_SV > 0 )); then
      DNM_ARGS=(--dnm-regions "${PREFIX}.child_only_svs.tsv")
      echo "  ${N_SV} child-only SVs pass the filters; discovery will report on each."
    else
      echo "  No child-only SV passes the filters; running without --dnm-regions."
    fi
    # Inherited SVs (child and one parent carry them) are the negative
    # controls: their k-mers are in a parent, so discovery should not
    # support them; where it does, that is the false-positive rate.
    N_NEG=$(grep -cv '^#' "${PREFIX}.inherited_svs.tsv" || true)
    if (( N_NEG > 0 )); then
      DNM_ARGS+=(--dnm-negative-regions "${PREFIX}.inherited_svs.tsv")
      echo "  ${N_NEG} inherited SVs go in as negative controls."
    fi
  fi
else
  echo "[1/3] SV call set disabled (USE_SV_CALLSET=0)."
fi

# ── 2. kmer-discovery ────────────────────────────────────────────────────────
CMD=(
  kmer-discovery
  --child "${CHILD_CRAM}" --father "${FATHER_CRAM}" --mother "${MOTHER_CRAM}"
  --ref-fasta "${REF_FASTA}" --ref-jf "${REF_JF}"
  --out-prefix "${PREFIX}"
  --kmer-size "${KMER_SIZE}" --min-child-count "${MIN_CHILD_COUNT}"
  --threads "${THREADS}" --memory "${MEMORY_GB}"
  --tmp-dir "${TMP_DIR}"
  --report "${PREFIX}.report.html"
)
(( ${#DNM_ARGS[@]} )) && CMD+=("${DNM_ARGS[@]}")
if [[ -n "${KRAKEN2_DB}" ]]; then
  # The non-human fraction screen. A database kmer-discovery cannot use is
  # recorded as a failed screen in the run's metrics; the run still succeeds.
  CMD+=(--kraken2-db "${KRAKEN2_DB}")
  [[ -d "${KRAKEN2_DB}" ]] || echo "WARNING: KRAKEN2_DB=${KRAKEN2_DB} is not a directory; the screen will be recorded as failed."
fi
if [[ -n "${DISCOVERY_EXTRA_ARGS}" ]]; then
  # shellcheck disable=SC2206
  CMD+=(${DISCOVERY_EXTRA_ARGS})
fi
echo "[2/3] Running: ${RUNTIME} exec --bind ${BINDS} ${SIF_IMAGE} ${CMD[*]}"
T0=$(date +%s)
set +e
"${RUNTIME}" exec --bind "${BINDS}" "${SIF_IMAGE}" "${CMD[@]}" 2>&1 | tee "${PREFIX}.log"
RC=${PIPESTATUS[0]}
set -e
T1=$(date +%s)
rm -rf "${TMP_DIR:?}"
if (( RC != 0 )); then
  echo "ERROR: kmer-discovery exited ${RC} after $(( T1 - T0 )) s; see ${PREFIX}.log. CRAMs kept for a retry."
  exit "${RC}"
fi
for ext in bed sv.bedpe metrics.json summary.txt; do
  if [[ ! -s "${PREFIX}.${ext}" ]]; then
    echo "ERROR: expected output missing: ${PREFIX}.${ext}"
    exit 1
  fi
done

# ── 3. Bookkeeping ───────────────────────────────────────────────────────────
N_REGIONS=$(grep -cv '^#' "${PREFIX}.bed" || true)
printf 'trio\tchild\twall_s\tstarted\tfinished\tnode\tthreads\tmemory_gb\tjob_id\timage\tregions\n' > "${PREFIX}.timing.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${TRIO}" "${CHILD}" "$(( T1 - T0 ))" "${T0}" "${T1}" "$(hostname)" "${THREADS}" "${MEMORY_GB}" \
  "${SLURM_JOB_ID:-}" "${CONTAINER_URI}" "${N_REGIONS}" >> "${PREFIX}.timing.tsv"
echo "[3/3] Done in $(( T1 - T0 )) s: ${N_REGIONS} candidate regions."
date '+%Y-%m-%d %H:%M:%S' > "${PREFIX}.done"
cleanup_crams
echo "============================================================"
echo " ${TRIO} complete. Finished: $(date)"
echo "============================================================"
