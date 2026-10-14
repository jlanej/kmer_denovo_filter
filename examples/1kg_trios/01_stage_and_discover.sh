#!/usr/bin/env bash
# =============================================================================
# 01_stage_and_discover.sh — stage each trio's CRAMs, then run discovery on it
# =============================================================================
#
# One long-running manager job stages the cohort a few trios at a time -
# TRIO_SLOTS concurrent trios, each downloading its three CRAMs one after the
# other over the S3 mirror or ENA (aria2c, else parallel curl, else wget) and
# verifying every CRAM against the index MD5 - and the moment a trio is
# complete on disk, its discovery run is submitted as an independent SLURM
# job (01b_discover_trio.sh). The WAN is the scarce resource, so it gets few
# fast transfers; discovery concurrency is whatever the scheduler grants.
#
# The sweep is idempotent: finished trios are skipped, trios whose discovery
# job is queued or running are skipped, and CRAMs already on disk (a dead
# run, or files staged with Globus) are verified and used without
# downloading. Re-running the script continues where the last sweep stopped.
# Run one manager at a time.
#
# Usage:
#   bash 01_stage_and_discover.sh                 # submits the manager job
#   MANAGER_LOCAL=1 bash 01_stage_and_discover.sh # run in place (DTN/login node, tmux)
#   TRIO_LIMIT=2 bash 01_stage_and_discover.sh    # smoke test: two trios
# =============================================================================

#SBATCH --job-name=1kgtrio_manager
#SBATCH --output=logs/manager_%j.out
#SBATCH --error=logs/manager_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=96:00:00

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
source "${CONFIG_FILE}"

mkdir -p "${LOG_DIR}" "${CRAM_DIR}" "${RESULTS_DIR}" "${STATE_DIR}" "${TMP_ROOT}"
DISCOVERY_JOB_SCRIPT="${CONFIG_DIR}/01b_discover_trio.sh"

# ── Preconditions ────────────────────────────────────────────────────────────
for f in "${TRIOS_TSV}" "${TRIO_FILES_TSV}"; do
  [[ -s "${f}" ]] || { echo "ERROR: manifest missing or empty: ${f} - run 00_setup.sh first."; exit 1; }
done
[[ -s "${SIF_IMAGE}" ]] || { echo "ERROR: container image missing: ${SIF_IMAGE} - run 00_setup.sh first."; exit 1; }
[[ -s "${REF_FASTA}" && -s "${REF_FASTA}.fai" ]] || { echo "ERROR: reference FASTA or .fai missing: ${REF_FASTA} - run 00_setup.sh first."; exit 1; }
if [[ ! -s "${REF_JF}" ]]; then
  echo "ERROR: reference Jellyfish index missing: ${REF_JF}"
  echo "       00_setup.sh submits its build; wait for that job (squeue -n ${JOB_PREFIX}_refidx) or run 00b_build_ref_index.sh."
  exit 1
fi
TOTAL_TRIOS=$(tail -n +2 "${TRIOS_TSV}" | wc -l | tr -d '[:space:]')

if [[ -z "${SLURM_JOB_ID:-}" && "${MANAGER_LOCAL:-0}" != "1" ]]; then
  # ── Submission mode ────────────────────────────────────────────────────────
  echo "Submitting the stage-and-discover manager (${TRIO_SLOTS} trios at a time)..."
  echo "  Manager log:        ${LOG_DIR}/manager_<jobid>.out"
  echo "  Per-trio stage log: ${LOG_DIR}/stage_<trio>.log"
  # shellcheck disable=SC2086
  exec sbatch ${SBATCH_EXTRA_ARGS} \
    --output="${LOG_DIR}/manager_%j.out" \
    --error="${LOG_DIR}/manager_%j.err" \
    "${BASH_SOURCE[0]}"
fi

# ── Transports ───────────────────────────────────────────────────────────────

# ENA serves the FTP tree over HTTPS with range support.
https_url_for() {
  echo "${1/ftp:\/\/ftp.sra.ebi.ac.uk/${ENA_HTTPS_BASE}}"
}

# The AWS Open Data mirror of the NYGC CRAMs (byte-identical to ENA's; every
# file still passes the index MD5). Its layout drops ENA's run-prefix level
# and splits the two release batches:
#   vol1/run/ERR324/ERR3239480/NA12718.final.cram (batch 2504)
#     -> 1000G_2504_high_coverage/data/ERR3239480/NA12718.final.cram
#   batch 698 -> 1000G_2504_high_coverage/additional_698_related/data/...
s3_url_for() {
  local ftp_url="$1" batch="$2"
  [[ -n "${S3_HTTPS_BASE}" ]] || { echo ""; return 0; }
  if [[ "${ftp_url}" != ftp://ftp.sra.ebi.ac.uk/vol1/run/*/*/* ]]; then
    echo ""
    return 0
  fi
  local rest="${ftp_url#ftp://ftp.sra.ebi.ac.uk/vol1/run/}"
  local run_and_file="${rest#*/}"
  case "${batch}" in
    2504) echo "${S3_HTTPS_BASE}/1000G_2504_high_coverage/data/${run_and_file}" ;;
    698)  echo "${S3_HTTPS_BASE}/1000G_2504_high_coverage/additional_698_related/data/${run_and_file}" ;;
    *)    echo "" ;;
  esac
}

remote_size() {
  curl -sS -I -L --max-time 60 "$1" 2>/dev/null \
    | awk '{ if (tolower($0) ~ /^content-length:/) v = $2 }
           END { gsub(/\r/, "", v); print v }'
}

# aria2c as a command array: the host's when present, else the container's.
ARIA2C=()
ARIA2C_HOW=""
resolve_aria2c() {
  if command -v aria2c >/dev/null 2>&1; then
    ARIA2C=(aria2c)
    ARIA2C_HOW="host"
    return 0
  fi
  local runtime
  if runtime="$(find_container_runtime)" \
     && "${runtime}" exec "${SIF_IMAGE}" aria2c --version >/dev/null 2>&1; then
    ARIA2C=("${runtime}" exec --bind "${CRAM_DIR}" "${SIF_IMAGE}" aria2c)
    ARIA2C_HOW="via ${SIF_IMAGE}"
    return 0
  fi
  return 1
}

download_aria2() {
  local url="$1" dest="$2"
  (( ${#ARIA2C[@]} )) || return 1
  if ! "${ARIA2C[@]}" \
    --quiet=true --continue=true \
    --max-connection-per-server="${DOWNLOAD_CONNECTIONS}" \
    --split="${DOWNLOAD_CONNECTIONS}" --min-split-size=1M \
    --max-tries=5 --retry-wait=10 \
    --auto-file-renaming=false --allow-overwrite=true \
    --dir "$(dirname "${dest}")" --out "$(basename "${dest}")" \
    "${url}"; then
    rm -f "${dest:?}" "${dest:?}.aria2"
    return 1
  fi
  rm -f "${dest:?}.aria2"
  # The S3 mirror has served a 0-byte CRAI with HTTP 200; nothing fetched
  # here is legitimately empty, so an empty result is a failed transport.
  if [[ ! -s "${dest}" ]]; then
    echo "  WARNING: aria2c returned an empty file for ${url}; trying the next source."
    rm -f "${dest:?}"
    return 1
  fi
}

# Parallel byte-range download with curl alone, for hosts without aria2c.
download_curl_parallel() {
  local url="$1" dest="$2"
  local total chunk parts_dir rc=0 i pid
  command -v curl >/dev/null 2>&1 || return 1
  total="$(remote_size "${url}")"
  [[ "${total}" =~ ^[0-9]+$ ]] || return 1
  (( total > 0 )) || return 1
  local n="${DOWNLOAD_CONNECTIONS}"
  (( n < 1 )) && n=1
  chunk=$(( (total + n - 1) / n ))
  parts_dir="${dest}.parts.$$"
  rm -rf "${parts_dir:?}"
  mkdir -p "${parts_dir}" || return 1
  local pids=()
  for (( i = 0; i < n; i++ )); do
    local start=$(( i * chunk ))
    local end=$(( start + chunk - 1 ))
    (( start >= total )) && break
    (( end >= total )) && end=$(( total - 1 ))
    curl -sS -L --fail --retry 5 --retry-delay 10 \
      -r "${start}-${end}" -o "${parts_dir}/part.$(printf '%04d' "${i}")" "${url}" &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "${pid}" || rc=1
  done
  if (( rc != 0 )) || ! cat "${parts_dir}"/part.* > "${dest}"; then
    rm -rf "${parts_dir:?}"; rm -f "${dest:?}"
    return 1
  fi
  rm -rf "${parts_dir:?}"
  if [[ "$(file_size "${dest}")" != "${total}" ]]; then
    echo "  ERROR: size mismatch after parallel download (expected ${total}, got $(file_size "${dest}"))"
    rm -f "${dest:?}"
    return 1
  fi
}

download_wget() {
  local url="$1" dest="$2"
  command -v wget >/dev/null 2>&1 || return 1
  if ! wget -q -c -O "${dest}" "${url}"; then
    rm -f "${dest:?}"
    return 1
  fi
  [[ -s "${dest}" ]] || { rm -f "${dest:?}"; return 1; }
}

# Try each source (S3 mirror, then ENA) with each transport in descending
# order of speed; wget on the FTP URL is the last resort.
download_file() {
  local url="$1" dest="$2" label="$3" batch="${4:-}"
  local sources=() names=() https i
  https="$(s3_url_for "${url}" "${batch}")"
  if [[ -n "${https}" ]]; then
    sources+=("${https}"); names+=("S3")
  fi
  sources+=("$(https_url_for "${url}")"); names+=("ENA")
  for i in "${!sources[@]}"; do
    if download_aria2 "${sources[${i}]}" "${dest}"; then
      echo "  ${label}: aria2c download complete (${names[${i}]}, ${DOWNLOAD_CONNECTIONS} streams)."
      return 0
    fi
    if download_curl_parallel "${sources[${i}]}" "${dest}"; then
      echo "  ${label}: parallel curl download complete (${names[${i}]}, ${DOWNLOAD_CONNECTIONS} ranges)."
      return 0
    fi
  done
  if download_wget "${url}" "${dest}"; then
    echo "  ${label}: wget download complete (single stream)."
    return 0
  fi
  echo "  ERROR: all download methods failed for ${url}"
  return 1
}

# ── Staging ──────────────────────────────────────────────────────────────────

LOCK_DIR="${STATE_DIR}/locks"

# A parent shared by two selected trios must not be downloaded twice at
# once; the per-sample lock is a directory (mkdir is atomic).
lock_sample() {
  local sample="${1:?}"
  until mkdir "${LOCK_DIR:?}/${sample}" 2>/dev/null; do
    sleep 15
  done
}
unlock_sample() {
  local sample="${1:?}"
  rmdir "${LOCK_DIR:?}/${sample}" 2>/dev/null || true
}

# stage_sample SAMPLE CRAM_URL CRAI_URL MD5 BATCH: the CRAM and CRAI are on
# disk and the CRAM matches its MD5 when this returns 0.
stage_sample() {
  local sample="$1" cram_url="$2" crai_url="$3" cram_md5="$4" batch="$5"
  local cram="${CRAM_DIR}/${sample}.cram" crai="${CRAM_DIR}/${sample}.cram.crai"
  local attempt actual
  lock_sample "${sample}"
  # shellcheck disable=SC2034
  for attempt in 1 2; do
    if [[ ! -s "${cram}" ]]; then
      echo "  ${sample}: downloading CRAM..."
      download_file "${cram_url}" "${cram}" "${sample} CRAM" "${batch}" || { unlock_sample "${sample}"; return 1; }
    else
      echo "  ${sample}: CRAM already present."
    fi
    if [[ ! -s "${crai}" ]]; then
      echo "  ${sample}: downloading CRAI..."
      download_file "${crai_url}" "${crai}" "${sample} CRAI" "${batch}" || { unlock_sample "${sample}"; return 1; }
    fi
    echo "  ${sample}: verifying MD5..."
    actual="$(file_md5 "${cram}")"
    if [[ "${actual}" == "${cram_md5}" ]]; then
      echo "  ${sample}: MD5 verified."
      unlock_sample "${sample}"
      return 0
    fi
    echo "  ${sample}: MD5 mismatch (expected ${cram_md5}, got ${actual}); removing and re-downloading."
    rm -f "${cram:?}" "${crai:?}"
  done
  echo "  ERROR: ${sample} failed MD5 after a fresh download."
  unlock_sample "${sample}"
  return 1
}

# stage_trio TRIO: stage its three samples, then submit the discovery job.
stage_trio() {
  local trio="$1" role sample cram_url crai_url md5 batch jobid
  echo "============================================================"
  echo " Trio ${trio} - started $(date)"
  echo "============================================================"
  while IFS=$'\t' read -r role sample cram_url crai_url md5 batch; do
    [[ -n "${sample}" ]] || continue
    echo "[${role}] ${sample}"
    stage_sample "${sample}" "${cram_url}" "${crai_url}" "${md5}" "${batch}" || return 1
  done < <(trio_rows "${trio}")
  echo "Submitting discovery for ${trio} (${DISCOVERY_CPUS} CPUs, ${DISCOVERY_MEM}, ${DISCOVERY_TIME})..."
  if jobid=$(submit_job \
      --job-name="${JOB_PREFIX}_disc_${trio}" \
      --output="${LOG_DIR}/discover_${trio}_%j.out" \
      --error="${LOG_DIR}/discover_${trio}_%j.err" \
      --cpus-per-task="${DISCOVERY_CPUS}" --mem="${DISCOVERY_MEM}" --time="${DISCOVERY_TIME}" \
      "${DISCOVERY_JOB_SCRIPT}" "${trio}"); then
    echo "  Job ${jobid} submitted for ${trio}."
    return 0
  fi
  echo "  ERROR: could not submit the discovery job for ${trio}; the verified CRAMs are kept"
  echo "         and the next sweep dispatches them without re-downloading."
  return 1
}

# Trios with all three CRAMs on disk that are not finished: what the
# backpressure gate counts (about 45 GB each).
count_staged_trios() {
  local trio n=0 sample role _u1 _u2 _u3 _u4 complete
  while IFS=$'\t' read -r trio _u1; do
    [[ -n "${trio}" ]] || continue
    trio_done "${trio}" && continue
    complete=1
    while IFS=$'\t' read -r role sample _u2 _u3 _u4; do
      [[ -s "${CRAM_DIR}/${sample}.cram" ]] || { complete=0; break; }
    done < <(trio_rows "${trio}")
    (( complete )) && n=$(( n + 1 ))
  done < <(tail -n +2 "${TRIOS_TSV}")
  echo "${n}"
}

trio_needs_download() {
  local role sample _u2 _u3 _u4
  while IFS=$'\t' read -r role sample _u2 _u3 _u4; do
    [[ -s "${CRAM_DIR}/${sample}.cram" && -s "${CRAM_DIR}/${sample}.cram.crai" ]] || return 0
  done < <(trio_rows "$1")
  return 1
}

# ── Manager ──────────────────────────────────────────────────────────────────

RESULT_DIR="${STATE_DIR}/results"
rm -rf "${RESULT_DIR:?}" "${LOCK_DIR:?}"
mkdir -p "${RESULT_DIR}" "${LOCK_DIR}"

IN_FLIGHT_FILE="${STATE_DIR}/in_flight.$$"
: > "${IN_FLIGHT_FILE}"
if command -v squeue >/dev/null 2>&1; then
  squeue -h -u "${USER:-$(id -un)}" -o "%j" 2>/dev/null \
    | sed -n "s/^${JOB_PREFIX}_disc_//p" > "${IN_FLIGHT_FILE}" || true
fi

NEEDED=()
SKIPPED_DONE=0
SKIPPED_IN_FLIGHT=0
while IFS=$'\t' read -r TRIO _rest; do
  [[ -n "${TRIO}" ]] || continue
  if trio_done "${TRIO}"; then
    SKIPPED_DONE=$(( SKIPPED_DONE + 1 ))
    continue
  fi
  if grep -qxF "${TRIO}" "${IN_FLIGHT_FILE}"; then
    SKIPPED_IN_FLIGHT=$(( SKIPPED_IN_FLIGHT + 1 ))
    continue
  fi
  NEEDED+=("${TRIO}")
done < <(tail -n +2 "${TRIOS_TSV}")
rm -f "${IN_FLIGHT_FILE:?}"

echo "Stage-and-discover manager: ${TOTAL_TRIOS} trios in ${TRIOS_TSV}."
resolve_aria2c || true
if (( ${#ARIA2C[@]} )); then
  echo "  aria2c: ${ARIA2C_HOW}"
else
  echo "  aria2c: not available - HTTPS transfers use parallel curl (${DOWNLOAD_CONNECTIONS} ranges)"
fi
echo "  Finished:            ${SKIPPED_DONE}"
echo "  Discovery in flight: ${SKIPPED_IN_FLIGHT}"
echo "  To stage:            ${#NEEDED[@]} (${TRIO_SLOTS} at a time)"
(( TRIO_LIMIT > 0 )) && echo "  TRIO_LIMIT=${TRIO_LIMIT}: stopping after that many trios."

if (( ${#NEEDED[@]} == 0 )); then
  echo "Nothing to stage."
  exit 0
fi

DISPATCHED=0
for TRIO in "${NEEDED[@]}"; do
  if (( TRIO_LIMIT > 0 && DISPATCHED >= TRIO_LIMIT )); then
    break
  fi
  # Polling rather than wait -n: transfers run for an hour, and bash older
  # than 4.3 has no wait -n.
  while (( $(jobs -rp | wc -l) >= TRIO_SLOTS )); do
    sleep 10
  done
  # Disk backpressure, only for trios that would download: one already on
  # disk consumes no new space, and dispatching it is what frees some.
  if trio_needs_download "${TRIO}"; then
    while (( $(count_staged_trios) >= MAX_STAGED_TRIOS )); do
      echo "Backpressure: ${MAX_STAGED_TRIOS} trios on disk await discovery; pausing downloads for 60 s..."
      sleep 60
    done
  fi
  echo "[$(date '+%H:%M:%S')] slot -> ${TRIO} (log: ${LOG_DIR}/stage_${TRIO}.log)"
  (
    if stage_trio "${TRIO}"; then
      echo "ok" > "${RESULT_DIR}/${TRIO}"
    else
      echo "fail" > "${RESULT_DIR}/${TRIO}"
    fi
  ) > "${LOG_DIR}/stage_${TRIO}.log" 2>&1 &
  DISPATCHED=$(( DISPATCHED + 1 ))
done
wait || true

OK_COUNT=$(grep -lx "ok" "${RESULT_DIR}"/* 2>/dev/null | wc -l | tr -d '[:space:]') || true
FAIL_COUNT=$(grep -lx "fail" "${RESULT_DIR}"/* 2>/dev/null | wc -l | tr -d '[:space:]') || true
echo ""
echo "============================================================"
echo " Staging sweep finished: ${OK_COUNT} trios dispatched to discovery, ${FAIL_COUNT} failed."
if (( FAIL_COUNT > 0 )); then
  echo " Failed trios (see ${LOG_DIR}/stage_<trio>.log):"
  grep -lx "fail" "${RESULT_DIR}"/* | while IFS= read -r f; do echo "   $(basename "${f}")"; done
  echo " Re-run this script to retry them."
fi
echo " Discovery jobs run independently: squeue -u \$USER | grep ${JOB_PREFIX}_disc"
echo " Summarise finished trios any time: bash 02_summarize.sh"
echo "============================================================"
(( FAIL_COUNT == 0 )) || exit 1
