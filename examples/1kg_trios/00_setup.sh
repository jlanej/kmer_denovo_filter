#!/usr/bin/env bash
# =============================================================================
# 00_setup.sh — one-time setup for the 1000 Genomes trio screen
# =============================================================================
#
#   1. Creates the working tree under WORK_DIR.
#   2. Pulls the kmer_denovo_filter container (or uses a local .sif).
#   3. Downloads the GRCh38 reference the NYGC CRAMs were aligned to, and
#      its .fai index.
#   4. Downloads the NYGC sequence indexes and the 3,202-sample pedigree,
#      and selects the trios to run (select_trios.py -> trios.tsv,
#      trio_files.tsv).
#   5. Downloads the 3,202-sample SV call set and subsets it to the cohort's
#      samples (optional: USE_SV_CALLSET=0 skips it).
#   6. Submits the reference Jellyfish index build (00b_build_ref_index.sh);
#      SETUP_LOCAL=1 builds it in place instead.
#
# Every step is idempotent: existing files are kept. Delete trios.tsv and
# trio_files.tsv to select a different cohort.
#
# Usage:
#   export WORK_DIR=/scratch/$USER/1kg_trios      # and SBATCH_EXTRA_ARGS
#   N_TRIOS=10 bash 00_setup.sh
#   bash 00_setup.sh --work-dir /scratch/$USER/1kg_trios --n-trios 10   # same, as flags
#   bash 00_setup.sh --config my_site.sh          # or a site config file
#   bash 00_setup.sh --help                       # flags and effective values
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  CONFIG_FILE="${SCRIPT_DIR}/config.sh"
elif [[ -f "$(pwd)/config.sh" ]]; then
  CONFIG_FILE="$(pwd)/config.sh"
else
  echo "ERROR: Could not find config.sh (checked the script directory and the current directory)."
  exit 1
fi
source "${CONFIG_FILE}" "$@"
set -- ${CONFIG_REMAINING_ARGS[@]+"${CONFIG_REMAINING_ARGS[@]}"}
load_modules
reject_extra_args "$@"

# fetch URL DEST — resumable download with curl, else wget; atomic rename.
fetch() {
  local url="$1" dest="$2"
  mkdir -p "$(dirname "${dest}")"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 5 --retry-delay 10 -C - -o "${dest}.part" "${url}" \
      || wget -q -c -O "${dest}.part" "${url}"
  else
    wget -q -c -O "${dest}.part" "${url}"
  fi
  [[ -s "${dest}.part" ]] || { echo "ERROR: empty download: ${url}"; rm -f "${dest}.part"; return 1; }
  mv "${dest}.part" "${dest}"
}

validate_manifests() {
  local header
  header="$(head -n1 "${TRIOS_TSV}")"
  if [[ "${header}" != $'TRIO_ID\tCHILD\tFATHER\tMOTHER\tFAMILY_ID\tCHILD_SEX\tPOPULATION\tSUPERPOPULATION' ]]; then
    echo "ERROR: unexpected header in ${TRIOS_TSV}: ${header}"
    exit 1
  fi
  header="$(head -n1 "${TRIO_FILES_TSV}")"
  if [[ "${header}" != $'TRIO_ID\tROLE\tSAMPLE_ID\tCRAM_FTP_URL\tCRAI_FTP_URL\tCRAM_MD5\tRELEASE_BATCH' ]]; then
    echo "ERROR: unexpected header in ${TRIO_FILES_TSV}: ${header}"
    exit 1
  fi
  local n_trios n_rows bad
  n_trios=$(tail -n +2 "${TRIOS_TSV}" | wc -l | tr -d '[:space:]')
  n_rows=$(tail -n +2 "${TRIO_FILES_TSV}" | wc -l | tr -d '[:space:]')
  if (( n_trios == 0 || n_rows != 3 * n_trios )); then
    echo "ERROR: ${TRIOS_TSV} lists ${n_trios} trios but ${TRIO_FILES_TSV} has ${n_rows} rows (expected $(( 3 * n_trios )))."
    exit 1
  fi
  bad="$(awk -F'\t' -v prefix="${EXPECTED_FTP_PREFIX}" '
    NR == 1 { next }
    NF != 7 || $1 == "" || $3 == "" || index($4, prefix) != 1 || index($5, prefix) != 1 || $6 !~ /^[0-9a-f]{32}$/ { print NR; exit }
    $1 ~ /[ \r]/ || $3 ~ /[ \r]/ || $4 ~ /[ \r]/ { print NR; exit }
  ' "${TRIO_FILES_TSV}")"
  if [[ -n "${bad}" ]]; then
    echo "ERROR: ${TRIO_FILES_TSV} has an invalid row at line ${bad}. Delete both manifests and re-run setup."
    exit 1
  fi
}

echo "============================================================"
echo " 1000 Genomes trio screen — setup"
echo "============================================================"
echo " WORK_DIR:    ${WORK_DIR}"
echo " Container:   ${CONTAINER_URI}"
echo " Reference:   ${REF_FASTA}"
echo " Trios:       ${TRIOS_TSV} (N_TRIOS=${N_TRIOS}, seed=${TRIO_SEED:-none})"
echo " SV call set: $([[ "${USE_SV_CALLSET}" == "1" ]] && echo "${SV_VCF}" || echo "(disabled)")"
echo ""

# ── 1. Directories ───────────────────────────────────────────────────────────
echo "[1/6] Creating the working tree..."
mkdir -p "${WORK_DIR}" "${CRAM_DIR}" "${RESULTS_DIR}" "${LOG_DIR}" "${REF_DIR}" "${STATE_DIR}" "${TMP_ROOT}"
mkdir -p "$(dirname "${SV_VCF}")"

# ── 2. Container ─────────────────────────────────────────────────────────────
if [[ "${CONTAINER_URI}" == *.sif ]]; then
  [[ -s "${CONTAINER_URI}" ]] || { echo "ERROR: container .sif not found: ${CONTAINER_URI}"; exit 1; }
  if [[ "$(cd "$(dirname "${CONTAINER_URI}")" && pwd)/$(basename "${CONTAINER_URI}")" != "${SIF_IMAGE}" ]]; then
    ln -sf "${CONTAINER_URI}" "${SIF_IMAGE}"
  fi
  echo "[2/6] Using local container: ${CONTAINER_URI}"
elif [[ -s "${SIF_IMAGE}" ]]; then
  echo "[2/6] Container already present: ${SIF_IMAGE}"
else
  RUNTIME="$(find_container_runtime)" || { echo "ERROR: apptainer/singularity not found on PATH (module load apptainer?)."; exit 1; }
  echo "[2/6] Pulling ${CONTAINER_URI} -> ${SIF_IMAGE}"
  "${RUNTIME}" pull "${SIF_IMAGE}" "${CONTAINER_URI}"
fi

# ── 3. Reference ─────────────────────────────────────────────────────────────
if [[ -s "${REF_FASTA}" ]]; then
  echo "[3/6] Reference already present: ${REF_FASTA}"
else
  echo "[3/6] Downloading the reference (~3 GB): ${REF_BASE_URL}/${REF_NAME}"
  fetch "${REF_BASE_URL}/${REF_NAME}" "${REF_FASTA}"
fi
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo "  Fetching the FASTA index..."
  if ! fetch "${REF_BASE_URL}/${REF_NAME}.fai" "${REF_FASTA}.fai" 2>/dev/null; then
    RUNTIME="${RUNTIME:-$(find_container_runtime)}"
    echo "  Not available for download; building it with samtools faidx."
    "${RUNTIME}" exec --bind "$(bind_list "$(dirname "${REF_FASTA}")")" "${SIF_IMAGE}" samtools faidx "${REF_FASTA}"
  fi
fi

# ── 4. Indexes, pedigree, trio selection ─────────────────────────────────────
echo "[4/6] Sequence indexes and pedigree..."
[[ -s "${INDEX_FILE_2504}" ]] || fetch "${INDEX_URL_2504}" "${INDEX_FILE_2504}"
[[ -s "${INDEX_FILE_698}" ]]  || fetch "${INDEX_URL_698}" "${INDEX_FILE_698}"
[[ -s "${PED_FILE}" ]]        || fetch "${PED_URL}" "${PED_FILE}"
if [[ -s "${TRIOS_TSV}" && -s "${TRIO_FILES_TSV}" ]]; then
  echo "  Manifests already exist (delete them to select a different cohort)."
else
  SELECT_ARGS=(
    --ped "${PED_FILE}"
    --index "${INDEX_FILE_2504}" --index "${INDEX_FILE_698}"
    --out-trios "${TRIOS_TSV}" --out-files "${TRIO_FILES_TSV}"
    --n "${N_TRIOS}"
    --expected-complete "${EXPECTED_COMPLETE_TRIOS}"
  )
  [[ -n "${TRIO_SEED}" ]]             && SELECT_ARGS+=(--seed "${TRIO_SEED}")
  [[ -n "${TRIO_POPULATIONS}" ]]      && SELECT_ARGS+=(--populations "${TRIO_POPULATIONS}")
  [[ -n "${TRIO_SUPERPOPULATIONS}" ]] && SELECT_ARGS+=(--superpopulations "${TRIO_SUPERPOPULATIONS}")
  [[ -n "${TRIO_CHILDREN}" ]]         && SELECT_ARGS+=(--children "${TRIO_CHILDREN}")
  [[ "${TRIO_ALL_SIBLINGS}" == "1" ]] && SELECT_ARGS+=(--all-siblings)
  if command -v python3 >/dev/null 2>&1; then
    python3 "${SCRIPT_DIR}/select_trios.py" "${SELECT_ARGS[@]}"
  else
    RUNTIME="${RUNTIME:-$(find_container_runtime)}"
    "${RUNTIME}" exec --bind "$(bind_list "${WORK_DIR}")" "${SIF_IMAGE}" \
      python3 "${CONTAINER_EXAMPLE_DIR}/select_trios.py" "${SELECT_ARGS[@]}"
  fi
fi
validate_manifests
N_SELECTED=$(tail -n +2 "${TRIOS_TSV}" | wc -l | tr -d '[:space:]')
echo "  ${N_SELECTED} trios selected:"
head -n 6 "${TRIOS_TSV}" | cut -f1-4,7,8 | { column -t -s$'\t' 2>/dev/null || cat; } | sed 's/^/    /'
(( N_SELECTED > 5 )) && echo "    ..."

# ── 5. SV call set ───────────────────────────────────────────────────────────
if [[ "${USE_SV_CALLSET}" != "1" ]]; then
  echo "[5/6] SV call set disabled (USE_SV_CALLSET=0); trios run without --dnm-regions."
else
  if [[ -s "${SV_VCF}" && -s "${SV_VCF}.tbi" ]]; then
    echo "[5/6] SV call set already present: ${SV_VCF}"
  else
    echo "[5/6] Downloading the 3,202-sample SV call set (~1.7 GB): ${SV_VCF_URL}"
    fetch "${SV_VCF_URL}" "${SV_VCF}"
    fetch "${SV_VCF_URL}.tbi" "${SV_VCF}.tbi"
  fi
  # Subset once to the cohort's samples, so each trio job reads a small file.
  # Redone whenever the manifest is newer than the subset.
  if [[ -s "${SV_COHORT_VCF}" && ! "${TRIO_FILES_TSV}" -nt "${SV_COHORT_VCF}" ]]; then
    echo "  Cohort subset up to date: ${SV_COHORT_VCF}"
  else
    echo "  Subsetting the call set to the cohort's $(( 3 * N_SELECTED )) samples..."
    SAMPLES_FILE="${STATE_DIR}/cohort_samples.txt"
    tail -n +2 "${TRIO_FILES_TSV}" | cut -f3 | sort -u > "${SAMPLES_FILE}"
    if command -v bcftools >/dev/null 2>&1; then
      BCFTOOLS=(bcftools)
    else
      RUNTIME="${RUNTIME:-$(find_container_runtime)}"
      BCFTOOLS=("${RUNTIME}" exec --bind "$(bind_list "${WORK_DIR}" "$(dirname "${SV_VCF}")")" "${SIF_IMAGE}" bcftools)
    fi
    # --force-samples: a sample absent from the call set is a warning, not a
    # stop; -c 1 drops sites where no cohort member carries the SV.
    "${BCFTOOLS[@]}" view -S "${SAMPLES_FILE}" --force-samples -c 1 -Oz -o "${SV_COHORT_VCF}.part" "${SV_VCF}"
    "${BCFTOOLS[@]}" index -t -f "${SV_COHORT_VCF}.part"
    mv "${SV_COHORT_VCF}.part.tbi" "${SV_COHORT_VCF}.tbi"
    mv "${SV_COHORT_VCF}.part" "${SV_COHORT_VCF}"
    echo "  Wrote ${SV_COHORT_VCF}"
  fi
fi

# ── 6. Reference Jellyfish index ─────────────────────────────────────────────
if [[ -s "${REF_JF}" ]]; then
  echo "[6/6] Reference index already present: ${REF_JF}"
elif [[ "${SETUP_LOCAL:-0}" == "1" ]]; then
  echo "[6/6] Building the reference index in place (SETUP_LOCAL=1; ~23 GB RAM, under an hour)..."
  write_config_snapshot "${CONFIG_SNAPSHOT}"
  bash "${SCRIPT_DIR}/00b_build_ref_index.sh" --config "${CONFIG_SNAPSHOT}"
elif command -v sbatch >/dev/null 2>&1; then
  echo "[6/6] Submitting the reference index build (${REF_INDEX_CPUS} CPUs, ${REF_INDEX_MEM}, ${REF_INDEX_TIME})..."
  # The job reads the effective configuration back from the snapshot, so it
  # sees the same paths whether or not the scheduler exports the environment.
  write_config_snapshot "${CONFIG_SNAPSHOT}"
  JOBID=$(submit_job \
    --job-name="${JOB_PREFIX}_refidx" \
    --output="${LOG_DIR}/refidx_%j.out" --error="${LOG_DIR}/refidx_%j.err" \
    --cpus-per-task="${REF_INDEX_CPUS}" --mem="${REF_INDEX_MEM}" --time="${REF_INDEX_TIME}" \
    "${SCRIPT_DIR}/00b_build_ref_index.sh" --config "${CONFIG_SNAPSHOT}")
  echo "  Job ${JOBID}; the manager refuses to start until ${REF_JF} exists."
else
  echo "[6/6] sbatch not found: build the reference index with SETUP_LOCAL=1 bash 00_setup.sh, or run 00b_build_ref_index.sh on a node with ~32 GB RAM."
  write_config_snapshot "${CONFIG_SNAPSHOT}"
fi

echo ""
echo "============================================================"
echo " Setup complete."
echo ""
echo " Next:"
echo "   1. Wait for the reference index: ls -l ${REF_JF}"
echo "   2. Stage the trios and run discovery:   bash 01_stage_and_discover.sh"
echo "   3. Summarise across trios (any time):    bash 02_summarize.sh"
echo "============================================================"
