#!/usr/bin/env bash
# =============================================================================
# 00b_build_ref_index.sh — Jellyfish index of the reference, built once
# =============================================================================
#
# kmer-discovery subtracts reference k-mers from the child's, using a
# canonical Jellyfish index of the reference FASTA. Building it takes about
# 23 GB of RAM and under an hour; every trio job then reuses it through
# --ref-jf instead of rebuilding it. Submitted by 00_setup.sh; safe to run
# by hand, and a no-op when the index already exists.
#
# Usage:
#   sbatch 00b_build_ref_index.sh [--config FILE | --work-dir DIR ...]
#   bash 00b_build_ref_index.sh          # in place, on a node with ~32 GB RAM
# =============================================================================

#SBATCH --job-name=1kgtrio_refidx
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00

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
reject_extra_args "$@"

if [[ -s "${REF_JF}" ]]; then
  echo "Reference index already exists: ${REF_JF}"
  exit 0
fi
if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: reference FASTA missing: ${REF_FASTA} - run 00_setup.sh first."
  exit 1
fi
if [[ ! -s "${SIF_IMAGE}" ]]; then
  echo "ERROR: container image missing: ${SIF_IMAGE} - run 00_setup.sh first."
  exit 1
fi
RUNTIME="$(find_container_runtime)" || { echo "ERROR: apptainer/singularity not found on PATH."; exit 1; }

THREADS="${SLURM_CPUS_PER_TASK:-${REF_INDEX_CPUS}}"
TMP_JF="${REF_JF}.tmp$$"
echo "Building ${REF_JF} (k=${KMER_SIZE}, -s ${REF_JF_HASH_SIZE}, ${THREADS} threads) from ${REF_FASTA}"
echo "Started: $(date)"
t0=$(date +%s)
# -C (canonical) is what kmer-discovery's own builder uses, and what its
# start-up check of an existing index requires.
"${RUNTIME}" exec \
  --bind "$(bind_list "$(dirname "${REF_FASTA}")" "$(dirname "${REF_JF}")")" \
  "${SIF_IMAGE}" \
  jellyfish count -m "${KMER_SIZE}" -s "${REF_JF_HASH_SIZE}" -t "${THREADS}" -C \
    -o "${TMP_JF}" "${REF_FASTA}"
if [[ ! -s "${TMP_JF}" ]]; then
  echo "ERROR: jellyfish wrote no index."
  rm -f "${TMP_JF}"*
  exit 1
fi
# A table that filled would have been dumped as several files; the single
# -o path is only complete when none of them exists.
if ls "${TMP_JF}"_* >/dev/null 2>&1; then
  echo "ERROR: jellyfish split the index (${TMP_JF}_*): the hash filled. Raise REF_JF_HASH_SIZE and retry."
  rm -f "${TMP_JF}"*
  exit 1
fi
mv "${TMP_JF}" "${REF_JF}"
echo "Done in $(( $(date +%s) - t0 )) s: ${REF_JF} ($(( $(file_size "${REF_JF}") / 1024 / 1024 / 1024 )) GB)"
