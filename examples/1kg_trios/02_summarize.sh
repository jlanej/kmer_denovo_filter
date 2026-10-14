#!/usr/bin/env bash
# =============================================================================
# 02_summarize.sh — cross-trio summary of the discovery runs
# =============================================================================
#
# Reads every finished trio under RESULTS_DIR and writes the cohort tables
# (see summarize_trios.py) to COHORT_DIR. Pure Python, standard library only;
# runs on a login node in seconds and can be re-run while trios are still
# finishing.
#
# Usage:
#   bash 02_summarize.sh
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  CONFIG_FILE="${SCRIPT_DIR}/config.sh"
elif [[ -f "$(pwd)/config.sh" ]]; then
  CONFIG_FILE="$(pwd)/config.sh"
else
  echo "ERROR: Could not find config.sh." >&2
  exit 1
fi
source "${CONFIG_FILE}"

if [[ ! -s "${TRIOS_TSV}" ]]; then
  echo "ERROR: manifest missing: ${TRIOS_TSV} - run 00_setup.sh first." >&2
  exit 1
fi

if command -v python3 >/dev/null 2>&1; then
  PYTHON=(python3)
else
  RUNTIME="$(find_container_runtime)" || { echo "ERROR: neither python3 nor apptainer found." >&2; exit 1; }
  PYTHON=("${RUNTIME}" exec --bind "$(bind_list "${WORK_DIR}" "${RESULTS_DIR}")" "${SIF_IMAGE}" python3)
  SCRIPT_DIR="${CONTAINER_EXAMPLE_DIR}"
fi

mkdir -p "${COHORT_DIR}"
"${PYTHON[@]}" "${SCRIPT_DIR}/summarize_trios.py" \
  --trios "${TRIOS_TSV}" \
  --results-dir "${RESULTS_DIR}" \
  --out-dir "${COHORT_DIR}" \
  --recurrent-min-trios "${RECURRENT_MIN_TRIOS}" \
  --slack "${RECURRENT_SLACK_BP}"
