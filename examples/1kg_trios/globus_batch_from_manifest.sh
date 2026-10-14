#!/usr/bin/env bash
# =============================================================================
# globus_batch_from_manifest.sh — Globus batch files for the selected trios
# =============================================================================
#
# When direct downloads crawl from your site, Globus is usually the fastest
# path: EBI exposes ftp.sra.ebi.ac.uk through the "EMBL-EBI Public Data"
# collection with the FTP tree mounted under /vol1, so a manifest URL becomes
# a Globus source path by stripping the host. Each emitted line renames the
# file on arrival to the <SAMPLE>.cram / <SAMPLE>.cram.crai names the
# pipeline expects. Indexes and CRAMs go in separate batch files because
# Globus orders files within a task itself.
#
# Files staged this way are treated by 01_stage_and_discover.sh as already
# present: it verifies the index MD5 and dispatches the trio without
# downloading. Run the manager after the transfer completes.
#
# Usage:
#   bash globus_batch_from_manifest.sh [output_prefix]     # default: globus_batch
#   globus transfer 47772002-3e5b-4fd3-b97c-18cee38d6df2 <YOUR_COLLECTION_UUID> \
#     --batch globus_batch_crai.txt --sync-level checksum --label "1kg trios crai"
#   globus transfer 47772002-3e5b-4fd3-b97c-18cee38d6df2 <YOUR_COLLECTION_UUID> \
#     --batch globus_batch_cram.txt --sync-level checksum --label "1kg trios cram"
#
# Find your site's collection with: globus endpoint search "<your site>"
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

OUT_PREFIX="${1:-globus_batch}"
CRAI_BATCH="${OUT_PREFIX}_crai.txt"
CRAM_BATCH="${OUT_PREFIX}_cram.txt"

if [[ ! -s "${TRIO_FILES_TSV}" ]]; then
  echo "ERROR: manifest missing or empty: ${TRIO_FILES_TSV} - run 00_setup.sh first." >&2
  exit 1
fi

EMITTED=0
SKIPPED=0
SEEN_FILE="$(mktemp)"
trap 'rm -f "${SEEN_FILE}"' EXIT
: > "${CRAI_BATCH}"
: > "${CRAM_BATCH}"
while IFS=$'\t' read -r TRIO ROLE SAMPLE CRAM_URL CRAI_URL _MD5 _BATCH; do
  [[ -z "${TRIO}" ]] && continue
  if trio_done "${TRIO}"; then
    SKIPPED=$(( SKIPPED + 1 ))
    continue
  fi
  # A parent shared by two selected trios is transferred once.
  if grep -qxF "${SAMPLE}" "${SEEN_FILE}"; then
    continue
  fi
  echo "${SAMPLE}" >> "${SEEN_FILE}"
  if [[ "${CRAM_URL}" != "${EXPECTED_FTP_PREFIX}"* ]]; then
    echo "WARNING: skipping ${SAMPLE} (${ROLE} of ${TRIO}) - not an ENA URL: ${CRAM_URL}" >&2
    SKIPPED=$(( SKIPPED + 1 ))
    continue
  fi
  printf '"%s" "%s"\n' "${CRAI_URL#ftp://ftp.sra.ebi.ac.uk}" "${CRAM_DIR}/${SAMPLE}.cram.crai" >> "${CRAI_BATCH}"
  printf '"%s" "%s"\n' "${CRAM_URL#ftp://ftp.sra.ebi.ac.uk}" "${CRAM_DIR}/${SAMPLE}.cram" >> "${CRAM_BATCH}"
  EMITTED=$(( EMITTED + 1 ))
done < <(tail -n +2 "${TRIO_FILES_TSV}")

echo "Emitted ${EMITTED} samples: ${CRAI_BATCH} (indexes) and ${CRAM_BATCH} (CRAMs); skipped ${SKIPPED} rows of finished trios or bad URLs." >&2
echo "Submit the CRAI batch first, both with --sync-level checksum (see the header of this script)." >&2
