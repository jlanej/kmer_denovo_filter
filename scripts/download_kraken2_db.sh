#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./scripts/download_kraken2_db.sh --db /path/to/kraken_db [--threads N]

Downloads/builds a Kraken2 standard database (taxonomy + standard libraries).

Requires the modern ``k2`` wrapper (Kraken2 ≥ 2.17).  ``k2`` downloads via
HTTP with built-in retry/resume and does not require rsync or --use-ftp.
Older Kraken2 installations that only ship ``kraken2-build`` are not
supported; please upgrade to Kraken2 ≥ 2.17.

Options:
  --db PATH       Target Kraken2 database directory (required)
  --threads N     Threads for the build (default: nproc or 4)
  -h, --help      Show this help

Examples:
  ./scripts/download_kraken2_db.sh --db /data/kraken2/standard --threads 16

  docker run --rm \
    -v "$PWD:/work" \
    ghcr.io/jlanej/kmer_denovo_filter:latest \
    bash /work/scripts/download_kraken2_db.sh --db /work/kraken2_db --threads 16
EOF
}

DB_PATH=""
ENV_THREADS="${THREADS:-}"
THREADS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db)
            DB_PATH="${2:-}"
            shift 2
            ;;
        --threads)
            THREADS="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 2
            ;;
    esac
done

if [[ -z "$DB_PATH" ]]; then
    echo "Error: --db is required" >&2
    usage
    exit 2
fi

if [[ -z "$THREADS" ]]; then
    if [[ -n "$ENV_THREADS" ]]; then
        THREADS="$ENV_THREADS"
    else
        THREADS="$(nproc 2>/dev/null || echo 4)"
    fi
fi

# Require the modern k2 wrapper (Kraken2 >= 2.17).
if ! command -v k2 >/dev/null 2>&1; then
    echo "Error: k2 not found on PATH. Kraken2 >= 2.17 is required." >&2
    echo "Please upgrade Kraken2: https://github.com/DerrickWood/kraken2" >&2
    exit 1
fi

mkdir -p "$DB_PATH"

echo "[kraken2-db] Building standard Kraken2 database at: $DB_PATH"
echo "[kraken2-db] Threads: $THREADS"

# Modern Kraken2 (>= 2.17): k2 downloads via HTTP with built-in
# retry and resume.  No rsync or --use-ftp workaround needed.
echo "[kraken2-db] Using k2 wrapper (Kraken2 >= 2.17, HTTP downloads)."
k2 build --standard --db "$DB_PATH" --threads "$THREADS"

# Validate key files expected by Kraken2Runner lineage-aware matching.
for req in "hash.k2d" "opts.k2d" "taxo.k2d" "taxonomy/nodes.dmp"; do
    if [[ ! -f "$DB_PATH/$req" ]]; then
        echo "Error: missing required database file: $DB_PATH/$req" >&2
        exit 1
    fi
done

echo "[kraken2-db] Complete."
