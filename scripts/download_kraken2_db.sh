#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./scripts/download_kraken2_db.sh --db /path/to/kraken_db [--threads N]

Downloads/builds a Kraken2 standard database (taxonomy + standard libraries)
using kraken2-build best-practice defaults.

Options:
  --db PATH       Target Kraken2 database directory (required)
  --threads N     Threads for kraken2-build (default: nproc or 4)
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
THREADS="${THREADS:-}"

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
    if command -v nproc >/dev/null 2>&1; then
        THREADS="$(nproc)"
    else
        THREADS="4"
    fi
fi

if ! command -v kraken2-build >/dev/null 2>&1; then
    echo "Error: kraken2-build not found on PATH" >&2
    exit 1
fi

mkdir -p "$DB_PATH"

echo "[kraken2-db] Building standard Kraken2 database at: $DB_PATH"
echo "[kraken2-db] Threads: $THREADS"
kraken2-build --standard --db "$DB_PATH" --threads "$THREADS"

# Validate key files expected by Kraken2Runner lineage-aware matching.
for req in "hash.k2d" "opts.k2d" "taxo.k2d" "taxonomy/nodes.dmp"; do
    if [[ ! -f "$DB_PATH/$req" ]]; then
        echo "Error: expected database file missing: $DB_PATH/$req" >&2
        exit 1
    fi
done

echo "[kraken2-db] Complete."
