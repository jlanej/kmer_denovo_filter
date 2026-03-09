#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./scripts/download_kraken2_db.sh --db /path/to/kraken_db [--threads N]

Downloads/builds a Kraken2 standard database (taxonomy + standard libraries).

Prefers the modern ``k2`` wrapper (Kraken2 ≥ 2.17) when available; falls
back to ``kraken2-build`` for older installations.  In both cases
KRAKEN2_USE_FTP=1 is exported so that Kraken2's internal download scripts
use wget/FTP instead of rsync (rsync is not required).

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

# Export KRAKEN2_USE_FTP=1 so that Kraken2's internal download scripts
# (download_taxonomy.sh, download_genomic_library.sh, etc.) use
# wget/FTP instead of rsync.  This is the belt-and-suspenders fix for
# environments where rsync is not available.
export KRAKEN2_USE_FTP=1

mkdir -p "$DB_PATH"

echo "[kraken2-db] Building standard Kraken2 database at: $DB_PATH"
echo "[kraken2-db] Threads: $THREADS"

if command -v k2 >/dev/null 2>&1; then
    # Modern Kraken2 (>= 2.17): k2 downloads via HTTP with built-in
    # retry and resume.  KRAKEN2_USE_FTP=1 is still exported above as
    # an extra safeguard in case any internal helper falls back to rsync.
    echo "[kraken2-db] Using k2 wrapper (Kraken2 >= 2.17, HTTP downloads)."
    k2 build --standard --db "$DB_PATH" --threads "$THREADS"
elif command -v kraken2-build >/dev/null 2>&1; then
    # Legacy Kraken2 (< 2.17): use --use-ftp so wget is preferred over
    # rsync for all downloads.  KRAKEN2_USE_FTP=1 (exported above) is
    # the belt-and-suspenders fix honoured by internal helper scripts.
    echo "[kraken2-db] k2 not found; falling back to kraken2-build --use-ftp."
    kraken2-build --standard --db "$DB_PATH" --threads "$THREADS" --use-ftp
else
    echo "Error: neither k2 nor kraken2-build found on PATH." >&2
    echo "Please install Kraken2: https://github.com/DerrickWood/kraken2" >&2
    exit 1
fi

# Validate key files expected by Kraken2Runner lineage-aware matching.
for req in "hash.k2d" "opts.k2d" "taxo.k2d" "taxonomy/nodes.dmp"; do
    if [[ ! -f "$DB_PATH/$req" ]]; then
        echo "Error: missing required database file: $DB_PATH/$req" >&2
        exit 1
    fi
done

echo "[kraken2-db] Complete."
