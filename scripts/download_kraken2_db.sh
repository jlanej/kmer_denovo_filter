#!/usr/bin/env bash
set -euo pipefail

# PrackenDB — curated Kraken2 database with one NCBI reference genome per
# species.  Contains bacteria, archaea, protists, fungi, the human genome,
# RefSeq viral genomes, and UniVec Core (as of October 7 2025).
# https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads
DEFAULT_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_NCBI_reference_20251007.tar.gz"

usage() {
    cat <<EOF
Usage:
  ./scripts/download_kraken2_db.sh --db /path/to/kraken_db [--url URL]

Downloads the pre-built PrackenDB Kraken2 database (NCBI reference
assemblies, one genome per species) from the Kraken2 project.

The database contains all NCBI reference assemblies of bacteria, archaea,
protists, and fungi (Oct 2025), plus the human genome, RefSeq viral
genomes, and UniVec Core.  Because it keeps a single reference genome per
species it is well suited for methods that count k-mers per species.

See: https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads

Options:
  --db PATH       Target Kraken2 database directory (required)
  --url URL       Override the download URL
                  (default: $DEFAULT_URL)
  -h, --help      Show this help

Examples:
  ./scripts/download_kraken2_db.sh --db /data/kraken2/prackendb

  docker run --rm \\
    -v "\$PWD:/work" \\
    ghcr.io/jlanej/kmer_denovo_filter:latest \\
    bash /work/scripts/download_kraken2_db.sh --db /work/kraken2_db
EOF
}

DB_PATH=""
URL="$DEFAULT_URL"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db)
            DB_PATH="${2:-}"
            shift 2
            ;;
        --url)
            URL="${2:-}"
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

# wget is used for the HTTP download.
if ! command -v wget >/dev/null 2>&1; then
    echo "Error: wget not found on PATH." >&2
    echo "Install wget (e.g. apt-get install wget) and retry." >&2
    exit 1
fi

mkdir -p "$DB_PATH"

TARBALL="$DB_PATH/kraken2_db.tar.gz"

echo "[kraken2-db] Downloading PrackenDB to: $DB_PATH"
echo "[kraken2-db] URL: $URL"

wget -O "$TARBALL" "$URL"

echo "[kraken2-db] Extracting database..."
tar -xzf "$TARBALL" -C "$DB_PATH"
rm -f "$TARBALL"

# Validate key database files expected by kraken2.
for req in "hash.k2d" "opts.k2d" "taxo.k2d"; do
    if [[ ! -f "$DB_PATH/$req" ]]; then
        echo "Error: missing required database file: $DB_PATH/$req" >&2
        exit 1
    fi
done

# taxonomy/nodes.dmp is used for lineage-aware bacterial classification.
# Pre-built databases may omit it; warn but do not fail—Kraken2Runner
# falls back to exact taxid==2 matching when it is absent.
if [[ ! -f "$DB_PATH/taxonomy/nodes.dmp" ]]; then
    echo "[kraken2-db] Warning: taxonomy/nodes.dmp not found." >&2
    echo "[kraken2-db] Lineage-aware bacterial classification will" >&2
    echo "[kraken2-db] fall back to exact taxid==2 matching." >&2
fi

echo "[kraken2-db] Complete."
