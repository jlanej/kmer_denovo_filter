#!/usr/bin/env bash
# =============================================================================
# config.local.example.sh — a site configuration for the 1000 Genomes trio screen
# =============================================================================
#
# Copy this file to config.local.sh beside config.sh (it is gitignored) and
# edit it; every script sources it before applying the defaults in
# config.sh. Point KMER_TRIOS_CONFIG at it, or pass --config FILE, to keep
# it somewhere else. Use the  VAR="${VAR:-value}"  form so that flags and
# exported variables still take precedence over this file.
#
# Only the settings you change belong here; see config.sh for the full list.
# =============================================================================

# Where everything lives: scratch with a few hundred GB free.
WORK_DIR="${WORK_DIR:-/scratch/${USER:-$(id -un)}/1kg_trios}"

# Reuse a reference and index that already exist on the site, if any.
# REF_FASTA="${REF_FASTA:-/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa}"
# REF_JF="${REF_JF:-${REF_FASTA}.k31.jf}"

# A container image pulled earlier, instead of pulling from ghcr.io.
# SIF_IMAGE="${SIF_IMAGE:-/path/to/kmer_denovo_filter.sif}"

# Scheduler specifics, applied to every sbatch: partition, account, QoS.
SBATCH_EXTRA_ARGS="${SBATCH_EXTRA_ARGS:---partition=compute --account=mylab}"

# Environment modules the scripts need on this site (space-separated).
MODULES_TO_LOAD="${MODULES_TO_LOAD:-apptainer}"

# Resources for the per-trio discovery job, tuned from observed runs.
# DISCOVERY_CPUS="${DISCOVERY_CPUS:-16}"
# DISCOVERY_MEM="${DISCOVERY_MEM:-128G}"
# DISCOVERY_TIME="${DISCOVERY_TIME:-36:00:00}"

# Node-local scratch for discovery intermediates, only if it is large enough
# (about 40 GB per running 30x trio: the child's ~29 GB k-mer index while it
# is dumped, plus the non-reference parts and their concatenation) and not a
# RAM-backed filesystem.
# DISCOVERY_TMP_DIR="${DISCOVERY_TMP_DIR:-/local/scratch/${USER:-$(id -un)}}"

# A Kraken2 database for the non-human fraction screen of each run's
# informative reads (scripts/download_kraken2_db.sh fetches PrackenDB).
# Empty = no screen. The database loads into the job's memory.
# KRAKEN2_DB="${KRAKEN2_DB:-/path/to/kraken2_db}"
