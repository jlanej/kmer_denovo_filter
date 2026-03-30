"""Jellyfish subprocess wrappers.

All functions shell out to the ``jellyfish`` binary.  None of them
accept an ``argparse.Namespace``; callers pass explicit, typed
parameters instead.
"""

import glob
import logging
import os
import subprocess
import tempfile
import time

from kmer_denovo_filter.core.memory_utils import (
    _log_memory,
    _log_subprocess_memory,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Formatting helpers (local copies to avoid circular imports)
# ---------------------------------------------------------------------------

def _format_elapsed(seconds):
    """Format elapsed seconds as a human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = seconds % 60
    return f"{hours}h {minutes}m {secs:.0f}s"


def _format_file_size(path):
    """Return human-readable file size, or '?' if unavailable."""
    try:
        size = os.path.getsize(path)
    except OSError:
        return "?"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"


# ---------------------------------------------------------------------------
# Jellyfish file helpers
# ---------------------------------------------------------------------------


def _find_jf_files(base_path):
    """Find all Jellyfish output files for a given base path.

    When Jellyfish's hash table overflows, it produces numbered files
    (``base.jf_0``, ``base.jf_1``, …) rather than a single ``base.jf``.
    Returns a sorted list of existing ``.jf`` / ``.jf_N`` files.
    """
    files = []
    if os.path.exists(base_path):
        files.append(base_path)
    files.extend(sorted(glob.glob(base_path + "_[0-9]*")))
    return files


def _estimate_jf_hash_size(bam_path, kmer_size, default="1G"):
    """Estimate a reasonable ``-s`` hash size for Jellyfish count.

    Uses the BAM file size as a rough proxy for the number of k-mers
    (approximately ``file_size / 4`` for compressed BAM → approximate
    bases, divided by reasonable overhead).  The result is capped to
    keep initial memory reasonable.

    Returns a string suitable for Jellyfish's ``-s`` argument
    (e.g. ``"2G"``, ``"500M"``).
    """
    try:
        file_size = os.path.getsize(bam_path)
    except OSError:
        return default

    # Rough heuristic: BAM compresses ~3:1, so uncompressed bases ≈
    # file_size * 3.  Distinct k-mers are typically much less than
    # total bases.  Jellyfish hash entry is ~(k+8) bytes.
    # Aim for hash that can hold ~75% of distinct k-mers without
    # resizing while staying within memory budget.
    #
    # Conservative estimate: 1 entry per 10 uncompressed bases,
    # capped between 100M and 4G.
    estimated_bases = file_size * 3
    estimated_distinct = estimated_bases // 10

    # Cap between 100M and 4G entries (not bytes)
    min_entries = 100_000_000   # 100M
    max_entries = 4_000_000_000  # 4G
    entries = max(min_entries, min(estimated_distinct, max_entries))

    if entries >= 1_000_000_000:
        return f"{entries // 1_000_000_000}G"
    return f"{entries // 1_000_000}M"


# ---------------------------------------------------------------------------
# Jellyfish subprocess wrappers
# ---------------------------------------------------------------------------


def _scan_parent_jellyfish(
    parent_bam, ref_fasta, kmer_fasta, kmer_size, parent_dir, threads=4,
    n_filter_kmers=None,
):
    """Scan a parent BAM and find which child k-mers are present.

    Uses ``jellyfish count`` with the ``--if`` filter so only child k-mers
    are tracked while the parent BAM is streamed.  The dump output is
    streamed line-by-line and the parent Jellyfish index is removed after
    the dump to free disk space.

    The hash size is set to ``2 × n_filter_kmers`` (clamped to a
    10 M minimum) so that Jellyfish has enough room for all filtered
    k-mers.  When *n_filter_kmers* is not provided, the entries in
    *kmer_fasta* are counted.  If Jellyfish still produces multiple
    chunk files (hash overflow), they are merged automatically.

    Args:
        parent_bam: Path to parent BAM/CRAM.
        ref_fasta: Path to reference FASTA (or None).
        kmer_fasta: FASTA file of k-mers to track (``--if`` filter).
        kmer_size: K-mer length.
        parent_dir: Working directory for the index.
        threads: Number of threads for jellyfish count.
        n_filter_kmers: Number of k-mers in *kmer_fasta*.  When ``None``
            this is estimated by counting entries in the file.

    Returns:
        Dict mapping canonical k-mer string to its count in the parent.
    """
    os.makedirs(parent_dir, exist_ok=True)
    jf_output = os.path.join(parent_dir, "parent.jf")

    # Size hash to fit the filter k-mers without overflow.
    if n_filter_kmers is None:
        n_filter_kmers = 0
        with open(kmer_fasta) as fh:
            for line in fh:
                if line.rstrip() and not line.startswith(">"):
                    n_filter_kmers += 1
    hash_size = max(n_filter_kmers * 2, 10_000_000)
    hash_size_str = str(hash_size)

    samtools_threads = max(1, threads // 4)
    samtools_cmd = [
        "samtools", "fasta", "-F", "0xD00",
        "-@", str(samtools_threads),
        parent_bam,
    ]
    if ref_fasta:
        samtools_cmd.extend(["--reference", ref_fasta])

    jellyfish_cmd = [
        "jellyfish", "count",
        "-m", str(kmer_size),
        "-s", hash_size_str,
        "-t", str(threads),
        "-C",
        "--if", kmer_fasta,
        "-o", jf_output,
        "/dev/fd/0",
    ]

    bam_size = _format_file_size(parent_bam)
    logger.info(
        "Scanning parent BAM (%s): %s", bam_size, parent_bam,
    )
    logger.info(
        "  samtools fasta → jellyfish count (k=%d, threads=%d, hash=%s)",
        kmer_size, threads, hash_size_str,
    )

    scan_start = time.monotonic()

    p_samtools = subprocess.Popen(
        samtools_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    p_jellyfish = subprocess.Popen(
        jellyfish_cmd,
        stdin=p_samtools.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    p_samtools.stdout.close()

    # Poll for completion, logging periodic progress
    poll_interval = 30  # seconds between progress updates
    last_log = scan_start
    while True:
        try:
            p_jellyfish.wait(timeout=poll_interval)
            break  # process finished
        except subprocess.TimeoutExpired:
            now = time.monotonic()
            elapsed = now - scan_start
            # Report jellyfish output file size as a proxy for progress
            jf_files_progress = _find_jf_files(jf_output)
            if jf_files_progress:
                total_size = sum(
                    os.path.getsize(f) for f in jf_files_progress
                    if os.path.exists(f)
                )
                jf_size = f"{total_size / (1024**3):.1f} GB" if total_size >= 1024**3 else (
                    _format_file_size(jf_files_progress[0])
                    if len(jf_files_progress) == 1 else
                    f"{total_size / (1024**2):.1f} MB"
                )
            else:
                jf_size = "pending"
            logger.info(
                "  … still scanning (%s elapsed, jf index: %s)",
                _format_elapsed(elapsed), jf_size,
            )
            _log_memory("parent scanning")
            _log_subprocess_memory(p_jellyfish, "jellyfish-count")
            _log_subprocess_memory(p_samtools, "samtools-fasta")
            last_log = now

    p_samtools.communicate()
    jf_stderr = p_jellyfish.stderr.read()

    if p_jellyfish.returncode != 0:
        raise RuntimeError(
            f"jellyfish count failed: {jf_stderr.decode()}"
        )

    # Handle multi-file output (hash overflow) — merge if needed
    jf_files = _find_jf_files(jf_output)
    if len(jf_files) > 1:
        merged_path = os.path.join(parent_dir, "parent_merged.jf")
        jf_output = _merge_jf_files(jf_files, merged_path)
    elif jf_files and jf_files[0] != jf_output:
        # Single numbered file (e.g. parent.jf_0) — rename to expected name
        os.rename(jf_files[0], jf_output)

    scan_elapsed = time.monotonic() - scan_start
    logger.info(
        "  Jellyfish counting complete (%s)", _format_elapsed(scan_elapsed),
    )

    found_kmers = {}
    if os.path.exists(jf_output):
        jf_size = _format_file_size(jf_output)
        logger.info("  Dumping jellyfish results (%s index)…", jf_size)
        dump_cmd = ["jellyfish", "dump", "-c", "-L", "1", jf_output]
        with tempfile.TemporaryFile(mode="w+") as stderr_f:
            p_dump = subprocess.Popen(
                dump_cmd, stdout=subprocess.PIPE, stderr=stderr_f,
                text=True,
            )
            for line in p_dump.stdout:
                line = line.rstrip("\n")
                if line:
                    parts = line.split()
                    found_kmers[parts[0]] = int(parts[1])
            p_dump.wait()
            if p_dump.returncode != 0:
                stderr_f.seek(0)
                raise RuntimeError(
                    f"jellyfish dump (parent) failed: {stderr_f.read()}"
                )

        # Remove the parent jellyfish index to free disk/cache.
        os.remove(jf_output)

    return found_kmers


def _ensure_ref_jf(ref_fasta, kmer_size, threads, ref_jf=None):
    """Ensure a Jellyfish reference index exists, building it if necessary.

    Args:
        ref_fasta: Path to the reference FASTA file.
        kmer_size: K-mer size.
        threads: Number of threads for jellyfish.
        ref_jf: Explicit path to the Jellyfish index; when *None*,
            defaults to ``{ref_fasta}.k{kmer_size}.jf``.

    Returns:
        Path to the Jellyfish reference index.
    """
    if ref_jf is None:
        ref_jf = f"{ref_fasta}.k{kmer_size}.jf"

    if os.path.isfile(ref_jf):
        logger.info("Reference Jellyfish index found: %s", ref_jf)
        return ref_jf

    logger.info(
        "Building reference Jellyfish index: %s (k=%d, threads=%d)",
        ref_jf, kmer_size, threads,
    )
    ref_hash_size = _estimate_jf_hash_size(ref_fasta, kmer_size, default="3G")
    logger.info("  Reference JF hash size: %s", ref_hash_size)
    build_start = time.monotonic()
    cmd = [
        "jellyfish", "count",
        "-m", str(kmer_size),
        "-s", ref_hash_size,
        "-t", str(threads),
        "-C",
        ref_fasta,
        "-o", ref_jf,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"jellyfish count (reference) failed: {result.stderr}"
        )
    logger.info(
        "Reference index built in %s (%s)",
        _format_elapsed(time.monotonic() - build_start),
        _format_file_size(ref_jf),
    )
    return ref_jf


def _merge_jf_files(jf_files, merged_path, threads=4):
    """Merge multiple Jellyfish chunk files into one.

    When Jellyfish count produces multiple output files (hash overflow),
    they must be merged before querying.  Uses ``jellyfish merge`` which
    streams chunks and requires memory proportional to one chunk at a
    time.
    """
    if len(jf_files) <= 1:
        return jf_files[0] if jf_files else None

    logger.info(
        "Merging %d Jellyfish chunks into %s…",
        len(jf_files), merged_path,
    )
    merge_start = time.monotonic()
    cmd = ["jellyfish", "merge", "-o", merged_path] + jf_files
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"jellyfish merge failed: {result.stderr}")

    # Remove chunk files to free disk
    for f in jf_files:
        if f != merged_path and os.path.exists(f):
            os.remove(f)

    logger.info(
        "Jellyfish merge complete (%s, merged: %s)",
        _format_elapsed(time.monotonic() - merge_start),
        _format_file_size(merged_path),
    )
    return merged_path


def _build_proband_jf_index(proband_unique_fa, kmer_size, tmpdir,
                            n_proband_unique=None):
    """Build a Jellyfish index from proband-unique k-mers.

    Creates a ``.jf`` hash file that can be queried via
    ``jellyfish query`` or :class:`JellyfishKmerQuery` to check k-mer
    membership without loading the set into Python memory.

    For WGS-scale data with hundreds of millions of proband-unique k-mers,
    the resulting index is typically 2–10 GB on disk and memory-mapped by
    each ``jellyfish query`` subprocess.  Multiple workers share the same
    OS page-cache mapping, so N workers ≈ 1× the hash file memory.

    Args:
        proband_unique_fa: Path to FASTA of proband-unique k-mers.
        kmer_size: K-mer length.
        tmpdir: Working directory for the index.
        n_proband_unique: Number of k-mers (for hash size estimation).

    Returns:
        Path to the Jellyfish index file.
    """
    if n_proband_unique is None:
        n_proband_unique = 0
        with open(proband_unique_fa) as fh:
            for line in fh:
                if line.rstrip() and not line.startswith(">"):
                    n_proband_unique += 1

    # Set hash size slightly larger than the k-mer count to avoid
    # overflow / multi-file output.
    hash_size = max(n_proband_unique * 2, 1_000_000)
    hash_size_str = f"{hash_size}"

    proband_jf = os.path.join(tmpdir, "proband_unique.jf")

    logger.info(
        "Building Jellyfish index from %d proband-unique k-mers "
        "(hash size: %s)…",
        n_proband_unique, hash_size_str,
    )

    jf_cmd = [
        "jellyfish", "count",
        "-m", str(kmer_size),
        "-s", hash_size_str,
        "-t", "1",
        "-C",
        "-o", proband_jf,
        proband_unique_fa,
    ]

    build_start = time.monotonic()
    result = subprocess.run(
        jf_cmd, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"jellyfish count (proband index) failed: {result.stderr}"
        )

    logger.info(
        "Proband Jellyfish index built (%s, index: %s)",
        _format_elapsed(time.monotonic() - build_start),
        _format_file_size(proband_jf),
    )
    _log_memory("after proband index build")
    return proband_jf
