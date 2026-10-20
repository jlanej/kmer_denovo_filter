"""Jellyfish subprocess wrappers.

All functions shell out to the ``jellyfish`` binary.  None of them
accept an ``argparse.Namespace``; callers pass explicit, typed
parameters instead.
"""

import glob
import json
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


# ---------------------------------------------------------------------------
# Sizing the child's k-mer count (Module 1)
# ---------------------------------------------------------------------------
#
# Jellyfish 2 facts the sizing rests on (sub_commands/count_main.cc):
#   * ``-s`` is the INITIAL hash size. Without ``--disk`` a full table is
#     doubled in memory - the old and the new table coexist during the copy -
#     so a 4G table that needs 32G passes through 23, 68, 45, 133, 88 and then
#     264 GB of RSS: that is how every WGS child used to be OOM-killed.  With
#     ``--disk`` (``do_size_doubling(false)``) a full table is dumped to
#     ``<out>_N`` files and jellyfish merges them into ``<out>`` at the end.
#   * ``--bf-size N`` puts a Bloom filter in front of the table: a k-mer is
#     inserted only once the filter has seen it, so k-mers seen once never
#     enter the hash and every count in the hash is ONE LESS than the true
#     count.  The caller dumps with ``-L (min_count - 1)`` to compensate.
#
# Memory (observed on the jobs that died, k=31: 23.1 GB at 2^32 entries, 45.2
# at 2^33, 88.0 at 2^34): a hash entry takes about (2k - log2(size) + 16) / 8
# bytes; a Bloom filter at the default 1% false-positive rate about 1.2 bytes
# per element.

GENOME_DISTINCT_KMERS = 3.1e9   # distinct k-mers of a human genome, errors aside
READ_ERROR_RATE = 0.004         # a sequencing error makes up to k novel k-mers
BAM_COMPRESSION = 3.0           # uncompressed bases per BAM byte (rough)
CRAM_COMPRESSION = 4.5          # ... per CRAM byte (reference-based compression)
HASH_LOAD = 0.8                 # the share of a jellyfish table used before it fills
HASH_MIN_ENTRIES = 100_000_000
BLOOM_MIN_ELEMENTS = 10_000_000
BLOOM_FP = 0.01
MEMORY_SHARE = 0.6              # the share of the memory budget the count may use (filter + table)


FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def _is_fasta(path):
    return str(path).lower().endswith(FASTA_SUFFIXES)


def _estimate_bases(bam_path):
    """Bases in the file, from its size and a compression ratio (0 when the size is unknown): reads in a BAM or CRAM,
    the sequence itself in a FASTA (a reference, uncompressed or gzipped)."""
    try:
        size = os.path.getsize(bam_path)
    except OSError:
        return 0
    name = str(bam_path).lower()
    if _is_fasta(name):
        ratio = 4.0 if name.endswith(".gz") else 1.0
    else:
        ratio = CRAM_COMPRESSION if name.endswith(".cram") else BAM_COMPRESSION
    return int(size * ratio)


def _genome_kmers(bases):
    """The error-free distinct k-mers the reads can hold: the genome's, or fewer bases than that for a small input."""
    return min(GENOME_DISTINCT_KMERS, float(bases))


def _estimate_distinct_kmers(bam_path, kmer_size):
    """Distinct k-mers the child's reads hold, errors included: the genome's (or the input's, when smaller) plus one
    novel k-mer per error per position of the window (bases x error rate x k)."""
    bases = _estimate_bases(bam_path)
    if _is_fasta(bam_path):                                            # a reference: no sequencing errors
        return int(_genome_kmers(bases))
    return int(_genome_kmers(bases) + bases * READ_ERROR_RATE * kmer_size)


def _jf_bytes_per_entry(kmer_size, entries):
    """Bytes per hash entry: the key bits jellyfish stores (2k minus the bits the position encodes), plus the
    reprobe and counter fields (about 14 bits)."""
    import math
    bits = 2 * kmer_size - math.log2(max(entries, 2)) + 16
    return max(bits, 16) / 8.0


def _bloom_bytes_per_element(fp=BLOOM_FP):
    import math
    return -math.log(fp) / (math.log(2) ** 2) / 8.0


def _parse_size(text):
    """'4G', '500M', '1k' or a plain integer -> entries (binary suffixes, as jellyfish reads them)."""
    t = str(text).strip()
    mult = {"k": 2**10, "m": 2**20, "g": 2**30, "t": 2**40}
    if t and t[-1].lower() in mult:
        return int(float(t[:-1]) * mult[t[-1].lower()])
    return int(float(t))


def _pow2_at_most(n):
    import math
    return 1 << max(int(math.floor(math.log2(max(n, 1)))), 0)


def _pow2_at_least(n):
    import math
    return 1 << max(int(math.ceil(math.log2(max(n, 1)))), 0)


def size_child_count(bam_path, kmer_size, memory_gb=None, bloom_filter=True, jf_hash_size=None):
    """Size the child count against a memory budget: the Bloom filter from the input (what the reads hold),
    the hash from what is left of the budget, never more than the input needs.

    Returns a dict: hash_entries, bloom_elements (0 when the filter is off), the estimated GB of each, the
    budget used, est_distinct, and a one-line ``note`` for the log.  memory_gb None: the system's total
    memory (``_get_available_memory_gb``), else 64 GB when that is unknown too.  jf_hash_size (a string
    such as '4G') overrides the hash size; the filter is still sized from the input."""
    from kmer_denovo_filter.core.memory_utils import _get_available_memory_gb
    import math
    if memory_gb is None:
        total, _avail = _get_available_memory_gb()
        memory_gb = total if total else 64.0
    budget = memory_gb * MEMORY_SHARE
    bases = _estimate_bases(bam_path)
    genome = _genome_kmers(bases)
    est_distinct = _estimate_distinct_kmers(bam_path, kmer_size)
    bloom_elements = 0
    bloom_gb = 0.0
    if bloom_filter:
        bloom_elements = max(BLOOM_MIN_ELEMENTS, int(est_distinct * 1.2))
        bloom_gb = bloom_elements * _bloom_bytes_per_element() / 1024**3
        # what reaches the hash: the genome's k-mers (and the variants), the errors that recur, the filter's false positives
        need = genome * 1.3 + (est_distinct - genome) * 0.05 + est_distinct * BLOOM_FP
    else:
        need = float(est_distinct)
    if jf_hash_size:
        hash_entries = _parse_size(jf_hash_size)
    else:
        # jellyfish rounds -s up to a power of two: ask for the power of two that holds the need at HASH_LOAD, unless
        # that exceeds what the budget leaves beside the filter - then the largest power of two the budget holds
        left_gb = max(budget - bloom_gb, 0.0)
        wanted = _pow2_at_least(max(int(need / HASH_LOAD), HASH_MIN_ENTRIES))
        cap = int(left_gb * 1024**3 / _jf_bytes_per_entry(kmer_size, wanted))
        hash_entries = wanted if wanted <= cap else max(_pow2_at_most(cap), HASH_MIN_ENTRIES)
    hash_gb = hash_entries * _jf_bytes_per_entry(kmer_size, hash_entries) / 1024**3
    overflow = need > hash_entries * HASH_LOAD
    note = ("hash %s entries (~%.1f GB)%s; %s; budget %.0f GB of %.0f; input ~%.2g distinct k-mers%s" % (
        f"{hash_entries:,}", hash_gb, " (given)" if jf_hash_size else "",
        ("Bloom filter %s elements (~%.1f GB)" % (f"{bloom_elements:,}", bloom_gb)) if bloom_filter else "no Bloom filter",
        budget, memory_gb, est_distinct, "; the table will fill and dump to disk (merged at the end)" if overflow else ""))
    return dict(hash_entries=int(hash_entries), hash_gb=hash_gb, bloom_elements=int(bloom_elements), bloom_gb=bloom_gb, budget_gb=budget,
                memory_gb=memory_gb, est_distinct=est_distinct, overflow=overflow, note=note)


def _estimate_jf_hash_size(bam_path, kmer_size, default="1G", memory_gb=None):
    """The ``-s`` for the child count as a string of entries, sized by :func:`size_child_count` (kept for callers
    of the old name; the file-size heuristic with its 4G cap is gone - it sent every WGS child into jellyfish's
    in-memory doubling)."""
    try:
        return str(size_child_count(bam_path, kmer_size, memory_gb=memory_gb, bloom_filter=False)["hash_entries"])
    except OSError:
        return default


# ---------------------------------------------------------------------------
# Jellyfish subprocess wrappers
# ---------------------------------------------------------------------------


def _stderr_tail(fh, max_lines=20):
    """Return the last *max_lines* lines written to a stderr temp file."""
    fh.seek(0)
    lines = fh.read().decode(errors="replace").strip().splitlines()
    return "\n".join(lines[-max_lines:])


def _run_samtools_jellyfish(samtools_cmd, jellyfish_cmd, label,
                            poll_interval=30, on_poll=None):
    """Run ``samtools fasta | jellyfish count``; raise if either side fails.

    Checking jellyfish's exit status alone is not enough: when samtools
    dies part-way through (truncated BAM, CRAM reference problems),
    jellyfish sees an early EOF and exits cleanly having counted only the
    reads it received.  samtools' own exit status is not enough either:
    a file cut at a block boundary reads cleanly, and some versions (e.g.
    1.16 with ``-@ 1``) exit 0 part-way through a truncated file.  htslib
    does warn that the EOF marker is absent in all of these cases, so
    that warning is treated as fatal too.

    Both stderr streams go to temporary files rather than pipes, so a
    verbose process cannot fill an unread pipe buffer and stall the
    pipeline.

    Args:
        samtools_cmd: ``samtools fasta`` command writing reads to stdout.
        jellyfish_cmd: ``jellyfish count`` command reading from stdin.
        label: Input description for error messages
            (e.g. ``"Mother: mother.bam"``).
        poll_interval: Seconds between *on_poll* calls.
        on_poll: Optional progress callback, called as
            ``on_poll(elapsed_seconds, p_samtools, p_jellyfish)`` every
            *poll_interval* seconds while jellyfish is running.

    Raises:
        RuntimeError: If samtools or jellyfish exits with a non-zero
            status, or samtools reports that its input is truncated.
    """
    start = time.monotonic()
    with tempfile.TemporaryFile() as samtools_err, \
            tempfile.TemporaryFile() as jellyfish_err:
        p_samtools = subprocess.Popen(
            samtools_cmd, stdout=subprocess.PIPE, stderr=samtools_err,
        )
        p_jellyfish = subprocess.Popen(
            jellyfish_cmd, stdin=p_samtools.stdout,
            stdout=subprocess.DEVNULL, stderr=jellyfish_err,
        )
        # Only jellyfish should hold the read end, so samtools gets
        # SIGPIPE instead of blocking if jellyfish exits early.
        p_samtools.stdout.close()

        while True:
            try:
                p_jellyfish.wait(timeout=poll_interval)
                break
            except subprocess.TimeoutExpired:
                if on_poll is not None:
                    on_poll(time.monotonic() - start, p_samtools, p_jellyfish)
        p_samtools.wait()

        if p_jellyfish.returncode != 0:
            raise RuntimeError(
                f"jellyfish count ({label}) failed with exit code "
                f"{p_jellyfish.returncode}: {_stderr_tail(jellyfish_err)}"
            )
        if p_samtools.returncode != 0:
            raise RuntimeError(
                f"samtools fasta ({label}) failed with exit code "
                f"{p_samtools.returncode}, so jellyfish counted only part "
                f"of the input: {_stderr_tail(samtools_err)}"
            )
        samtools_err.seek(0)
        if b"EOF marker is absent" in samtools_err.read():
            raise RuntimeError(
                f"samtools fasta ({label}): the input appears truncated "
                f"(EOF marker is absent), so jellyfish may have counted "
                f"only part of it: {_stderr_tail(samtools_err)}"
            )


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

    def _log_progress(elapsed, p_samtools, p_jellyfish):
        # Report jellyfish output file size as a proxy for progress
        jf_files_progress = _find_jf_files(jf_output)
        if jf_files_progress:
            total_size = sum(
                os.path.getsize(f) for f in jf_files_progress
                if os.path.exists(f)
            )
            if total_size >= 1024**3:
                jf_size = f"{total_size / (1024**3):.1f} GB"
            elif len(jf_files_progress) == 1:
                jf_size = _format_file_size(jf_files_progress[0])
            else:
                jf_size = f"{total_size / (1024**2):.1f} MB"
        else:
            jf_size = "pending"
        logger.info(
            "  … still scanning (%s elapsed, jf index: %s)",
            _format_elapsed(elapsed), jf_size,
        )
        _log_memory("parent scanning")
        _log_subprocess_memory(p_jellyfish, "jellyfish-count")
        _log_subprocess_memory(p_samtools, "samtools-fasta")

    scan_start = time.monotonic()
    _run_samtools_jellyfish(
        samtools_cmd, jellyfish_cmd, f"parent: {parent_bam}",
        on_poll=_log_progress,
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


def _check_ref_jf(ref_jf, kmer_size):
    """Check that a reference index suits this run; reads only its header.

    Reference subtraction queries canonical child k-mers of length
    *kmer_size*.  ``jellyfish query`` does not fail on a mismatched index:
    with a different k it re-splits every query into k-mers of the
    index's length, and without canonical mode it misses reference
    k-mers that occur only on the reverse strand.

    Raises:
        RuntimeError: If the header can't be read, or the index holds
            k-mers of another length or was built without ``-C``.
    """
    result = subprocess.run(
        ["jellyfish", "info", "-j", ref_jf], capture_output=True, text=True,
    )
    try:
        header = json.loads(result.stdout) if result.returncode == 0 else None
    except ValueError:
        header = None
    # jellyfish info exits 0 on a non-Jellyfish file, printing a stub
    # header without key_len.
    if not isinstance(header, dict) or not header.get("key_len"):
        raise RuntimeError(
            f"Cannot read Jellyfish index {ref_jf}: "
            f"{result.stderr.strip() or 'no Jellyfish header found'}"
        )
    index_k = header["key_len"] // 2
    if index_k != kmer_size:
        raise RuntimeError(
            f"Reference Jellyfish index {ref_jf} holds {index_k}-mers but "
            f"--kmer-size is {kmer_size}; pass an index built with "
            f"-m {kmer_size}, or delete this one so it is rebuilt"
        )
    if not header.get("canonical"):
        raise RuntimeError(
            f"Reference Jellyfish index {ref_jf} was built without -C, so "
            f"reference k-mers seen only on the reverse strand would not be "
            f"subtracted; rebuild it with jellyfish count -C"
        )


def _ensure_ref_jf(ref_fasta, kmer_size, threads, ref_jf=None):
    """Ensure a Jellyfish reference index exists, building it if necessary.

    An existing index is checked with :func:`_check_ref_jf`.  A new one
    is written under a temporary name and renamed when complete, so an
    interrupted build never leaves a partial index for later runs to use.

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
        _check_ref_jf(ref_jf, kmer_size)
        logger.info("Reference Jellyfish index found: %s", ref_jf)
        return ref_jf

    logger.info(
        "Building reference Jellyfish index: %s (k=%d, threads=%d)",
        ref_jf, kmer_size, threads,
    )
    ref_hash_size = _estimate_jf_hash_size(ref_fasta, kmer_size, default="3G")
    logger.info("  Reference JF hash size: %s", ref_hash_size)
    build_start = time.monotonic()
    tmp_jf = f"{ref_jf}.tmp{os.getpid()}"
    cmd = [
        "jellyfish", "count",
        "-m", str(kmer_size),
        "-s", ref_hash_size,
        "-t", str(threads),
        "-C",
        ref_fasta,
        "-o", tmp_jf,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        for f in _find_jf_files(tmp_jf):
            os.remove(f)
        raise RuntimeError(
            f"jellyfish count (reference) failed: {result.stderr}"
        )
    os.replace(tmp_jf, ref_jf)
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

    Creates a ``.jf`` hash file of the proband-unique k-mers: the index
    discovery mode keeps beside its output (``[out-prefix].proband_unique.jf``,
    tens of MB for a genome) for cohort-wide use, where the k-mers
    recurring in unrelated probands are errors.  Module 3 no longer reads
    it: the scan loads the FASTA into a sorted array itself (see
    :mod:`kmer_denovo_filter.core.kmer_set`).

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
