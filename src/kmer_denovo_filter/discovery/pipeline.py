"""Discovery-mode pipeline for VCF-free de novo k-mer analysis."""

import bisect
import collections
import concurrent.futures
import contextlib
import json
import logging
import math
import multiprocessing
import os
import random
import select
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

import pysam

try:
    import fcntl
except ImportError:  # not a Unix: no pipe sizing
    fcntl = None

from kmer_denovo_filter.discovery import nhf
from kmer_denovo_filter.core.bam_scanner import (
    _clip_side,
    _collect_read_alignment_metadata,
    _init_scan_worker,
    _parse_cigar,
    _query_span,
    _reset_scan_worker,
    _scan_contig_for_hits,
)
from kmer_denovo_filter.core.jellyfish_wrappers import (
    _build_proband_jf_index,
    _ensure_ref_jf,
    _merge_jf_files,
    _run_samtools_jellyfish,
    _scan_parent_jellyfish,
)
from kmer_denovo_filter.core.kmer_set import KmerSet, MAX_CODED_K
from kmer_denovo_filter.core.memory_utils import (
    _get_available_memory_gb,
    _log_children_memory,
    _log_dir_size,
    _log_disk_usage,
    _log_memory,
    _log_subprocess_memory,
    _log_subprocesses_memory,
)
from kmer_denovo_filter.kmer_utils import (
    canonicalize,
    read_supports_alt,
)
from kmer_denovo_filter.utils import (
    _check_tool,
    _collect_kmer_ref_positions,
    _estimate_fasta_sequence_count,
    _estimate_jf_hash_size,
    size_child_count,
    _find_jf_files,
    _format_elapsed,
    _format_file_size,
    _infer_sv_type,
    _is_tmpfs,
    _resolve_tmp_dir,
    _validate_inputs,
    _write_kmer_fasta,
)

logger = logging.getLogger(__name__)

#: How many times one scan task (a contig, or a window of one) is attempted
#: in Module 3 before the run fails (see :func:`_anchor_and_cluster`).
_SCAN_TASK_ATTEMPTS = 3

#: Module 3 scans a contig in windows of about this many bases, one task
#: each, so that the last tasks of a genome are short and the workers stay
#: busy to the end (chr1 alone took the last 45 minutes of a 16-worker run
#: when a contig was one task).  A read belongs to the window its alignment
#: starts in.
_SCAN_WINDOW_BP = 20_000_000

#: What a scan worker costs beyond the k-mer set: the interpreter with
#: pysam and numpy and the reads of one batch, and the transient arrays of
#: the batch (a few tens of bytes per base of it; see
#: ``core.bam_scanner._SCAN_BATCH_BASES``).  Measured in the image on a
#: 3 M-read genome with a 6 M-k-mer set: 0.39 GB per worker at most.
_SCAN_WORKER_BASE_GB = 0.25
_SCAN_BYTES_PER_BATCH_BASE = 56

#: The reference subtraction's splitter moves the dump in blocks of this
#: size, each cut at a record boundary, and asks for this much pipe buffer
#: per query process (Linux; the default pipe size elsewhere).
_SUBTRACT_BLOCK_BYTES = 1 << 20
_SUBTRACT_PIPE_BYTES = 1 << 20
_F_SETPIPE_SZ = 1031  # fcntl(2) on Linux

#: The filter behind each ``jellyfish query`` of the reference subtraction:
#: keeps the k-mers whose reference count is 0 as FASTA records (headers
#: unique across the parts), and leaves the number of k-mers seen and kept
#: in a counts file.
_NON_REF_AWK = (
    '$2 == 0 { n++; print ">" p "." n; print $1 } '
    'END { print NR, n + 0 > counts }'
)


# ── Discovery-mode helpers ────────────────────────────────────────


def _child_count_command(child_jf, kmer_size, threads, hash_entries, bloom_elements, lower_count):
    """``jellyfish count`` for the child, reading the FASTA stream on stdin: canonical k-mers, the hash size given,
    ``--disk`` (a full table is dumped and merged, never doubled in memory), the Bloom filter when *bloom_elements*
    is positive, and ``-L`` so that only k-mers at or above *lower_count* are written to the index."""
    cmd = [
        "jellyfish", "count",
        "-m", str(kmer_size),
        "-s", str(hash_entries),
        "-t", str(threads),
        "-C",
        "--disk",
    ]
    if bloom_elements:
        cmd.extend(["--bf-size", str(bloom_elements)])
    cmd.extend(["-L", str(lower_count), "-o", child_jf, "/dev/fd/0"])
    return cmd


def _extract_child_kmers_discovery(child_bam, ref_fasta, kmer_size,
                                   min_child_count, threads, tmpdir,
                                   jf_hash_size=None, memory_gb=None,
                                   bloom_filter=True):
    """Module 1a: count the child's k-mers.

    Child k-mers are counted with ``jellyfish count -C`` (canonical mode)
    to match the reference index built by :func:`_ensure_ref_jf`, so that
    the reference subtraction (:func:`_subtract_reference_kmers`) compares
    k-mers in the same canonical orientation.

    The count is bounded in memory (``jellyfish_wrappers.size_child_count``):
    ``--disk`` stops jellyfish from doubling a full table in RAM (the old
    behaviour: a 4G table for a 70x child doubled to 8G, 16G and died
    asking for 32G, 264 GB with the copy in flight) - a full table is dumped
    to ``child.jf_N`` files that jellyfish merges into ``child.jf`` at the end;
    the hash size comes from *memory_gb* (the ``--memory`` budget) and from
    what the input needs, never from the file size alone; and, with
    *bloom_filter* (the default when *min_child_count* is 2 or more), a Bloom
    filter in front of the table keeps the k-mers seen once - most of a
    WGS child's distinct k-mers are sequencing errors seen once - out of
    the hash.  The filter admits a k-mer from its second occurrence, so every
    count in the hash is one less than the truth: the dump's lower bound is
    ``min_child_count - 1``.  *jf_hash_size* (a string such as '4G') overrides
    the hash size.

    When Jellyfish leaves multiple chunk files (``--no-merge``, an older
    jellyfish), they are merged in a streaming fashion.

    Returns:
        child_jf: path of the child's index (``None`` when jellyfish wrote
            nothing), to be dumped and subtracted by
            :func:`_subtract_reference_kmers`.
        dump_lower: the ``-L`` bound that selects the k-mers with count
            >= *min_child_count* from it.
    """
    child_jf = os.path.join(tmpdir, "child.jf")

    # Step 1: Count all child k-mers, within the memory budget
    use_bf = bool(bloom_filter) and min_child_count >= 2
    if bloom_filter and not use_bf:
        logger.info(
            "Bloom filter off: --min-child-count %d keeps k-mers seen once, which the filter would drop",
            min_child_count,
        )
    sizing = size_child_count(child_bam, kmer_size, memory_gb=memory_gb, bloom_filter=use_bf, jf_hash_size=jf_hash_size)
    jf_hash_size = str(sizing["hash_entries"])
    dump_lower = min_child_count - 1 if use_bf else min_child_count
    logger.info(
        "Extracting child k-mers from BAM (k=%d; %s; --disk: a full table is dumped and merged, never doubled in memory)…",
        kmer_size, sizing["note"],
    )
    if use_bf:
        logger.info(
            "  Bloom filter: a k-mer enters the hash from its second occurrence, so the hash counts are one less than "
            "the truth; k-mers with count >= %d are dumped with -L %d", min_child_count, dump_lower,
        )
    samtools_threads = max(1, threads // 4)
    samtools_cmd = [
        "samtools", "fasta", "-F", "0xD00",
        "-@", str(samtools_threads),
        child_bam,
    ]
    if ref_fasta:
        samtools_cmd.extend(["--reference", ref_fasta])

    jellyfish_cmd = _child_count_command(child_jf, kmer_size, threads, jf_hash_size, sizing["bloom_elements"] if use_bf else 0, dump_lower)

    def _log_progress(elapsed, p_samtools, p_jellyfish):
        jf_files = _find_jf_files(child_jf)
        if jf_files:
            total_size = sum(
                os.path.getsize(f) for f in jf_files
                if os.path.exists(f)
            )
            jf_size = f"{total_size / (1024**3):.1f} GB"
            n_chunks = len(jf_files)
        else:
            jf_size = "pending"
            n_chunks = 0
        logger.info(
            "  … child k-mer counting (%s elapsed, jf index: %s, "
            "files: %d - more than one: the table filled and was dumped; jellyfish merges them at the end)",
            _format_elapsed(elapsed), jf_size, n_chunks,
        )
        _log_memory("child k-mer counting")
        _log_subprocess_memory(p_jellyfish, "jellyfish-count")
        _log_subprocess_memory(p_samtools, "samtools-fasta")
        _log_disk_usage(tmpdir, "tmpdir during counting")

    extract_start = time.monotonic()
    _run_samtools_jellyfish(
        samtools_cmd, jellyfish_cmd, f"child: {child_bam}",
        poll_interval=60, on_poll=_log_progress,
    )

    # Check for multi-file output (hash overflow)
    jf_files = _find_jf_files(child_jf)
    total_jf_size = sum(
        os.path.getsize(f) for f in jf_files if os.path.exists(f)
    )
    logger.info(
        "Child k-mer counting complete (%s, index: %.1f GB, files: %d)",
        _format_elapsed(time.monotonic() - extract_start),
        total_jf_size / (1024**3), len(jf_files),
    )
    _log_memory("after child k-mer counting")

    # Merge chunk files if Jellyfish produced multiple outputs
    if len(jf_files) > 1:
        merged_jf = os.path.join(tmpdir, "child_merged.jf")
        return _merge_jf_files(jf_files, merged_jf, threads), dump_lower
    if jf_files:
        return jf_files[0], dump_lower
    logger.warning("No Jellyfish output files found")
    return None, dump_lower


def _stream_records_to_pipes(source_fd, sink_fds, block_bytes=_SUBTRACT_BLOCK_BYTES,
                             on_block=None):
    """Move a FASTA stream from *source_fd* into the pipes *sink_fds*,
    whole records at a time.

    Blocks of about *block_bytes* are read from the source and cut at the
    last record start (``\n>``), so every sink receives complete records
    (a k-mer's header and sequence line together).  Each block goes to an
    idle sink, chosen round-robin among those that have nothing pending,
    and the pending block of a sink is written as the pipe drains
    (``select``): a slow consumer simply gets fewer blocks, and no block
    waits in Python for a full pipe while another pipe is empty.  Nothing
    is buffered beyond one block per sink.

    *on_block* is called after every block with ``(bytes_moved,
    records_moved)`` so far.  Returns ``(bytes_moved, records_moved)``;
    a record is a line starting with ``>``.

    Raises:
        BrokenPipeError: when a sink's process is gone.
    """
    for fd in sink_fds:
        os.set_blocking(fd, False)
        if fcntl is not None:
            try:
                fcntl.fcntl(fd, _F_SETPIPE_SZ, _SUBTRACT_PIPE_BYTES)
            except OSError:
                pass  # not Linux, or the size is not allowed: the default holds
    pending = {fd: memoryview(b"") for fd in sink_fds}
    carry = b""
    n_bytes = 0
    n_records = 0
    eof = False
    next_sink = 0

    def _assign(chunk):
        nonlocal next_sink, n_bytes, n_records
        for _ in range(len(sink_fds)):
            fd = sink_fds[next_sink]
            next_sink = (next_sink + 1) % len(sink_fds)
            if not pending[fd]:
                pending[fd] = memoryview(chunk)
                break
        else:
            raise AssertionError("no idle sink")
        n_bytes += len(chunk)
        n_records += chunk.count(b">")
        if on_block is not None:
            on_block(n_bytes, n_records)

    while True:
        busy = [fd for fd in sink_fds if pending[fd]]
        idle = len(busy) < len(sink_fds)
        if eof and not busy:
            break
        readable = [source_fd] if (idle and not eof) else []
        ready_r, ready_w, _ = select.select(readable, busy, [])
        for fd in ready_w:
            try:
                written = os.write(fd, pending[fd])
            except BlockingIOError:
                continue
            pending[fd] = pending[fd][written:]
        if ready_r:
            block = os.read(source_fd, block_bytes)
            if not block:
                eof = True
                if carry:
                    _assign(carry)
                    carry = b""
                continue
            block = carry + block if carry else block
            cut = block.rfind(b"\n>")
            if cut < 0:
                carry = block
                continue
            carry = block[cut + 1:]
            _assign(block[:cut + 1])
    return n_bytes, n_records


def _stderr_tail(fh, max_lines=20):
    fh.seek(0)
    lines = fh.read().decode(errors="replace").strip().splitlines()
    return "\n".join(lines[-max_lines:])


def _subtract_reference_kmers(ref_jf, child_jf, dump_lower, tmpdir, threads,
                              min_child_count=None, poll_interval=60):
    """Module 1b: the child's candidate k-mers that are absent from the
    reference, streamed from the child's index through parallel queries.

    ``jellyfish dump -L dump_lower child_jf`` writes the candidate k-mers
    (count >= *min_child_count*; *dump_lower* carries the Bloom filter's
    offset) as FASTA, one record each.  That stream is cut at record
    boundaries and fed, block by block, into the stdin pipes of *threads*
    ``jellyfish query ref_jf -s /dev/stdin`` processes
    (:func:`_stream_records_to_pipes`).  Every query process memory-maps
    the one reference index, which the page cache holds once, so this
    costs no more memory than the single query it replaces; each one's
    output (``k-mer count`` lines) goes through ``awk``, which keeps the
    k-mers whose reference count is 0 as FASTA records in a part file.
    The parts are concatenated into ``child_non_ref_kmers.fa``.  Nothing
    of the dump is written to disk, and nothing of it is held in Python
    beyond one block per pipe (the old way wrote the whole dump, 100 GB
    for a 30x genome, to a FASTA and ran one query over it).

    Both the reference index and the child's are built with ``-C``, so
    the k-mers compare in the same canonical orientation.  The child's
    index is removed once the dump has been read.

    Returns:
        child_non_ref_fa: path to the FASTA of non-reference child k-mers.
        n_candidates: number of candidate k-mers dumped.
        n_non_ref: number of them absent from the reference.
    """
    n_query = max(1, int(threads))
    if min_child_count is None:
        min_child_count = dump_lower
    parts_dir = os.path.join(tmpdir, "non_ref_parts")
    os.makedirs(parts_dir, exist_ok=True)
    child_non_ref_fa = os.path.join(tmpdir, "child_non_ref_kmers.fa")
    logger.info(
        "Dumping child k-mers with count >= %d (-L %d) from %s (%s) into %d "
        "parallel jellyfish query processes against %s (memory-mapped, held "
        "once by the page cache); nothing of the dump touches disk…",
        min_child_count, dump_lower, child_jf, _format_file_size(child_jf),
        n_query, ref_jf,
    )
    start = time.monotonic()
    last_log = start
    dump_cmd = ["jellyfish", "dump", "-L", str(dump_lower), child_jf]
    query_cmd = ["jellyfish", "query", ref_jf, "-s", "/dev/stdin"]
    procs = []        # (label, Popen, stderr file)
    queries = []
    part_paths = []
    count_paths = []
    n_bytes = n_records = 0
    with contextlib.ExitStack() as stack:
        dump_err = stack.enter_context(tempfile.TemporaryFile())
        p_dump = subprocess.Popen(
            dump_cmd, stdout=subprocess.PIPE, stderr=dump_err, bufsize=0,
        )
        procs.append(("jellyfish dump (child)", p_dump, dump_err))
        for j in range(n_query):
            part = os.path.join(parts_dir, f"part_{j:03d}.fa")
            counts = os.path.join(parts_dir, f"part_{j:03d}.counts")
            part_paths.append(part)
            count_paths.append(counts)
            q_err = stack.enter_context(tempfile.TemporaryFile())
            p_query = subprocess.Popen(
                query_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=q_err, bufsize=0,
            )
            procs.append((f"jellyfish query (ref subtraction, pipe {j})",
                          p_query, q_err))
            queries.append(p_query)
            f_err = stack.enter_context(tempfile.TemporaryFile())
            part_fh = stack.enter_context(open(part, "wb"))
            p_filter = subprocess.Popen(
                ["awk", "-v", f"p={j}", "-v", f"counts={counts}",
                 _NON_REF_AWK],
                stdin=p_query.stdout, stdout=part_fh, stderr=f_err,
            )
            p_query.stdout.close()  # only awk holds the read end
            procs.append((f"awk (ref subtraction, pipe {j})", p_filter,
                          f_err))

        def _on_block(moved_bytes, moved_records):
            nonlocal last_log
            now = time.monotonic()
            if now - last_log >= poll_interval:
                elapsed = now - start
                logger.info(
                    "  … %d k-mers dumped and queried so far (%s, %.1f GB "
                    "streamed, %.2g k-mers/s)",
                    moved_records, _format_elapsed(elapsed),
                    moved_bytes / 1024**3,
                    moved_records / elapsed if elapsed else 0.0,
                )
                _log_memory("during reference subtraction")
                _log_subprocess_memory(p_dump, "jellyfish-dump")
                _log_subprocesses_memory(queries, "jellyfish-query")
                last_log = now

        failure = None
        try:
            n_bytes, n_records = _stream_records_to_pipes(
                p_dump.stdout.fileno(), [q.stdin.fileno() for q in queries],
                on_block=_on_block,
            )
        except BrokenPipeError as exc:
            failure = exc
        finally:
            for q in queries:
                try:
                    q.stdin.close()
                except OSError:
                    pass
            p_dump.stdout.close()
        if failure is not None:
            for _label, proc, _err in procs:
                if proc.poll() is None:
                    proc.kill()
        for _label, proc, _err in procs:
            proc.wait()
        failed = [(label, proc, err) for label, proc, err in procs
                  if proc.returncode != 0]
        if failed:
            label, proc, err = failed[0]
            raise RuntimeError(
                f"{label} failed with exit code {proc.returncode}: "
                f"{_stderr_tail(err)}"
            )
        if failure is not None:
            raise RuntimeError(
                "reference subtraction: a query pipe closed early although "
                "every process exited 0"
            ) from failure

    dump_elapsed = time.monotonic() - start
    logger.info(
        "Child k-mer dump complete (%s, %d candidates, %.1f GB streamed at "
        "%.2g k-mers/s)",
        _format_elapsed(dump_elapsed), n_records, n_bytes / 1024**3,
        n_records / dump_elapsed if dump_elapsed else 0.0,
    )

    # Remove the child jellyfish index immediately – it is no longer
    # needed and can be very large (30 GB for a 30x genome).
    for f in set(_find_jf_files(os.path.join(tmpdir, "child.jf")) + [child_jf]):
        if os.path.exists(f):
            os.remove(f)
    logger.info("Removed child Jellyfish index to free disk/cache")
    _log_memory("after child index removal")
    logger.info(
        "Child candidate k-mers (count >= %d): %d",
        min_child_count, n_records,
    )

    # Every k-mer the splitter sent must have been queried and counted.
    n_queried = 0
    n_non_ref = 0
    for counts in count_paths:
        with open(counts) as fh:
            seen, kept = fh.read().split()
        n_queried += int(seen)
        n_non_ref += int(kept)
    if n_queried != n_records:
        raise RuntimeError(
            f"reference subtraction: {n_records} k-mers were dumped but "
            f"the query processes reported {n_queried}"
        )
    concat_start = time.monotonic()
    with open(child_non_ref_fa, "wb") as out:
        for part in part_paths:
            with open(part, "rb") as fh:
                shutil.copyfileobj(fh, out, 16 << 20)
            os.remove(part)
    shutil.rmtree(parts_dir, ignore_errors=True)

    logger.info(
        "Non-reference child k-mers after subtraction: %d (%s; %s to "
        "collect the %d parts into %s)",
        n_non_ref, _format_elapsed(time.monotonic() - start),
        _format_elapsed(time.monotonic() - concat_start), n_query,
        _format_file_size(child_non_ref_fa),
    )
    return child_non_ref_fa, n_records, n_non_ref


def _count_parent_jellyfish(parent_bam, ref_fasta, kmer_fasta, kmer_size,
                            parent_dir, threads, label="Parent",
                            n_filter_kmers=None):
    """Count child k-mers in a parent BAM using jellyfish.

    Runs ``samtools fasta | jellyfish count --if`` to produce a Jellyfish
    index that tracks only the child candidate k-mers.  Returns the path
    to the ``.jf`` file — the caller is responsible for querying and
    deleting it.

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
        label: Human-readable label for log messages.
        n_filter_kmers: Number of k-mers in *kmer_fasta*.  When ``None``
            this is estimated from a sampled FASTA prefix and file size.

    Returns:
        Path to the Jellyfish index file.
    """
    os.makedirs(parent_dir, exist_ok=True)
    jf_output = os.path.join(parent_dir, "parent.jf")

    # Size hash to fit the filter k-mers without overflow.
    if n_filter_kmers is None:
        n_filter_kmers, n_filter_kmers_is_extrapolated = _estimate_fasta_sequence_count(
            kmer_fasta
        )
        if n_filter_kmers_is_extrapolated:
            logger.info(
                "  estimated filter_kmers from FASTA size/sample: ~%d",
                n_filter_kmers,
            )
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
        "%s: scanning BAM (%s): %s", label, bam_size, parent_bam,
    )
    logger.info(
        "  samtools fasta → jellyfish count (k=%d, threads=%d, "
        "hash=%s, filter_kmers=%d)",
        kmer_size, threads, hash_size_str, n_filter_kmers,
    )

    def _log_progress(elapsed, p_samtools, p_jellyfish):
        jf_files = _find_jf_files(jf_output)
        if jf_files:
            total_size = 0
            for f in jf_files:
                try:
                    total_size += os.path.getsize(f)
                except FileNotFoundError:
                    pass
            jf_size = f"{total_size / (1024**3):.1f} GB"
        elif os.path.exists(jf_output):
            jf_size = _format_file_size(jf_output)
        else:
            jf_size = "pending"
        logger.info(
            "  … %s still scanning (%s elapsed, jf index: %s)",
            label, _format_elapsed(elapsed), jf_size,
        )
        _log_memory(f"{label} counting")
        _log_subprocess_memory(p_jellyfish, f"jellyfish-count ({label})")
        _log_subprocess_memory(p_samtools, f"samtools-fasta ({label})")

    scan_start = time.monotonic()
    _run_samtools_jellyfish(
        samtools_cmd, jellyfish_cmd, f"{label}: {parent_bam}",
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

    logger.info(
        "  %s jellyfish counting complete (%s, index: %s)",
        label, _format_elapsed(time.monotonic() - scan_start),
        _format_file_size(jf_output),
    )
    return jf_output


def _filter_parents_discovery(mother_bam, father_bam, ref_fasta,
                              child_non_ref_fa, kmer_size, threads, tmpdir,
                              parent_max_count=0):
    """Module 2: Filter non-reference child k-mers against both parents.

    Uses streaming jellyfish queries to avoid loading all non-reference
    k-mers into Python memory simultaneously.  Each parent's Jellyfish
    index is built, queried, and deleted before proceeding to the next,
    keeping peak memory proportional to disk I/O rather than in-memory
    data structures.

    After the mother scan, only the surviving k-mers are written to a
    reduced FASTA used as the ``--if`` filter for the father scan, which
    further reduces jellyfish memory.

    Returns a tuple of:
        n_proband_unique: Number of proband-unique k-mers (int).
        proband_unique_fa: Path to FASTA file of proband-unique k-mers,
            or None when no k-mers survive filtering.

    Args:
        parent_max_count: Maximum k-mer count allowed in a parent before
            the k-mer is considered parental.  K-mers with count >
            parent_max_count in either parent are removed.
    """
    n_input, n_input_is_extrapolated = _estimate_fasta_sequence_count(
        child_non_ref_fa
    )

    if n_input == 0:
        return 0, None

    if n_input_is_extrapolated:
        logger.info(
            "Filtering ~%d non-reference k-mers against parents…", n_input,
        )
    else:
        logger.info(
            "Filtering %d non-reference k-mers against parents…", n_input,
        )
    _log_memory("before parent filtering")

    # ── Mother scan ────────────────────────────────────────────────
    mother_jf = _count_parent_jellyfish(
        mother_bam, ref_fasta, child_non_ref_fa, kmer_size,
        os.path.join(tmpdir, "mother"), threads, label="Mother",
        n_filter_kmers=n_input,
    )

    # Stream-filter: keep k-mers with count <= parent_max_count
    after_mother_fa = os.path.join(tmpdir, "after_mother.fa")
    n_surviving = 0
    n_removed_mother = 0
    query_cmd = [
        "jellyfish", "query", mother_jf, "-s", child_non_ref_fa,
    ]
    with tempfile.TemporaryFile(mode="w+") as stderr_f:
        p_query = subprocess.Popen(
            query_cmd, stdout=subprocess.PIPE, stderr=stderr_f, text=True,
        )
        with open(after_mother_fa, "w") as fh:
            for line in p_query.stdout:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 2 and int(parts[1]) <= parent_max_count:
                    fh.write(f">{n_surviving}\n{parts[0]}\n")
                    n_surviving += 1
                else:
                    n_removed_mother += 1
        p_query.wait()
        if p_query.returncode != 0:
            stderr_f.seek(0)
            raise RuntimeError(
                f"jellyfish query (mother filter) failed: {stderr_f.read()}"
            )

    # Remove mother index immediately
    if os.path.exists(mother_jf):
        os.remove(mother_jf)

    logger.info(
        "Mother: %d / %d non-ref k-mers found (count > %d), %d surviving",
        n_removed_mother, n_input, parent_max_count, n_surviving,
    )
    _log_memory("after mother filtering")

    if n_surviving == 0:
        return 0, None

    # ── Father scan (uses reduced k-mer set from mother filter) ────
    father_jf = _count_parent_jellyfish(
        father_bam, ref_fasta, after_mother_fa, kmer_size,
        os.path.join(tmpdir, "father"), threads, label="Father",
        n_filter_kmers=n_surviving,
    )

    # Stream-filter: write surviving k-mers to FASTA (no in-memory set
    # to avoid holding billions of k-mer strings in Python memory).
    proband_unique_fa = os.path.join(tmpdir, "proband_unique.fa")
    n_removed_father = 0
    n_proband = 0
    query_cmd = [
        "jellyfish", "query", father_jf, "-s", after_mother_fa,
    ]
    with tempfile.TemporaryFile(mode="w+") as stderr_f:
        p_query = subprocess.Popen(
            query_cmd, stdout=subprocess.PIPE, stderr=stderr_f, text=True,
        )
        with open(proband_unique_fa, "w") as fh:
            for line in p_query.stdout:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 2 and int(parts[1]) <= parent_max_count:
                    kmer = parts[0]
                    fh.write(f">{n_proband}\n{kmer}\n")
                    n_proband += 1
                else:
                    n_removed_father += 1
        p_query.wait()
        if p_query.returncode != 0:
            stderr_f.seek(0)
            raise RuntimeError(
                f"jellyfish query (father filter) failed: {stderr_f.read()}"
            )

    # Remove father index and intermediate file
    if os.path.exists(father_jf):
        os.remove(father_jf)
    if os.path.exists(after_mother_fa):
        os.remove(after_mother_fa)

    logger.info(
        "Father: %d / %d surviving k-mers found (count > %d), "
        "%d proband-unique",
        n_removed_father, n_surviving, parent_max_count,
        n_proband,
    )
    logger.info(
        "Proband-unique k-mers (absent from both parents): %d / %d",
        n_proband, n_input,
    )
    logger.info(
        "Proband-unique FASTA: %s (%s)",
        proband_unique_fa, _format_file_size(proband_unique_fa),
    )
    _log_memory("after parent filtering")
    return n_proband, proband_unique_fa


def _scan_windows(contigs, lengths, window_size):
    """The scan tasks of a genome: ``(contig, start, end)`` per window of
    about *window_size* bases (a contig shorter than that is one window,
    scanned whole: ``(contig, None, None)``), in header order, then
    ``(None, None, None)`` for the unplaced unmapped reads."""
    windows = []
    for contig, length in zip(contigs, lengths):
        n_windows = (max(1, math.ceil(length / window_size))
                     if window_size and length else 1)
        if n_windows == 1:
            windows.append((contig, None, None))
            continue
        size = math.ceil(length / n_windows)
        for w in range(n_windows):
            windows.append((contig, w * size, min(length, (w + 1) * size)))
    windows.append((None, None, None))
    return windows


def _window_label(contig, start, end):
    if contig is None:
        return "(unmapped)"
    if start is None:
        return contig
    return f"{contig}:{start + 1:,}-{end:,}"


def _scan_pool_context():
    """The start method of the scan workers: ``fork`` where it exists, so
    the workers share the parent's k-mer set copy-on-write (and are its
    direct children, which the memory log reads); the platform's default
    elsewhere (spawn on macOS, where each worker gets a copy)."""
    if (sys.platform.startswith("linux")
            and "fork" in multiprocessing.get_all_start_methods()):
        return multiprocessing.get_context("fork")
    return None


def _anchor_and_cluster(child_bam, ref_fasta, proband_unique_kmers,
                        kmer_size, merge_distance=500, threads=1,
                        min_distinct_kmers_per_read=1,
                        proband_unique_fa=None,
                        n_proband_unique=None,
                        tmpdir=None,
                        memory_limit_gb=None,
                        informative_bam=None,
                        window_size=_SCAN_WINDOW_BP):
    """Module 3: Find reads containing proband-unique k-mers and cluster regions.

    Scans **all** primary, non-duplicate child reads — including unmapped
    and low-MAPQ reads — so that no proband-unique k-mers are missed at
    this stage.  The proband-unique k-mers are loaded once, here, as a
    :class:`KmerSet` (a sorted array of 2-bit codes, about 8 bytes per
    k-mer) that the workers share copy-on-write; each worker tests its
    reads against it in vectorized batches (see
    :mod:`kmer_denovo_filter.core.kmer_set`).  The scan runs one task per
    window of about *window_size* bases of each contig (a shorter contig
    is one task) plus one for unplaced unmapped reads, in a process pool
    when more than one worker is available and in-process otherwise.
    Results are combined in task order, so the output does not depend on
    which worker finishes first, nor on the window size.

    Args:
        child_bam: Path to child BAM file.
        ref_fasta: Path to reference FASTA (or None).
        proband_unique_kmers: Iterable of k-mer strings, or None when
            *proband_unique_fa* is provided.
        kmer_size: K-mer size.
        merge_distance: Maximum gap for merging adjacent regions.
        threads: Number of parallel workers (default 1).
        min_distinct_kmers_per_read: Minimum distinct proband-unique
            k-mers a read must carry to be retained (default 1).
        proband_unique_fa: Optional FASTA path of the proband-unique
            k-mers (one per record), loaded into the set.
        n_proband_unique: Number of proband-unique k-mers (for logging).
        tmpdir: Writable directory for temporary files (the per-task
            informative-read shards).
        memory_limit_gb: Explicit memory limit in GB.  When provided,
            overrides auto-detected system memory for worker-count
            planning.  Useful on HPC where SLURM allocations differ
            from total node memory.
        informative_bam: Optional output path.  When given, every read
            retained here (i.e. carrying at least
            *min_distinct_kmers_per_read* distinct proband-unique k-mers)
            is written to this BAM with a ``dk:i:1`` tag during the same
            scan; the BAM is then sorted and indexed.
        window_size: Bases per scan task; ``None`` or 0 scans each contig
            as one task.

    Returns:
        Tuple of (regions, region_reads, total_informative, region_kmers,
        unmapped_informative, read_sv_meta, kmer_coverage, read_coverage).
    """
    anchor_start = time.monotonic()

    # ── The k-mer set, loaded once and shared with the workers ──────
    load_start = time.monotonic()
    if proband_unique_fa:
        kmer_set = KmerSet.from_fasta(proband_unique_fa, kmer_size)
    else:
        kmer_set = KmerSet.from_kmers(proband_unique_kmers or (), kmer_size)
    n_kmer_count = len(kmer_set)
    mp_context = _scan_pool_context()
    if kmer_size > MAX_CODED_K:
        logger.warning(
            "  k=%d does not fit a 64-bit code (k <= %d does): the scan "
            "uses the per-k-mer fallback, which is much slower",
            kmer_size, MAX_CODED_K,
        )
    logger.info(
        "  Proband-unique k-mer set: %d k-mers, %.1f MB as sorted 2-bit "
        "codes (%s; loaded in %s)",
        n_kmer_count, kmer_set.nbytes / 1024**2,
        "shared copy-on-write by the workers" if mp_context is not None
        else "copied to each worker",
        _format_elapsed(time.monotonic() - load_start),
    )
    _log_memory("after loading the k-mer set")

    # ── Memory planning ────────────────────────────────────────────
    total_mem_gb, avail_mem_gb = _get_available_memory_gb()
    # When the caller supplies an explicit memory limit (e.g. from
    # --memory on HPC), treat it as both total and available so that
    # worker-count planning uses the user's allocation rather than
    # system-reported values.
    if memory_limit_gb is not None:
        total_mem_gb = memory_limit_gb
        avail_mem_gb = memory_limit_gb
        logger.info(
            "  [MemoryPlanning] Using explicit memory limit: %.1f GB",
            memory_limit_gb,
        )
    from kmer_denovo_filter.core.bam_scanner import _SCAN_BATCH_BASES
    batch_gb = _SCAN_BATCH_BASES * _SCAN_BYTES_PER_BATCH_BASE / 1024**3
    set_gb = kmer_set.nbytes / 1024**3
    est_per_worker_gb = _SCAN_WORKER_BASE_GB + batch_gb + (
        0.0 if mp_context is not None else set_gb)
    logger.info(
        "  [MemoryPlanning] k-mer set %.3f GB (one copy%s); about %.2f GB "
        "per worker: %.2f GB of transient arrays per batch of %d bases "
        "plus the interpreter",
        set_gb, "" if mp_context is not None else " per worker",
        est_per_worker_gb, batch_gb, _SCAN_BATCH_BASES,
    )
    if total_mem_gb is not None:
        logger.info(
            "  [MemoryPlanning] System memory: %.1f GB total, %s available",
            total_mem_gb,
            f"{avail_mem_gb:.1f} GB" if avail_mem_gb is not None
            else "(unknown)",
        )

    bam = pysam.AlignmentFile(
        child_bam,
        reference_filename=ref_fasta if ref_fasta else None,
    )
    contigs = list(bam.references)
    lengths = list(bam.lengths)
    bam.close()

    # One task per window of a contig + one for unmapped reads.  With an
    # informative BAM requested, each task writes its reads to its own
    # shard, and the shards are combined in task (coordinate) order
    # afterwards.
    windows = _scan_windows(contigs, lengths, window_size)
    shard_dir = None
    shard_paths = [None] * len(windows)
    if informative_bam:
        shard_dir = tempfile.mkdtemp(prefix="informative_reads_", dir=tmpdir)
        shard_paths = [
            os.path.join(shard_dir, f"{i:06d}.bam")
            for i in range(len(windows))
        ]
    tasks = [
        (child_bam, ref_fasta, contig, shard, start, end)
        for (contig, start, end), shard in zip(windows, shard_paths)
    ]
    labels = [_window_label(contig, start, end)
              for contig, start, end in windows]

    # ── Cap workers by memory ───────────────────────────────────────
    max_workers_by_mem = threads
    if avail_mem_gb is not None and est_per_worker_gb > 0:
        usable_gb = avail_mem_gb * 0.8 - set_gb
        max_workers_by_mem = max(1, int(usable_gb / est_per_worker_gb))
    elif total_mem_gb is not None and est_per_worker_gb > 0:
        usable_gb = total_mem_gb * 0.7 - set_gb
        max_workers_by_mem = max(1, int(usable_gb / est_per_worker_gb))
    n_workers = max(1, min(threads, len(tasks), max_workers_by_mem))

    logger.info(
        "  Anchoring: %d contigs in %d scan tasks (windows of %s bases, "
        "plus the unmapped reads), %d workers (requested=%d)",
        len(contigs), len(tasks),
        f"{window_size:,}" if window_size else "the whole contig",
        n_workers, threads,
    )
    init_args = (kmer_set, kmer_size, min_distinct_kmers_per_read)

    read_hits = []
    reads_seen = set()
    read_sv_meta = {}
    kmer_coverage = collections.defaultdict(collections.Counter)
    read_coverage = collections.defaultdict(collections.Counter)
    unmapped_informative = 0
    total_reads_scanned = 0
    completed_tasks = 0
    informative_received = 0  # before cross-task dedup; for progress logs
    # Hits and SV metadata of finished tasks that are waiting for an
    # earlier task to finish, keyed by task index (see _merge_result).
    held = {}
    next_to_fold = 0
    # The keys reported by the folded windows of each contig, for the
    # coverage a window could not settle on its own (see _merge_result).
    contig_seen = collections.defaultdict(set)
    last_progress_log = time.monotonic()
    progress_interval = 300  # seconds (5 minutes)

    def _merge_result(index, result):
        """Fold one task's scan results into the running totals.

        Coverage and counts are added straight away.  Hits and SV metadata
        are de-duplicated across tasks by (read name, is_supplementary),
        where the first task to report a key wins — e.g. when both mates
        of a pair are informative on different chromosomes, or in
        different windows of one.  Folding them strictly in task order
        keeps that choice independent of which worker finishes first, and
        identical to a single-worker run over whole contigs.

        A scan of a whole contig skips a read whose key it has seen (the
        second informative mate of a pair, say) entirely, coverage
        included, while reads on different contigs each count.  A window
        cannot know what the earlier windows of its contig saw, so it
        defers the coverage of the reads whose key they may have reported
        (``deferred``), and it is added here only when none of them did.
        """
        nonlocal unmapped_informative, total_reads_scanned
        nonlocal completed_tasks, informative_received, next_to_fold
        (hits, seen, unmapped, scanned,
         sv_meta, worker_cov, worker_read_cov, deferred) = result
        total_reads_scanned += scanned
        unmapped_informative += unmapped
        completed_tasks += 1
        informative_received += len(hits) + unmapped
        # Merge k-mer coverage (update is faster than += for Counters)
        for chrom, cov in worker_cov.items():
            kmer_coverage[chrom].update(cov)
        # Merge read coverage
        for chrom, cov in worker_read_cov.items():
            read_coverage[chrom].update(cov)

        held[index] = (hits, seen, sv_meta, deferred, windows[index][0])
        while next_to_fold in held:
            (task_hits, task_seen, task_meta, task_deferred,
             task_contig) = held.pop(next_to_fold)
            next_to_fold += 1
            for hit in task_hits:
                # hit: (ref_name, start, end, qname, kmers, is_supp)
                dedup_key = (hit[3], hit[5])  # (qname, is_supplementary)
                if dedup_key not in reads_seen:
                    reads_seen.add(dedup_key)
                    read_hits.append(hit)
            # Merge SV metadata (only for keys not already seen)
            for key, meta in task_meta.items():
                if key not in read_sv_meta:
                    read_sv_meta[key] = meta
            # Track all seen keys for cross-task dedup
            reads_seen.update(task_seen)
            if task_deferred:
                earlier = contig_seen[task_contig]
                for key, cov in task_deferred:
                    if key in earlier:
                        continue
                    kmer_coverage[task_contig].update(cov)
                    for pos in cov:
                        read_coverage[task_contig][pos] += 1
            if task_contig is not None:
                contig_seen[task_contig].update(task_seen)

        # Log individual task completion for large ones
        if scanned >= 1_000_000:
            logger.info(
                "  [Anchoring] %s complete: %d reads "
                "scanned, %d informative hits (%s elapsed)",
                labels[index], scanned, len(hits) + unmapped,
                _format_elapsed(time.monotonic() - anchor_start),
            )

    def _log_progress():
        """Log progress every 5 minutes or every 100 tasks."""
        nonlocal last_progress_log
        now = time.monotonic()
        at_milestone = (
            completed_tasks % 100 == 0
            or completed_tasks == len(tasks)
        )
        if now - last_progress_log >= progress_interval or at_milestone:
            logger.info(
                "  [Anchoring] Progress: %d/%d tasks complete "
                "(%.0f%%), %d reads scanned, %d informative (%s)",
                completed_tasks, len(tasks),
                100 * completed_tasks / len(tasks),
                total_reads_scanned,
                informative_received,
                _format_elapsed(now - anchor_start),
            )
            _log_memory("during anchoring")
            _log_children_memory("during anchoring")
            last_progress_log = now

    # A task that fails is run again, up to _SCAN_TASK_ATTEMPTS times in
    # all, before the run gives up: on a network filesystem a worker can
    # lose the file it reads or a shard it is writing ("Stale file
    # handle") for reasons that are gone a moment later, and the scan of
    # a whole genome should not be lost to that.  A rerun starts the
    # task's shard afresh and reports its hits in task order like any
    # other, so the result is the same as an untroubled run.
    attempts = collections.Counter()

    def _retry_or_raise(index, exc):
        attempts[index] += 1
        if attempts[index] >= _SCAN_TASK_ATTEMPTS:
            logger.error(
                "Worker failed for %s (%d attempts)",
                labels[index], attempts[index],
            )
            raise exc
        logger.warning(
            "Worker failed for %s (attempt %d of %d): %s; running "
            "the task again",
            labels[index], attempts[index], _SCAN_TASK_ATTEMPTS,
            str(exc).splitlines()[0] if str(exc) else type(exc).__name__,
        )

    # Submission order: the unmapped reads first (their number is unknown),
    # then the longest windows, so that the last tasks to finish are short.
    # Results are still folded in task (index) order.
    def _span(index):
        contig, start, end = windows[index]
        if contig is None:
            return float("inf")
        if start is None:
            return lengths[contigs.index(contig)]
        return end - start
    submit_order = sorted(range(len(tasks)), key=_span, reverse=True)

    try:
        if n_workers > 1:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=n_workers,
                initializer=_init_scan_worker,
                initargs=init_args,
                mp_context=mp_context,
            ) as executor:
                futures = {
                    executor.submit(_scan_contig_for_hits, *tasks[i]): i
                    for i in submit_order
                }

                # Use wait() with timeout for time-based progress
                # reporting so users get feedback even when large
                # tasks take long.
                pending = set(futures.keys())
                while pending:
                    done, pending = concurrent.futures.wait(
                        pending, timeout=progress_interval,
                        return_when=concurrent.futures.FIRST_COMPLETED,
                    )

                    if not done:
                        # Timeout — no task completed; log heartbeat
                        logger.info(
                            "  [Anchoring] Heartbeat: %d/%d tasks "
                            "complete, %d reads scanned, %d informative "
                            "(%s elapsed)",
                            completed_tasks, len(tasks),
                            total_reads_scanned,
                            informative_received,
                            _format_elapsed(
                                time.monotonic() - anchor_start,
                            ),
                        )
                        _log_memory("during anchoring")
                        _log_children_memory("during anchoring")
                        last_progress_log = time.monotonic()
                        continue

                    for future in done:
                        index = futures.pop(future)
                        try:
                            result = future.result()
                        except Exception as exc:
                            _retry_or_raise(index, exc)
                            retry = executor.submit(
                                _scan_contig_for_hits, *tasks[index],
                            )
                            futures[retry] = index
                            pending.add(retry)
                            continue
                        _merge_result(index, result)
                    _log_progress()
        else:
            # A single worker: run the same per-task scan in-process.
            _init_scan_worker(*init_args)
            try:
                for index, task in enumerate(tasks):
                    while True:
                        try:
                            result = _scan_contig_for_hits(*task)
                        except Exception as exc:
                            _retry_or_raise(index, exc)
                            continue
                        break
                    _merge_result(index, result)
                    _log_progress()
            finally:
                _reset_scan_worker()

        if informative_bam:
            n_written = _merge_informative_shards(
                shard_paths, child_bam, ref_fasta, informative_bam,
                os.path.join(shard_dir, "merged.unsorted.bam"),
            )
            logger.info(
                "Informative reads BAM written: %s (%d reads)",
                informative_bam, n_written,
            )
    finally:
        if shard_dir is not None:
            shutil.rmtree(shard_dir, ignore_errors=True)

    logger.info(
        "  Anchoring: %d reads scanned, %d informative (%s) [%d workers]",
        total_reads_scanned,
        len(read_hits) + unmapped_informative,
        _format_elapsed(time.monotonic() - anchor_start),
        n_workers,
    )

    _log_memory("after anchoring complete")

    total_informative = len(read_hits) + unmapped_informative
    logger.info(
        "Anchoring complete: %d informative reads (%d mapped, %d unmapped) "
        "from %d scanned (%s)",
        total_informative, len(read_hits), unmapped_informative,
        total_reads_scanned,
        _format_elapsed(time.monotonic() - anchor_start),
    )

    if not read_hits:
        return ([], {}, total_informative, {}, unmapped_informative,
                read_sv_meta, kmer_coverage, read_coverage)

    # Sort by chrom, start
    read_hits.sort(key=lambda x: (x[0], x[1]))

    # Cluster into regions
    regions = []
    region_reads = {}
    region_kmers = {}
    current_chrom = read_hits[0][0]
    current_start = read_hits[0][1]
    current_end = read_hits[0][2]
    current_names = {read_hits[0][3]}
    current_kmers = set(read_hits[0][4])

    for chrom, start, end, name, unique_in_read, _is_supp in read_hits[1:]:
        if chrom == current_chrom and start <= current_end + merge_distance:
            current_end = max(current_end, end)
            current_names.add(name)
            current_kmers.update(unique_in_read)
        else:
            region_key = (current_chrom, current_start, current_end)
            regions.append(region_key)
            region_reads[region_key] = current_names
            region_kmers[region_key] = current_kmers
            current_chrom = chrom
            current_start = start
            current_end = end
            current_names = {name}
            current_kmers = set(unique_in_read)

    # Don't forget the last region
    region_key = (current_chrom, current_start, current_end)
    regions.append(region_key)
    region_reads[region_key] = current_names
    region_kmers[region_key] = current_kmers

    logger.info(
        "Clustered %d mapped informative reads into %d regions",
        len(read_hits), len(regions),
    )

    return (regions, region_reads, total_informative, region_kmers,
            unmapped_informative, read_sv_meta, kmer_coverage,
            read_coverage)


def _merge_informative_shards(shard_paths, child_bam, ref_fasta, output_bam,
                              unsorted_path):
    """Combine per-task informative-read shards into a sorted, indexed BAM.

    Shards are read in task order (the windows of the contigs in header
    order, then unplaced unmapped reads).  Each ``(query_name, is_supplementary)`` key — the key
    used to de-duplicate reads while anchoring — is written once, so when
    both mates of a pair are informative on different contigs only the
    first is kept.

    Args:
        shard_paths: Shard BAM paths in task order.  Missing files (tasks
            without informative reads) are skipped.  A read that two
            windows of one contig both report (both mates of a pair
            straddling a window boundary share the key) is written once.
        child_bam: Child BAM/CRAM whose header the output reuses.
        ref_fasta: Reference FASTA for CRAM input (or None).
        output_bam: Final BAM path; a ``.bai`` index is written beside it.
        unsorted_path: Scratch path for the BAM before sorting.

    Returns:
        Number of reads written.
    """
    written = set()
    with pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    ) as src, pysam.AlignmentFile(
        unsorted_path, "wb", header=src.header,
    ) as out:
        for shard in shard_paths:
            if not os.path.exists(shard):
                continue
            with pysam.AlignmentFile(shard) as fh:
                for read in fh.fetch(until_eof=True):
                    key = (read.query_name, read.is_supplementary)
                    if key not in written:
                        written.add(key)
                        out.write(read)
    pysam.sort("-o", output_bam, unsorted_path)
    pysam.index(output_bam)
    os.remove(unsorted_path)
    return len(written)


def _region_mapq_stats(regions, region_reads, read_sv_meta):
    """Per region, the mapping quality of the reads that made it.

    Discovery anchors every primary read whatever its mapping quality,
    because a proband-unique k-mer is novel sequence wherever the read
    landed.  Where it landed still matters to the reader: a region made
    only of MAPQ-0 reads sits in sequence present at several places in
    the reference, where a k-mer that is a sequencing error reaches the
    child-count floor by being seen across the copies.  Nothing is
    filtered on it here; the numbers are written beside each region so
    that a cohort summary can rank or drop such regions after the fact.

    Returns:
        Dict mapping each region tuple to ``(max_mapq, mapq0_reads)``:
        the highest mapping quality among the region's informative reads
        and how many of them have mapping quality 0.  A read with both a
        primary and a supplementary record counts once, by its primary.
    """
    stats = {}
    for region_key in regions:
        best = 0
        n_zero = 0
        for name in region_reads.get(region_key, ()):
            meta = read_sv_meta.get((name, False))
            if meta is None:
                meta = read_sv_meta.get((name, True))
            if meta is None or "mapq" not in meta:
                continue
            mapq = meta["mapq"]
            if mapq > best:
                best = mapq
            if mapq == 0:
                n_zero += 1
        stats[region_key] = (best, n_zero)
    return stats


def _write_bed(regions, region_reads, region_kmers, bed_path,
               region_annotations=None, filters=None, region_mapq=None):
    """Write clustered regions to a BED file with read/k-mer counts and SV annotations.

    Args:
        regions: List of (chrom, start, end) region tuples.
        region_reads: Dict mapping region tuple to set of read names.
        region_kmers: Dict mapping region tuple to set of k-mer strings.
        bed_path: Output BED file path.
        region_annotations: Optional dict of SV annotations per region.
        filters: Optional dict of applied filter parameters to record
            in the file header (e.g. min_supporting_reads,
            min_distinct_kmers, min_distinct_kmers_per_read).
        region_mapq: Optional dict from :func:`_region_mapq_stats`; the
            two mapping-quality columns are 0 when it is missing.

    The last two columns, ``nhf_reads`` and ``nhf``, are written as ``.``
    (unscreened); :func:`nhf.annotate_bed` fills them in from the
    non-human fraction screen, during the run or afterwards.
    """
    with open(bed_path, "w") as fh:
        if filters:
            parts = " ".join(f"{k}={v}" for k, v in sorted(filters.items()))
            fh.write(f"#filters: {parts}\n")
        fh.write(
            "#chrom\tstart\tend\treads\tunique_kmers"
            "\tsplit_reads\tdiscordant_pairs"
            "\tmax_clip_len\tunmapped_mates\tclass"
            "\tbreakpoint_reads\tlarge_indel_reads\tsv_type"
            "\tmax_mapq\tmapq0_reads"
            "\t" + "\t".join(nhf.NHF_COLUMNS) + "\n"
        )
        for chrom, start, end in regions:
            region_key = (chrom, start, end)
            n_reads = len(region_reads.get(region_key, set()))
            n_kmers = len(region_kmers.get(region_key, set()))
            ann = (region_annotations or {}).get(region_key, {})
            split_reads = ann.get("split_reads", 0)
            discordant_pairs = ann.get("discordant_pairs", 0)
            max_clip_len = ann.get("max_clip_len", 0)
            unmapped_mates = ann.get("unmapped_mates", 0)
            region_class = ann.get("class", "SMALL")
            breakpoint_reads = ann.get("breakpoint_reads", 0)
            large_indel_reads = ann.get("large_indel_reads", 0)
            sv_type = ann.get("sv_type", ".")
            max_mapq, mapq0_reads = (region_mapq or {}).get(region_key, (0, 0))
            fh.write(
                f"{chrom}\t{start}\t{end}\t{n_reads}\t{n_kmers}"
                f"\t{split_reads}\t{discordant_pairs}"
                f"\t{max_clip_len}\t{unmapped_mates}\t{region_class}"
                f"\t{breakpoint_reads}\t{large_indel_reads}\t{sv_type}"
                f"\t{max_mapq}\t{mapq0_reads}"
                f"\t{nhf.UNSCREENED}\t{nhf.UNSCREENED}\n"
            )
    logger.info("BED file written: %s (%d regions)", bed_path, len(regions))


def _write_bedgraph(kmer_coverage, bedgraph_path, read_coverage=None,
                    min_reads=3):
    """Write a 4-column bedGraph of novel k-mer reference coverage.

    Adjacent positions with identical coverage are merged into single
    intervals to minimise file size.  When *read_coverage* is provided,
    only positions supported by at least *min_reads* distinct reads are
    included, which dramatically reduces output size at WGS scale.

    Positions are sorted once per chromosome and filtered inline to
    avoid intermediate dict copies at WGS scale.

    Args:
        kmer_coverage: Dict mapping chrom to Counter of reference
            positions with coverage counts (total k-mer base overlaps).
        bedgraph_path: Output file path.
        read_coverage: Optional dict mapping chrom to Counter of reference
            positions with the number of distinct reads touching each
            position.  When provided, positions with fewer than
            *min_reads* reads are filtered out.
        min_reads: Minimum number of distinct reads with at least one
            de novo k-mer required at a position for it to be included
            in the bedGraph (default: 3).
    """
    total_intervals = 0
    total_filtered = 0
    with open(bedgraph_path, "w") as fh:
        fh.write(
            f"#track type=bedGraph "
            f"description=\"De novo k-mer coverage (unique k-mer base "
            f"overlaps per position, min_reads>={min_reads})\"\n"
        )
        for chrom in sorted(kmer_coverage):
            positions = kmer_coverage[chrom]
            if not positions:
                continue
            rc = read_coverage.get(chrom, {}) if read_coverage else None
            sorted_pos = sorted(positions)
            run_start = None
            run_val = None
            run_end = None
            for pos in sorted_pos:
                if rc is not None and rc.get(pos, 0) < min_reads:
                    total_filtered += 1
                    if run_start is not None:
                        fh.write(
                            f"{chrom}\t{run_start}\t{run_end}\t{run_val}\n"
                        )
                        total_intervals += 1
                        run_start = None
                    continue
                val = positions[pos]
                if run_start is None:
                    run_start = pos
                    run_val = val
                    run_end = pos + 1
                elif pos == run_end and val == run_val:
                    run_end = pos + 1
                else:
                    fh.write(
                        f"{chrom}\t{run_start}\t{run_end}\t{run_val}\n"
                    )
                    total_intervals += 1
                    run_start = pos
                    run_val = val
                    run_end = pos + 1
            if run_start is not None:
                fh.write(
                    f"{chrom}\t{run_start}\t{run_end}\t{run_val}\n"
                )
                total_intervals += 1
    if total_filtered:
        logger.info(
            "bedGraph file written: %s (%d intervals, %d positions "
            "filtered by min_reads=%d)",
            bedgraph_path, total_intervals, total_filtered, min_reads,
        )
    else:
        logger.info(
            "bedGraph file written: %s (%d intervals)",
            bedgraph_path, total_intervals,
        )


def _write_read_coverage_bed(kmer_coverage, read_coverage, bed_path,
                             min_reads=3):
    """Write a BED file with per-position read count and average k-mers per read.

    For every reference position where at least *min_reads* distinct reads
    carry a de novo k-mer, writes a BED interval with two value columns:

    - **read_count**: number of distinct reads touching the position with
      at least one de novo k-mer.
    - **avg_kmers_per_read**: ``kmer_coverage / read_count`` rounded to
      one decimal — the average number of unique k-mer overlaps per read
      at this position, measuring per-read k-mer signal density.

    Adjacent positions with identical (read_count, avg_kmers) are merged
    into intervals.

    Args:
        kmer_coverage: Dict mapping chrom to Counter of reference
            positions with k-mer base overlap counts.
        read_coverage: Dict mapping chrom to Counter of reference
            positions with distinct read counts.
        bed_path: Output BED file path.
        min_reads: Minimum distinct reads at a position (default: 3).
    """
    total_intervals = 0
    with open(bed_path, "w") as fh:
        fh.write(
            f"#track description=\"De novo k-mer read support "
            f"(min_reads>={min_reads})\"\n"
            f"#chrom\tstart\tend\tread_count\tavg_kmers_per_read\n"
        )
        for chrom in sorted(read_coverage):
            rc = read_coverage[chrom]
            kc = kmer_coverage.get(chrom, {})
            # Build filtered positions
            filtered = {}
            for pos, n_reads in rc.items():
                if n_reads >= min_reads:
                    avg_k = round(kc.get(pos, 0) / n_reads, 1)
                    filtered[pos] = (n_reads, avg_k)
            if not filtered:
                continue
            sorted_pos = sorted(filtered)
            run_start = sorted_pos[0]
            run_val = filtered[run_start]
            run_end = run_start + 1
            for pos in sorted_pos[1:]:
                val = filtered[pos]
                if pos == run_end and val == run_val:
                    run_end = pos + 1
                else:
                    fh.write(
                        f"{chrom}\t{run_start}\t{run_end}"
                        f"\t{run_val[0]}\t{run_val[1]}\n"
                    )
                    total_intervals += 1
                    run_start = pos
                    run_val = val
                    run_end = pos + 1
            fh.write(
                f"{chrom}\t{run_start}\t{run_end}"
                f"\t{run_val[0]}\t{run_val[1]}\n"
            )
            total_intervals += 1
    logger.info(
        "Read coverage BED written: %s (%d intervals)",
        bed_path, total_intervals,
    )


#: Minimum mapping quality for an alignment to count as SV evidence or to
#: link two regions: the informative read, its supplementary alignment
#: (SA tag) and, when the MQ tag records it, its mate.
_MIN_LINK_MAPQ = 20

#: Soft clips starting within this many bp of one another mark one
#: breakpoint.
_SV_CLIP_TOLERANCE = 5

#: Reads clipped on their right may end up to this many bp after reads
#: clipped on their left start and still flank one insertion: a
#: target-site duplication, or bases matching the insertion by chance.
_INS_MAX_OVERLAP = 50

#: Per-region evidence counts; two molecules of any one kind make an SV.
_SV_EVIDENCE = (
    "split_reads", "discordant_pairs", "unmapped_mates",
    "breakpoint_reads", "large_indel_reads",
)


def _unique_majority(values):
    """Return the most common value, or None if absent or tied."""
    ranked = collections.Counter(values).most_common(2)
    if not ranked or (len(ranked) == 2 and ranked[0][1] == ranked[1][1]):
        return None
    return ranked[0][0]


def _largest_clip_cluster(clips, tolerance=_SV_CLIP_TOLERANCE):
    """Return the most molecules whose clips start within *tolerance* bp.

    Args:
        clips: List of (reference position, read name) pairs.

    Returns:
        (molecules, position): the count, and the median clip position of
        that cluster (None without clips).
    """
    clips = sorted(clips)
    in_window = collections.Counter()
    best, best_pos, lo = 0, None, 0
    for hi, (pos, qname) in enumerate(clips):
        in_window[qname] += 1
        while pos - clips[lo][0] > tolerance:
            dropped = clips[lo][1]
            in_window[dropped] -= 1
            if not in_window[dropped]:
                del in_window[dropped]
            lo += 1
        if len(in_window) > best:
            best, best_pos = len(in_window), clips[(lo + hi) // 2][0]
    return best, best_pos


def _clipped_insertion(clips):
    """Whether reads are clipped from both sides at one point.

    Reads reaching into an insertion from its left are clipped on their
    right ("+") where it starts, and reads from its right on their left
    ("-") at the same point, or earlier by a target-site duplication.
    Needs two molecules on each side.

    Args:
        clips: List of (reference position, side, read name).
    """
    n_right, right = _largest_clip_cluster(
        [(pos, qname) for pos, side, qname in clips if side == "+"])
    n_left, left = _largest_clip_cluster(
        [(pos, qname) for pos, side, qname in clips if side == "-"])
    return (n_right >= 2 and n_left >= 2
            and -_SV_CLIP_TOLERANCE <= right - left <= _INS_MAX_OVERLAP)


def _split_is_insertion(part_a, part_b):
    """Whether two collinear parts of a split read flank an insertion.

    Each part is (reference start, reference end, read start, read end),
    in reference orientation; *part_a* is clipped on its right and
    *part_b* on its left.  More read than reference lies between them for
    an insertion; more reference for a deletion; and for a tandem
    duplication, part_b starts on the reference before part_a ends, by
    more than the read between them.
    """
    ref_gap = part_b[0] - part_a[1]
    read_gap = part_b[2] - part_a[3]
    return read_gap > abs(ref_gap)


def _annotate_and_link_from_metadata(regions, read_sv_meta, link_slack=0):
    """Annotate regions and link breakpoints using pre-collected metadata.

    Uses per-read SV metadata collected during the anchoring scan
    (Module 3), so no additional BAM I/O is needed.

    Evidence is counted per molecule (read name), at most once per count
    per region, and only from alignments with MAPQ >= ``_MIN_LINK_MAPQ``:
    reads of sequence missing from the reference that resembles a repeat
    (a new mobile-element copy, say) align to the repeat's reference
    copies with low MAPQ, and would pile up false evidence there.
    Pair-level evidence (split alignment, unmapped mate, discordant pair)
    counts in every region where the molecule has such an informative
    alignment.  Soft clips and CIGAR indels of at least 50 bp count in
    the region of the alignment that shows them; ``breakpoint_reads`` is
    the most molecules clipped at one breakpoint, or 0 when that is fewer
    than two.  An unmapped informative read counts as an unmapped mate
    in the region where it is placed (unless the MQ tag gives its mate a
    low MAPQ).  ``max_clip_len`` is the longest clip of any alignment.

    Two regions are linked by a molecule with informative alignments in
    both, by a supplementary alignment (SA tag) of a read in one that
    falls in the other, or by a discordant pair with one end in each.
    Linking alignments need MAPQ >= ``_MIN_LINK_MAPQ``, and SA and mate
    positions may lie up to *link_slack* bp outside the target region.

    Each joining molecule also gives the breakpoint orientation: from
    which side of each split-read segment is clipped, or from the read
    strands of a discordant pair (FR libraries).  A forward-reverse pair
    with both ends in one region gives none, as its insert may be too
    long (a deletion) or too short (an insertion).  The molecules joining
    two regions vote for an SV type (DEL, DUP, INV; BND across
    chromosomes), and each orientation of the majority type is one link:
    a junction, with the molecules showing it.  Both junctions of a
    balanced inversion or reciprocal translocation are thus reported.
    With no majority, the regions get one link with unknown orientation
    (INTRA, or BND across chromosomes).  A split read whose parts are
    collinear is an insertion (INS) rather than a deletion or tandem
    duplication when more of the read than of the reference lies between
    them.  A region's ``sv_type`` is the majority over its molecules,
    including CIGAR indels (DEL, INS) and joins within the region; with
    no votes, INS when reads are clipped from both sides at one point,
    otherwise ``.``.

    Args:
        regions: List of (chrom, start, end) tuples.
        read_sv_meta: Dict mapping (query_name, is_supplementary) to the
            per-read metadata recorded by ``_process_informative_read``.
        link_slack: How far (bp) an SA or mate position may fall outside
            a region and still link to it.

    Returns:
        (annotations, links) where:
        - annotations: Dict mapping region tuple to a dict with the counts
          in ``_SV_EVIDENCE``, ``max_clip_len`` and ``sv_type``.
        - links: List of dicts, one per junction, with keys: region_a,
          region_b, supporting_reads (read names), strands (a (strand1,
          strand2) tuple, or None when unknown), sv_type_hint.
    """
    regions_by_chrom = collections.defaultdict(list)
    for region in sorted(regions):
        regions_by_chrom[region[0]].append(region)
    starts_by_chrom = {
        chrom: [r[1] for r in rlist]
        for chrom, rlist in regions_by_chrom.items()
    }

    def region_at(chrom, pos, slack=0):
        """The region containing *pos*, else the nearest within *slack* bp."""
        rlist = regions_by_chrom.get(chrom)
        if not rlist:
            return None
        idx = bisect.bisect_right(starts_by_chrom[chrom], pos) - 1
        best, best_dist = None, slack + 1
        if idx >= 0:
            before = rlist[idx]
            dist = 0 if pos < before[2] else pos - before[2] + 1
            if dist < best_dist:
                best, best_dist = before, dist
        if idx + 1 < len(rlist):
            after = rlist[idx + 1]
            if after[1] - pos < best_dist:
                best = after
        return best

    # The regions where each molecule has a confidently placed
    # informative alignment; only there does it count as evidence
    placed = collections.defaultdict(set)
    for (qname, _is_supp), meta in read_sv_meta.items():
        if "pos" in meta and meta["mapq"] >= _MIN_LINK_MAPQ:
            region = region_at(*meta["pos"])
            if region is not None:
                placed[qname].add(region)

    annotations = {
        r: {**{name: 0 for name in _SV_EVIDENCE}, "max_clip_len": 0}
        for r in regions
    }

    # ── Annotation from metadata ──
    # Each molecule adds at most one to each count per region, even when
    # both its primary and supplementary alignments are informative and
    # carry the same flags; otherwise one molecule could reach the SV
    # threshold of two on its own.
    counted = set()  # (qname, region, count name)

    def count(qname, region, name):
        if (qname, region, name) not in counted:
            counted.add((qname, region, name))
            annotations[region][name] += 1

    clips_by_region = collections.defaultdict(list)  # (pos, side, qname)
    type_votes = collections.defaultdict(dict)  # region -> {qname: SV type}
    for (qname, _is_supp), meta in read_sv_meta.items():
        if "placed_unmapped" in meta:
            region = region_at(*meta["placed_unmapped"])
            mate_mapq = meta.get("mate_mapq")
            if region is not None and (mate_mapq is None
                                       or mate_mapq >= _MIN_LINK_MAPQ):
                count(qname, region, "unmapped_mates")
            continue

        own = region_at(*meta["pos"]) if "pos" in meta else None
        if own is not None:
            ann = annotations[own]
            ann["max_clip_len"] = max(ann["max_clip_len"], meta["max_clip"])
        if meta["mapq"] < _MIN_LINK_MAPQ:
            continue

        evidence = []
        if meta["has_sa"]:
            evidence.append("split_reads")
        if meta["is_paired"]:
            if meta["mate_is_unmapped"]:
                evidence.append("unmapped_mates")
            elif not meta["is_proper_pair"]:
                evidence.append("discordant_pairs")
        for region in placed.get(qname, ()):
            for name in evidence:
                count(qname, region, name)

        if own is None:
            continue
        for pos, side in meta["clips"]:
            clips_by_region[own].append((pos, side, qname))
        if meta["large_indel"]:
            count(qname, own, "large_indel_reads")
            type_votes[own].setdefault(qname, meta["large_indel"])

    # One long clip is often an adapter or a low-quality tail; a
    # breakpoint needs at least two molecules clipped at it.
    for region, clips in clips_by_region.items():
        n_clipped, _pos = _largest_clip_cluster(
            [(pos, qname) for pos, _side, qname in clips])
        if n_clipped >= 2:
            annotations[region]["breakpoint_reads"] = n_clipped

    # ── Linking and SV type ──
    # Each molecule joining two breakends (in two regions, or both in one)
    # supports the link between their regions and votes for an SV type
    # from the breakpoint orientation.  A breakend is (chrom, position,
    # side): the breakpoint follows a "+" alignment and precedes a "-" one.
    bridges = {}
    link_votes = collections.defaultdict(dict)  # link -> {qname: (type, strands)}

    def breakend(chrom, start, end, side):
        return (chrom, end if side == "+" else start, side)

    def join(end_a, region_a, end_b, region_b, qname, pair=False,
             insertion=False):
        if region_a is None or region_b is None:
            return
        # By position; on a tie, a "+" end comes first
        if (end_b[:2], end_b[2] != "+") < (end_a[:2], end_a[2] != "+"):
            end_a, end_b = end_b, end_a
            region_a, region_b = region_b, region_a
        strands = (end_a[2], end_b[2]) if end_a[2] and end_b[2] else None
        sv_type = "INS" if insertion else _infer_sv_type(
            region_a, region_b, strands)
        if pair and sv_type == "DEL" and region_a == region_b:
            # A forward-reverse pair within one region is either too far
            # apart (a deletion) or too close (an insertion, e.g. one
            # longer than the reads); the insert size alone can't say which
            return
        if region_a != region_b:
            key = tuple(sorted((region_a, region_b)))
            bridges.setdefault(key, set()).add(qname)
            if strands:
                if key[0] != region_a:
                    strands = strands[::-1]
                link_votes[key].setdefault(qname, (sv_type, strands))
        if sv_type != "INTRA":
            for region in (region_a, region_b):
                type_votes[region].setdefault(qname, sv_type)

    molecule_ends = collections.defaultdict(list)
    for (qname, _is_supp), meta in read_sv_meta.items():
        if "pos" not in meta or meta["mapq"] < _MIN_LINK_MAPQ:
            continue
        chrom, start = meta["pos"]
        own = region_at(chrom, start)
        if own is None:
            continue
        own_end = breakend(chrom, start, meta["end"], meta["clip_side"])
        molecule_ends[qname].append((own_end, own))

        for sa_entry in (meta["sa_str"] or "").rstrip(";").split(";"):
            parts = sa_entry.split(",")
            if len(parts) < 5:
                continue
            try:
                sa_start = int(parts[1]) - 1  # 1-based to 0-based
                sa_mapq = int(parts[4])
            except ValueError:
                continue
            if sa_mapq < _MIN_LINK_MAPQ:
                continue
            sa_cigar = _parse_cigar(parts[3])
            sa_end = sa_start + sum(
                length for op, length in sa_cigar if op in (0, 2, 3, 7, 8)
            )
            sa_side = _clip_side(sa_cigar)
            # Collinear parts (one chromosome and strand, one clipped on
            # each side) flank a deletion, a tandem duplication or an
            # insertion; the read between them tells which
            insertion = False
            if parts[0] == chrom and (parts[2] == "-") == meta["is_reverse"]:
                own_part = (start, meta["end"], *meta["qspan"])
                sa_part = (sa_start, sa_end, *_query_span(sa_cigar))
                sides = (meta["clip_side"], sa_side)
                if sides == ("+", "-"):
                    insertion = _split_is_insertion(own_part, sa_part)
                elif sides == ("-", "+"):
                    insertion = _split_is_insertion(sa_part, own_part)
            join(own_end, own,
                 breakend(parts[0], sa_start, sa_end, sa_side),
                 region_at(parts[0], sa_start, link_slack), qname,
                 insertion=insertion)

        # Paired-end orientation (FR libraries): the breakpoint lies after
        # a forward read and before a reverse one.  The two reads are
        # ordered by where they start, so a forward read overlapping its
        # reverse mate still reads as forward-reverse.
        if meta["mate"] is not None:
            mate_chrom, mate_pos, mate_mapq, mate_reverse = meta["mate"]
            if mate_mapq is None or mate_mapq >= _MIN_LINK_MAPQ:
                join((chrom, start, "-" if meta["is_reverse"] else "+"), own,
                     (mate_chrom, mate_pos, "-" if mate_reverse else "+"),
                     region_at(mate_chrom, mate_pos, link_slack), qname,
                     pair=True)

    # A molecule with informative alignments in two regions (e.g. both
    # parts of a split read) joins them.
    for qname, ends in molecule_ends.items():
        for i, (end_a, region_a) in enumerate(ends):
            for end_b, region_b in ends[i + 1:]:
                if region_a != region_b:
                    join(end_a, region_a, end_b, region_b, qname)

    for region in regions:
        votes = type_votes[region].values()
        sv_type = _unique_majority(votes)
        if not votes and _clipped_insertion(clips_by_region[region]):
            # No alignment says where the clipped sequence belongs: new
            # sequence, such as a mobile element or one longer than a read
            sv_type = "INS"
        annotations[region]["sv_type"] = sv_type or "."

    # Build links, one per junction: each breakpoint orientation of the
    # majority type among the molecules joining two regions.  Both
    # junctions of a balanced inversion (+ +, - -) or reciprocal
    # translocation join the same two regions, so each gets its own link;
    # molecules of a minority type are outvoted.
    links = []
    for key in sorted(bridges):
        votes = link_votes.get(key, {})
        sv_type = _unique_majority(t for t, _ in votes.values())
        junctions = collections.defaultdict(set)
        for qname, (vote_type, strands) in votes.items():
            if vote_type == sv_type:
                junctions[strands].add(qname)
        if not junctions:
            # No orientation known, or the molecules disagree on the type
            links.append({
                "region_a": key[0],
                "region_b": key[1],
                "supporting_reads": bridges[key],
                "strands": None,
                "sv_type_hint": _infer_sv_type(*key),
            })
        for strands, qnames in sorted(junctions.items(),
                                      key=lambda j: (-len(j[1]), j[0])):
            links.append({
                "region_a": key[0],
                "region_b": key[1],
                "supporting_reads": qnames,
                "strands": strands,
                "sv_type_hint": sv_type,
            })

    return annotations, links


def _write_bedpe(links, bedpe_path):
    """Write linked SV breakpoint pairs to a BEDPE file, one per junction.

    Uses the standard BEDPE layout, so tools such as bedtools can read it:
    name (``SV_n``) and score (supporting reads) in columns 7–8, the
    breakpoint orientation in columns 9–10 (``.`` when unknown), and the
    SV type as an extra column 11.

    Args:
        links: List of link dicts from ``_annotate_and_link_from_metadata()``.
        bedpe_path: Output BEDPE file path.
    """
    with open(bedpe_path, "w") as fh:
        fh.write(
            "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2"
            "\tsv_id\tsupporting_reads\tstrand1\tstrand2\tsv_type\n"
        )
        for idx, link in enumerate(links, 1):
            ra = link["region_a"]
            rb = link["region_b"]
            n_support = len(link["supporting_reads"])
            strand1, strand2 = link.get("strands") or (".", ".")
            sv_type = link["sv_type_hint"]
            fh.write(
                f"{ra[0]}\t{ra[1]}\t{ra[2]}"
                f"\t{rb[0]}\t{rb[1]}\t{rb[2]}"
                f"\tSV_{idx}\t{n_support}\t{strand1}\t{strand2}"
                f"\t{sv_type}\n"
            )
    logger.info("BEDPE file written: %s (%d links)", bedpe_path, len(links))


def _classify_regions(regions, region_annotations, sv_links):
    """Assign SV classification to each region.

    - ``SV``: at least two molecules show one kind of evidence: one of
      ``_SV_EVIDENCE`` (split alignments, discordant pairs, unmapped
      mates, soft clips at one breakpoint, CIGAR indels of 50 bp or
      more), or links to one other region (over the pair's junctions)
    - ``SMALL``: none of that evidence and no links
    - ``AMBIGUOUS``: otherwise (evidence from a single molecule)

    Updates region_annotations in place with a ``class`` key.
    """
    # Molecules joining each pair of regions, over the pair's junctions
    pair_reads = collections.defaultdict(set)
    for link in sv_links:
        pair_reads[link["region_a"], link["region_b"]].update(
            link["supporting_reads"])
    linking_reads = collections.Counter()
    for pair, qnames in pair_reads.items():
        for region in pair:
            linking_reads[region] = max(linking_reads[region], len(qnames))

    for region_key in regions:
        ann = region_annotations.get(region_key, {})
        strongest = max(
            [ann.get(name, 0) for name in _SV_EVIDENCE]
            + [linking_reads[region_key]]
        )
        if strongest >= 2:
            ann["class"] = "SV"
        elif strongest == 0:
            ann["class"] = "SMALL"
        else:
            ann["class"] = "AMBIGUOUS"
        region_annotations[region_key] = ann


def _parse_candidate_summary(summary_path, dka_dkt_min=0.25, dka_min=10):
    """Parse a VCF-mode summary.txt and return high-quality de novo candidates.

    Filters the Per-Variant Results table for candidates meeting both
    ``DKA_DKT > dka_dkt_min`` and ``DKA > dka_min``.

    Args:
        summary_path: Path to a VCF-mode summary.txt file.
        dka_dkt_min: Minimum DKA_DKT proportion (exclusive).
        dka_min: Minimum DKA count (exclusive).

    Returns:
        List of dicts with keys: chrom, pos (1-based), ref, alt, dka,
        dka_dkt, call.
    """
    candidates = []
    in_table = False
    with open(summary_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.strip().startswith("Variant") and "DKU" in line:
                in_table = True
                continue
            if in_table and line.strip().startswith("-------"):
                continue
            if in_table and line.strip() == "":
                break
            if in_table and line.strip().startswith("="):
                break
            if in_table:
                parts = line.split()
                if len(parts) < 12:
                    continue
                # Columns: Variant R>A DKU DKT DKA DKU_DKT DKA_DKT ...
                # e.g. "chr11:55003995" "T>C" "21" "46" "21" "0.4565" "0.4565" ...
                variant = parts[0]  # chr:pos
                ref_alt = parts[1]  # R>A
                dku = int(parts[2])
                dkt = int(parts[3])
                dka = int(parts[4])
                dku_dkt = float(parts[5])
                dka_dkt = float(parts[6])
                call = parts[-1]
                chrom, pos_str = variant.rsplit(":", 1)
                pos = int(pos_str)
                ref, alt = ref_alt.split(">")

                if dka_dkt > dka_dkt_min and dka > dka_min:
                    candidates.append({
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "dka": dka,
                        "dka_dkt": dka_dkt,
                        "call": call,
                    })
    return candidates


def _compare_candidates_to_regions(candidates, regions):
    """Compare high-quality VCF candidates to discovery regions.

    For each candidate, checks whether its 1-based position falls within
    any discovery region (0-based half-open BED coordinates).

    Args:
        candidates: List of dicts from ``_parse_candidate_summary()``.
        regions: List of (chrom, start, end) tuples (0-based, half-open).

    Returns:
        List of dicts with original candidate fields plus:
        - captured (bool): Whether the candidate falls in a region.
        - region (str or None): The matching region label, if captured.
    """
    results = []
    for cand in candidates:
        captured = False
        match_region = None
        for chrom, start, end in regions:
            if cand["chrom"] == chrom and start < cand["pos"] <= end:
                captured = True
                match_region = f"{chrom}:{start + 1}-{end}"
                break
        results.append({**cand, "captured": captured, "region": match_region})
    return results


def _load_dnm_regions(path):
    """Load known de novo events for :func:`_evaluate_dnm_regions`.

    The file is tab-separated, one event per line: chrom, 1-based
    position of the event start, size in bp (``.`` or ``0`` when
    unknown) and an event-type label.  Blank lines and lines starting
    with ``#`` are skipped.

    Returns:
        List of ``(chrom, pos, size_or_None, event_type)`` tuples.

    Raises:
        ValueError: If a line is malformed (the message gives its line
            number) or the file lists no events.
    """
    regions = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 4:
                raise ValueError(
                    f"{path}:{lineno}: expected 4 tab-separated columns "
                    f"(chrom, pos, size, event_type), found {len(parts)}"
                )
            chrom, pos, size, event_type = parts
            try:
                pos = int(pos)
                size = None if size == "." else int(size)
            except ValueError:
                raise ValueError(
                    f"{path}:{lineno}: pos and size must be integers "
                    f"(size may be '.')"
                ) from None
            if pos < 1 or (size is not None and size < 0):
                raise ValueError(
                    f"{path}:{lineno}: pos must be >= 1 and size >= 0"
                )
            regions.append((chrom, pos, size or None, event_type))
    if not regions:
        raise ValueError(f"{path}: no de novo events listed")
    return regions


class _RegionIndex:
    """The discovery regions of one run, indexed for interval queries.

    Regions are grouped by chromosome and sorted by start; a running
    maximum of their ends lets a query find every region overlapping an
    interval with two bisections, so that evaluating thousands of loci
    against hundreds of thousands of regions stays cheap.
    """

    def __init__(self, discovery_regions):
        buckets = collections.defaultdict(list)
        for key in discovery_regions:
            buckets[key[0]].append(key)
        self._by_chrom = {}
        for chrom, keys in buckets.items():
            keys.sort(key=lambda k: (k[1], k[2]))
            max_ends = []
            running = -1
            for key in keys:
                running = max(running, key[2])
                max_ends.append(running)
            self._by_chrom[chrom] = ([k[1] for k in keys], max_ends, keys)

    def overlapping(self, chrom, start, end):
        """Region keys overlapping ``[start, end)`` (0-based, half-open)."""
        entry = self._by_chrom.get(chrom)
        if entry is None:
            return []
        starts, max_ends, keys = entry
        hi = bisect.bisect_left(starts, end)
        lo = bisect.bisect_right(max_ends, start)
        return [k for k in keys[lo:hi] if k[2] > start]


def _evaluate_dnm_regions(discovery_regions, region_detail, dnm_regions,
                          slack=0):
    """Evaluate the discovered regions against a list of known events.

    Two questions are answered for each event (from ``--dnm-regions`` or
    ``--dnm-negative-regions``), and they are kept apart on purpose:

    * **Is there a region anywhere in the event?**  (``any_region``)  A
      discovery region overlapping the event's span, or lying within
      *slack* bp of it.  On a real genome this is weak evidence: with
      a region every few kilobases, a span of some kilobases holds one
      by chance, and a long event always does.
    * **Is a breakpoint supported?**  (``breakpoint_support``, the
      count of supported ends, and ``detected`` when it is at least 1)
      Split reads, discordant pairs or breakpoint-clipped reads in the
      regions within *slack* bp of an END of the event: its start, and
      for a sized event its end.  This is the evidence a junction leaves
      and the one that should be compared with the random-locus null
      the pipeline reports beside it.

    The slack matters for deletions: their breakpoint regions flank the
    deleted interval rather than overlap it, and curated breakpoints can
    be off by a few bases.

    Args:
        discovery_regions: List of (chrom, start, end) tuples (0-based,
            half-open) from the discovery BED.
        region_detail: List of dicts with per-region metrics (the
            ``regions`` array from the discovery metrics JSON).
        dnm_regions: List of (chrom, pos, size_or_None, event_type)
            tuples, as returned by :func:`_load_dnm_regions`.  *pos* is
            1-based; each event covers *size* bp from *pos* (1 bp when
            the size is unknown).
        slack: How far (bp) outside an event or an end a region may lie
            and still count (the pipeline uses ``--cluster-distance``).

    Returns:
        List of dicts, one per event, with keys:

        - locus (str): ``chrom:pos``
        - event_type (str): Event type label.
        - event_size (int or None): Expected event size in bp.
        - any_region (bool): A region overlaps the span or lies within
          the slack of it.
        - ends (int): Ends the event has: 2 for a sized event, else 1.
        - breakpoint_support (int): Ends with split, discordant or
          breakpoint-clipped reads within the slack.
        - breakpoints (dict): Per end (``start``, ``end``): the regions
          within the slack and their split, discordant and clipped read
          counts, and whether the end is supported.
        - detected (bool): ``breakpoint_support >= 1``.
        - discovery_regions (list of str): Labels of the regions in the
          span (the ``any_region`` matches).
        - total_reads, total_unique_kmers, max_clip_len, unmapped_mates,
          discordant_pairs, split_reads: Sums over those regions.
        - sv_class (str): Most severe class among them, ``NONE`` if none.
        - kmer_signal (float): ``total_unique_kmers / span_bp``.
        - assessment (str): ``SUPPORTED_BOTH_ENDS``, ``SUPPORTED``,
          ``REGION_ONLY`` or ``NOT_DETECTED``.
    """
    detail_by_key = {}
    for rd in region_detail:
        detail_by_key[(rd["chrom"], rd["start"], rd["end"])] = rd
    index = _RegionIndex(discovery_regions)
    class_priority = {"SV": 3, "AMBIGUOUS": 2, "SMALL": 1}

    results = []
    for chrom, pos, size, event_type in dnm_regions:
        # 1-based pos → 0-based half-open, like the discovery regions
        dnm_start = pos - 1
        dnm_end = dnm_start + (size if size else 1)  # point if no size

        # ── A region anywhere in the span (or within the slack of it) ──
        matches = index.overlapping(chrom, dnm_start - slack, dnm_end + slack)
        total_reads = total_kmers = max_clip = 0
        total_unmapped = total_discordant = total_split = 0
        region_labels = []
        sv_classes = []
        span_start, span_end = dnm_start, dnm_end
        for m_key in matches:
            rd = detail_by_key.get(m_key, {})
            total_reads += rd.get("reads", 0)
            total_kmers += rd.get("unique_kmers", 0)
            max_clip = max(max_clip, rd.get("max_clip_len", 0))
            total_unmapped += rd.get("unmapped_mates", 0)
            total_discordant += rd.get("discordant_pairs", 0)
            total_split += rd.get("split_reads", 0)
            sv_classes.append(rd.get("class", "SMALL"))
            region_labels.append(f"{m_key[0]}:{m_key[1] + 1}-{m_key[2]}")
            span_start = min(span_start, m_key[1])
            span_end = max(span_end, m_key[2])
        any_region = bool(matches)
        span_bp = max(span_end - span_start, 1)
        kmer_signal = total_kmers / span_bp if any_region else 0.0
        sv_class = (max(sv_classes, key=lambda c: class_priority.get(c, 0))
                    if sv_classes else "NONE")

        # ── Breakpoint evidence at each end ────────────────────────────
        ends = [("start", dnm_start)]
        if size:
            ends.append(("end", dnm_end))
        breakpoints = {}
        n_supported = 0
        for label, bp in ends:
            # The junction lies between bases bp - 1 and bp: a region that
            # ends on the one or starts on the other touches it before any
            # slack is added (a deletion's flanks never overlap it).
            near = index.overlapping(chrom, bp - 1 - slack, bp + 1 + slack)
            split = discordant = clipped = 0
            for key in near:
                rd = detail_by_key.get(key, {})
                split += rd.get("split_reads", 0)
                discordant += rd.get("discordant_pairs", 0)
                clipped += rd.get("breakpoint_reads", 0)
            supported = (split + discordant + clipped) > 0
            n_supported += supported
            breakpoints[label] = {
                "regions": len(near),
                "split_reads": split,
                "discordant_pairs": discordant,
                "breakpoint_reads": clipped,
                "supported": supported,
            }
        detected = n_supported > 0
        if detected and size and n_supported == len(ends):
            assessment = "SUPPORTED_BOTH_ENDS"
        elif detected:
            assessment = "SUPPORTED"
        elif any_region:
            assessment = "REGION_ONLY"
        else:
            assessment = "NOT_DETECTED"

        results.append({
            "locus": f"{chrom}:{pos}",
            "event_type": event_type,
            "event_size": size,
            "any_region": any_region,
            "ends": len(ends),
            "breakpoint_support": n_supported,
            "breakpoints": breakpoints,
            "detected": detected,
            "discovery_regions": region_labels,
            "total_reads": total_reads,
            "total_unique_kmers": total_kmers,
            "max_clip_len": max_clip,
            "unmapped_mates": total_unmapped,
            "discordant_pairs": total_discordant,
            "split_reads": total_split,
            "sv_class": sv_class,
            "kmer_signal": round(kmer_signal, 4),
            "assessment": assessment,
        })

    return results


#: Seed of the random loci the known-event null draws (see
#: :func:`_null_dnm_regions`): the same run gives the same null.
_DNM_NULL_SEED = 1


def _null_dnm_regions(dnm_regions, contig_lengths, per_event, seed=_DNM_NULL_SEED):
    """Random loci matched to a list of events: the null of the evaluation.

    For each event, *per_event* loci of the same size and type are drawn
    uniformly along the same chromosome (a chromosome absent from
    *contig_lengths*, or shorter than the event, draws none).  Evaluating
    them like the events says how often chance alone puts a region, or
    breakpoint evidence, at a locus of that size on this genome.
    """
    rng = random.Random(seed)
    null = []
    for chrom, pos, size, event_type in dnm_regions:
        length = contig_lengths.get(chrom, 0)
        span = size or 1
        if length <= span + 1:
            continue
        for _ in range(per_event):
            null.append((chrom, rng.randrange(1, length - span + 1), size,
                         event_type))
    return null


def _summarise_dnm_evaluation(loci):
    """The headline counts of an evaluation: what the metrics report."""
    n = len(loci)
    supported = sum(1 for e in loci if e["detected"])
    any_region = sum(1 for e in loci if e["any_region"])
    both = sum(1 for e in loci if e["ends"] == 2 and e["breakpoint_support"] == 2)
    return {
        "total_loci": n,
        "detected": supported,
        "detection_rate": supported / n if n else 0.0,
        "both_ends": both,
        "any_region": any_region,
        "any_region_rate": any_region / n if n else 0.0,
    }


def _contig_lengths(bam_path, ref_fasta=None):
    """Contig lengths from an alignment file's header."""
    with pysam.AlignmentFile(
        bam_path, reference_filename=ref_fasta if ref_fasta else None,
    ) as bam:
        return dict(zip(bam.references, bam.lengths))


def _write_discovery_summary(summary_path, regions, region_reads,
                             region_kmers, metrics,
                             candidate_comparison=None,
                             region_annotations=None,
                             dnm_evaluation=None, nhf_block=None):
    """Write a human-readable summary for the discovery pipeline.

    Analogous to ``_write_summary()`` in VCF mode, but reports
    per-region statistics instead of per-variant annotations.

    Args:
        summary_path: Output file path for the summary text.
        regions: List of (chrom, start, end) tuples (0-based, half-open).
        region_reads: Dict mapping region tuple to set of read names.
        region_kmers: Dict mapping region tuple to set of proband-unique
            k-mer strings observed in that region's reads.
        metrics: Dict with overall discovery pipeline statistics.
        candidate_comparison: Optional list of comparison dicts from
            ``_compare_candidates_to_regions()``.
        region_annotations: Optional dict mapping region tuple to SV
            annotation dict.
        dnm_evaluation: Optional list of dicts from
            ``_evaluate_dnm_regions()``.
        nhf_block: Optional non-human fraction screen block (see
            :mod:`kmer_denovo_filter.discovery.nhf`); its section comes
            last, where :func:`nhf.update_summary` can rewrite it.
    """
    n_regions = metrics["candidate_regions"]
    n_reads_total = metrics["informative_reads"]
    n_unmapped = metrics.get("unmapped_informative_reads", 0)
    n_unique_kmers = metrics["proband_unique_kmers"]
    n_candidates = metrics["child_candidate_kmers"]
    n_non_ref = metrics["non_ref_kmers"]

    lines = []
    lines.append("=" * 60)
    lines.append("  kmer-denovo  —  Discovery Mode Summary")
    lines.append("=" * 60)
    lines.append("")
    lines.append("K-mer Filtering")
    lines.append("-" * 40)
    lines.append(f"  Child candidate k-mers:      {n_candidates:>8}")
    lines.append(f"  Non-reference k-mers:        {n_non_ref:>8}")
    lines.append(f"  Proband-unique k-mers:       {n_unique_kmers:>8}")
    lines.append("")
    lines.append("Region Counts")
    lines.append("-" * 40)
    lines.append(f"  Candidate regions:           {n_regions:>8}")
    lines.append(f"  Total informative reads:     {n_reads_total:>8}")
    if n_unmapped > 0:
        lines.append(f"    (unmapped informative):     {n_unmapped:>8}")
    lines.append("")

    if regions:
        reads_per_region = [
            len(region_reads.get(r, set())) for r in regions
        ]
        kmers_per_region = [
            len(region_kmers.get(r, set())) for r in regions
        ]
        sizes = [end - start for _, start, end in regions]

        lines.append("Region Statistics")
        lines.append("-" * 40)
        lines.append(
            f"  Reads/region   mean: {sum(reads_per_region) / len(reads_per_region):>6.1f}"
            f"   median: {statistics.median(reads_per_region):>4}"
            f"   max: {max(reads_per_region):>4}"
        )
        lines.append(
            f"  K-mers/region  mean: {sum(kmers_per_region) / len(kmers_per_region):>6.1f}"
            f"   median: {statistics.median(kmers_per_region):>4}"
            f"   max: {max(kmers_per_region):>4}"
        )
        lines.append(
            f"  Region size    mean: {sum(sizes) / len(sizes):>6.0f} bp"
            f"   median: {statistics.median(sizes):>4} bp"
            f"   max: {max(sizes):>4} bp"
        )
        lines.append("")

    if regions:
        lines.append("Per-Region Results")
        lines.append("-" * 140)
        lines.append(
            f"  {'Region':<35s} {'Size':>8s} {'Reads':>6s}"
            f" {'Unique K-mers':>14s}"
            f" {'Split':>6s} {'Disc':>5s} {'MaxClip':>8s}"
            f" {'UnmapMate':>10s} {'BkptClip':>9s} {'Indel50':>8s}"
            f" {'Class':>10s} {'Type':>5s}"
        )
        lines.append(
            f"  {'------':<35s} {'----':>8s} {'-----':>6s}"
            f" {'-------------':>14s}"
            f" {'-----':>6s} {'----':>5s} {'-------':>8s}"
            f" {'---------':>10s} {'--------':>9s} {'-------':>8s}"
            f" {'-----':>10s} {'----':>5s}"
        )

        for chrom, start, end in regions:
            region_key = (chrom, start, end)
            n_reads = len(region_reads.get(region_key, set()))
            n_kmers = len(region_kmers.get(region_key, set()))
            ann = (region_annotations or {}).get(region_key, {})
            # Display as 1-based coordinates for human readability
            label = f"{chrom}:{start + 1}-{end}"
            size = end - start
            lines.append(
                f"  {label:<35s} {size:>7d}bp {n_reads:>6d}"
                f" {n_kmers:>14d}"
                f" {ann.get('split_reads', 0):>6d}"
                f" {ann.get('discordant_pairs', 0):>5d}"
                f" {ann.get('max_clip_len', 0):>8d}"
                f" {ann.get('unmapped_mates', 0):>10d}"
                f" {ann.get('breakpoint_reads', 0):>9d}"
                f" {ann.get('large_indel_reads', 0):>8d}"
                f" {ann.get('class', 'SMALL'):>10s}"
                f" {ann.get('sv_type', '.'):>5s}"
            )

    if candidate_comparison:
        n_total = len(candidate_comparison)
        n_captured = sum(1 for c in candidate_comparison if c["captured"])
        pct = (n_captured / n_total * 100) if n_total else 0.0

        lines.append("Candidate Comparison (DKA_DKT > 0.25, DKA > 10)")
        lines.append("-" * 80)
        lines.append(f"  High-quality candidates:     {n_total:>8}")
        lines.append(
            f"  Captured by discovery:       {n_captured:>8}"
            f" / {n_total} ({pct:.1f}%)"
        )
        lines.append("")
        lines.append(
            f"  {'Candidate':<30s}  {'DKA':>4s}  {'DKA_DKT':>8s}"
            f"  {'Region':>35s}"
        )
        lines.append(
            f"  {'---------':<30s}  {'---':>4s}  {'-------':>8s}"
            f"  {'------':>35s}"
        )
        for c in candidate_comparison:
            var_label = f"{c['chrom']}:{c['pos']} {c['ref']}>{c['alt']}"
            region_label = c["region"] if c["captured"] else "NOT CAPTURED"
            lines.append(
                f"  {var_label:<30s}  {c['dka']:>4d}  {c['dka_dkt']:>8.4f}"
                f"  {region_label:>35s}"
            )
        lines.append("")

    evaluation = metrics.get("dnm_evaluation") or {}
    if dnm_evaluation or evaluation.get("negatives"):
        slack = evaluation.get("slack_bp")
        source = evaluation.get("source")
        lines.append(
            f"Known Event Evaluation ({source})" if source
            else "Known Event Evaluation"
        )
        lines.append("-" * 80)
        if slack is not None:
            lines.append(
                f"  An end is supported by split, discordant or"
                f" breakpoint-clipped reads within {slack} bp of it;"
            )
            lines.append(
                f"  a region anywhere in the span (within {slack} bp) is"
                f" reported, not counted: on a genome with a"
            )
            lines.append(
                "  region every few kb it is chance, which the random loci"
                " below measure."
            )
        if dnm_evaluation:
            n_total = evaluation.get("total_loci", len(dnm_evaluation))
            n_detected = evaluation.get("detected", 0)
            pct = (n_detected / n_total * 100) if n_total else 0.0
            lines.append(f"  Known events:                  {n_total:>8}")
            lines.append(
                f"  Supported at a breakpoint:     {n_detected:>8}"
                f" / {n_total} ({pct:.1f}%)"
            )
            lines.append(
                f"  Supported at both ends:        "
                f"{evaluation.get('both_ends', 0):>8}"
            )
            lines.append(
                f"  A region anywhere in the span: "
                f"{evaluation.get('any_region', 0):>8} / {n_total}"
            )
        null = evaluation.get("null")
        if null:
            lines.append(
                f"  Random loci of the same sizes: "
                f"{null['total_loci']:>8} ({null['loci_per_event']} per event):"
                f" {100 * null['detection_rate']:.1f}% supported at a"
                f" breakpoint, {100 * null['any_region_rate']:.1f}% with a"
                f" region in the span"
            )
        negatives = evaluation.get("negatives")
        if negatives:
            n_neg = negatives["total_loci"]
            n_neg_sup = negatives["detected"]
            neg_pct = (n_neg_sup / n_neg * 100) if n_neg else 0.0
            lines.append(
                f"  Negative controls ({negatives.get('source')}):"
                f" {n_neg_sup} / {n_neg} supported at a breakpoint"
                f" ({neg_pct:.1f}%; inherited events, whose k-mers are in a"
                f" parent: support here is a false positive)"
            )
        for title, loci in (("", dnm_evaluation or []),
                            ("Negative controls", (negatives or {}).get("loci", []))):
            if not loci:
                continue
            lines.append("")
            if title:
                lines.append(f"  {title}:")
            lines.append(
                f"  {'Locus':<20s} {'Event':>14s} {'Size':>8s} {'Ends':>5s}"
                f" {'Reads':>6s} {'Split':>5s} {'Disc':>5s} {'Clip':>5s}"
                f" {'Class':>10s} {'Status':>20s}"
            )
            lines.append(
                f"  {'-----':<20s} {'-----':>14s} {'----':>8s} {'----':>5s}"
                f" {'-----':>6s} {'-----':>5s} {'----':>5s} {'----':>5s}"
                f" {'-----':>10s} {'------':>20s}"
            )
            for e in loci:
                size_str = (f"{e['event_size']}bp"
                            if e["event_size"] else "–")
                ends_str = f"{e['breakpoint_support']}/{e['ends']}"
                bp = e.get("breakpoints", {})
                split = sum(v["split_reads"] for v in bp.values())
                disc = sum(v["discordant_pairs"] for v in bp.values())
                clip = sum(v["breakpoint_reads"] for v in bp.values())
                lines.append(
                    f"  {e['locus']:<20s}"
                    f" {e['event_type'][:14]:>14s}"
                    f" {size_str:>8s}"
                    f" {ends_str:>5s}"
                    f" {e['total_reads']:>6d}"
                    f" {split:>5d}"
                    f" {disc:>5d}"
                    f" {clip:>5d}"
                    f" {e['sv_class']:>10s}"
                    f" {e['assessment']:>20s}"
                )
        lines.append("")

    if nhf_block is not None:
        if lines[-1] != "":
            lines.append("")
        lines.extend(nhf.summary_lines(nhf_block))

    lines.append("=" * 60)
    lines.append("")

    text = "\n".join(lines)

    with open(summary_path, "w") as fh:
        fh.write(text)

    return text


def _write_empty_discovery_outputs(bed_path, metrics_path, summary_path,
                                   metrics, bedpe_path=None, nhf_block=None):
    """Write empty discovery outputs for early-exit cases."""
    _write_bed([], {}, {}, bed_path)
    if bedpe_path:
        _write_bedpe([], bedpe_path)
    if nhf_block is not None:
        nhf.annotate_metrics(metrics, nhf_block)
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    _write_discovery_summary(summary_path, [], {}, {}, metrics,
                             nhf_block=nhf_block)


def run_discovery_pipeline(args):
    """Run the VCF-free de novo k-mer discovery pipeline."""
    pipeline_start = time.monotonic()

    logging.basicConfig(
        level=logging.DEBUG if args.debug_kmers else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    # ── Pre-flight checks ──────────────────────────────────────────
    for tool in ("samtools", "jellyfish"):
        if not _check_tool(tool):
            logger.error("%s not found in PATH", tool)
            sys.exit(1)

    _validate_inputs(args)

    # The non-human fraction screen (see nhf.py) is optional and never
    # fatal: a database that cannot be used is reported now, the run goes
    # on, and the screen is recorded as failed at the end; kmer-nhf can
    # add it to the finished run afterwards.
    kraken2_db = getattr(args, "kraken2_db", None)
    nhf_options = {
        "cluster_distance": getattr(args, "cluster_distance", 500),
        "sample_reads": (getattr(args, "nhf_sample_reads", None)
                         or nhf.DEFAULT_SAMPLE_READS),
        "confidence": getattr(args, "kraken2_confidence", 0.0) or 0.0,
        "memory_mapping": bool(getattr(args, "kraken2_memory_mapping", False)),
        "threads": args.threads,
    }
    nhf_preflight = None
    if kraken2_db is not None:
        nhf_preflight = nhf.kraken2_db_problems(kraken2_db) or None
        if nhf_preflight:
            logger.warning(
                "[NHF] The non-human fraction screen cannot run: %s. The run "
                "continues; the screen is recorded as failed, and kmer-nhf "
                "can add it to the finished run.", "; ".join(nhf_preflight),
            )

    def _early_nhf_block():
        """The screen's block for a run that stops before anchoring."""
        params = nhf.screen_params(**nhf_options)
        if kraken2_db is None:
            return nhf.disabled_block(params)
        if nhf_preflight:
            return nhf.failed_block("; ".join(nhf_preflight), db=kraken2_db,
                                    params=params,
                                    engine=nhf.engine_info(kraken2_db))
        return nhf.empty_run_block(kraken2_db, params,
                                   nhf.engine_info(kraken2_db))

    # Parse the known de novo events up front so a malformed file fails
    # before hours of counting rather than after.
    dnm_regions_path = getattr(args, "dnm_regions", None)
    dnm_regions = None
    if dnm_regions_path:
        try:
            dnm_regions = _load_dnm_regions(dnm_regions_path)
        except ValueError as exc:
            logger.error("Validation error: %s", exc)
            sys.exit(1)
    negative_regions_path = getattr(args, "dnm_negative_regions", None)
    negative_regions = None
    if negative_regions_path:
        try:
            negative_regions = _load_dnm_regions(negative_regions_path)
        except ValueError as exc:
            logger.error("Validation error: %s", exc)
            sys.exit(1)

    out_prefix = args.out_prefix
    bed_path = f"{out_prefix}.bed"
    info_bam_path = f"{out_prefix}.informative.bam"
    metrics_path = f"{out_prefix}.metrics.json"
    summary_path = f"{out_prefix}.summary.txt"
    bedpe_path = getattr(args, "sv_bedpe", None) or f"{out_prefix}.sv.bedpe"
    bedgraph_path = f"{out_prefix}.kmer_coverage.bedgraph"
    read_cov_bed_path = f"{out_prefix}.read_coverage.bed"
    min_bedgraph_reads = getattr(args, "min_bedgraph_reads", 3)
    min_dk_per_read = getattr(args, "min_distinct_kmers_per_read", None)
    if min_dk_per_read is None:
        min_dk_per_read = max(1, args.kmer_size // 4)
    jf_hash_size = getattr(args, "jf_hash_size", None)
    memory_limit_gb = getattr(args, "memory", None)

    # ── Configuration summary ──────────────────────────────────────
    logger.info("=" * 60)
    logger.info("  kmer-denovo  —  discovery pipeline starting")
    logger.info("=" * 60)
    logger.info(
        "  Child BAM/CRAM:    %s (%s)", args.child,
        _format_file_size(args.child),
    )
    logger.info(
        "  Mother BAM/CRAM:   %s (%s)", args.mother,
        _format_file_size(args.mother),
    )
    logger.info(
        "  Father BAM/CRAM:   %s (%s)", args.father,
        _format_file_size(args.father),
    )
    logger.info("  Reference FASTA:   %s", args.ref_fasta or "(not set)")
    logger.info("  Reference JF:      %s", getattr(args, 'ref_jf', None) or "(auto)")
    logger.info("  Output prefix:     %s", out_prefix)
    logger.info("  k-mer size:        %d", args.kmer_size)
    logger.info("  Min child count:   %d", args.min_child_count)
    logger.info("  Min distinct kmers/read: %d", min_dk_per_read)
    logger.info("  JF hash size:      %s", jf_hash_size or "(from the memory budget and the input)")
    logger.info("  Child count:       --disk (no in-memory doubling); Bloom filter %s",
                "off (--no-child-bloom-filter)" if getattr(args, "no_child_bloom_filter", False) else "on (k-mers seen once never enter the hash)")
    logger.info("  Threads:           %d", args.threads)
    logger.info(
        "  Memory limit:      %s",
        f"{memory_limit_gb:.1f} GB" if memory_limit_gb is not None
        else "(auto-detect)",
    )
    logger.info("  Tmp dir:           %s", getattr(args, 'tmp_dir', None) or "(auto)")
    logger.info(
        "  DNM regions:       %s",
        f"{dnm_regions_path} ({len(dnm_regions)} events)" if dnm_regions
        else "(none)",
    )
    logger.info(
        "  Kraken2 DB:        %s",
        kraken2_db if kraken2_db is not None
        else "(none: no non-human fraction screen)",
    )
    total_mem_gb, avail_mem_gb = _get_available_memory_gb()
    if total_mem_gb is not None:
        logger.info(
            "  System memory:     %.1f GB total, %s available",
            total_mem_gb,
            f"{avail_mem_gb:.1f} GB" if avail_mem_gb is not None
            else "(unknown)",
        )
    logger.info("=" * 60)
    _log_memory("pipeline start")

    # Resolve temp directory — avoid RAM-backed /tmp on HPC systems
    out_dir = os.path.dirname(os.path.abspath(out_prefix)) or "."
    tmp_root = _resolve_tmp_dir(args.tmp_dir, out_dir)
    logger.info("  Temp directory root: %s", tmp_root)
    if _is_tmpfs(tmp_root):
        logger.warning(
            "  ⚠ Temp directory %s appears to be on tmpfs (RAM-backed)! "
            "Large intermediate files (100+ GB for WGS) will consume RAM. "
            "Consider using --tmp-dir to point to a disk-backed filesystem.",
            tmp_root,
        )
    _log_disk_usage(tmp_root, "tmpdir filesystem")

    with tempfile.TemporaryDirectory(prefix="kmer_denovo_disc_",
                                     dir=tmp_root) as tmpdir:
        logger.info("  Working temp directory: %s", tmpdir)

        # ── Module 0: Reference K-mer Indexing ─────────────────────
        step_start = time.monotonic()
        logger.info("[Module 0] Ensuring reference Jellyfish index")
        ref_jf = _ensure_ref_jf(
            args.ref_fasta, args.kmer_size, args.threads,
            getattr(args, 'ref_jf', None),
        )
        logger.info(
            "[Module 0] Reference index ready (%s)",
            _format_elapsed(time.monotonic() - step_start),
        )
        _log_memory("after Module 0")

        # ── Module 1: Child K-merization & Reference Subtraction ───
        step_start = time.monotonic()
        logger.info("[Module 1] Child k-mer extraction & reference subtraction")
        _log_dir_size(tmpdir, "before Module 1")
        child_jf, dump_lower = _extract_child_kmers_discovery(
            args.child, args.ref_fasta, args.kmer_size,
            args.min_child_count, args.threads, tmpdir,
            jf_hash_size=jf_hash_size, memory_gb=memory_limit_gb,
            bloom_filter=not getattr(args, "no_child_bloom_filter", False),
        )
        if child_jf is None:
            n_candidates, n_non_ref, child_non_ref_fa = 0, 0, None
        else:
            child_non_ref_fa, n_candidates, n_non_ref = (
                _subtract_reference_kmers(
                    ref_jf, child_jf, dump_lower, tmpdir, args.threads,
                    min_child_count=args.min_child_count,
                )
            )

        if n_candidates == 0:
            logger.warning("No child candidate k-mers found; writing empty outputs")
            empty_metrics = {
                "mode": "discovery",
                "child_candidate_kmers": 0,
                "non_ref_kmers": 0,
                "proband_unique_kmers": 0,
                "informative_reads": 0,
                "unmapped_informative_reads": 0,
                "candidate_regions": 0,
            }
            _write_empty_discovery_outputs(
                bed_path, metrics_path, summary_path, empty_metrics,
                bedpe_path=bedpe_path, nhf_block=_early_nhf_block(),
            )
            logger.info(
                "Pipeline finished in %s",
                _format_elapsed(time.monotonic() - pipeline_start),
            )
            return

        logger.info(
            "[Module 1] Complete (%s)",
            _format_elapsed(time.monotonic() - step_start),
        )
        _log_memory("after Module 1")
        _log_dir_size(tmpdir, "after Module 1")
        _log_disk_usage(tmpdir, "tmpdir filesystem after Module 1")

        if n_non_ref == 0:
            logger.warning(
                "All child k-mers are in the reference; writing empty outputs"
            )
            empty_metrics = {
                "mode": "discovery",
                "child_candidate_kmers": n_candidates,
                "non_ref_kmers": 0,
                "proband_unique_kmers": 0,
                "informative_reads": 0,
                "unmapped_informative_reads": 0,
                "candidate_regions": 0,
            }
            _write_empty_discovery_outputs(
                bed_path, metrics_path, summary_path, empty_metrics,
                bedpe_path=bedpe_path, nhf_block=_early_nhf_block(),
            )
            logger.info(
                "Pipeline finished in %s",
                _format_elapsed(time.monotonic() - pipeline_start),
            )
            return

        # ── Module 2: Parent Filtering ─────────────────────────────
        step_start = time.monotonic()
        logger.info("[Module 2] Parent filtering")
        _log_dir_size(tmpdir, "before Module 2")
        n_proband_unique, proband_unique_fa = _filter_parents_discovery(
            args.mother, args.father, args.ref_fasta,
            child_non_ref_fa, args.kmer_size, args.threads, tmpdir,
            parent_max_count=args.parent_max_count,
        )
        logger.info(
            "[Module 2] Complete (%s)",
            _format_elapsed(time.monotonic() - step_start),
        )
        _log_memory("after Module 2")
        _log_dir_size(tmpdir, "after Module 2")
        _log_disk_usage(tmpdir, "tmpdir filesystem after Module 2")

        if n_proband_unique == 0:
            logger.warning(
                "No proband-unique k-mers after parent filtering; "
                "writing empty outputs"
            )
            empty_metrics = {
                "mode": "discovery",
                "child_candidate_kmers": n_candidates,
                "non_ref_kmers": n_non_ref,
                "proband_unique_kmers": 0,
                "informative_reads": 0,
                "unmapped_informative_reads": 0,
                "candidate_regions": 0,
            }
            _write_empty_discovery_outputs(
                bed_path, metrics_path, summary_path, empty_metrics,
                bedpe_path=bedpe_path, nhf_block=_early_nhf_block(),
            )
            logger.info(
                "Pipeline finished in %s",
                _format_elapsed(time.monotonic() - pipeline_start),
            )
            return

        # ── Module 2b: the kept Jellyfish index of proband-unique k-mers ──
        #
        # The proband-unique k-mer index is small (tens of MB for a 30x
        # genome) and is the one intermediate worth keeping: across a
        # cohort, k-mers recurring in unrelated probands are errors, and
        # that set can only be built from these files.  It is copied out
        # before the temporary directory goes.  Module 3 does not read it:
        # it loads the FASTA itself.
        step_start = time.monotonic()
        if getattr(args, "discard_proband_kmers", False):
            logger.info("[Module 2b] Proband-unique k-mer index not built "
                        "(--discard-proband-kmers)")
            logger.info(
                "[Module 2b] Complete (%s, index: not kept)",
                _format_elapsed(time.monotonic() - step_start),
            )
        else:
            logger.info(
                "[Module 2b] Building Jellyfish index of %d proband-unique "
                "k-mers", n_proband_unique,
            )
            _log_dir_size(tmpdir, "before proband index build")
            proband_jf = _build_proband_jf_index(
                proband_unique_fa, args.kmer_size, tmpdir,
                n_proband_unique=n_proband_unique,
            )
            kept_jf = f"{out_prefix}.proband_unique.jf"
            shutil.copyfile(proband_jf, kept_jf)
            os.remove(proband_jf)
            logger.info(
                "[Module 2b] Complete (%s, index: %s)",
                _format_elapsed(time.monotonic() - step_start),
                _format_file_size(kept_jf),
            )
            logger.info("  Proband-unique k-mer index kept: %s (%s)",
                        kept_jf, _format_file_size(kept_jf))
            _log_memory("after proband index build")

        # ── Module 3: Anchoring & Region Clustering ────────────────
        # The informative reads BAM is written during this same scan.
        step_start = time.monotonic()
        logger.info(
            "[Module 3] Anchoring %d proband-unique k-mers to child reads "
            "(from %s, %s)",
            n_proband_unique, proband_unique_fa,
            _format_file_size(proband_unique_fa),
        )
        _log_memory("before Module 3")
        (regions, region_reads, total_informative, region_kmers,
         unmapped_informative, read_sv_meta, kmer_coverage,
         read_coverage) = (
            _anchor_and_cluster(
                args.child, args.ref_fasta, None,
                args.kmer_size, merge_distance=args.cluster_distance,
                threads=args.threads,
                min_distinct_kmers_per_read=min_dk_per_read,
                proband_unique_fa=proband_unique_fa,
                n_proband_unique=n_proband_unique,
                tmpdir=tmpdir,
                memory_limit_gb=memory_limit_gb,
                informative_bam=info_bam_path,
            )
        )
        logger.info(
            "[Module 3] Complete (%s)",
            _format_elapsed(time.monotonic() - step_start),
        )
        _log_memory("after Module 3")

    # ── tmpdir cleaned up — all temp files removed ──────────────────
    logger.info("Temporary directory cleaned up")
    _log_memory("after tmpdir cleanup")

    # Remove the tmp_root if it was auto-created and is now empty
    try:
        if not getattr(args, "tmp_dir", None) and os.path.isdir(tmp_root):
            os.rmdir(tmp_root)  # only succeeds if empty
    except OSError:
        pass

    # ── Region filtering ───────────────────────────────────────────
    min_reads = args.min_supporting_reads
    min_kmers = args.min_distinct_kmers
    if min_reads > 1 or min_kmers > 1:
        pre_filter = len(regions)
        filtered_regions = []
        for region_key in regions:
            n_reads = len(region_reads.get(region_key, set()))
            n_kmers = len(region_kmers.get(region_key, set()))
            if n_reads >= min_reads and n_kmers >= min_kmers:
                filtered_regions.append(region_key)
            else:
                region_reads.pop(region_key, None)
                region_kmers.pop(region_key, None)
        regions = filtered_regions
        logger.info(
            "Region filtering: %d → %d regions "
            "(min-supporting-reads=%d, min-distinct-kmers=%d)",
            pre_filter, len(regions), min_reads, min_kmers,
        )

    # ── Module 4: Output ───────────────────────────────────────────
    step_start = time.monotonic()
    logger.info("[Module 4] Writing output files")

    # SV annotation and linking (from metadata — no extra BAM scan)
    logger.info("[Module 4] Annotating regions and linking breakpoints")
    region_annotations, sv_links = _annotate_and_link_from_metadata(
        regions, read_sv_meta, link_slack=args.cluster_distance,
    )
    _classify_regions(regions, region_annotations, sv_links)
    # Mapping quality beside each region: annotated, never filtered on.
    region_mapq = _region_mapq_stats(regions, region_reads, read_sv_meta)

    bed_filters = {
        "min_distinct_kmers_per_read": min_dk_per_read,
        "min_supporting_reads": min_reads,
        "min_distinct_kmers": min_kmers,
    }
    _write_bed(regions, region_reads, region_kmers, bed_path,
               region_annotations=region_annotations,
               filters=bed_filters, region_mapq=region_mapq)

    _write_bedgraph(kmer_coverage, bedgraph_path,
                    read_coverage=read_coverage,
                    min_reads=min_bedgraph_reads)

    _write_read_coverage_bed(kmer_coverage, read_coverage,
                             read_cov_bed_path,
                             min_reads=min_bedgraph_reads)

    # Free coverage data after writing — can be large at WGS scale
    logger.info(
        "  Coverage data: kmer_coverage=%d chroms, read_coverage=%d chroms",
        len(kmer_coverage), len(read_coverage),
    )
    total_positions = sum(len(v) for v in kmer_coverage.values())
    logger.info("  Total tracked positions: %d", total_positions)
    del kmer_coverage
    del read_coverage
    _log_memory("after freeing coverage data")

    _write_bedpe(sv_links, bedpe_path)

    # ── Non-human fraction screen ──────────────────────────────────
    # With the regions and junctions on disk: the informative reads of
    # every SV and AMBIGUOUS region (and of every junction end), and a
    # sample of all of them, classified with kraken2 (see nhf.py).  The
    # BED's two columns are filled in from the result; "." is unscreened.
    # kmer-nhf does the same to a finished run from its files.
    nhf_start = time.monotonic()
    region_dicts = [
        {
            "chrom": chrom, "start": start, "end": end,
            "class": region_annotations.get(
                (chrom, start, end), {}).get("class", "SMALL"),
            "sv_type": region_annotations.get(
                (chrom, start, end), {}).get("sv_type", "."),
        }
        for chrom, start, end in regions
    ]
    junction_ends = ({link["region_a"] for link in sv_links}
                     | {link["region_b"] for link in sv_links})
    nhf_tmp = None
    if kraken2_db is not None and not nhf_preflight:
        logger.info("[Module 4] Non-human fraction screen against %s",
                    kraken2_db)
        nhf_tmp = _resolve_tmp_dir(getattr(args, "tmp_dir", None), out_dir)
    nhf_block = nhf.run_screen(
        info_bam_path, region_dicts, junction_ends, kraken2_db,
        preflight_problems=nhf_preflight, tmpdir=nhf_tmp, **nhf_options,
    )
    nhf.annotate_bed(bed_path, nhf_block)
    if nhf_block["status"] == nhf.STATUS_FAILED:
        logger.warning("[NHF] %s", nhf.describe(nhf_block))
    elif nhf_block["status"] == nhf.STATUS_SCREENED:
        logger.info("[Module 4] %s (%s)", nhf.describe(nhf_block),
                    _format_elapsed(time.monotonic() - nhf_start))
    try:
        if nhf_tmp and not getattr(args, "tmp_dir", None):
            os.rmdir(nhf_tmp)  # only succeeds if empty
    except OSError:
        pass

    # ── Optional candidate comparison ──────────────────────────────
    candidate_comparison = None
    candidate_summary = getattr(args, "candidate_summary", None)
    if candidate_summary and os.path.isfile(candidate_summary):
        logger.info("[Module 4] Comparing to candidate summary: %s",
                     candidate_summary)
        hq_candidates = _parse_candidate_summary(candidate_summary)
        candidate_comparison = _compare_candidates_to_regions(
            hq_candidates, regions,
        )
        n_captured = sum(1 for c in candidate_comparison if c["captured"])
        logger.info(
            "[Module 4] High-quality candidates: %d, captured: %d",
            len(candidate_comparison), n_captured,
        )

    metrics = {
        "mode": "discovery",
        "child_candidate_kmers": n_candidates,
        "non_ref_kmers": n_non_ref,
        "proband_unique_kmers": n_proband_unique,
        "informative_reads": total_informative,
        "unmapped_informative_reads": unmapped_informative,
        "candidate_regions": len(regions),
        "filters": {
            "min_distinct_kmers_per_read": min_dk_per_read,
            "min_supporting_reads": min_reads,
            "min_distinct_kmers": min_kmers,
            "min_bedgraph_reads": min_bedgraph_reads,
        },
        "regions": [
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "size": end - start,
                "reads": len(region_reads.get((chrom, start, end), set())),
                "unique_kmers": len(
                    region_kmers.get((chrom, start, end), set())
                ),
                "split_reads": region_annotations.get(
                    (chrom, start, end), {},
                ).get("split_reads", 0),
                "discordant_pairs": region_annotations.get(
                    (chrom, start, end), {},
                ).get("discordant_pairs", 0),
                "max_clip_len": region_annotations.get(
                    (chrom, start, end), {},
                ).get("max_clip_len", 0),
                "unmapped_mates": region_annotations.get(
                    (chrom, start, end), {},
                ).get("unmapped_mates", 0),
                "breakpoint_reads": region_annotations.get(
                    (chrom, start, end), {},
                ).get("breakpoint_reads", 0),
                "large_indel_reads": region_annotations.get(
                    (chrom, start, end), {},
                ).get("large_indel_reads", 0),
                "sv_type": region_annotations.get(
                    (chrom, start, end), {},
                ).get("sv_type", "."),
                "class": region_annotations.get(
                    (chrom, start, end), {},
                ).get("class", "SMALL"),
                "max_mapq": region_mapq.get((chrom, start, end), (0, 0))[0],
                "mapq0_reads": region_mapq.get(
                    (chrom, start, end), (0, 0),
                )[1],
            }
            for chrom, start, end in regions
        ],
    }
    nhf.annotate_metrics(metrics, nhf_block)
    if candidate_comparison is not None:
        n_total = len(candidate_comparison)
        n_captured = sum(1 for c in candidate_comparison if c["captured"])
        metrics["candidate_comparison"] = {
            "hq_candidates": n_total,
            "captured": n_captured,
            "capture_rate": (n_captured / n_total) if n_total else 0.0,
            "candidates": [
                {
                    "variant": (f"{c['chrom']}:{c['pos']}"
                                f" {c['ref']}>{c['alt']}"),
                    "dka": c["dka"],
                    "dka_dkt": c["dka_dkt"],
                    "captured": c["captured"],
                    "region": c["region"],
                }
                for c in candidate_comparison
            ],
        }

    # ── Optional evaluation against known events ───────────────────
    # Breakpoint evidence within the slack of an event's ends is what
    # counts as support; a region anywhere in the span is reported but
    # is only chance on a genome with a region every few kilobases.  A
    # random-locus null says how much chance, and the negative controls
    # (inherited events, whose k-mers are in a parent) how often support
    # appears where there should be none.
    dnm_evaluation = None
    if dnm_regions:
        dnm_evaluation = _evaluate_dnm_regions(
            regions, metrics["regions"], dnm_regions,
            slack=args.cluster_distance,
        )
        headline = _summarise_dnm_evaluation(dnm_evaluation)
        logger.info(
            "[Module 4] Known events: %d / %d supported at a breakpoint "
            "(%d at both ends; %d with a region anywhere in the span)",
            headline["detected"], headline["total_loci"],
            headline["both_ends"], headline["any_region"],
        )
        metrics["dnm_evaluation"] = {
            "source": os.path.basename(dnm_regions_path),
            "slack_bp": args.cluster_distance,
            **headline,
            "loci": dnm_evaluation,
        }
        n_null = getattr(args, "dnm_null", 0) or 0
        if n_null > 0:
            null_loci = _null_dnm_regions(
                dnm_regions, _contig_lengths(args.child, args.ref_fasta),
                n_null,
            )
            null_headline = _summarise_dnm_evaluation(_evaluate_dnm_regions(
                regions, metrics["regions"], null_loci,
                slack=args.cluster_distance,
            ))
            logger.info(
                "[Module 4] Random loci of the same sizes (%d per event): "
                "%.1f%% supported at a breakpoint, %.1f%% with a region "
                "in the span",
                n_null, 100 * null_headline["detection_rate"],
                100 * null_headline["any_region_rate"],
            )
            metrics["dnm_evaluation"]["null"] = {
                "loci_per_event": n_null,
                "seed": _DNM_NULL_SEED,
                **null_headline,
            }
    if negative_regions:
        negative_evaluation = _evaluate_dnm_regions(
            regions, metrics["regions"], negative_regions,
            slack=args.cluster_distance,
        )
        negative_headline = _summarise_dnm_evaluation(negative_evaluation)
        logger.info(
            "[Module 4] Negative controls: %d / %d supported at a "
            "breakpoint (%d with a region in the span)",
            negative_headline["detected"], negative_headline["total_loci"],
            negative_headline["any_region"],
        )
        negatives = {
            "source": os.path.basename(negative_regions_path),
            **negative_headline,
            "loci": negative_evaluation,
        }
        if "dnm_evaluation" in metrics:
            metrics["dnm_evaluation"]["negatives"] = negatives
        else:
            metrics["dnm_evaluation"] = {
                "source": None,
                "slack_bp": args.cluster_distance,
                **_summarise_dnm_evaluation([]),
                "loci": [],
                "negatives": negatives,
            }

    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    logger.info("[Module 4] Metrics written to: %s", metrics_path)

    logger.info("[Module 4] Writing summary: %s", summary_path)
    _write_discovery_summary(
        summary_path, regions, region_reads, region_kmers, metrics,
        candidate_comparison=candidate_comparison,
        region_annotations=region_annotations,
        dnm_evaluation=dnm_evaluation,
        nhf_block=nhf_block,
    )

    logger.info(
        "[Module 4] Output complete (%s)",
        _format_elapsed(time.monotonic() - step_start),
    )

    # ── Optional interactive HTML report ───────────────────────────
    report_path = getattr(args, "report", None)
    if report_path:
        logger.info("[Report] Generating interactive HTML report: %s", report_path)
        from kmer_denovo_filter.report import generate_report
        generate_report(
            output_path=report_path,
            discovery_metrics_path=metrics_path,
            discovery_summary_path=summary_path,
        )

    # ── User guidance ──────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("  Discovery pipeline complete!")
    logger.info("=" * 60)
    logger.info("  Candidate regions: %s", bed_path)
    logger.info("  K-mer coverage:    %s", bedgraph_path)
    logger.info("  Read coverage:     %s", read_cov_bed_path)
    logger.info("  Informative BAM:   %s", info_bam_path)
    logger.info("  SV breakpoints:    %s", bedpe_path)
    logger.info("  Metrics:           %s", metrics_path)
    logger.info("  Summary:           %s", summary_path)
    logger.info("")
    logger.info(
        "  Next step: pass %s to a genotyper such as", bed_path,
    )
    logger.info(
        "  GATK HaplotypeCaller (--intervals) or DeepVariant for"
    )
    logger.info("  robust VCF generation.")
    logger.info("=" * 60)

    logger.info(
        "Pipeline finished successfully in %s",
        _format_elapsed(time.monotonic() - pipeline_start),
    )

