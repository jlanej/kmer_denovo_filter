"""Discovery-mode pipeline for VCF-free de novo k-mer analysis."""

import bisect
import collections
import concurrent.futures
import json
import logging
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

import pysam

from kmer_denovo_filter.core.bam_scanner import (
    _clip_side,
    _collect_read_alignment_metadata,
    _init_scan_worker,
    _parse_cigar,
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
from kmer_denovo_filter.core.memory_utils import (
    _get_available_memory_gb,
    _log_children_memory,
    _log_dir_size,
    _log_disk_usage,
    _log_memory,
    _log_subprocess_memory,
)
from kmer_denovo_filter.kmer_utils import (
    canonicalize,
    estimate_automaton_memory_gb,
    read_supports_alt,
)
from kmer_denovo_filter.utils import (
    _check_tool,
    _collect_kmer_ref_positions,
    _estimate_fasta_sequence_count,
    _estimate_jf_hash_size,
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


# ── Discovery-mode helpers ────────────────────────────────────────


def _extract_child_kmers_discovery(child_bam, ref_fasta, kmer_size,
                                   min_child_count, threads, tmpdir,
                                   jf_hash_size=None):
    """Module 1: Extract all child k-mers and filter by minimum count.

    Child k-mers are counted with ``jellyfish count -C`` (canonical mode)
    to match the reference index built by :func:`_ensure_ref_jf`.  This
    ensures that ``jellyfish query`` during reference subtraction compares
    k-mers in the same canonical orientation.

    The hash size for Jellyfish is estimated dynamically from the BAM
    file size (unless overridden via *jf_hash_size*) to avoid hash
    overflow that can produce multi-hundred-GB monolithic index files.

    When Jellyfish produces multiple chunk files (hash overflow), they
    are merged in a streaming fashion before dumping.

    The ``jellyfish dump`` output is streamed line-by-line so the full
    dump text is never held in memory.  Each Jellyfish chunk file is
    removed as soon as it has been processed to free disk space and OS
    page cache.

    Returns:
        child_candidates_fa: path to FASTA of candidate child k-mers
            (count >= min_child_count).
        n_candidates: number of candidate k-mers.
    """
    child_jf = os.path.join(tmpdir, "child.jf")

    # Step 1: Count all child k-mers
    if jf_hash_size is None:
        jf_hash_size = _estimate_jf_hash_size(child_bam, kmer_size, default="1G")
    logger.info(
        "Extracting child k-mers from BAM (k=%d, jf hash size=%s)…",
        kmer_size, jf_hash_size,
    )
    samtools_threads = max(1, threads // 4)
    samtools_cmd = [
        "samtools", "fasta", "-F", "0xD00",
        "-@", str(samtools_threads),
        child_bam,
    ]
    if ref_fasta:
        samtools_cmd.extend(["--reference", ref_fasta])

    jellyfish_cmd = [
        "jellyfish", "count",
        "-m", str(kmer_size),
        "-s", jf_hash_size,
        "-t", str(threads),
        "-C",
        "-o", child_jf,
        "/dev/fd/0",
    ]

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
            "chunks: %d)",
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
        child_jf_final = _merge_jf_files(jf_files, merged_jf, threads)
    elif jf_files:
        child_jf_final = jf_files[0]
    else:
        logger.warning("No Jellyfish output files found")
        child_candidates_fa = os.path.join(tmpdir, "child_candidates.fa")
        with open(child_candidates_fa, "w"):
            pass
        return child_candidates_fa, 0

    # Step 2: Dump k-mers with count >= min_child_count (streamed to
    # avoid holding the full dump text in memory).
    logger.info(
        "Dumping child k-mers with count >= %d from %s (%s)…",
        min_child_count, child_jf_final,
        _format_file_size(child_jf_final),
    )
    dump_start = time.monotonic()
    dump_cmd = [
        "jellyfish", "dump", "-c",
        "-L", str(min_child_count),
        child_jf_final,
    ]
    child_candidates_fa = os.path.join(tmpdir, "child_candidates.fa")
    n_candidates = 0
    with tempfile.TemporaryFile(mode="w+") as stderr_f:
        p_dump = subprocess.Popen(
            dump_cmd, stdout=subprocess.PIPE, stderr=stderr_f, text=True,
        )
        # Monitor dump subprocess memory during streaming
        last_dump_log = dump_start
        with open(child_candidates_fa, "w") as fh:
            for line in p_dump.stdout:
                line = line.rstrip("\n")
                if line:
                    kmer = line.split()[0]
                    fh.write(f">{n_candidates}\n{kmer}\n")
                    n_candidates += 1
                    if n_candidates % 5_000_000 == 0:
                        now = time.monotonic()
                        if now - last_dump_log >= 30:
                            logger.info(
                                "  … dumped %d k-mers so far (%s, "
                                "FASTA: %s)",
                                n_candidates,
                                _format_elapsed(now - dump_start),
                                _format_file_size(child_candidates_fa),
                            )
                            _log_subprocess_memory(p_dump, "jellyfish-dump")
                            _log_memory("during child dump")
                            _log_disk_usage(tmpdir, "tmpdir during dump")
                            last_dump_log = now
        p_dump.wait()
        if p_dump.returncode != 0:
            stderr_f.seek(0)
            raise RuntimeError(
                f"jellyfish dump (child) failed: {stderr_f.read()}"
            )

    logger.info(
        "Child k-mer dump complete (%s, %d candidates, FASTA: %s)",
        _format_elapsed(time.monotonic() - dump_start),
        n_candidates, _format_file_size(child_candidates_fa),
    )

    # Remove the child jellyfish index immediately – it is no longer
    # needed and can be very large (100+ GB for WGS).
    for f in _find_jf_files(child_jf):
        if os.path.exists(f):
            os.remove(f)
    if child_jf_final != child_jf and os.path.exists(child_jf_final):
        os.remove(child_jf_final)
    logger.info("Removed child Jellyfish index to free disk/cache")
    _log_memory("after child index removal")

    logger.info(
        "Child candidate k-mers (count >= %d): %d",
        min_child_count, n_candidates,
    )
    return child_candidates_fa, n_candidates


def _subtract_reference_kmers(ref_jf, child_candidates_fa, tmpdir):
    """Subtract reference genome k-mers from child candidates.

    Both the reference index (*ref_jf*) and the child candidates FASTA are
    produced with ``jellyfish -C`` (canonical mode), so ``jellyfish query``
    correctly matches k-mers regardless of strand orientation.

    The query output is streamed line-by-line to avoid holding the full
    result text in memory.  The *child_candidates_fa* file is removed
    after the query completes because it is no longer needed.

    Returns:
        child_non_ref_fa: path to FASTA of non-reference child k-mers.
        n_non_ref: number of surviving k-mers.
    """
    query_cmd = [
        "jellyfish", "query", ref_jf, "-s", child_candidates_fa,
    ]

    child_non_ref_fa = os.path.join(tmpdir, "child_non_ref_kmers.fa")
    n_non_ref = 0
    with tempfile.TemporaryFile(mode="w+") as stderr_f:
        p_query = subprocess.Popen(
            query_cmd, stdout=subprocess.PIPE, stderr=stderr_f, text=True,
        )
        with open(child_non_ref_fa, "w") as fh:
            for line in p_query.stdout:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[1] == "0":
                    fh.write(f">{n_non_ref}\n{parts[0]}\n")
                    n_non_ref += 1
        p_query.wait()
        if p_query.returncode != 0:
            stderr_f.seek(0)
            raise RuntimeError(
                f"jellyfish query (ref subtraction) failed: {stderr_f.read()}"
            )

    # Remove the candidates FASTA – the non-ref subset is now on disk.
    if os.path.exists(child_candidates_fa):
        os.remove(child_candidates_fa)

    logger.info(
        "Non-reference child k-mers after subtraction: %d", n_non_ref,
    )
    return child_non_ref_fa, n_non_ref


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


def _anchor_and_cluster(child_bam, ref_fasta, proband_unique_kmers,
                        kmer_size, merge_distance=500, threads=1,
                        min_distinct_kmers_per_read=1,
                        proband_unique_fa=None,
                        proband_jf=None,
                        n_proband_unique=None,
                        tmpdir=None,
                        memory_limit_gb=None,
                        informative_bam=None):
    """Module 3: Find reads containing proband-unique k-mers and cluster regions.

    Scans **all** primary, non-duplicate child reads — including unmapped
    and low-MAPQ reads — so that no proband-unique k-mers are missed at
    this stage.  The scan runs one task per contig plus one for unplaced
    unmapped reads, in a process pool when more than one worker is
    available and in-process otherwise.  Results are combined in task
    order, so the output does not depend on which worker finishes first.

    Supports two scanning backends:

    - **Jellyfish-backed** (when *proband_jf* is provided) — each worker
      opens a ``jellyfish query`` subprocess that memory-maps the hash
      file.  The OS page cache is shared across workers, so N workers
      ≈ 1× the hash file in memory.  This is the path for WGS discovery
      with hundreds of millions of proband-unique k-mers.
    - **Aho-Corasick** (when *proband_unique_fa* or *proband_unique_kmers*
      is provided) — builds an in-memory automaton.  Fast for small k-mer
      sets (VCF mode or unit tests).

    Args:
        child_bam: Path to child BAM file.
        ref_fasta: Path to reference FASTA (or None).
        proband_unique_kmers: Set of canonical k-mer strings, or None
            when *proband_unique_fa* or *proband_jf* is provided.
        kmer_size: K-mer size.
        merge_distance: Maximum gap for merging adjacent regions.
        threads: Number of parallel workers (default 1).
        min_distinct_kmers_per_read: Minimum distinct proband-unique
            k-mers a read must carry to be retained (default 1).
        proband_unique_fa: Optional FASTA path.  Workers load k-mers
            from disk to build an Aho-Corasick automaton.
        proband_jf: Optional Jellyfish index path (.jf).  Workers query
            the index via ``jellyfish query`` — low memory, suited for
            large k-mer sets.
        n_proband_unique: Number of proband-unique k-mers (for logging).
        tmpdir: Writable directory for temporary files (the per-contig
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

    Returns:
        Tuple of (regions, region_reads, total_informative, region_kmers,
        unmapped_informative, read_sv_meta, kmer_coverage, read_coverage).
    """
    anchor_start = time.monotonic()

    # Determine scanning mode
    use_jellyfish = proband_jf is not None

    # ── Count k-mers if not provided ────────────────────────────────
    if n_proband_unique is None and proband_unique_fa:
        n_proband_unique = 0
        with open(proband_unique_fa) as fh:
            for line in fh:
                if line and not line.startswith(">"):
                    n_proband_unique += 1

    n_kmer_count = n_proband_unique or (
        len(proband_unique_kmers) if proband_unique_kmers else 0
    )

    # ── Log memory planning ─────────────────────────────────────────
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
    if use_jellyfish:
        # Jellyfish-backed: workers share OS page cache, memory ≈ 1×
        # the .jf file size regardless of worker count.
        jf_size_gb = 0.0
        try:
            jf_size_gb = os.path.getsize(proband_jf) / (1024**3)
        except OSError:
            pass
        logger.info(
            "  [MemoryPlanning] jellyfish-backed mode: %d k-mers, "
            "index: %.1f GB (shared via page cache)",
            n_kmer_count, jf_size_gb,
        )
    else:
        est_per_worker_gb = estimate_automaton_memory_gb(n_kmer_count)
        logger.info(
            "  [MemoryPlanning] Aho-Corasick mode: %d k-mers, "
            "estimated %.1f GB per worker",
            n_kmer_count, est_per_worker_gb,
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
    bam.close()

    # One task per contig + one for unmapped reads.  With an informative
    # BAM requested, each task writes its reads to its own shard, and the
    # shards are combined in task (coordinate) order afterwards.
    shard_dir = None
    shard_paths = [None] * (len(contigs) + 1)
    if informative_bam:
        shard_dir = tempfile.mkdtemp(prefix="informative_reads_", dir=tmpdir)
        shard_paths = [
            os.path.join(shard_dir, f"{i:06d}.bam")
            for i in range(len(contigs) + 1)
        ]
    tasks = [
        (child_bam, ref_fasta, contig, shard)
        for contig, shard in zip(contigs + [None], shard_paths)
    ]

    # ── Dynamically cap workers based on available memory ──────────
    if use_jellyfish:
        # Jellyfish workers share page cache; no per-worker penalty.
        # All workers share the same memory-mapped .jf file, so
        # worker count is bounded only by threads and task count.
        n_workers = min(threads, len(tasks))
    else:
        max_workers_by_mem = threads
        if avail_mem_gb is not None and est_per_worker_gb > 0:
            usable_gb = avail_mem_gb * 0.8
            max_workers_by_mem = max(1, int(usable_gb / est_per_worker_gb))
        elif total_mem_gb is not None and est_per_worker_gb > 0:
            usable_gb = total_mem_gb * 0.7
            max_workers_by_mem = max(1, int(usable_gb / est_per_worker_gb))
        n_workers = min(threads, len(tasks), max_workers_by_mem)
    n_workers = max(n_workers, 1)

    logger.info(
        "  Anchoring: %d contigs, %d workers (requested=%d, mode=%s)",
        len(contigs), n_workers, threads,
        "jellyfish" if use_jellyfish else "aho-corasick",
    )

    # Worker init: choose data source based on mode
    if use_jellyfish:
        init_args = (proband_jf, kmer_size,
                     min_distinct_kmers_per_read)
        init_mode = "jellyfish-query"
    elif proband_unique_fa:
        init_args = (proband_unique_fa, kmer_size,
                     min_distinct_kmers_per_read)
        init_mode = "fasta-file"
    else:
        init_args = (proband_unique_kmers or set(), kmer_size,
                     min_distinct_kmers_per_read)
        init_mode = "in-memory-set"

    logger.info(
        "  Worker init mode: %s", init_mode,
    )

    read_hits = []
    reads_seen = set()
    read_sv_meta = {}
    kmer_coverage = collections.defaultdict(collections.Counter)
    read_coverage = collections.defaultdict(collections.Counter)
    unmapped_informative = 0
    total_reads_scanned = 0
    completed_contigs = 0
    informative_received = 0  # before cross-contig dedup; for progress logs
    # Hits and SV metadata of finished tasks that are waiting for an
    # earlier task to finish, keyed by task index (see _merge_result).
    held = {}
    next_to_fold = 0
    last_progress_log = time.monotonic()
    progress_interval = 300  # seconds (5 minutes)

    def _merge_result(index, result):
        """Fold one task's scan results into the running totals.

        Coverage and counts are added straight away.  Hits and SV metadata
        are de-duplicated across contigs by (read name, is_supplementary),
        where the first contig to report a key wins — e.g. when both mates
        of a pair are informative on different chromosomes.  Folding them
        strictly in task order keeps that choice independent of which
        worker finishes first, and identical to a single-worker run.
        """
        nonlocal unmapped_informative, total_reads_scanned
        nonlocal completed_contigs, informative_received, next_to_fold
        (hits, seen, unmapped, scanned,
         sv_meta, worker_cov, worker_read_cov) = result
        total_reads_scanned += scanned
        unmapped_informative += unmapped
        completed_contigs += 1
        informative_received += len(hits) + unmapped
        # Merge k-mer coverage (update is faster than += for Counters)
        for chrom, cov in worker_cov.items():
            kmer_coverage[chrom].update(cov)
        # Merge read coverage
        for chrom, cov in worker_read_cov.items():
            read_coverage[chrom].update(cov)

        held[index] = (hits, seen, sv_meta)
        while next_to_fold in held:
            task_hits, task_seen, task_meta = held.pop(next_to_fold)
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
            # Track all seen keys for cross-contig dedup
            reads_seen.update(task_seen)

        # Log individual contig completion for large contigs
        if scanned >= 1_000_000:
            logger.info(
                "  [Anchoring] Contig %s complete: %d reads "
                "scanned, %d informative hits (%s elapsed)",
                tasks[index][2] or "(unmapped)", scanned,
                len(hits) + unmapped,
                _format_elapsed(time.monotonic() - anchor_start),
            )

    def _log_progress():
        """Log progress every 5 minutes or every 100 contigs."""
        nonlocal last_progress_log
        now = time.monotonic()
        at_milestone = (
            completed_contigs % 100 == 0
            or completed_contigs == len(tasks)
        )
        if now - last_progress_log >= progress_interval or at_milestone:
            logger.info(
                "  [Anchoring] Progress: %d/%d contigs complete "
                "(%.0f%%), %d reads scanned, %d informative (%s)",
                completed_contigs, len(tasks),
                100 * completed_contigs / len(tasks),
                total_reads_scanned,
                informative_received,
                _format_elapsed(now - anchor_start),
            )
            _log_memory("during anchoring")
            _log_children_memory("during anchoring")
            last_progress_log = now

    try:
        if n_workers > 1:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=n_workers,
                initializer=_init_scan_worker,
                initargs=init_args,
            ) as executor:
                futures = {
                    executor.submit(_scan_contig_for_hits, *t): i
                    for i, t in enumerate(tasks)
                }

                # Use wait() with timeout for time-based progress
                # reporting so users get feedback even when large
                # contigs take hours.
                pending = set(futures.keys())
                while pending:
                    done, pending = concurrent.futures.wait(
                        pending, timeout=progress_interval,
                        return_when=concurrent.futures.FIRST_COMPLETED,
                    )

                    if not done:
                        # Timeout — no contig completed; log heartbeat
                        logger.info(
                            "  [Anchoring] Heartbeat: %d/%d contigs "
                            "complete, %d reads scanned, %d informative "
                            "(%s elapsed)",
                            completed_contigs, len(tasks),
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
                        index = futures[future]
                        try:
                            result = future.result()
                        except Exception:
                            logger.error(
                                "Worker failed for contig=%s",
                                tasks[index][2],
                            )
                            raise
                        _merge_result(index, result)
                    _log_progress()
        else:
            # A single worker: run the same per-contig scan in-process.
            _init_scan_worker(*init_args)
            try:
                for index, task in enumerate(tasks):
                    _merge_result(index, _scan_contig_for_hits(*task))
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
    """Combine per-contig informative-read shards into a sorted, indexed BAM.

    Shards are read in task order (contigs in header order, then unplaced
    unmapped reads).  Each ``(query_name, is_supplementary)`` key — the key
    used to de-duplicate reads while anchoring — is written once, so when
    both mates of a pair are informative on different contigs only the
    first is kept.

    Args:
        shard_paths: Shard BAM paths in task order.  Missing files (tasks
            without informative reads) are skipped.
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


def _write_bed(regions, region_reads, region_kmers, bed_path,
               region_annotations=None, filters=None):
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
    """
    with open(bed_path, "w") as fh:
        if filters:
            parts = " ".join(f"{k}={v}" for k, v in sorted(filters.items()))
            fh.write(f"#filters: {parts}\n")
        fh.write(
            "#chrom\tstart\tend\treads\tunique_kmers"
            "\tsplit_reads\tdiscordant_pairs"
            "\tmax_clip_len\tunmapped_mates\tclass"
            "\tbreakpoint_reads\tlarge_indel_reads\tsv_type\n"
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
            fh.write(
                f"{chrom}\t{start}\t{end}\t{n_reads}\t{n_kmers}"
                f"\t{split_reads}\t{discordant_pairs}"
                f"\t{max_clip_len}\t{unmapped_mates}\t{region_class}"
                f"\t{breakpoint_reads}\t{large_indel_reads}\t{sv_type}\n"
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


#: Minimum mapping quality for an alignment to link two regions: the
#: informative read, its supplementary alignment (SA tag) and, when the
#: MQ tag records it, its mate.
_MIN_LINK_MAPQ = 20

#: Soft clips starting within this many bp of one another mark one
#: breakpoint.
_SV_CLIP_TOLERANCE = 5

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
    """
    clips = sorted(clips)
    in_window = collections.Counter()
    best = lo = 0
    for pos, qname in clips:
        in_window[qname] += 1
        while pos - clips[lo][0] > tolerance:
            dropped = clips[lo][1]
            in_window[dropped] -= 1
            if not in_window[dropped]:
                del in_window[dropped]
            lo += 1
        best = max(best, len(in_window))
    return best


def _annotate_and_link_from_metadata(regions, region_reads, read_sv_meta,
                                     link_slack=0):
    """Annotate regions and link breakpoints using pre-collected metadata.

    Uses per-read SV metadata collected during the anchoring scan
    (Module 3), so no additional BAM I/O is needed.

    Evidence is counted per molecule (read name), at most once per count
    per region.  Pair-level evidence (split alignment, unmapped mate,
    discordant pair) counts in every region where the molecule has an
    informative alignment.  Soft clips, CIGAR indels of at least 50 bp
    and ``max_clip_len`` count in the region of the alignment that shows
    them; ``breakpoint_reads`` is the most molecules clipped at one
    breakpoint, or 0 when that is fewer than two.  An unmapped
    informative read counts as an unmapped mate in the region where it
    is placed.

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
    (INTRA, or BND across chromosomes).  A region's ``sv_type`` is the
    majority over its molecules, including CIGAR indels (DEL, INS) and
    joins within the region, or ``.`` when there is none.

    Args:
        regions: List of (chrom, start, end) tuples.
        region_reads: Dict mapping region tuple to set of read names.
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

    # Build lookup from read name to regions it belongs to
    read_to_regions = {}
    for region_key in regions:
        for qname in region_reads.get(region_key, set()):
            read_to_regions.setdefault(qname, set()).add(region_key)

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

    clips_by_region = collections.defaultdict(list)
    type_votes = collections.defaultdict(dict)  # region -> {qname: SV type}
    for (qname, _is_supp), meta in read_sv_meta.items():
        placed = meta.get("placed_unmapped")
        if placed is not None:
            region = region_at(*placed)
            if region is not None:
                count(qname, region, "unmapped_mates")
            continue

        evidence = []
        if meta["has_sa"]:
            evidence.append("split_reads")
        if meta["is_paired"]:
            if meta["mate_is_unmapped"]:
                evidence.append("unmapped_mates")
            elif not meta["is_proper_pair"]:
                evidence.append("discordant_pairs")
        for region in read_to_regions.get(qname, ()):
            for name in evidence:
                count(qname, region, name)

        own = region_at(*meta["pos"]) if "pos" in meta else None
        if own is None:
            continue
        ann = annotations[own]
        ann["max_clip_len"] = max(ann["max_clip_len"], meta["max_clip"])
        for pos in meta["clip_positions"]:
            clips_by_region[own].append((pos, qname))
        if meta["large_indel"]:
            count(qname, own, "large_indel_reads")
            type_votes[own].setdefault(qname, meta["large_indel"])

    # One long clip is often an adapter or a low-quality tail; a
    # breakpoint needs at least two molecules clipped at it.
    for region, clips in clips_by_region.items():
        n_clipped = _largest_clip_cluster(clips)
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

    def join(end_a, region_a, end_b, region_b, qname, pair=False):
        if region_a is None or region_b is None:
            return
        # By position; on a tie, a "+" end comes first
        if (end_b[:2], end_b[2] != "+") < (end_a[:2], end_a[2] != "+"):
            end_a, end_b = end_b, end_a
            region_a, region_b = region_b, region_a
        strands = (end_a[2], end_b[2]) if end_a[2] and end_b[2] else None
        sv_type = _infer_sv_type(region_a, region_b, strands)
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
            join(own_end, own,
                 breakend(parts[0], sa_start, sa_end, _clip_side(sa_cigar)),
                 region_at(parts[0], sa_start, link_slack), qname)

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
        annotations[region]["sv_type"] = (
            _unique_majority(type_votes[region].values()) or "."
        )

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

    - ``SV``: at least two molecules show one kind of evidence in
      ``_SV_EVIDENCE`` (split alignments, discordant pairs, unmapped
      mates, soft clips at one breakpoint, CIGAR indels of 50 bp or
      more), or the region is linked to another region
    - ``SMALL``: none of that evidence and not linked
    - ``AMBIGUOUS``: otherwise (evidence from a single molecule)

    Updates region_annotations in place with a ``class`` key.
    """
    linked_regions = set()
    for link in sv_links:
        linked_regions.add(link["region_a"])
        linked_regions.add(link["region_b"])

    for region_key in regions:
        ann = region_annotations.get(region_key, {})
        strongest = max(ann.get(name, 0) for name in _SV_EVIDENCE)
        if strongest >= 2 or region_key in linked_regions:
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


def _evaluate_dnm_regions(discovery_regions, region_detail, dnm_regions,
                          slack=0):
    """Evaluate how well VCF-free discovery captures known DNM regions.

    For each known de novo event (from ``--dnm-regions``), determines
    whether it was nominated by the discovery pipeline and collects
    quantitative k-mer and SV-signal evidence from the discovery
    region(s) overlapping it or within *slack* bp of it.  The slack
    matters for deletions: their breakpoint regions flank the deleted
    interval rather than overlap it, and curated breakpoints can be off
    by a few bases.

    The evaluation provides a simple genotype-like assessment per region:

    - **DETECTED**: ≥1 discovery region matches the curated locus with
      informative reads carrying proband-unique k-mers.
    - **NOT_DETECTED**: No matching discovery region found.

    For detected regions, a *k-mer signal score* summarises evidence
    strength as ``unique_kmers / region_size_bp``.  Higher density
    indicates more child-specific sequence variation concentrated in the
    region, consistent with a de novo SV.

    Args:
        discovery_regions: List of (chrom, start, end) tuples (0-based,
            half-open) from the discovery BED.
        region_detail: List of dicts with per-region metrics (the
            ``regions`` array from the discovery metrics JSON).
        dnm_regions: List of (chrom, pos, size_or_None, event_type)
            tuples, as returned by :func:`_load_dnm_regions`.  *pos* is
            1-based; each event covers *size* bp from *pos* (1 bp when
            the size is unknown).
        slack: How far (bp) outside an event a region may lie and still
            count toward it (the pipeline uses ``--cluster-distance``).

    Returns:
        List of dicts, one per known event, with keys:

        - locus (str): ``chrom:pos``
        - event_type (str): Event type label.
        - event_size (int or None): Expected event size in bp.
        - detected (bool): Whether ≥1 discovery region overlaps.
        - discovery_regions (list of str): Matching region labels.
        - total_reads (int): Sum of informative reads across matches.
        - total_unique_kmers (int): Sum of proband-unique k-mers.
        - max_clip_len (int): Max soft-clip length across matches.
        - unmapped_mates (int): Sum of unmapped-mate counts.
        - discordant_pairs (int): Sum of discordant-pair counts.
        - split_reads (int): Sum of split-read counts.
        - sv_class (str): Most severe SV class across matches.
        - kmer_signal (float): ``total_unique_kmers / span_bp``.
        - assessment (str): ``DETECTED`` or ``NOT_DETECTED``.
    """
    # Build an index from region tuple → detail dict
    detail_by_key = {}
    for rd in region_detail:
        key = (rd["chrom"], rd["start"], rd["end"])
        detail_by_key[key] = rd

    results = []
    for chrom, pos, size, event_type in dnm_regions:
        # 1-based pos → 0-based half-open, like the discovery regions
        dnm_start = pos - 1
        dnm_end = dnm_start + (size if size else 1)  # point if no size

        # Find discovery regions overlapping the event or within slack
        matches = []
        for dr_key in discovery_regions:
            dr_chrom, dr_start, dr_end = dr_key
            if dr_chrom != chrom:
                continue
            # Both 0-based half-open; a gap shorter than slack counts
            if dr_start < dnm_end + slack and dnm_start - slack < dr_end:
                matches.append(dr_key)

        detected = len(matches) > 0

        # Aggregate evidence across matching regions
        total_reads = 0
        total_kmers = 0
        max_clip = 0
        total_unmapped = 0
        total_discordant = 0
        total_split = 0
        region_labels = []
        sv_classes = []
        span_start = dnm_start
        span_end = dnm_end

        for m_key in matches:
            rd = detail_by_key.get(m_key, {})
            total_reads += rd.get("reads", 0)
            total_kmers += rd.get("unique_kmers", 0)
            clip = rd.get("max_clip_len", 0)
            if clip > max_clip:
                max_clip = clip
            total_unmapped += rd.get("unmapped_mates", 0)
            total_discordant += rd.get("discordant_pairs", 0)
            total_split += rd.get("split_reads", 0)
            sv_classes.append(rd.get("class", "SMALL"))
            label = f"{m_key[0]}:{m_key[1] + 1}-{m_key[2]}"
            region_labels.append(label)
            # Expand span to cover all matching regions
            if m_key[1] < span_start:
                span_start = m_key[1]
            if m_key[2] > span_end:
                span_end = m_key[2]

        span_bp = max(span_end - span_start, 1)
        kmer_signal = total_kmers / span_bp if detected else 0.0

        # Most severe SV class
        class_priority = {"SV": 3, "AMBIGUOUS": 2, "SMALL": 1}
        sv_class = max(sv_classes, key=lambda c: class_priority.get(c, 0)) \
            if sv_classes else "NONE"

        assessment = "DETECTED" if detected else "NOT_DETECTED"

        results.append({
            "locus": f"{chrom}:{pos}",
            "event_type": event_type,
            "event_size": size,
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


def _write_discovery_summary(summary_path, regions, region_reads,
                             region_kmers, metrics,
                             candidate_comparison=None,
                             region_annotations=None,
                             dnm_evaluation=None):
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

    if dnm_evaluation:
        n_total = len(dnm_evaluation)
        n_detected = sum(1 for e in dnm_evaluation if e["detected"])
        pct = (n_detected / n_total * 100) if n_total else 0.0
        source = metrics.get("dnm_evaluation", {}).get("source")
        lines.append(
            f"Curated DNM Region Evaluation ({source})" if source
            else "Curated DNM Region Evaluation"
        )
        lines.append("-" * 80)
        lines.append(f"  Curated DNM loci:            {n_total:>8}")
        lines.append(
            f"  Detected by discovery:       {n_detected:>8}"
            f" / {n_total} ({pct:.1f}%)"
        )
        slack = metrics.get("dnm_evaluation", {}).get("slack_bp")
        if slack is not None:
            lines.append(
                f"  Regions counted within:      {slack:>8} bp of an event"
            )
        lines.append("")
        lines.append(
            f"  {'Locus':<20s} {'Event':>25s} {'Size':>8s}"
            f" {'Reads':>6s} {'Kmers':>6s} {'Signal':>7s}"
            f" {'MaxClip':>8s} {'Class':>10s} {'Status':>14s}"
        )
        lines.append(
            f"  {'-----':<20s} {'-----':>25s} {'----':>8s}"
            f" {'-----':>6s} {'-----':>6s} {'------':>7s}"
            f" {'-------':>8s} {'-----':>10s} {'------':>14s}"
        )
        for e in dnm_evaluation:
            size_str = (f"{e['event_size']}bp"
                        if e["event_size"] else "–")
            lines.append(
                f"  {e['locus']:<20s}"
                f" {e['event_type']:>25s}"
                f" {size_str:>8s}"
                f" {e['total_reads']:>6d}"
                f" {e['total_unique_kmers']:>6d}"
                f" {e['kmer_signal']:>7.4f}"
                f" {e['max_clip_len']:>8d}"
                f" {e['sv_class']:>10s}"
                f" {e['assessment']:>14s}"
            )
        lines.append("")

    lines.append("=" * 60)
    lines.append("")

    text = "\n".join(lines)

    with open(summary_path, "w") as fh:
        fh.write(text)

    return text


def _write_empty_discovery_outputs(bed_path, metrics_path, summary_path,
                                   metrics, bedpe_path=None):
    """Write empty discovery outputs for early-exit cases."""
    _write_bed([], {}, {}, bed_path)
    if bedpe_path:
        _write_bedpe([], bedpe_path)
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    _write_discovery_summary(summary_path, [], {}, {}, metrics)


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
    logger.info("  JF hash size:      %s", jf_hash_size or "(auto)")
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
        child_candidates_fa, n_candidates = _extract_child_kmers_discovery(
            args.child, args.ref_fasta, args.kmer_size,
            args.min_child_count, args.threads, tmpdir,
            jf_hash_size=jf_hash_size,
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
                bedpe_path=bedpe_path,
            )
            logger.info(
                "Pipeline finished in %s",
                _format_elapsed(time.monotonic() - pipeline_start),
            )
            return

        child_non_ref_fa, n_non_ref = _subtract_reference_kmers(
            ref_jf, child_candidates_fa, tmpdir,
        )
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
                bedpe_path=bedpe_path,
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
                bedpe_path=bedpe_path,
            )
            logger.info(
                "Pipeline finished in %s",
                _format_elapsed(time.monotonic() - pipeline_start),
            )
            return

        # ── Module 2b: Build Jellyfish index of proband-unique k-mers ──
        #
        # Build a .jf hash index from the proband-unique FASTA so that
        # Module 3 workers can query k-mer membership via disk-backed
        # jellyfish query instead of loading them into Python memory.
        # The .jf is memory-mapped and shared across workers via the
        # OS page cache (N workers ≈ 1× hash file in RAM).
        step_start = time.monotonic()
        logger.info(
            "[Module 2b] Building Jellyfish index of %d proband-unique k-mers",
            n_proband_unique,
        )
        _log_dir_size(tmpdir, "before proband index build")
        proband_jf = _build_proband_jf_index(
            proband_unique_fa, args.kmer_size, tmpdir,
            n_proband_unique=n_proband_unique,
        )
        logger.info(
            "[Module 2b] Complete (%s, index: %s)",
            _format_elapsed(time.monotonic() - step_start),
            _format_file_size(proband_jf),
        )
        _log_memory("after proband index build")
        _log_dir_size(tmpdir, "after proband index build")

        # ── Module 3: Anchoring & Region Clustering ────────────────
        # The informative reads BAM is written during this same scan.
        step_start = time.monotonic()
        logger.info(
            "[Module 3] Anchoring %d proband-unique k-mers to child reads "
            "(jellyfish index: %s)",
            n_proband_unique,
            _format_file_size(proband_jf),
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
                proband_jf=proband_jf,
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
        regions, region_reads, read_sv_meta,
        link_slack=args.cluster_distance,
    )
    _classify_regions(regions, region_annotations, sv_links)

    bed_filters = {
        "min_distinct_kmers_per_read": min_dk_per_read,
        "min_supporting_reads": min_reads,
        "min_distinct_kmers": min_kmers,
    }
    _write_bed(regions, region_reads, region_kmers, bed_path,
               region_annotations=region_annotations,
               filters=bed_filters)

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
            }
            for chrom, start, end in regions
        ],
    }
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

    # ── Optional evaluation against known de novo events ──────────
    dnm_evaluation = None
    if dnm_regions:
        dnm_evaluation = _evaluate_dnm_regions(
            regions, metrics["regions"], dnm_regions,
            slack=args.cluster_distance,
        )
        n_dnm_detected = sum(1 for e in dnm_evaluation if e["detected"])
        logger.info(
            "[Module 4] Curated DNM evaluation: %d / %d detected",
            n_dnm_detected, len(dnm_evaluation),
        )
        metrics["dnm_evaluation"] = {
            "source": os.path.basename(dnm_regions_path),
            "slack_bp": args.cluster_distance,
            "total_loci": len(dnm_evaluation),
            "detected": n_dnm_detected,
            "detection_rate": n_dnm_detected / len(dnm_evaluation),
            "loci": dnm_evaluation,
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

