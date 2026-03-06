"""Shared utility helpers for kmer_denovo_filter.

This module contains self-contained utility functions used across the
package.  Functions are organised into sections:

* **Formatting** — human-readable elapsed time and file sizes.
* **Tool / system checks** — PATH look-ups, tmpfs detection.
* **Filesystem helpers** — temp-dir resolution, Jellyfish file discovery.
* **FASTA I/O** — reading and writing simple k-mer FASTA files.
* **Logging / monitoring** — disk, memory, and subprocess diagnostics.
* **Estimation** — heuristic sizing for Jellyfish hash tables and FASTA
  entry counts.
* **Alignment helpers** — mapping k-mer query hits to reference
  coordinates, SV-type inference.
"""

import collections
import glob
import logging
import os
import shutil

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Formatting helpers
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
# Tool / system checks
# ---------------------------------------------------------------------------


def _check_tool(name):
    """Check if an external tool is available on PATH."""
    return shutil.which(name) is not None


def _is_tmpfs(path):
    """Check if *path* resides on a tmpfs (RAM-backed) filesystem.

    Returns True on Linux when the filesystem type is tmpfs.
    Returns False on other platforms or when detection fails.
    """
    try:
        # Linux: check /proc/mounts
        real = os.path.realpath(path)
        best_mount = ""
        best_fstype = ""
        with open("/proc/mounts") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    mount_point, fstype = parts[1], parts[2]
                    if real.startswith(mount_point) and len(mount_point) > len(best_mount):
                        best_mount = mount_point
                        best_fstype = fstype
        return best_fstype == "tmpfs"
    except (FileNotFoundError, PermissionError, OSError):
        return False


# ---------------------------------------------------------------------------
# Filesystem helpers
# ---------------------------------------------------------------------------


def _resolve_tmp_dir(args, fallback_dir):
    """Resolve the temporary directory for intermediate files.

    Uses ``args.tmp_dir`` when provided.  Otherwise creates a
    subdirectory next to the output files to avoid RAM-backed ``/tmp``
    (tmpfs) on HPC systems.

    Args:
        args: Parsed CLI arguments (may have ``tmp_dir`` attribute).
        fallback_dir: Directory to use when ``--tmp-dir`` is not set
            (typically the parent directory of the output prefix or
            output file).

    Returns:
        Absolute path to the temporary directory root (created if needed).
    """
    tmp_dir = getattr(args, "tmp_dir", None)
    if tmp_dir:
        os.makedirs(tmp_dir, exist_ok=True)
        return os.path.abspath(tmp_dir)

    # Default: create a subdirectory next to the output
    tmp_root = os.path.join(fallback_dir, "kmer_denovo_tmp")
    os.makedirs(tmp_root, exist_ok=True)
    return os.path.abspath(tmp_root)


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
# FASTA I/O helpers
# ---------------------------------------------------------------------------


def _write_kmer_fasta(kmers, filepath):
    """Write k-mers to a FASTA file for jellyfish ``--if``."""
    with open(filepath, "w") as fh:
        for i, kmer in enumerate(kmers):
            fh.write(f">{i}\n{kmer}\n")


def _load_kmers_from_fasta(fasta_path):
    """Load k-mer strings from a FASTA file.

    Reads the file line-by-line, yielding only sequence lines (not
    header lines starting with ``>``).  This avoids building a large
    intermediate list.
    """
    kmers = set()
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line and not line.startswith(">"):
                kmers.add(line)
    return kmers


def _estimate_fasta_sequence_count(fasta_path, sample_lines=1000):
    """Estimate FASTA sequence-entry count from a sampled prefix.

    Returns a tuple of ``(count, extrapolated)`` where ``extrapolated`` is
    ``True`` when the count is estimated from file size and sampled bytes, and
    ``False`` when the full file was read (small files) or empty.

    Args:
        fasta_path: Path to FASTA file to sample.
        sample_lines: Number of leading lines to sample before extrapolating.
            Defaults to 1000.
    """
    if sample_lines <= 0:
        raise ValueError("sample_lines must be > 0")

    try:
        file_size = os.path.getsize(fasta_path)
    except OSError:
        return 0, False

    if file_size == 0:
        return 0, False

    sampled_bytes = 0
    sampled_entries = 0
    lines_read = 0
    hit_eof = False

    with open(fasta_path, "rb") as fh:
        while lines_read < sample_lines:
            line = fh.readline()
            if not line:
                hit_eof = True
                break
            sampled_bytes += len(line)
            lines_read += 1
            stripped = line.strip()
            if stripped and stripped.startswith(b">"):
                sampled_entries += 1

    if sampled_bytes == 0:
        return 0, False
    if sampled_entries == 0:
        return 0, False

    if hit_eof:
        return sampled_entries, False

    estimated = int(round((sampled_entries / sampled_bytes) * file_size))
    return max(estimated, 1), True


# ---------------------------------------------------------------------------
# Logging / monitoring helpers
# ---------------------------------------------------------------------------


def _log_disk_usage(path, label=""):
    """Log disk usage and available space for the filesystem containing *path*."""
    try:
        stat = os.statvfs(path)
        total_gb = (stat.f_blocks * stat.f_frsize) / (1024**3)
        avail_gb = (stat.f_bavail * stat.f_frsize) / (1024**3)
        used_gb = total_gb - avail_gb
        logger.info(
            "  [Disk] %s — %.1f GB used / %.1f GB total (%.1f GB available) — %s",
            label, used_gb, total_gb, avail_gb, path,
        )
    except OSError:
        pass


def _log_dir_size(path, label=""):
    """Log the total size of files in a directory."""
    try:
        total = 0
        for entry in os.scandir(path):
            if entry.is_file(follow_symlinks=False):
                total += entry.stat().st_size
        logger.info(
            "  [TmpDir] %s — %.2f GB in %s",
            label, total / (1024**3), path,
        )
    except OSError:
        pass


def _log_memory(label=""):
    """Log current and peak process memory usage.

    Works on Linux (``/proc/self/status``) and macOS/BSD
    (``resource.getrusage``).  Falls back silently when neither source
    is available.
    """
    try:
        info = {}
        # Linux: read from /proc
        try:
            with open("/proc/self/status") as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        info["RSS"] = int(line.split()[1]) / (1024 * 1024)
                    elif line.startswith("VmPeak:"):
                        info["Peak"] = int(line.split()[1]) / (1024 * 1024)
        except FileNotFoundError:
            pass
        # Fallback: resource module (works on macOS & Linux)
        if not info:
            import resource
            rusage = resource.getrusage(resource.RUSAGE_SELF)
            import platform
            if platform.system() == "Darwin":
                # macOS reports maxrss in bytes
                info["Peak_RSS"] = rusage.ru_maxrss / (1024**3)
            else:
                # Linux reports maxrss in KB
                info["Peak_RSS"] = rusage.ru_maxrss / (1024 * 1024)
        if info:
            parts = [f"{k}={v:.2f} GB" for k, v in sorted(info.items())]
            logger.info("  [Memory] %s — %s", label, ", ".join(parts))
    except Exception:
        pass


def _log_subprocess_memory(proc, label=""):
    """Log memory usage of a subprocess (Linux only).

    Reads ``VmRSS`` from ``/proc/{pid}/status`` for the given
    :class:`subprocess.Popen` object.  Silently skipped on non-Linux
    or when the process has already exited.
    """
    if proc is None or proc.poll() is not None:
        return
    try:
        rss_kb = 0
        with open(f"/proc/{proc.pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    rss_kb = int(line.split()[1])
                    break
        if rss_kb:
            logger.info(
                "  [SubprocessMem] %s (pid=%d) — RSS=%.2f GB",
                label, proc.pid, rss_kb / (1024 * 1024),
            )
    except Exception:
        pass


def _get_available_memory_gb():
    """Return total system memory in GB, or None if unavailable.

    On Linux reads ``/proc/meminfo`` for ``MemTotal`` and ``MemAvailable``.
    On macOS/BSD uses ``os.sysconf`` for total memory.

    Returns:
        Tuple of (total_gb, available_gb).  *available_gb* may be None
        on macOS where ``MemAvailable`` is not reported.
    """
    total_gb = None
    available_gb = None

    # Linux: parse /proc/meminfo
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    total_gb = int(line.split()[1]) / (1024 * 1024)
                elif line.startswith("MemAvailable:"):
                    available_gb = int(line.split()[1]) / (1024 * 1024)
        if total_gb is not None:
            return total_gb, available_gb
    except (FileNotFoundError, PermissionError, OSError):
        pass

    # macOS / POSIX fallback
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        if pages > 0 and page_size > 0:
            total_gb = (pages * page_size) / (1024**3)
    except (ValueError, OSError, AttributeError):
        pass

    return total_gb, available_gb


def _log_children_memory(label=""):
    """Log aggregate memory of all child processes (Linux only).

    Reads ``/proc/{pid}/status`` for each child process to report total
    RSS across all subprocesses.
    """
    try:
        my_pid = os.getpid()
        total_rss_kb = 0
        n_children = 0

        # Walk /proc for children
        proc_path = f"/proc/{my_pid}/task/{my_pid}/children"
        try:
            with open(proc_path) as f:
                child_pids = f.read().split()
        except (FileNotFoundError, PermissionError):
            # Alternative: use /proc/*/stat
            child_pids = []
            for entry in os.listdir("/proc"):
                if not entry.isdigit():
                    continue
                try:
                    with open(f"/proc/{entry}/stat") as f:
                        stat = f.read().split()
                    ppid = stat[3]
                    if int(ppid) == my_pid:
                        child_pids.append(entry)
                except (FileNotFoundError, PermissionError, IndexError):
                    continue

        for cpid in child_pids:
            try:
                with open(f"/proc/{cpid}/status") as f:
                    for line in f:
                        if line.startswith("VmRSS:"):
                            total_rss_kb += int(line.split()[1])
                            n_children += 1
                            break
            except (FileNotFoundError, PermissionError):
                continue

        if n_children > 0:
            logger.info(
                "  [ChildProcessMem] %s — %d children, total RSS=%.2f GB",
                label, n_children, total_rss_kb / (1024 * 1024),
            )
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Estimation helpers
# ---------------------------------------------------------------------------


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
# Alignment / SV helpers
# ---------------------------------------------------------------------------


def _collect_kmer_ref_positions(read, kmer_hit_indices, kmer_size):
    """Map query-level k-mer hit positions to reference coordinates.

    Args:
        read: A pysam.AlignedSegment (must be mapped).
        kmer_hit_indices: Set of query start indices where novel k-mers
            were found.
        kmer_size: Length of k-mers.

    Returns:
        Counter keyed by reference position with coverage counts.
    """
    cov = collections.Counter()
    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    query_to_ref = {qpos: rpos for qpos, rpos in aligned_pairs}
    for start_idx in kmer_hit_indices:
        for qpos in range(start_idx, start_idx + kmer_size):
            rpos = query_to_ref.get(qpos)
            if rpos is not None:
                cov[rpos] += 1
    return cov


def _infer_sv_type(region_a, region_b):
    """Infer SV type from two linked regions.

    Returns one of: INTRA (same chromosome) or BND (translocation).
    Distinguishing DEL/DUP/INV would require SA tag strand information
    which is not carried in the region data structure.
    """
    if region_a[0] != region_b[0]:
        return "BND"
    return "INTRA"
