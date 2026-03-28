"""Memory, disk, and subprocess monitoring utilities.

Functions in this module log process-level and system-level resource
usage.  They are safe to call from any platform — unsupported metrics
are silently skipped.
"""

import logging
import os

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Disk monitoring
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


# ---------------------------------------------------------------------------
# Process / system memory
# ---------------------------------------------------------------------------


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
