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

Several functions have been moved to :mod:`kmer_denovo_filter.core`
sub-modules and are re-exported here for backward compatibility.
"""

import logging
import os
import shutil
import sys

from kmer_denovo_filter.core.bam_scanner import (  # noqa: F401
    _collect_kmer_ref_positions,
    _infer_sv_type,
)
from kmer_denovo_filter.core.jellyfish_wrappers import (  # noqa: F401
    _estimate_jf_hash_size,
    _find_jf_files,
)
from kmer_denovo_filter.core.memory_utils import (  # noqa: F401
    _get_available_memory_gb,
    _log_children_memory,
    _log_dir_size,
    _log_disk_usage,
    _log_memory,
    _log_subprocess_memory,
)

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


def _resolve_tmp_dir(tmp_dir, fallback_dir):
    """Resolve the temporary directory for intermediate files.

    Uses *tmp_dir* when provided.  Otherwise creates a subdirectory
    next to the output files to avoid RAM-backed ``/tmp`` (tmpfs) on
    HPC systems.

    Args:
        tmp_dir: Explicit temporary directory path, or ``None``.
            For backward compatibility this may also be an
            ``argparse.Namespace`` with a ``tmp_dir`` attribute.
        fallback_dir: Directory to use when *tmp_dir* is not set
            (typically the parent directory of the output prefix or
            output file).

    Returns:
        Absolute path to the temporary directory root (created if needed).
    """
    # Backward compatibility: accept an argparse.Namespace transparently.
    resolved = getattr(tmp_dir, "tmp_dir", tmp_dir)
    if resolved:
        os.makedirs(resolved, exist_ok=True)
        return os.path.abspath(resolved)

    # Default: create a subdirectory next to the output
    tmp_root = os.path.join(fallback_dir, "kmer_denovo_tmp")
    os.makedirs(tmp_root, exist_ok=True)
    return os.path.abspath(tmp_root)


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
# Input validation
# ---------------------------------------------------------------------------


def _validate_inputs(args):
    """Validate pipeline inputs before starting computation.

    Args:
        args: ``argparse.Namespace`` produced by :func:`~kmer_denovo_filter.cli.parse_args`.
            Expected attributes include *child*, *mother*, *father*,
            *ref_fasta*, *vcf* (``None`` in discovery mode), *kmer_size*,
            *min_baseq*, *threads*, and mode-specific options such as
            *ref_jf*, *min_child_count* (discovery) or *min_mapq* (VCF).

    Raises:
        SystemExit: with a clear error message for any invalid input.
    """
    errors = []

    # Check required input files exist
    required_files = [
        ("Child BAM/CRAM (--child)", args.child),
        ("Mother BAM/CRAM (--mother)", args.mother),
        ("Father BAM/CRAM (--father)", args.father),
    ]
    if args.vcf is not None:
        required_files.append(("Input VCF (--vcf)", args.vcf))
    for label, path in required_files:
        if not os.path.isfile(path):
            errors.append(f"{label}: file not found: {path}")

    # Check ref-fasta exists when provided
    if args.ref_fasta is not None:
        if not os.path.isfile(args.ref_fasta):
            errors.append(
                f"Reference FASTA (--ref-fasta): file not found: "
                f"{args.ref_fasta}"
            )

    # CRAM files require a reference FASTA
    for label, path in [
        ("--child", args.child),
        ("--mother", args.mother),
        ("--father", args.father),
    ]:
        if path.endswith(".cram") and args.ref_fasta is None:
            errors.append(
                f"{label} is a CRAM file but --ref-fasta was not provided"
            )

    # Check BAM/CRAM indexes exist
    for label, path in [
        ("--child", args.child),
        ("--mother", args.mother),
        ("--father", args.father),
    ]:
        if os.path.isfile(path):
            bai = path + ".bai"
            csi = path + ".csi"
            crai = path + ".crai"
            # .bam.bai or .bai (some tools drop the .bam prefix)
            alt_bai = path.rsplit(".", 1)[0] + ".bai" if "." in path else None
            if not any(
                os.path.isfile(p) for p in [bai, csi, crai]
                if p is not None
            ) and not (alt_bai and os.path.isfile(alt_bai)):
                errors.append(
                    f"{label}: no index found for {path} "
                    f"(expected .bai, .csi, or .crai)"
                )

    # Validate parameter ranges
    if args.kmer_size < 3:
        errors.append(
            f"--kmer-size must be >= 3, got {args.kmer_size}"
        )
    if args.kmer_size > 201:
        errors.append(
            f"--kmer-size must be <= 201, got {args.kmer_size}"
        )
    if args.kmer_size % 2 == 0:
        errors.append(
            f"--kmer-size should be odd for canonical k-mer symmetry, "
            f"got {args.kmer_size}"
        )
    if args.min_baseq < 0:
        errors.append(
            f"--min-baseq must be >= 0, got {args.min_baseq}"
        )
    if args.threads < 1:
        errors.append(
            f"--threads must be >= 1, got {args.threads}"
        )

    # Discovery-mode-specific validation
    if args.vcf is None:
        if args.ref_fasta is None and getattr(args, 'ref_jf', None) is None:
            errors.append(
                "Discovery mode requires --ref-fasta (or --ref-jf) "
                "to subtract reference k-mers"
            )
        ref_jf = getattr(args, 'ref_jf', None)
        if ref_jf is not None and not os.path.isfile(ref_jf):
            errors.append(
                f"Reference Jellyfish index (--ref-jf): file not found: "
                f"{ref_jf}"
            )
        min_child_count = getattr(args, 'min_child_count', 3)
        if min_child_count < 1:
            errors.append(
                f"--min-child-count must be >= 1, got {min_child_count}"
            )

    # VCF-mode-specific validation
    if args.vcf is not None:
        if args.min_mapq < 0:
            errors.append(
                f"--min-mapq must be >= 0, got {args.min_mapq}"
            )

    if errors:
        for err in errors:
            logger.error("Validation error: %s", err)
        sys.exit(1)
