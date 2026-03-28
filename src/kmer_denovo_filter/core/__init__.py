"""Shared infrastructure for kmer_denovo_filter.

This package contains the generic, domain-agnostic building blocks used
by both VCF-mode and discovery-mode pipelines.  Nothing in this package
references VCF schemas, argparse namespaces, or pipeline-specific
business logic.

Sub-modules
-----------
* :mod:`memory_utils` — process / system memory and disk diagnostics.
* :mod:`jellyfish_wrappers` — Jellyfish subprocess helpers (count,
  merge, dump, query).
* :mod:`bam_scanner` — BAM/CRAM scanning, CIGAR parsing, and the
  multi-backend k-mer hit engine.
"""

from kmer_denovo_filter.core.bam_scanner import (  # noqa: F401
    _collect_kmer_ref_positions,
    _collect_read_alignment_metadata,
    _extract_softclips,
    _infer_sv_type,
    _init_scan_worker,
    _process_informative_read,
    _scan_contig_for_hits,
)
from kmer_denovo_filter.core.jellyfish_wrappers import (  # noqa: F401
    _build_proband_jf_index,
    _ensure_ref_jf,
    _estimate_jf_hash_size,
    _find_jf_files,
    _merge_jf_files,
    _scan_parent_jellyfish,
)
from kmer_denovo_filter.core.memory_utils import (  # noqa: F401
    _get_available_memory_gb,
    _log_children_memory,
    _log_dir_size,
    _log_disk_usage,
    _log_memory,
    _log_subprocess_memory,
)
