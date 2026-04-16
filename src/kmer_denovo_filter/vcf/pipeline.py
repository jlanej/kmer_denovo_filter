"""VCF-mode pipeline for de novo variant k-mer analysis."""

import collections
import json
import logging
import os
import statistics
import subprocess
import sys
import tempfile
import time

import pysam

from kmer_denovo_filter.core.bam_scanner import (
    _collect_read_alignment_metadata,
    _extract_softclips,
    _init_scan_worker,
    _scan_contig_for_hits,
)
from kmer_denovo_filter.core.jellyfish_wrappers import (
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
    JellyfishKmerQuery,
    Kraken2Runner,
    _HUMAN_TAXID,
    _extract_read_kmers,
    _is_symbolic,
    build_kmer_automaton,
    canonicalize,
    estimate_automaton_memory_gb,
    extract_variant_spanning_kmers,
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
    _load_kmers_from_fasta,
    _resolve_tmp_dir,
    _validate_inputs,
    _write_kmer_fasta,
)

logger = logging.getLogger(__name__)
_FRACTION_PRECISION = 4


def _run_kraken2_on_reads(
    child_bam, ref_fasta, read_names, kraken2_db,
    confidence=0.0, threads=1, tmpdir=None,
    informative_reads_by_variant=None, memory_mapping=False,
):
    """Classify child reads with kraken2 and return a result summary.

    Extracts sequences for *read_names* from *child_bam*, classifies them
    with kraken2, and returns a :class:`Kraken2Runner.Result` with tallied
    bacterial / archaeal / fungal / protist / viral / univec_core /
    non-human / human / root counts.

    Args:
        child_bam: Path to child BAM/CRAM.
        ref_fasta: Path to reference FASTA (may be None for BAM).
        read_names: Set of read names to classify.
        kraken2_db: Path to the kraken2 database directory.
        confidence: Kraken2 confidence threshold.
        threads: Number of threads for kraken2.
        tmpdir: Optional directory for temporary files.
        informative_reads_by_variant: Optional dict mapping internal
            variant keys (``chrom:pos`` with 0-based ``pos`` as produced
            by this pipeline) to informative read-name sets. When
            provided, only those loci are fetched from the BAM/CRAM to
            avoid a whole-file scan.
        memory_mapping: Whether to pass ``--memory-mapping`` to Kraken2
            to reduce RAM usage by memory-mapping DB files.

    Returns:
        A :class:`Kraken2Runner.Result`.
    """
    if not read_names:
        return Kraken2Runner.Result()

    # Collect sequences from BAM.
    # Prefer targeted locus fetches when variant→read mappings are
    # available; otherwise fall back to a whole-file scan.
    sequences = {}
    bam = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )
    used_targeted_fetch = False
    if informative_reads_by_variant:
        loci_to_names = {}
        for var_key, names in informative_reads_by_variant.items():
            if not names:
                continue
            parts = var_key.split(":")
            if len(parts) < 2:
                logger.warning(
                    "[Kraken2] Skipping malformed variant key (missing ':'): %s",
                    var_key,
                )
                continue
            chrom = parts[0]
            pos_str = parts[1]
            try:
                pos = int(pos_str)
            except ValueError:
                logger.warning(
                    "[Kraken2] Skipping malformed variant key (non-integer pos): %s",
                    var_key,
                )
                continue
            target_names = set(names).intersection(read_names)
            if not target_names:
                continue
            loci_to_names.setdefault((chrom, pos), set()).update(target_names)

        if loci_to_names:
            used_targeted_fetch = True
            for (chrom, pos), target_names in sorted(loci_to_names.items()):
                for read in bam.fetch(chrom, pos, pos + 1):
                    if (
                        read.query_name in target_names
                        and read.query_sequence
                        and read.query_name not in sequences
                    ):
                        sequences[read.query_name] = read.query_sequence

    if not used_targeted_fetch:
        for read in bam.fetch(until_eof=True):
            if read.query_name in read_names and read.query_sequence:
                if read.query_name not in sequences:
                    sequences[read.query_name] = read.query_sequence
    bam.close()

    if not sequences:
        return Kraken2Runner.Result()

    kr = Kraken2Runner(
        kraken2_db,
        confidence=confidence,
        threads=threads,
        memory_mapping=memory_mapping,
    )
    return kr.classify_sequences(sequences, tmpdir=tmpdir)


def _parse_kmer_votes(kmer_string, name_map=None, top_n=10):
    """Parse a Kraken2 kmer_detail_string into vote summaries.

    Args:
        kmer_string: Raw Kraken2 kmer detail string (space-separated
            ``taxid:count`` tokens; ``|:|`` separates paired-end mates).
        name_map: Optional ``{taxid: name}`` dict for named output.
        top_n: Maximum number of entries to include (sorted descending
            by count).

    Returns:
        Tuple ``(kmer_votes, kmer_votes_named, total_kmers,
        human_kmer_count)`` where ``kmer_votes`` is
        ``taxid1:count1;taxid2:count2;...`` and ``kmer_votes_named``
        replaces taxids with names.  Taxid ``0`` is rendered as
        ``unclassified`` in the named column.  ``A`` (ambiguous) tokens
        are excluded.
    """
    if not kmer_string:
        return ("", "", 0, 0)

    counts = {}  # taxid (int) -> total count
    for token in kmer_string.replace("|:|", " ").split():
        taxid_str, _, count_str = token.partition(":")
        if not taxid_str or not count_str:
            continue
        try:
            tid = int(taxid_str)
            cnt = int(count_str)
        except ValueError:
            continue
        counts[tid] = counts.get(tid, 0) + cnt

    total_kmers = sum(counts.values())
    human_kmer_count = counts.get(_HUMAN_TAXID, 0)

    # Sort descending by count, truncate
    sorted_votes = sorted(counts.items(), key=lambda x: (-x[1], x[0]))
    top_votes = sorted_votes[:top_n]

    kmer_votes = ";".join(f"{tid}:{cnt}" for tid, cnt in top_votes)

    def _name_for(tid):
        if tid == 0:
            return "unclassified"
        if name_map and tid in name_map:
            return name_map[tid]
        return str(tid)

    kmer_votes_named = ";".join(
        f"{_name_for(tid)}:{cnt}" for tid, cnt in top_votes
    )

    return (kmer_votes, kmer_votes_named, total_kmers, human_kmer_count)


def _write_kraken2_read_detail_bed(
    output_path,
    informative_reads_by_variant,
    informative_alt_reads_by_variant,
    kraken2_result,
    name_map,
):
    """Write per-read Kraken2 classification detail as a BED file.

    Produces a bgzipped, tabix-indexed BED file with one row per
    (variant, read) pair.  The first three columns are BED-standard
    ``chrom``, ``chromStart``, ``chromEnd``; the remaining columns
    capture read-level classification detail.

    The file includes a ``#``-prefixed header line with column names so
    that downstream tools can parse the schema.

    Args:
        output_path: Destination path (should end with ``.bed.gz``).
            A ``.tbi`` index is written alongside.
        informative_reads_by_variant: ``{variant_key: set_of_read_names}``
            for all DKU reads.
        informative_alt_reads_by_variant: ``{variant_key: set_of_read_names}``
            for DKA reads (subset of DKU).
        kraken2_result: :class:`Kraken2Runner.Result` with
            ``per_read_detail`` populated.
        name_map: ``{taxid: name}`` dict from
            :meth:`Kraken2Runner._load_name_map`, or ``None``.
    """
    _BED_COLUMNS = [
        "#chrom", "chromStart", "chromEnd", "variant", "read_name",
        "read_set", "kraken2_status", "assigned_taxid", "assigned_taxon",
        "domain", "guard_status", "is_nonhuman", "kmer_votes",
        "kmer_votes_named", "total_kmers", "human_kmer_count",
    ]

    # Collect and sort rows: (chrom, pos_int, variant_key, read_name)
    row_keys = []
    for var_key in informative_reads_by_variant:
        parts = var_key.split(":")
        if len(parts) < 4:
            continue
        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        ref = parts[2]
        for rname in informative_reads_by_variant[var_key]:
            row_keys.append((chrom, pos, ref, var_key, rname))

    # Sort by chrom (lexicographic), then pos (numeric), then read_name
    row_keys.sort(key=lambda x: (x[0], x[1], x[4]))

    # Write uncompressed BED to a temp file, then bgzip
    raw_path = output_path.replace(".bed.gz", ".bed")
    if raw_path == output_path:
        raw_path = output_path + ".tmp"

    with open(raw_path, "w") as fh:
        fh.write("\t".join(_BED_COLUMNS) + "\n")

        for chrom, pos, ref, var_key, rname in row_keys:
            detail = kraken2_result.per_read_detail.get(rname)
            if detail is None:
                continue

            dka_names = informative_alt_reads_by_variant.get(var_key, set())
            read_set = "DKA" if rname in dka_names else "DKU"

            taxid = detail["taxid"]
            status = detail["status"]
            domain = detail["domain"]
            guard_status = detail["guard_status"]
            is_nonhuman = detail["is_nonhuman"]
            kmer_string = detail["kmer_string"]

            # Taxon name
            if status == "U" or taxid == 0:
                assigned_taxon = "."
            elif name_map and taxid in name_map:
                assigned_taxon = name_map[taxid]
            else:
                assigned_taxon = str(taxid)

            # K-mer vote summaries
            kmer_votes, kmer_votes_named, total_kmers, human_kmer_count = (
                _parse_kmer_votes(kmer_string, name_map)
            )

            # BED coordinates: 0-based start, exclusive end
            chrom_start = pos
            chrom_end = pos + len(ref)

            fields = [
                chrom,
                str(chrom_start),
                str(chrom_end),
                var_key,
                rname,
                read_set,
                status,
                str(taxid),
                assigned_taxon,
                domain,
                guard_status,
                "true" if is_nonhuman else "false",
                kmer_votes,
                kmer_votes_named,
                str(total_kmers),
                str(human_kmer_count),
            ]
            fh.write("\t".join(fields) + "\n")

    # bgzip and tabix index
    pysam.tabix_compress(raw_path, output_path, force=True)
    try:
        os.unlink(raw_path)
    except OSError:
        pass
    pysam.tabix_index(
        output_path, force=True, preset="bed",
        meta_char="#",
    )


def _build_span_bed_rows(
    alignment_meta,
    informative_reads_by_variant,
    informative_alt_reads_by_variant,
    kraken2_result,
    name_map,
):
    """Build per-alignment row data for span BED output.

    Returns a list of tuples
    ``(chrom, start, read_name, is_supplementary, rec_dict, annotation_dict)``
    sorted by ``(chrom, start, read_name)``.

    ``rec_dict`` contains the alignment record fields (``start``, ``end``,
    ``softclip_left``, ``softclip_right``, etc.).  ``annotation_dict``
    contains read-level classification fields (``taxon_name``, ``domain``,
    ``guard_status``, ``is_nonhuman``, ``variant_str``, ``read_set``,
    ``is_split``).
    """
    # Build read_name → set of variant keys mapping
    read_to_variants = {}
    for var_key, names in informative_reads_by_variant.items():
        for rname in names:
            read_to_variants.setdefault(rname, set()).add(var_key)

    # Build read_name → read_set (DKA if in any DKA set, otherwise DKU)
    dka_reads = set()
    for names in informative_alt_reads_by_variant.values():
        dka_reads.update(names)

    rows = []
    for rname, records in alignment_meta.items():
        detail = kraken2_result.per_read_detail.get(rname)
        if detail is None:
            continue

        var_keys = read_to_variants.get(rname, set())
        if not var_keys:
            continue

        variant_str = ",".join(sorted(var_keys))
        read_set = "DKA" if rname in dka_reads else "DKU"

        taxid = detail["taxid"]
        status = detail["status"]
        domain = detail["domain"]
        guard_status = detail["guard_status"]
        is_nonhuman = detail["is_nonhuman"]

        # Taxon name
        if status == "U" or taxid == 0:
            taxon_name = "Unclassified"
        elif name_map and taxid in name_map:
            taxon_name = name_map[taxid]
        else:
            taxon_name = f"Unknown_taxid_{taxid}"

        is_split = any(r["has_sa"] for r in records)

        annotation = {
            "taxon_name": taxon_name,
            "domain": domain,
            "guard_status": guard_status,
            "is_nonhuman": is_nonhuman,
            "variant_str": variant_str,
            "read_set": read_set,
            "is_split": is_split,
            "rname": rname,
        }

        for rec in records:
            rows.append((
                rec["chrom"], rec["start"], rname,
                rec["is_supplementary"], rec, annotation,
            ))

    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    return rows


_SPAN_BED_COLUMNS = [
    "#chrom", "start", "end", "taxon_name", "domain",
    "guard_status", "is_nonhuman", "read_name", "variant",
    "read_set", "mapq", "softclip_left", "softclip_right",
    "is_split", "is_supplementary",
]

_EXPANDED_SPAN_BED_COLUMNS = _SPAN_BED_COLUMNS + [
    "aligned_start", "aligned_end",
]


def _format_span_row(rec, ann):
    """Format a single standard span BED row from record and annotation."""
    return [
        rec["chrom"],
        str(rec["start"]),
        str(rec["end"]),
        ann["taxon_name"],
        ann["domain"],
        ann["guard_status"],
        "true" if ann["is_nonhuman"] else "false",
        ann["rname"],
        ann["variant_str"],
        ann["read_set"],
        str(rec["mapq"]),
        str(rec["softclip_left"]),
        str(rec["softclip_right"]),
        "true" if ann["is_split"] else "false",
        "true" if rec["is_supplementary"] else "false",
    ]


def _format_expanded_span_row(rec, ann):
    """Format a single expanded span BED row from record and annotation.

    The expanded BED extends coordinates by the observed soft-clip
    lengths::

        expanded_start = max(0, reference_start - softclip_left)
        expanded_end   = reference_end + softclip_right

    Two extra columns (``aligned_start``, ``aligned_end``) preserve the
    original mapped coordinates.
    """
    expanded_start = max(0, rec["start"] - rec["softclip_left"])
    expanded_end = rec["end"] + rec["softclip_right"]
    fields = [
        rec["chrom"],
        str(expanded_start),
        str(expanded_end),
        ann["taxon_name"],
        ann["domain"],
        ann["guard_status"],
        "true" if ann["is_nonhuman"] else "false",
        ann["rname"],
        ann["variant_str"],
        ann["read_set"],
        str(rec["mapq"]),
        str(rec["softclip_left"]),
        str(rec["softclip_right"]),
        "true" if ann["is_split"] else "false",
        "true" if rec["is_supplementary"] else "false",
        str(rec["start"]),
        str(rec["end"]),
    ]
    return fields


def _write_bed_from_rows(output_path, columns, rows, format_fn):
    """Write a bgzipped, tabix-indexed BED file from pre-built rows.

    All rows are formatted first, then sorted by the actual BED
    coordinates (``chrom``, ``start``) to guarantee the output is
    position-sorted regardless of the coordinate transform applied by
    *format_fn*.  This is essential for the expanded span BED where
    soft-clip extension can reorder positions relative to aligned starts.

    Args:
        output_path: Destination ``.bed.gz`` path.
        columns: Header column names (first element should start with ``#``).
        rows: Row tuples from :func:`_build_span_bed_rows` (need not
            be pre-sorted; sorting by output coordinates is handled
            internally).
        format_fn: Callable ``(rec, ann) -> list[str]`` producing tab fields.
    """
    raw_path = output_path.replace(".bed.gz", ".bed")
    if raw_path == output_path:
        raw_path = output_path + ".tmp"

    # Format all rows, then sort by actual BED coordinates so that
    # tabix indexing never encounters out-of-order positions.
    formatted = [format_fn(rec, ann) for _, _, _, _, rec, ann in rows]
    formatted.sort(key=lambda f: (f[0], int(f[1])))

    with open(raw_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for fields in formatted:
            fh.write("\t".join(fields) + "\n")

    pysam.tabix_compress(raw_path, output_path, force=True)
    try:
        os.unlink(raw_path)
    except OSError:
        pass
    pysam.tabix_index(
        output_path, force=True, preset="bed",
        meta_char="#",
    )


def _write_kraken2_span_bed(
    output_path,
    alignment_meta,
    informative_reads_by_variant,
    informative_alt_reads_by_variant,
    kraken2_result,
    name_map,
):
    """Write species-annotated genomic span BED for classified reads.

    Produces a bgzipped, tabix-indexed BED file with one row per
    (variant, read, alignment record) combination.  The BED coordinates
    are the read's **aligned reference span** (``reference_start`` to
    ``reference_end``), and the species label applies to the **entire
    read** (Kraken2 classifies the full read, not sub-regions).

    For split reads (SA tag), both alignment segments are emitted as
    separate BED intervals with the same classification, linked by
    ``read_name``.

    Args:
        output_path: Destination path (should end with ``.bed.gz``).
            A ``.tbi`` index is written alongside.
        alignment_meta: Dict from :func:`_collect_read_alignment_metadata`
            mapping ``read_name`` → list of alignment record dicts.
        informative_reads_by_variant: ``{variant_key: set_of_read_names}``
            for all DKU reads.
        informative_alt_reads_by_variant: ``{variant_key: set_of_read_names}``
            for DKA reads (subset of DKU).
        kraken2_result: :class:`Kraken2Runner.Result` with
            ``per_read_detail`` populated.
        name_map: ``{taxid: name}`` dict from
            :meth:`Kraken2Runner._load_name_map`, or ``None``.
    """
    rows = _build_span_bed_rows(
        alignment_meta, informative_reads_by_variant,
        informative_alt_reads_by_variant, kraken2_result, name_map,
    )
    _write_bed_from_rows(output_path, _SPAN_BED_COLUMNS, rows,
                         _format_span_row)


def _write_kraken2_expanded_span_bed(
    output_path,
    alignment_meta,
    informative_reads_by_variant,
    informative_alt_reads_by_variant,
    kraken2_result,
    name_map,
):
    """Write soft-clip–expanded genomic span BED for classified reads.

    Like :func:`_write_kraken2_span_bed`, but the BED coordinates are
    naively extended by the observed soft-clip lengths::

        expanded_start = max(0, reference_start - softclip_left)
        expanded_end   = reference_end + softclip_right

    Two extra columns (``aligned_start``, ``aligned_end``) preserve the
    original mapped coordinates for reference.

    **The expanded spans are for visualization only** — they do not
    represent verified reference alignments.  Clusters of expanded
    non-human intervals at a locus are suggestive of contamination,
    integration breakpoints, or library-prep chimeras.

    Args:
        output_path: Destination path (should end with ``.bed.gz``).
            A ``.tbi`` index is written alongside.
        alignment_meta: See :func:`_write_kraken2_span_bed`.
        informative_reads_by_variant: See :func:`_write_kraken2_span_bed`.
        informative_alt_reads_by_variant: See :func:`_write_kraken2_span_bed`.
        kraken2_result: See :func:`_write_kraken2_span_bed`.
        name_map: See :func:`_write_kraken2_span_bed`.
    """
    rows = _build_span_bed_rows(
        alignment_meta, informative_reads_by_variant,
        informative_alt_reads_by_variant, kraken2_result, name_map,
    )
    _write_bed_from_rows(output_path, _EXPANDED_SPAN_BED_COLUMNS, rows,
                         _format_expanded_span_row)


def _collect_child_kmers(
    child_bam, ref_fasta, variants, kmer_size, min_baseq, min_mapq,
    debug_kmers, kmer_fasta,
    flush_threshold=500_000,
):
    """Extract child k-mers spanning each variant position.

    K-mers are written in batches to *kmer_fasta* (FASTA format) to
    avoid holding millions of strings in memory at once.  Each batch is
    partially deduplicated in-memory before being flushed.

    Returns:
        total_child_kmers: approximate number of unique child k-mers written
        variant_read_kmers: dict mapping variant key to list of
            (read_name, kmer_set, supports_alt) tuples
    """
    bam = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )
    batch = set()
    total_written = 0
    total_reads_scanned = 0
    fasta_fh = open(kmer_fasta, "w")
    variant_read_kmers = {}
    n_variants = len(variants)
    log_interval = max(1, n_variants // 10)
    extract_start = time.monotonic()

    def _flush_batch():
        nonlocal total_written
        for kmer in batch:
            fasta_fh.write(f">{total_written}\n{kmer}\n")
            total_written += 1
        batch.clear()

    for var_idx, var in enumerate(variants, 1):
        chrom = var["chrom"]
        pos = var["pos"]  # 0-based
        ref = var["ref"]
        alts = var["alts"]
        alt = alts[0] if alts else None
        if alts and len(alts) > 1:
            logger.warning(
                "Multiallelic variant %s:%d has %d ALT alleles; "
                "only the first ALT (%s) will be evaluated",
                chrom, pos + 1, len(alts), alt,
            )
        alt_str = alt if alt is not None else "."
        var_key = f"{chrom}:{pos}:{ref}:{alt_str}"
        if alt is not None and _is_symbolic(alt):
            logger.debug(
                "Skipping variant %s:%d with symbolic allele %s",
                chrom, pos, alt,
            )
            variant_read_kmers[var_key] = []
            continue
        read_kmers = []

        for read in bam.fetch(chrom, pos, pos + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.is_duplicate:
                continue
            if not (read.reference_start <= pos < read.reference_end):
                continue

            total_reads_scanned += 1
            aligned_pairs = read.get_aligned_pairs(matches_only=False)
            seq = read.query_sequence
            quals = read.query_qualities
            kmers = extract_variant_spanning_kmers(
                read, pos, kmer_size, min_baseq, ref=ref, alt=alt,
                aligned_pairs=aligned_pairs, seq=seq, quals=quals,
            )
            if kmers:
                supports = read_supports_alt(
                    read, pos, ref, alt, min_baseq=min_baseq,
                    aligned_pairs=aligned_pairs, seq=seq, quals=quals,
                )
                read_kmers.append((read.query_name, kmers, supports))
                batch.update(kmers)
                if len(batch) >= flush_threshold:
                    _flush_batch()

        variant_read_kmers[var_key] = read_kmers

        if debug_kmers:
            unique = (
                set().union(*(k for _, k, _ in read_kmers)) if read_kmers
                else set()
            )
            logger.info(
                "Variant %s: %d reads, %d unique k-mers",
                var_key, len(read_kmers), len(unique),
            )

        if var_idx % log_interval == 0 or var_idx == n_variants:
            elapsed = time.monotonic() - extract_start
            logger.info(
                "[Step 2/5]   Processed %d / %d variants (%.0f%%) — "
                "%d reads scanned, %d k-mers collected (%s)",
                var_idx, n_variants, 100 * var_idx / n_variants,
                total_reads_scanned, total_written + len(batch),
                _format_elapsed(elapsed),
            )

    # Flush remaining k-mers
    if batch:
        _flush_batch()

    fasta_fh.close()
    bam.close()
    return total_written, variant_read_kmers



def _parse_vcf_variants(vcf_path):
    """Parse VCF file and return a list of variant dicts.

    Each dict contains chrom, pos (0-based), ref, alts, and id.
    """
    vcf = pysam.VariantFile(vcf_path)
    variants = []
    for rec in vcf:
        variants.append({
            "chrom": rec.chrom,
            "pos": rec.start,  # 0-based
            "ref": rec.ref,
            "alts": rec.alts,
            "id": rec.id,
        })
    vcf.close()
    return variants


def _write_annotated_vcf(input_vcf, output_vcf, annotations, proband_id=None):
    """Write annotated VCF with de novo k-mer metrics.

    The output is always bgzipped and tabix-indexed.  If *output_vcf*
    does not already end with ``.gz``, ``.gz`` is appended.

    When *proband_id* matches a sample in the VCF, DKU and related fields are written
    as FORMAT fields on that sample.  Otherwise they are written as INFO
    fields.

    When non-human fraction annotations are present in *annotations*,
    DKU_BF/DKA_BF, DKU_AF/DKA_AF, DKU_FF/DKA_FF, DKU_PF/DKA_PF,
    DKU_VF/DKA_VF, DKU_UCF/DKA_UCF, DKU_NHF/DKA_NHF, DKU_UF/DKA_UF,
    and DKU_HLF/DKA_HLF are also added.

    Returns:
        The actual output path (with ``.gz`` suffix).
    """
    vcf_in = pysam.VariantFile(input_vcf)
    has_kraken_fractions = any(
        "dku_bacterial_fraction" in ann or "dku_nonhuman_fraction" in ann
        for ann in annotations.values()
    )

    use_format = proband_id is not None and proband_id in list(vcf_in.header.samples)

    if use_format:
        logger.info(
            "Proband '%s' found in VCF samples; annotating as FORMAT fields",
            proband_id,
        )
    elif proband_id is not None:
        logger.warning(
            "Proband '%s' not found in VCF samples (%s); "
            "falling back to INFO annotation",
            proband_id, list(vcf_in.header.samples),
        )

    category = "FORMAT" if use_format else "INFO"

    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKU"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Number of child fragments (unique read names) with at least one "
             "variant-spanning k-mer unique to child (absent from both parents)"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKT"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Total child fragments (unique read names) with variant-spanning k-mers"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKA"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Number of child fragments (unique read names) with at least one "
             "unique k-mer that also exactly supports the candidate allele"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKU_DKT"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description",
             "Proportion of child fragments with unique k-mers (DKU/DKT)"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKA_DKT"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description",
             "Proportion of child fragments with unique allele-supporting "
             "k-mers (DKA/DKT)"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "MAX_PKC"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Maximum k-mer count in parents for variant-spanning k-mers"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "AVG_PKC"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description",
             "Average k-mer count in parents for variant-spanning k-mers found in parents"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "MIN_PKC"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Minimum k-mer count in parents for variant-spanning k-mers"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "MAX_PKC_ALT"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Maximum k-mer count in parents for alt-allele-supporting k-mers"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "AVG_PKC_ALT"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description",
             "Average k-mer count in parents for alt-allele-supporting k-mers found in parents"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "MIN_PKC_ALT"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Minimum k-mer count in parents for alt-allele-supporting k-mers"),
        ],
    )
    if has_kraken_fractions:
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_BF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments classified as bacterial by "
                 "kraken2; denominator equals DKU (both are fragment-based)"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_BF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments classified as bacterial by "
                 "kraken2; DKA fragments are always a subset of DKU"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_AF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments classified as archaeal by "
                 "kraken2; denominator equals DKU (both are fragment-based)"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_AF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments classified as archaeal by "
                 "kraken2; DKA fragments are always a subset of DKU"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_FF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments classified as fungal by "
                 "kraken2; denominator equals DKU (both are fragment-based)"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_FF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments classified as fungal by "
                 "kraken2; DKA fragments are always a subset of DKU"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_PF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments classified as protist by "
                 "kraken2; denominator equals DKU (both are fragment-based)"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_PF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments classified as protist by "
                 "kraken2; DKA fragments are always a subset of DKU"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_VF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments classified as viral by "
                 "kraken2; denominator equals DKU (both are fragment-based). "
                 "Reads with any human k-mer evidence are excluded, which "
                 "conservatively handles viruses that integrate into human "
                 "DNA (e.g. endogenous retroviruses, HBV, HPV)"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_VF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments classified as viral by "
                 "kraken2; DKA fragments are always a subset of DKU"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_UCF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments classified as UniVec Core "
                 "(synthetic sequencing-vector/adapter sequences, taxid "
                 "81077) by kraken2; denominator equals DKU (both are "
                 "fragment-based). Reads with any human k-mer evidence "
                 "are excluded. UniVec Core reads are NOT included in "
                 "the non-human fraction (DKU_NHF)"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_UCF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments classified as UniVec Core "
                 "by kraken2; DKA fragments are always a subset of DKU"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_NHF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments classified as non-human by "
                 "kraken2; denominator equals DKU (both are fragment-based). "
                 "UniVec Core reads are excluded (see DKU_UCF)"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_NHF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments classified as non-human by "
                 "kraken2; DKA fragments are always a subset of DKU. "
                 "UniVec Core reads are excluded (see DKA_UCF)"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_UF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments that were unclassified by "
                 "kraken2 (no taxonomic assignment). Denominator equals "
                 "DKU (both are fragment-based). Together DKU_NHF + "
                 "DKU_UCF + DKU_HLF + DKU_UF = 1.0"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_UF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments that were unclassified by "
                 "kraken2; DKA fragments are always a subset of DKU. "
                 "Together DKA_NHF + DKA_UCF + DKA_HLF + DKA_UF = 1.0"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_HLF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKU fragments in the human lineage: "
                 "classified reads that are neither definitively "
                 "non-human (DKU_NHF) nor UniVec Core (DKU_UCF). "
                 "Includes reads directly classified as human, reads "
                 "cleared by the human homology guard (HHG), and reads "
                 "assigned to broad taxonomic ranks on the human-to-root "
                 "path (e.g. Eukaryota, Root). Together DKU_NHF + "
                 "DKU_UCF + DKU_HLF + DKU_UF = 1.0"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_HLF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of DKA fragments in the human lineage; "
                 "DKA fragments are always a subset of DKU. "
                 "Together DKA_NHF + DKA_UCF + DKA_HLF + DKA_UF = 1.0"),
            ],
        )

    if not output_vcf.endswith(".gz"):
        output_vcf = output_vcf + ".gz"

    vcf_out = pysam.VariantFile(output_vcf, "wz", header=vcf_in.header)

    for rec in vcf_in:
        alt_str = rec.alts[0] if rec.alts else "."
        var_key = f"{rec.chrom}:{rec.start}:{rec.ref}:{alt_str}"
        if var_key in annotations:
            ann = annotations[var_key]
            if use_format:
                rec.samples[proband_id]["DKU"] = ann["dku"]
                rec.samples[proband_id]["DKT"] = ann["dkt"]
                rec.samples[proband_id]["DKA"] = ann["dka"]
                rec.samples[proband_id]["DKU_DKT"] = ann["dku_dkt"]
                rec.samples[proband_id]["DKA_DKT"] = ann["dka_dkt"]
                rec.samples[proband_id]["MAX_PKC"] = ann["max_pkc"]
                rec.samples[proband_id]["AVG_PKC"] = ann["avg_pkc"]
                rec.samples[proband_id]["MIN_PKC"] = ann["min_pkc"]
                rec.samples[proband_id]["MAX_PKC_ALT"] = ann["max_pkc_alt"]
                rec.samples[proband_id]["AVG_PKC_ALT"] = ann["avg_pkc_alt"]
                rec.samples[proband_id]["MIN_PKC_ALT"] = ann["min_pkc_alt"]
                if has_kraken_fractions:
                    rec.samples[proband_id]["DKU_BF"] = ann.get(
                        "dku_bacterial_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_BF"] = ann.get(
                        "dka_bacterial_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKU_AF"] = ann.get(
                        "dku_archaeal_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_AF"] = ann.get(
                        "dka_archaeal_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKU_FF"] = ann.get(
                        "dku_fungal_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_FF"] = ann.get(
                        "dka_fungal_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKU_PF"] = ann.get(
                        "dku_protist_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_PF"] = ann.get(
                        "dka_protist_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKU_VF"] = ann.get(
                        "dku_viral_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_VF"] = ann.get(
                        "dka_viral_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKU_UCF"] = ann.get(
                        "dku_univec_core_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_UCF"] = ann.get(
                        "dka_univec_core_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKU_NHF"] = ann.get(
                        "dku_nonhuman_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_NHF"] = ann.get(
                        "dka_nonhuman_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKU_UF"] = ann.get(
                        "dku_unclassified_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_UF"] = ann.get(
                        "dka_unclassified_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKU_HLF"] = ann.get(
                        "dku_human_lineage_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_HLF"] = ann.get(
                        "dka_human_lineage_fraction", 0.0,
                    )
            else:
                rec.info["DKU"] = ann["dku"]
                rec.info["DKT"] = ann["dkt"]
                rec.info["DKA"] = ann["dka"]
                rec.info["DKU_DKT"] = ann["dku_dkt"]
                rec.info["DKA_DKT"] = ann["dka_dkt"]
                rec.info["MAX_PKC"] = ann["max_pkc"]
                rec.info["AVG_PKC"] = ann["avg_pkc"]
                rec.info["MIN_PKC"] = ann["min_pkc"]
                rec.info["MAX_PKC_ALT"] = ann["max_pkc_alt"]
                rec.info["AVG_PKC_ALT"] = ann["avg_pkc_alt"]
                rec.info["MIN_PKC_ALT"] = ann["min_pkc_alt"]
                if has_kraken_fractions:
                    rec.info["DKU_BF"] = ann.get("dku_bacterial_fraction", 0.0)
                    rec.info["DKA_BF"] = ann.get("dka_bacterial_fraction", 0.0)
                    rec.info["DKU_AF"] = ann.get("dku_archaeal_fraction", 0.0)
                    rec.info["DKA_AF"] = ann.get("dka_archaeal_fraction", 0.0)
                    rec.info["DKU_FF"] = ann.get("dku_fungal_fraction", 0.0)
                    rec.info["DKA_FF"] = ann.get("dka_fungal_fraction", 0.0)
                    rec.info["DKU_PF"] = ann.get("dku_protist_fraction", 0.0)
                    rec.info["DKA_PF"] = ann.get("dka_protist_fraction", 0.0)
                    rec.info["DKU_VF"] = ann.get("dku_viral_fraction", 0.0)
                    rec.info["DKA_VF"] = ann.get("dka_viral_fraction", 0.0)
                    rec.info["DKU_UCF"] = ann.get("dku_univec_core_fraction", 0.0)
                    rec.info["DKA_UCF"] = ann.get("dka_univec_core_fraction", 0.0)
                    rec.info["DKU_NHF"] = ann.get("dku_nonhuman_fraction", 0.0)
                    rec.info["DKA_NHF"] = ann.get("dka_nonhuman_fraction", 0.0)
                    rec.info["DKU_UF"] = ann.get("dku_unclassified_fraction", 0.0)
                    rec.info["DKA_UF"] = ann.get("dka_unclassified_fraction", 0.0)
                    rec.info["DKU_HLF"] = ann.get("dku_human_lineage_fraction", 0.0)
                    rec.info["DKA_HLF"] = ann.get("dka_human_lineage_fraction", 0.0)
        vcf_out.write(rec)

    vcf_out.close()
    vcf_in.close()

    pysam.tabix_index(output_vcf, preset="vcf", force=True)

    return output_vcf


def _write_informative_reads(
    child_bam, ref_fasta, informative_reads_by_variant, output_bam,
):
    """Write child reads carrying informative k-mers to a BAM file.

    Each output read is tagged with ``DV`` (the variant key it supports).
    Reads are sorted and indexed for IGV visualization.

    Args:
        child_bam: Path to the child BAM file.
        ref_fasta: Path to the reference FASTA.
        informative_reads_by_variant: dict mapping variant key
            (chrom:pos) to a set of read names.
        output_bam: Path for the output BAM file.
    """
    bam_in = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )

    unsorted_path = output_bam + ".unsorted.bam"
    bam_out = pysam.AlignmentFile(unsorted_path, "wb", header=bam_in.header)

    # Invert: read_name -> set of variant keys
    read_to_variants = {}
    for var_key, read_names in informative_reads_by_variant.items():
        for rname in read_names:
            read_to_variants.setdefault(rname, set()).add(var_key)

    # Collect unique regions to fetch
    regions = set()
    for var_key in informative_reads_by_variant:
        parts = var_key.split(":")
        chrom = parts[0]
        pos = int(parts[1])
        regions.add((chrom, pos))

    written = set()
    for chrom, pos in sorted(regions):
        for read in bam_in.fetch(chrom, pos, pos + 1):
            if read.query_name in read_to_variants and read.query_name not in written:
                var_keys = sorted(read_to_variants[read.query_name])
                read.set_tag("DV", ",".join(var_keys), value_type="Z")
                bam_out.write(read)
                written.add(read.query_name)

    bam_out.close()
    bam_in.close()

    pysam.sort("-o", output_bam, unsorted_path)
    pysam.index(output_bam)
    os.remove(unsorted_path)


def _write_summary(summary_path, variants, annotations):
    """Write a human-readable summary of variant stats and likely DNMs."""
    total = len(variants)
    likely_dnm = sum(1 for a in annotations.values() if a["dku"] > 0)
    inherited = total - likely_dnm

    dku_values = [a["dku"] for a in annotations.values()]
    dkt_values = [a["dkt"] for a in annotations.values()]
    dka_values = [a["dka"] for a in annotations.values()]
    dku_dkt_values = [a["dku_dkt"] for a in annotations.values()]
    dka_dkt_values = [a["dka_dkt"] for a in annotations.values()]
    max_pkc_values = [a["max_pkc"] for a in annotations.values()]
    avg_pkc_values = [a["avg_pkc"] for a in annotations.values()]
    min_pkc_values = [a["min_pkc"] for a in annotations.values()]
    max_pkc_alt_values = [a["max_pkc_alt"] for a in annotations.values()]
    avg_pkc_alt_values = [a["avg_pkc_alt"] for a in annotations.values()]
    min_pkc_alt_values = [a["min_pkc_alt"] for a in annotations.values()]
    dnm_dku = [a["dku"] for a in annotations.values() if a["dku"] > 0]

    lines = []
    lines.append("=" * 60)
    lines.append("  kmer-denovo  —  De Novo Variant Summary")
    lines.append("=" * 60)
    lines.append("")
    lines.append("Variant Counts")
    lines.append("-" * 40)
    lines.append(f"  Total candidates analyzed:   {total:>6}")
    lines.append(f"  Likely de novo (DKU > 0):    {likely_dnm:>6}")
    lines.append(f"  Inherited / unclear (DKU=0): {inherited:>6}")
    lines.append("")

    if dku_values:
        mean_dku = sum(dku_values) / len(dku_values)
        mean_dkt = sum(dkt_values) / len(dkt_values)
        mean_dka = sum(dka_values) / len(dka_values)
        mean_dku_dkt = sum(dku_dkt_values) / len(dku_dkt_values)
        mean_dka_dkt = sum(dka_dkt_values) / len(dka_dkt_values)
        median_dku = statistics.median(dku_values)
        mean_max_pkc = sum(max_pkc_values) / len(max_pkc_values)
        mean_avg_pkc = sum(avg_pkc_values) / len(avg_pkc_values)
        mean_min_pkc = sum(min_pkc_values) / len(min_pkc_values)
        mean_max_pkc_alt = sum(max_pkc_alt_values) / len(max_pkc_alt_values)
        mean_avg_pkc_alt = sum(avg_pkc_alt_values) / len(avg_pkc_alt_values)
        mean_min_pkc_alt = sum(min_pkc_alt_values) / len(min_pkc_alt_values)
        lines.append("Read Support Statistics")
        lines.append("-" * 40)
        lines.append(f"  DKU  mean:   {mean_dku:>6.1f}   median: {median_dku:>4}")
        lines.append(f"  DKT  mean:   {mean_dkt:>6.1f}")
        lines.append(f"  DKA  mean:   {mean_dka:>6.1f}")
        lines.append(f"  DKU_DKT  mean: {mean_dku_dkt:>6.4f}")
        lines.append(f"  DKA_DKT  mean: {mean_dka_dkt:>6.4f}")
        lines.append(f"  MAX_PKC  mean: {mean_max_pkc:>6.1f}")
        lines.append(f"  AVG_PKC  mean: {mean_avg_pkc:>6.1f}")
        lines.append(f"  MIN_PKC  mean: {mean_min_pkc:>6.1f}")
        lines.append(f"  MAX_PKC_ALT  mean: {mean_max_pkc_alt:>6.1f}")
        lines.append(f"  AVG_PKC_ALT  mean: {mean_avg_pkc_alt:>6.1f}")
        lines.append(f"  MIN_PKC_ALT  mean: {mean_min_pkc_alt:>6.1f}")
        lines.append("")

    if dnm_dku:
        mean_dnm_dku = sum(dnm_dku) / len(dnm_dku)
        lines.append(f"  Avg DKU among likely DNMs:   {mean_dnm_dku:>6.1f}")
        lines.append("")

    lines.append("Per-Variant Results")
    lines.append("-" * 120)
    lines.append(f"  {'Variant':<30s} {'DKU':>5s} {'DKT':>5s} {'DKA':>5s} {'DKU_DKT':>8s} {'DKA_DKT':>8s} {'MAX_PKC':>8s} {'AVG_PKC':>8s} {'MIN_PKC':>8s} {'MAX_PKC_ALT':>12s} {'AVG_PKC_ALT':>12s} {'MIN_PKC_ALT':>12s}  Call")
    lines.append(f"  {'-------':<30s} {'---':>5s} {'---':>5s} {'---':>5s} {'-------':>8s} {'-------':>8s} {'-------':>8s} {'-------':>8s} {'-------':>8s} {'-----------':>12s} {'-----------':>12s} {'-----------':>12s}  ----")

    for var in variants:
        ref = var["ref"]
        alts = var["alts"]
        alt = alts[0] if alts else "."
        var_key = f"{var['chrom']}:{var['pos']}:{ref}:{alt}"
        ann = annotations.get(var_key, {"dku": 0, "dkt": 0, "dka": 0, "dku_dkt": 0.0, "dka_dkt": 0.0, "max_pkc": 0, "avg_pkc": 0.0, "min_pkc": 0, "max_pkc_alt": 0, "avg_pkc_alt": 0.0, "min_pkc_alt": 0})
        label = f"{var['chrom']}:{var['pos'] + 1} {ref}>{alt}"
        call = "DE_NOVO" if ann["dku"] > 0 else "inherited"
        lines.append(f"  {label:<30s} {ann['dku']:>5d} {ann['dkt']:>5d} {ann['dka']:>5d} {ann['dku_dkt']:>8.4f} {ann['dka_dkt']:>8.4f} {ann['max_pkc']:>8d} {ann['avg_pkc']:>8.2f} {ann['min_pkc']:>8d} {ann['max_pkc_alt']:>12d} {ann['avg_pkc_alt']:>12.2f} {ann['min_pkc_alt']:>12d}  {call}")

    lines.append("")
    lines.append("=" * 60)
    lines.append("")

    text = "\n".join(lines)

    with open(summary_path, "w") as fh:
        fh.write(text)

    return text


def run_pipeline(args):
    """Run the de novo k-mer analysis pipeline."""
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

    kraken2_db = getattr(args, "kraken2_db", None)
    kraken2_confidence = getattr(args, "kraken2_confidence", 0.0)
    kraken2_memory_mapping = getattr(args, "kraken2_memory_mapping", False)
    if kraken2_db is not None:
        if not _check_tool("kraken2"):
            logger.error("kraken2 not found in PATH (required by --kraken2-db)")
            sys.exit(1)
        if not os.path.isdir(kraken2_db):
            logger.error("Kraken2 database not found: %s", kraken2_db)
            sys.exit(1)

    _validate_inputs(args)

    # ── Configuration summary ──────────────────────────────────────
    logger.info("=" * 60)
    logger.info("  kmer-denovo  —  pipeline starting")
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
    logger.info("  Input VCF:         %s", args.vcf)
    logger.info("  Output VCF:        %s", args.output)
    logger.info("  Reference FASTA:   %s", args.ref_fasta or "(not set)")
    logger.info("  k-mer size:        %d", args.kmer_size)
    logger.info("  Min base quality:  %d", args.min_baseq)
    logger.info("  Min mapping qual:  %d", args.min_mapq)
    logger.info("  Threads:           %d", args.threads)
    memory_limit_gb = getattr(args, "memory", None)
    logger.info(
        "  Memory limit:      %s",
        f"{memory_limit_gb:.1f} GB" if memory_limit_gb is not None
        else "(auto-detect)",
    )
    logger.info("  Proband ID:        %s", args.proband_id or "(not set)")
    logger.info("  Kraken2 DB:        %s", kraken2_db or "(disabled)")
    logger.info("=" * 60)

    # ── Step 1: Parse VCF ──────────────────────────────────────────
    step_start = time.monotonic()
    logger.info("[Step 1/5] Parsing VCF: %s", args.vcf)
    variants = _parse_vcf_variants(args.vcf)
    logger.info(
        "[Step 1/5] Found %d candidate variants (%s)",
        len(variants), _format_elapsed(time.monotonic() - step_start),
    )

    if not variants:
        logger.warning("No variants found in VCF; writing empty output")
        _write_annotated_vcf(args.vcf, args.output, {}, args.proband_id)
        if args.metrics:
            with open(args.metrics, "w") as fh:
                json.dump({"total_variants": 0}, fh, indent=2)
        logger.info(
            "Pipeline finished in %s",
            _format_elapsed(time.monotonic() - pipeline_start),
        )
        return

    # ── Step 2: Extract child k-mers ───────────────────────────────
    step_start = time.monotonic()
    logger.info(
        "[Step 2/5] Extracting child k-mers from %d variants (k=%d)",
        len(variants), args.kmer_size,
    )

    parent_found_kmers = collections.Counter()

    out_dir = os.path.dirname(os.path.abspath(args.output)) or "."
    tmp_root = _resolve_tmp_dir(args.tmp_dir, out_dir)
    logger.info("  Temp directory root: %s", tmp_root)
    if _is_tmpfs(tmp_root):
        logger.warning(
            "  ⚠ Temp directory %s appears to be on tmpfs (RAM-backed)! "
            "Consider using --tmp-dir to point to a disk-backed filesystem.",
            tmp_root,
        )
    _log_disk_usage(tmp_root, "tmpdir filesystem")

    with tempfile.TemporaryDirectory(prefix="kmer_denovo_",
                                     dir=tmp_root) as tmpdir:
        kmer_fasta = os.path.join(tmpdir, "child_kmers.fa")

        total_child_kmers, variant_read_kmers = _collect_child_kmers(
            args.child, args.ref_fasta, variants,
            args.kmer_size, args.min_baseq, args.min_mapq, args.debug_kmers,
            kmer_fasta,
        )
        logger.info(
            "[Step 2/5] Wrote %d child k-mers — partially deduplicated (%s)",
            total_child_kmers,
            _format_elapsed(time.monotonic() - step_start),
        )

        # ── Step 3: Scan parents ───────────────────────────────────
        step_start = time.monotonic()
        if total_child_kmers == 0:
            logger.info(
                "[Step 3/5] No child k-mers found; skipping parent scans"
            )
        else:
            logger.info(
                "[Step 3/5] Scanning parent BAMs for %d child k-mers",
                total_child_kmers,
            )

            parent_start = time.monotonic()
            logger.info(
                "[Step 3/5] ── Mother scan (1/2) ──",
            )
            mother_kmers = _scan_parent_jellyfish(
                args.mother, args.ref_fasta, kmer_fasta, args.kmer_size,
                os.path.join(tmpdir, "mother"), args.threads,
                n_filter_kmers=total_child_kmers,
            )
            parent_found_kmers.update(mother_kmers)
            logger.info(
                "[Step 3/5] Mother done — %d / %d child k-mers found in "
                "mother (%s)",
                len(mother_kmers), total_child_kmers,
                _format_elapsed(time.monotonic() - parent_start),
            )

            parent_start = time.monotonic()
            logger.info(
                "[Step 3/5] ── Father scan (2/2) ──",
            )
            father_kmers = _scan_parent_jellyfish(
                args.father, args.ref_fasta, kmer_fasta, args.kmer_size,
                os.path.join(tmpdir, "father"), args.threads,
                n_filter_kmers=total_child_kmers,
            )
            parent_found_kmers.update(father_kmers)
            logger.info(
                "[Step 3/5] Father done — %d / %d child k-mers found in "
                "father (%s)",
                len(father_kmers), total_child_kmers,
                _format_elapsed(time.monotonic() - parent_start),
            )

            logger.info(
                "[Step 3/5] Parent scanning complete — %d distinct "
                "child k-mers found across parents (%s)",
                len(parent_found_kmers),
                _format_elapsed(time.monotonic() - step_start),
            )

    child_unique_kmers = max(0, total_child_kmers - len(parent_found_kmers))

    # Remove the tmp_root if it was auto-created and is now empty
    try:
        if not getattr(args, "tmp_dir", None) and os.path.isdir(tmp_root):
            os.rmdir(tmp_root)  # only succeeds if empty
    except OSError:
        pass

    logger.info(
        "Child-unique k-mers (approx): %d / %d (%.1f%% unique)",
        child_unique_kmers,
        total_child_kmers,
        100 * child_unique_kmers / total_child_kmers if total_child_kmers else 0,
    )

    # ── Step 4: Annotate variants ──────────────────────────────────
    step_start = time.monotonic()
    logger.info(
        "[Step 4/5] Annotating %d variants with k-mer evidence",
        len(variants),
    )
    annotations = {}
    informative_reads_by_variant = {}
    informative_alt_reads_by_variant = {}
    n_variants = len(variants)
    log_interval = max(1, n_variants // 10)  # report ~10 times
    running_dnm = 0
    running_reads = 0

    # Materialise parent k-mer keys as a plain set once so that
    # per-read set operations (issubset / membership tests) are O(k)
    # instead of rebuilding a set from the Counter on every call.
    parent_kmer_set = set(parent_found_kmers)
    logger.info(
        "[Step 4/5] Parent k-mer lookup set: %d entries", len(parent_kmer_set),
    )

    for idx, var in enumerate(variants, 1):
        alt = var['alts'][0] if var['alts'] else "."
        var_key = f"{var['chrom']}:{var['pos']}:{var['ref']}:{alt}"
        read_kmers_list = variant_read_kmers.get(var_key, [])

        # Collect unique fragment names for all spanning reads.
        # When both mates of a paired-end read span the variant they
        # share the same query_name and are counted as one fragment.
        spanning_names = set()
        informative_names = set()
        informative_alt_names = set()
        all_variant_kmers = set()
        alt_variant_kmers = set()
        for read_name, kmers, supports_alt in read_kmers_list:
            spanning_names.add(read_name)
            all_variant_kmers.update(kmers)
            if supports_alt:
                alt_variant_kmers.update(kmers)
            # A fragment is informative if at least one of its
            # alignments has a variant-spanning k-mer absent from
            # both parents.
            if not kmers.issubset(parent_kmer_set):
                informative_names.add(read_name)
                if supports_alt:
                    informative_alt_names.add(read_name)

        # DKT, DKU, DKA are all fragment-based (unique read names)
        # so they are consistent with the BF denominators.
        dkt = len(spanning_names)
        dku = len(informative_names)
        dka = len(informative_alt_names)
        running_reads += dkt

        if dku > 0:
            running_dnm += 1

        # Compute parent k-mer count metrics for this variant
        parent_counts = [
            parent_found_kmers[k]
            for k in all_variant_kmers
            if k in parent_kmer_set
        ]
        max_pkc = max(parent_counts) if parent_counts else 0
        avg_pkc = round(statistics.mean(parent_counts), 2) if parent_counts else 0.0
        min_pkc = min(parent_counts) if parent_counts else 0

        # Compute parent k-mer count metrics for alt-allele-supporting k-mers
        alt_parent_counts = [
            parent_found_kmers[k]
            for k in alt_variant_kmers
            if k in parent_kmer_set
        ]
        max_pkc_alt = max(alt_parent_counts) if alt_parent_counts else 0
        avg_pkc_alt = round(statistics.mean(alt_parent_counts), 2) if alt_parent_counts else 0.0
        min_pkc_alt = min(alt_parent_counts) if alt_parent_counts else 0

        annotations[var_key] = {
            "dku": dku, "dkt": dkt, "dka": dka,
            "dku_dkt": round(dku / dkt, 4) if dkt > 0 else 0.0,
            "dka_dkt": round(dka / dkt, 4) if dkt > 0 else 0.0,
            "max_pkc": max_pkc, "avg_pkc": avg_pkc, "min_pkc": min_pkc,
            "max_pkc_alt": max_pkc_alt, "avg_pkc_alt": avg_pkc_alt, "min_pkc_alt": min_pkc_alt,
        }
        if informative_names:
            informative_reads_by_variant[var_key] = informative_names
        if informative_alt_names:
            informative_alt_reads_by_variant[var_key] = informative_alt_names

        if args.debug_kmers:
            logger.info("Variant %s: DKU=%d DKT=%d DKA=%d", var_key, dku, dkt, dka)

        if idx % log_interval == 0 or idx == n_variants:
            elapsed = time.monotonic() - step_start
            rate = idx / elapsed if elapsed > 0 else 0
            eta = (n_variants - idx) / rate if rate > 0 else 0
            logger.info(
                "[Step 4/5]   %d / %d variants (%.0f%%) — "
                "%d de novo so far, %d total reads "
                "(%.0f var/s, %s elapsed, ~%s remaining)",
                idx, n_variants, 100 * idx / n_variants,
                running_dnm, running_reads,
                rate, _format_elapsed(elapsed), _format_elapsed(eta),
            )

    likely_dnm = running_dnm
    logger.info(
        "[Step 4/5] Annotation complete — %d likely de novo, %d inherited (%s)",
        likely_dnm,
        n_variants - likely_dnm,
        _format_elapsed(time.monotonic() - step_start),
    )

    # ── Kraken2 non-human content flagging (VCF mode) ────────────────
    kraken2_result = None
    name_map = None
    if kraken2_db is not None:
        step_start = time.monotonic()
        all_informative_names = set()
        for names in informative_reads_by_variant.values():
            all_informative_names.update(names)
        logger.info(
            "[Kraken2] Classifying %d informative reads for non-human content",
            len(all_informative_names),
        )
        kraken2_result = _run_kraken2_on_reads(
            args.child, args.ref_fasta, all_informative_names,
            kraken2_db, confidence=kraken2_confidence,
            threads=args.threads,
            informative_reads_by_variant=informative_reads_by_variant,
            memory_mapping=kraken2_memory_mapping,
        )
        logger.info(
            "[Kraken2] %s (%s)",
            kraken2_result.summary(),
            _format_elapsed(time.monotonic() - step_start),
        )

        # Load taxon name map for per-read detail BED
        name_map = Kraken2Runner._load_name_map(kraken2_db)

        for var_key, ann in annotations.items():
            dku_names = informative_reads_by_variant.get(var_key, set())
            dka_names = informative_alt_reads_by_variant.get(var_key, set())

            for label, read_set in (
                ("bacterial", kraken2_result.bacterial_read_names),
                ("archaeal", kraken2_result.archaeal_read_names),
                ("fungal", kraken2_result.fungal_read_names),
                ("protist", kraken2_result.protist_read_names),
                ("viral", kraken2_result.viral_read_names),
                ("univec_core", kraken2_result.univec_core_read_names),
                ("nonhuman", kraken2_result.nonhuman_read_names),
                ("unclassified", kraken2_result.unclassified_read_names),
                ("human_lineage", kraken2_result.human_lineage_read_names),
            ):
                dku_count = len(dku_names.intersection(read_set))
                dka_count = len(dka_names.intersection(read_set))

                ann[f"dku_{label}_fraction"] = (
                    round(dku_count / len(dku_names), _FRACTION_PRECISION)
                    if dku_names else 0.0
                )
                ann[f"dka_{label}_fraction"] = (
                    round(dka_count / len(dka_names), _FRACTION_PRECISION)
                    if dka_names else 0.0
                )

    # ── Step 5: Write outputs ──────────────────────────────────────
    step_start = time.monotonic()
    logger.info("[Step 5/5] Writing output files")

    logger.info("[Step 5/5] Writing annotated VCF: %s", args.output)
    actual_output = _write_annotated_vcf(
        args.vcf, args.output, annotations, args.proband_id,
    )

    if args.informative_reads:
        logger.info(
            "[Step 5/5] Writing informative reads BAM: %s",
            args.informative_reads,
        )
        _write_informative_reads(
            args.child, args.ref_fasta,
            informative_reads_by_variant, args.informative_reads,
        )
        total_reads = sum(
            len(names)
            for names in informative_reads_by_variant.values()
        )
        logger.info(
            "[Step 5/5] Wrote %d informative reads across %d variants",
            total_reads, len(informative_reads_by_variant),
        )

    # Write per-read Kraken2 classification detail BED
    if kraken2_result is not None:
        kraken2_detail_path = getattr(args, "kraken2_read_detail", None)
        if kraken2_detail_path is None:
            # Auto-derive from --output
            base = args.output
            for ext in (".vcf.gz", ".vcf.bgz", ".vcf"):
                if base.endswith(ext):
                    base = base[: -len(ext)]
                    break
            kraken2_detail_path = base + ".kraken2_reads.bed.gz"
        logger.info(
            "[Step 5/5] Writing per-read Kraken2 detail BED: %s",
            kraken2_detail_path,
        )
        _write_kraken2_read_detail_bed(
            kraken2_detail_path,
            informative_reads_by_variant,
            informative_alt_reads_by_variant,
            kraken2_result,
            name_map,
        )
        logger.info(
            "[Step 5/5] Wrote per-read Kraken2 detail for %d variants",
            len(informative_reads_by_variant),
        )

        # Write species-annotated genomic span BED
        kraken2_span_path = getattr(args, "kraken2_span_bed", None)
        if kraken2_span_path is None:
            # Auto-derive from --output
            span_base = args.output
            for ext in (".vcf.gz", ".vcf.bgz", ".vcf"):
                if span_base.endswith(ext):
                    span_base = span_base[: -len(ext)]
                    break
            kraken2_span_path = span_base + ".kraken2_spans.bed.gz"
        logger.info(
            "[Step 5/5] Collecting alignment metadata for span BED",
        )
        alignment_meta = _collect_read_alignment_metadata(
            args.child, args.ref_fasta, all_informative_names,
            informative_reads_by_variant=informative_reads_by_variant,
        )
        logger.info(
            "[Step 5/5] Writing species-annotated span BED: %s",
            kraken2_span_path,
        )
        _write_kraken2_span_bed(
            kraken2_span_path,
            alignment_meta,
            informative_reads_by_variant,
            informative_alt_reads_by_variant,
            kraken2_result,
            name_map,
        )
        logger.info(
            "[Step 5/5] Wrote span BED for %d reads",
            len(alignment_meta),
        )

        # Write expanded span BED (unless --no-expanded-bed)
        if not getattr(args, "no_expanded_bed", False):
            expanded_path = kraken2_span_path.replace(
                ".kraken2_spans.bed.gz",
                ".kraken2_spans_expanded.bed.gz",
            )
            if expanded_path == kraken2_span_path:
                # Fallback if span path doesn't match expected pattern
                expanded_path = kraken2_span_path.replace(
                    ".bed.gz", "_expanded.bed.gz",
                )
            logger.info(
                "[Step 5/5] Writing expanded span BED: %s",
                expanded_path,
            )
            _write_kraken2_expanded_span_bed(
                expanded_path,
                alignment_meta,
                informative_reads_by_variant,
                informative_alt_reads_by_variant,
                kraken2_result,
                name_map,
            )
            logger.info(
                "[Step 5/5] Wrote expanded span BED for %d reads",
                len(alignment_meta),
            )

    if args.metrics:
        metrics = {
            "total_variants": len(variants),
            "total_child_kmers": total_child_kmers,
            "parent_found_kmers": len(parent_found_kmers),
            "child_unique_kmers": child_unique_kmers,
            "variants_with_unique_reads": likely_dnm,
        }
        if kraken2_result is not None:
            metrics["kraken2"] = {
                "total_reads_classified": kraken2_result.total,
                "classified": kraken2_result.classified,
                "unclassified": kraken2_result.unclassified,
                "bacterial_reads": kraken2_result.bacterial_count,
                "archaeal_reads": kraken2_result.archaeal_count,
                "fungal_reads": kraken2_result.fungal_count,
                "protist_reads": kraken2_result.protist_count,
                "viral_reads": kraken2_result.viral_count,
                "univec_core_reads": kraken2_result.univec_core_count,
                "nonhuman_reads": kraken2_result.nonhuman_count,
                "human_reads": kraken2_result.human_count,
                "root_reads": kraken2_result.root_count,
                "bacterial_fraction": kraken2_result.bacterial_fraction,
            }
        with open(args.metrics, "w") as fh:
            json.dump(metrics, fh, indent=2)
        logger.info("[Step 5/5] Metrics written to: %s", args.metrics)

    if args.summary:
        logger.info("[Step 5/5] Writing summary: %s", args.summary)
        _write_summary(args.summary, variants, annotations)

    # ── Optional interactive HTML report ───────────────────────────
    report_path = getattr(args, "report", None)
    if report_path:
        logger.info("[Step 5/5] Generating interactive HTML report: %s", report_path)
        from kmer_denovo_filter.report import generate_report
        generate_report(
            output_path=report_path,
            vcf_metrics_path=args.metrics,
            vcf_summary_path=args.summary,
            vcf_path=actual_output,
        )

    logger.info(
        "[Step 5/5] Output complete (%s)",
        _format_elapsed(time.monotonic() - step_start),
    )

    # ── Done ───────────────────────────────────────────────────────
    logger.info(
        "Pipeline finished successfully in %s",
        _format_elapsed(time.monotonic() - pipeline_start),
    )
