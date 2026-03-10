"""Main pipeline for de novo variant k-mer analysis."""

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

from kmer_denovo_filter.kmer_utils import (
    JellyfishKmerQuery,
    Kraken2Runner,
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
    _get_available_memory_gb,
    _infer_sv_type,
    _is_tmpfs,
    _load_kmers_from_fasta,
    _log_children_memory,
    _log_dir_size,
    _log_disk_usage,
    _log_memory,
    _log_subprocess_memory,
    _resolve_tmp_dir,
    _write_kmer_fasta,
)

logger = logging.getLogger(__name__)
_BACTERIAL_FRACTION_PRECISION = 4


def _run_kraken2_on_reads(
    child_bam, ref_fasta, read_names, kraken2_db,
    confidence=0.0, threads=1, tmpdir=None,
    informative_reads_by_variant=None,
):
    """Classify child reads with kraken2 and return a result summary.

    Extracts sequences for *read_names* from *child_bam*, classifies them
    with kraken2, and returns a :class:`Kraken2Runner.Result` with tallied
    bacterial / human / root counts.

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
        kraken2_db, confidence=confidence, threads=threads,
    )
    return kr.classify_sequences(sequences, tmpdir=tmpdir)


def _validate_inputs(args):
    """Validate pipeline inputs before starting computation.

    Raises SystemExit with a clear error message for any invalid input.
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


def _scan_parent_jellyfish(
    parent_bam, ref_fasta, kmer_fasta, kmer_size, parent_dir, threads=4,
):
    """Scan a parent BAM and find which child k-mers are present.

    Uses ``jellyfish count`` with the ``--if`` filter so only child k-mers
    are tracked while the parent BAM is streamed.  The dump output is
    streamed line-by-line and the parent Jellyfish index is removed after
    the dump to free disk space.

    Returns:
        Dict mapping canonical k-mer string to its count in the parent.
    """
    os.makedirs(parent_dir, exist_ok=True)
    jf_output = os.path.join(parent_dir, "parent.jf")

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
        "-s", "10M",
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
        "  samtools fasta → jellyfish count (k=%d, threads=%d)",
        kmer_size, threads,
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
            jf_size = _format_file_size(jf_output) if os.path.exists(jf_output) else "pending"
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

    When bacterial-fraction annotations are present in *annotations*,
    DKU_BF and DKA_BF are also added.

    Returns:
        The actual output path (with ``.gz`` suffix).
    """
    vcf_in = pysam.VariantFile(input_vcf)
    has_bacterial_fraction = any(
        "dku_bacterial_fraction" in ann or "dka_bacterial_fraction" in ann
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
             "Number of child reads with at least one variant-spanning "
             "k-mer unique to child (absent from both parents)"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKT"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Total child reads with variant-spanning k-mers"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKA"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description",
             "Number of child reads with at least one unique k-mer "
             "that also exactly support the candidate allele"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKU_DKT"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description",
             "Proportion of child reads with unique k-mers (DKU/DKT)"),
        ],
    )
    vcf_in.header.add_meta(
        category,
        items=[
            ("ID", "DKA_DKT"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description",
             "Proportion of child reads with unique allele-supporting "
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
    if has_bacterial_fraction:
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKU_BF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of unique DKU read names (fragments with >=1 "
                 "child-unique variant-spanning k-mer) classified as bacterial "
                 "by kraken2; denominator equals the number of unique fragments, "
                 "which may be less than DKU when both mates span the locus"),
            ],
        )
        vcf_in.header.add_meta(
            category,
            items=[
                ("ID", "DKA_BF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description",
                 "Fraction of unique DKA read names (DKU fragments that also "
                 "support the alternate allele) classified as bacterial by "
                 "kraken2; DKA names are always a subset of DKU names"),
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
                if has_bacterial_fraction:
                    rec.samples[proband_id]["DKU_BF"] = ann.get(
                        "dku_bacterial_fraction", 0.0,
                    )
                    rec.samples[proband_id]["DKA_BF"] = ann.get(
                        "dka_bacterial_fraction", 0.0,
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
                if has_bacterial_fraction:
                    rec.info["DKU_BF"] = ann.get("dku_bacterial_fraction", 0.0)
                    rec.info["DKA_BF"] = ann.get("dka_bacterial_fraction", 0.0)
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
        alt = alts[0] if alts else "."
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


# ── Discovery-mode helpers ────────────────────────────────────────


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

    extract_start = time.monotonic()
    p_samtools = subprocess.Popen(
        samtools_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    p_jellyfish = subprocess.Popen(
        jellyfish_cmd, stdin=p_samtools.stdout,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    p_samtools.stdout.close()

    # Poll for completion with periodic progress logging
    poll_interval = 60
    while True:
        try:
            p_jellyfish.wait(timeout=poll_interval)
            break
        except subprocess.TimeoutExpired:
            elapsed = time.monotonic() - extract_start
            jf_files = _find_jf_files(child_jf)
            if jf_files:
                total_size = sum(
                    os.path.getsize(f) for f in jf_files
                    if os.path.exists(f)
                )
                jf_size = _format_file_size.__wrapped__(total_size) \
                    if hasattr(_format_file_size, '__wrapped__') \
                    else f"{total_size / (1024**3):.1f} GB"
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

    jf_stderr = p_jellyfish.stderr.read() if p_jellyfish.stderr else b""
    p_samtools.communicate()

    if p_jellyfish.returncode != 0:
        raise RuntimeError(
            f"jellyfish count (child) failed: "
            f"{jf_stderr.decode() if jf_stderr else ''}"
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

    scan_start = time.monotonic()
    p_samtools = subprocess.Popen(
        samtools_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    p_jellyfish = subprocess.Popen(
        jellyfish_cmd, stdin=p_samtools.stdout,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    p_samtools.stdout.close()

    # Poll for completion with periodic progress logging
    poll_interval = 30
    while True:
        try:
            p_jellyfish.wait(timeout=poll_interval)
            break
        except subprocess.TimeoutExpired:
            elapsed = time.monotonic() - scan_start
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

    p_samtools.communicate()
    jf_stderr = p_jellyfish.stderr.read()

    if p_jellyfish.returncode != 0:
        raise RuntimeError(
            f"jellyfish count ({label}) failed: {jf_stderr.decode()}"
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


# ── Multiprocessing helpers for _anchor_and_cluster ────────────────────

_worker_automaton = None       # Aho-Corasick automaton (small k-mer sets)
_worker_jf_query = None        # JellyfishKmerQuery (large k-mer sets)
_worker_kmer_size = None
_worker_min_distinct_kmers_per_read = 1

# Number of reads to accumulate before issuing a single jellyfish
# subprocess call.  Larger batches amortize subprocess overhead
# but temporarily hold more read objects in memory.
_JF_READ_BATCH_SIZE = 5000


def _init_scan_worker(proband_data, kmer_size,
                      min_distinct_kmers_per_read=1):
    """Initializer for per-contig scan workers.

    *proband_data* may be:
    - A path ending in ``.jf`` — opens a :class:`JellyfishKmerQuery`
      that queries k-mers against a memory-mapped jellyfish hash.
      This is the low-memory path used for WGS discovery mode with
      hundreds of millions of proband-unique k-mers.
    - A FASTA file path — loads k-mers into a Python set and builds
      an Aho-Corasick automaton (fast but memory-intensive).
    - A Python set — builds an Aho-Corasick automaton directly.
    """
    global _worker_automaton, _worker_jf_query, _worker_kmer_size
    global _worker_min_distinct_kmers_per_read

    _worker_automaton = None
    _worker_jf_query = None

    if isinstance(proband_data, str) and proband_data.endswith(".jf"):
        # Jellyfish-backed mode: each query subprocess memory-maps
        # the same .jf file; the OS page cache is shared across workers.
        _worker_jf_query = JellyfishKmerQuery(proband_data)
    elif isinstance(proband_data, str):
        kmers = _load_kmers_from_fasta(proband_data)
        _worker_automaton = build_kmer_automaton(kmers)
        del kmers
    else:
        _worker_automaton = build_kmer_automaton(proband_data)

    _worker_kmer_size = kmer_size
    _worker_min_distinct_kmers_per_read = min_distinct_kmers_per_read


def _process_informative_read(read, unique_in_read, kmer_hit_indices,
                              kmer_size, reads_seen, read_hits,
                              read_sv_meta, kmer_coverage, read_coverage):
    """Record an informative read's hits, coverage, and SV metadata.

    Returns 1 if the read is unmapped-informative, 0 otherwise.
    Mutates *reads_seen*, *read_hits*, *read_sv_meta*, *kmer_coverage*,
    and *read_coverage* in place.
    """
    dedup_key = (read.query_name, read.is_supplementary)
    if dedup_key in reads_seen:
        return 0

    reads_seen.add(dedup_key)
    if read.is_unmapped:
        return 1

    read_hits.append((
        read.reference_name,
        read.reference_start,
        read.reference_end,
        read.query_name,
        unique_in_read,
        read.is_supplementary,
    ))
    # Map novel k-mer query positions to reference coords
    chrom = read.reference_name
    cov = _collect_kmer_ref_positions(
        read, kmer_hit_indices, kmer_size,
    )
    kmer_coverage[chrom] += cov
    # Count one read per touched position
    for pos in cov:
        read_coverage[chrom][pos] += 1

    # Collect SV metadata for this informative read
    max_clip = 0
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == 4 and length > max_clip:  # soft clip
                max_clip = length
    read_sv_meta[dedup_key] = {
        "has_sa": read.has_tag("SA"),
        "sa_str": read.get_tag("SA") if (
            read.has_tag("SA") and not read.is_supplementary
        ) else None,
        "is_paired": read.is_paired,
        "is_proper_pair": read.is_proper_pair,
        "mate_is_unmapped": (
            read.mate_is_unmapped if read.is_paired else False
        ),
        "max_clip": max_clip,
    }
    return 0


def _scan_contig_for_hits(child_bam, ref_fasta, contig):
    """Scan reads mapped to *contig* for proband-unique k-mers.

    When *contig* is ``None``, unmapped reads are scanned instead.

    Supports two scanning backends:
    - **Aho-Corasick automaton** (``_worker_automaton``) — fast C-level
      multi-pattern matching; used when the k-mer set fits in memory.
    - **JellyfishKmerQuery** (``_worker_jf_query``) — disk-backed
      queries via ``jellyfish query``; used for large k-mer sets in
      discovery mode where the set is too large for Aho-Corasick.

    Returns:
        (read_hits, reads_seen, unmapped_informative, total_reads_scanned,
         read_sv_meta, kmer_coverage, read_coverage)

    ``read_sv_meta`` is a dict keyed by ``(query_name, is_supplementary)``
    with per-read SV metadata (has_sa, sa_str, is_paired, is_proper_pair,
    mate_is_unmapped, max_clip) collected for each informative read so
    that annotation and linking can be done without re-scanning the BAM.

    ``kmer_coverage`` is a dict mapping chrom to a Counter of reference
    positions overlapped by novel k-mers (counts total k-mer base
    overlaps across all reads).

    ``read_coverage`` is a dict mapping chrom to a Counter of reference
    positions where at least one novel k-mer was found, counting the
    number of distinct reads touching each position.
    """
    automaton = _worker_automaton
    jf_query = _worker_jf_query
    kmer_size = _worker_kmer_size
    min_dk_per_read = _worker_min_distinct_kmers_per_read
    bam = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )

    read_hits = []
    reads_seen = set()
    read_sv_meta = {}
    kmer_coverage = collections.defaultdict(collections.Counter)
    read_coverage = collections.defaultdict(collections.Counter)
    unmapped_informative = 0
    total_reads_scanned = 0

    if contig is None:
        try:
            iterator = bam.fetch("*")
        except (ValueError, KeyError):
            bam.close()
            return (read_hits, reads_seen, unmapped_informative,
                    total_reads_scanned, read_sv_meta, kmer_coverage,
                    read_coverage)
    else:
        iterator = bam.fetch(contig=contig)

    if jf_query is not None:
        # ── Batched jellyfish path ─────────────────────────────────
        # Reads are accumulated in batches so that k-mers from many
        # reads are queried in a single jellyfish subprocess call.
        # This reduces subprocess overhead from O(n_reads) to
        # O(n_reads / batch_size).
        pending = []   # (read, canon_at_pos)
        pending_kmers = set()

        for read in iterator:
            if read.is_secondary:
                continue
            if read.is_duplicate:
                continue

            total_reads_scanned += 1
            seq = read.query_sequence
            if seq is None:
                continue

            canon_at_pos, unique_candidates = _extract_read_kmers(
                seq, kmer_size,
            )
            pending_kmers.update(unique_candidates)
            pending.append((read, canon_at_pos))

            if len(pending) < _JF_READ_BATCH_SIZE:
                continue

            # Query all unique k-mers from this batch in one subprocess.
            # Avoid unbounded per-worker cache growth by processing each
            # batch against this local hit set, then clearing cache.
            batch_hits = set()
            if pending_kmers:
                batch_hits = jf_query.query_batch(list(pending_kmers))
                pending_kmers = set()

            # Process each read against batch-level hits
            for read_obj, c_at_pos in pending:
                unique_in_read = set()
                kmer_hit_indices = set()
                for pos, canon in c_at_pos.items():
                    if canon in batch_hits:
                        unique_in_read.add(canon)
                        kmer_hit_indices.add(pos)

                if len(unique_in_read) < min_dk_per_read:
                    continue

                unmapped_informative += _process_informative_read(
                    read_obj, unique_in_read, kmer_hit_indices,
                    kmer_size, reads_seen, read_hits,
                    read_sv_meta, kmer_coverage, read_coverage,
                )
            pending = []
            jf_query.close()

        # Flush remaining reads
        if pending:
            batch_hits = set()
            if pending_kmers:
                batch_hits = jf_query.query_batch(list(pending_kmers))
            for read_obj, c_at_pos in pending:
                unique_in_read = set()
                kmer_hit_indices = set()
                for pos, canon in c_at_pos.items():
                    if canon in batch_hits:
                        unique_in_read.add(canon)
                        kmer_hit_indices.add(pos)

                if len(unique_in_read) < min_dk_per_read:
                    continue

                unmapped_informative += _process_informative_read(
                    read_obj, unique_in_read, kmer_hit_indices,
                    kmer_size, reads_seen, read_hits,
                    read_sv_meta, kmer_coverage, read_coverage,
                )
            jf_query.close()
    else:
        # ── Aho-Corasick path (or no backend) ──────────────────────
        for read in iterator:
            if read.is_secondary:
                continue
            if read.is_duplicate:
                continue

            total_reads_scanned += 1
            seq = read.query_sequence
            if seq is None:
                continue

            unique_in_read = set()
            kmer_hit_indices = set()
            if automaton is not None:
                for _end_idx, canonical_kmer in automaton.iter(seq):
                    unique_in_read.add(canonical_kmer)
                    kmer_hit_indices.add(_end_idx - kmer_size + 1)

            if len(unique_in_read) < min_dk_per_read:
                continue

            unmapped_informative += _process_informative_read(
                read, unique_in_read, kmer_hit_indices,
                kmer_size, reads_seen, read_hits,
                read_sv_meta, kmer_coverage, read_coverage,
            )

    bam.close()
    return (read_hits, reads_seen, unmapped_informative,
            total_reads_scanned, read_sv_meta, kmer_coverage,
            read_coverage)


def _anchor_and_cluster(child_bam, ref_fasta, proband_unique_kmers,
                        kmer_size, merge_distance=500, threads=1,
                        min_distinct_kmers_per_read=1,
                        proband_unique_fa=None,
                        proband_jf=None,
                        n_proband_unique=None,
                        tmpdir=None,
                        memory_limit_gb=None):
    """Module 3: Find reads containing proband-unique k-mers and cluster regions.

    Scans **all** primary, non-duplicate child reads — including unmapped
    and low-MAPQ reads — so that no proband-unique k-mers are missed at
    this stage.

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
        tmpdir: Writable directory for temporary files.
        memory_limit_gb: Explicit memory limit in GB.  When provided,
            overrides auto-detected system memory for worker-count
            planning.  Useful on HPC where SLURM allocations differ
            from total node memory.

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

    if threads > 1:
        # ── Parallel scanning by chromosome ────────────────────────
        bam = pysam.AlignmentFile(
            child_bam,
            reference_filename=ref_fasta if ref_fasta else None,
        )
        contigs = list(bam.references)
        bam.close()

        # One task per contig + one for unmapped reads
        tasks = [(child_bam, ref_fasta, c) for c in contigs]
        tasks.append((child_bam, ref_fasta, None))

        # ── Dynamically cap workers based on available memory ──────
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
            "  Parallel anchoring: %d contigs, %d workers "
            "(requested=%d, mode=%s)",
            len(contigs), n_workers, threads,
            "jellyfish" if use_jellyfish else "aho-corasick",
        )

        read_hits = []
        reads_seen = set()
        read_sv_meta = {}
        kmer_coverage = collections.defaultdict(collections.Counter)
        read_coverage = collections.defaultdict(collections.Counter)
        unmapped_informative = 0
        total_reads_scanned = 0
        completed_contigs = 0

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
            init_args = (proband_unique_kmers, kmer_size,
                         min_distinct_kmers_per_read)
            init_mode = "in-memory-set"

        logger.info(
            "  Worker init mode: %s", init_mode,
        )

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_init_scan_worker,
            initargs=init_args,
        ) as executor:
            futures = {
                executor.submit(_scan_contig_for_hits, *t): t[2]
                for t in tasks
            }

            # Use wait() with timeout for time-based progress reporting
            # so users get feedback even when large contigs take hours.
            pending = set(futures.keys())
            last_progress_log = time.monotonic()
            progress_interval = 300  # seconds (5 minutes)

            while pending:
                done, pending = concurrent.futures.wait(
                    pending, timeout=progress_interval,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )

                if not done:
                    # Timeout — no contig completed; log heartbeat
                    elapsed = time.monotonic() - anchor_start
                    logger.info(
                        "  [Anchoring] Heartbeat: %d/%d contigs complete, "
                        "%d reads scanned, %d informative (%s elapsed)",
                        completed_contigs, len(tasks),
                        total_reads_scanned,
                        len(read_hits) + unmapped_informative,
                        _format_elapsed(elapsed),
                    )
                    _log_memory("during anchoring")
                    _log_children_memory("during anchoring")
                    last_progress_log = time.monotonic()
                    continue

                for future in done:
                    contig = futures[future]
                    try:
                        (hits, seen, unmapped, scanned,
                         sv_meta, worker_cov,
                         worker_read_cov) = future.result()
                    except Exception:
                        logger.error(
                            "Worker failed for contig=%s", contig,
                        )
                        raise
                    total_reads_scanned += scanned
                    unmapped_informative += unmapped
                    completed_contigs += 1
                    for hit in hits:
                        # hit: (ref_name, start, end, qname, kmers, is_supp)
                        dedup_key = (hit[3], hit[5])  # (qname, is_supplementary)
                        if dedup_key not in reads_seen:
                            reads_seen.add(dedup_key)
                            read_hits.append(hit)
                    # Merge SV metadata (only for keys not already seen)
                    for key, meta in sv_meta.items():
                        if key not in read_sv_meta:
                            read_sv_meta[key] = meta
                    # Merge k-mer coverage (update is faster than += for Counters)
                    for chrom, cov in worker_cov.items():
                        kmer_coverage[chrom].update(cov)
                    # Merge read coverage
                    for chrom, cov in worker_read_cov.items():
                        read_coverage[chrom].update(cov)
                    # Track all seen keys for cross-contig dedup
                    reads_seen.update(seen)

                    # Log individual contig completion for large contigs
                    if scanned >= 1_000_000:
                        logger.info(
                            "  [Anchoring] Contig %s complete: %d reads "
                            "scanned, %d informative hits (%s elapsed)",
                            contig or "(unmapped)", scanned,
                            len(hits) + unmapped,
                            _format_elapsed(time.monotonic() - anchor_start),
                        )

                # Log progress: time-based (every 5 min) or milestone-based
                now = time.monotonic()
                time_since_log = now - last_progress_log
                at_milestone = (
                    completed_contigs % 100 == 0
                    or completed_contigs == len(tasks)
                )
                if time_since_log >= progress_interval or at_milestone:
                    pct = (
                        100 * completed_contigs / len(tasks)
                        if tasks else 0
                    )
                    logger.info(
                        "  [Anchoring] Progress: %d/%d contigs complete "
                        "(%.0f%%), %d reads scanned, %d informative (%s)",
                        completed_contigs, len(tasks),
                        pct,
                        total_reads_scanned,
                        len(read_hits) + unmapped_informative,
                        _format_elapsed(now - anchor_start),
                    )
                    _log_memory("during anchoring")
                    _log_children_memory("during anchoring")
                    last_progress_log = now

        logger.info(
            "  Anchoring: %d reads scanned, %d informative (%s) [%d workers]",
            total_reads_scanned,
            len(read_hits) + unmapped_informative,
            _format_elapsed(time.monotonic() - anchor_start),
            n_workers,
        )
    else:
        # ── Single-threaded scanning ────────────────────────────────
        jf_query_st = None
        automaton = None

        if use_jellyfish:
            jf_query_st = JellyfishKmerQuery(proband_jf)
            logger.info(
                "  Single-threaded scan using jellyfish query: %s",
                proband_jf,
            )
        elif proband_unique_fa:
            kmer_data = _load_kmers_from_fasta(proband_unique_fa)
            logger.info(
                "  Loaded %d k-mers from FASTA for single-threaded scan",
                len(kmer_data),
            )
            _log_memory("after k-mer load (single-threaded)")
            automaton = build_kmer_automaton(kmer_data)
            del kmer_data
        else:
            kmer_data = proband_unique_kmers or set()
            _log_memory("after k-mer load (single-threaded)")
            automaton = build_kmer_automaton(kmer_data)
            del kmer_data

        bam = pysam.AlignmentFile(
            child_bam,
            reference_filename=ref_fasta if ref_fasta else None,
        )

        read_hits = []
        reads_seen = set()
        read_sv_meta = {}
        kmer_coverage = collections.defaultdict(collections.Counter)
        read_coverage = collections.defaultdict(collections.Counter)
        unmapped_informative = 0
        total_reads_scanned = 0

        if jf_query_st is not None:
            # Batched jellyfish path (single-threaded)
            pending = []
            pending_kmers = set()

            for read in bam.fetch():
                if read.is_secondary:
                    continue
                if read.is_duplicate:
                    continue

                total_reads_scanned += 1
                seq = read.query_sequence
                if seq is None:
                    continue

                canon_at_pos, unique_candidates = _extract_read_kmers(
                    seq, kmer_size,
                )
                pending_kmers.update(unique_candidates)
                pending.append((read, canon_at_pos, unique_candidates))

                if len(pending) < _JF_READ_BATCH_SIZE:
                    continue

                if pending_kmers:
                    jf_query_st.query_batch(list(pending_kmers))
                    pending_kmers = set()

                for read_obj, c_at_pos, u_cands in pending:
                    hits = jf_query_st.query_batch(u_cands)
                    unique_in_read = set()
                    kmer_hit_indices = set()
                    for pos, canon in c_at_pos.items():
                        if canon in hits:
                            unique_in_read.add(canon)
                            kmer_hit_indices.add(pos)

                    if len(unique_in_read) < min_distinct_kmers_per_read:
                        continue

                    unmapped_informative += _process_informative_read(
                        read_obj, unique_in_read, kmer_hit_indices,
                        kmer_size, reads_seen, read_hits,
                        read_sv_meta, kmer_coverage, read_coverage,
                    )
                pending = []

                if total_reads_scanned % 1_000_000 == 0:
                    logger.info(
                        "  Anchoring: %d reads scanned, %d informative (%s)",
                        total_reads_scanned,
                        len(read_hits) + unmapped_informative,
                        _format_elapsed(time.monotonic() - anchor_start),
                    )

            # Flush remaining
            if pending:
                if pending_kmers:
                    jf_query_st.query_batch(list(pending_kmers))
                for read_obj, c_at_pos, u_cands in pending:
                    hits = jf_query_st.query_batch(u_cands)
                    unique_in_read = set()
                    kmer_hit_indices = set()
                    for pos, canon in c_at_pos.items():
                        if canon in hits:
                            unique_in_read.add(canon)
                            kmer_hit_indices.add(pos)

                    if len(unique_in_read) < min_distinct_kmers_per_read:
                        continue

                    unmapped_informative += _process_informative_read(
                        read_obj, unique_in_read, kmer_hit_indices,
                        kmer_size, reads_seen, read_hits,
                        read_sv_meta, kmer_coverage, read_coverage,
                    )
        else:
            for read in bam.fetch():
                if read.is_secondary:
                    continue
                if read.is_duplicate:
                    continue

                total_reads_scanned += 1
                seq = read.query_sequence
                if seq is None:
                    continue

                unique_in_read = set()
                kmer_hit_indices = set()
                if automaton is not None:
                    for _end_idx, canonical_kmer in automaton.iter(seq):
                        unique_in_read.add(canonical_kmer)
                        kmer_hit_indices.add(_end_idx - kmer_size + 1)

                # Per-read filter: require a minimum number of distinct kmers
                if len(unique_in_read) < min_distinct_kmers_per_read:
                    unique_in_read = set()
                    kmer_hit_indices = set()

                dedup_key = (read.query_name, read.is_supplementary)
                if unique_in_read and dedup_key not in reads_seen:
                    reads_seen.add(dedup_key)
                    if read.is_unmapped:
                        unmapped_informative += 1
                    else:
                        read_hits.append((
                            read.reference_name,
                            read.reference_start,
                            read.reference_end,
                            read.query_name,
                            unique_in_read,
                            read.is_supplementary,
                        ))
                        # Map novel k-mer query positions to reference coords
                        chrom = read.reference_name
                        cov = _collect_kmer_ref_positions(
                            read, kmer_hit_indices, kmer_size,
                        )
                        kmer_coverage[chrom] += cov
                        # Count one read per touched position
                        for pos in cov:
                            read_coverage[chrom][pos] += 1

                    # Collect SV metadata for this informative read
                    max_clip = 0
                    if read.cigartuples:
                        for op, length in read.cigartuples:
                            if op == 4 and length > max_clip:
                                max_clip = length
                    read_sv_meta[dedup_key] = {
                        "has_sa": read.has_tag("SA"),
                        "sa_str": read.get_tag("SA") if (
                            read.has_tag("SA") and not read.is_supplementary
                        ) else None,
                        "is_paired": read.is_paired,
                        "is_proper_pair": read.is_proper_pair,
                        "mate_is_unmapped": (
                            read.mate_is_unmapped if read.is_paired else False
                        ),
                        "max_clip": max_clip,
                    }

                if total_reads_scanned % 1_000_000 == 0:
                    logger.info(
                        "  Anchoring: %d reads scanned, %d informative (%s)",
                        total_reads_scanned,
                        len(read_hits) + unmapped_informative,
                        _format_elapsed(time.monotonic() - anchor_start),
                    )

        bam.close()
        if jf_query_st is not None:
            jf_query_st.close()

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
            "\tmax_clip_len\tunmapped_mates\tclass\n"
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
            fh.write(
                f"{chrom}\t{start}\t{end}\t{n_reads}\t{n_kmers}"
                f"\t{split_reads}\t{discordant_pairs}"
                f"\t{max_clip_len}\t{unmapped_mates}\t{region_class}\n"
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


def _annotate_and_link_from_metadata(regions, region_reads, read_sv_meta):
    """Annotate regions and link breakpoints using pre-collected metadata.

    Uses per-read SV metadata collected during the anchoring scan
    (Module 3), so no additional BAM I/O is needed.

    Args:
        regions: List of (chrom, start, end) tuples.
        region_reads: Dict mapping region tuple to set of read names.
        read_sv_meta: Dict mapping (query_name, is_supplementary) to
            per-read SV metadata dict with keys: has_sa, sa_str,
            is_paired, is_proper_pair, mate_is_unmapped, max_clip.

    Returns:
        (annotations, links) where:
        - annotations: Dict mapping region tuple to annotation dict with
          keys: split_reads, discordant_pairs, max_clip_len, unmapped_mates.
        - links: List of dicts with keys: region_a, region_b,
          supporting_reads, sv_type_hint.
    """
    # Build lookup from read name to regions it belongs to
    read_to_regions = {}
    for region_key in regions:
        for qname in region_reads.get(region_key, set()):
            read_to_regions.setdefault(qname, set()).add(region_key)

    # Initialize annotations
    annotations = {
        r: {"split_reads": 0, "discordant_pairs": 0,
            "max_clip_len": 0, "unmapped_mates": 0}
        for r in regions
    }

    if not read_to_regions:
        return annotations, []

    # ── Annotation from metadata ──
    # Track which (qname, region) pairs have already been counted for
    # split_reads so each molecule is counted at most once per region,
    # even when both primary and supplementary alignments are informative.
    split_read_counted = set()

    for dedup_key, meta in read_sv_meta.items():
        qname = dedup_key[0]
        if qname not in read_to_regions:
            continue

        for region_key in read_to_regions[qname]:
            ann = annotations[region_key]

            if meta["has_sa"]:
                sr_key = (qname, region_key)
                if sr_key not in split_read_counted:
                    ann["split_reads"] += 1
                    split_read_counted.add(sr_key)

            if meta["is_paired"]:
                if meta["mate_is_unmapped"]:
                    ann["unmapped_mates"] += 1
                elif not meta["is_proper_pair"]:
                    ann["discordant_pairs"] += 1

            if meta["max_clip"] > ann["max_clip_len"]:
                ann["max_clip_len"] = meta["max_clip"]

    # ── Linking via SA tags ──
    # Build interval index per chromosome for O(log n) SA target lookup
    region_by_chrom = {}
    for r in regions:
        region_by_chrom.setdefault(r[0], []).append(r)
    chrom_starts = {}
    chrom_regions_sorted = {}
    for chrom, rlist in region_by_chrom.items():
        rlist.sort(key=lambda x: x[1])
        chrom_starts[chrom] = [r[1] for r in rlist]
        chrom_regions_sorted[chrom] = rlist

    sa_bridges = {}

    for dedup_key, meta in read_sv_meta.items():
        qname = dedup_key[0]
        sa_str = meta.get("sa_str")
        if not sa_str:
            continue
        if qname not in read_to_regions:
            continue

        primary_regions = read_to_regions[qname]

        for sa_entry in sa_str.rstrip(";").split(";"):
            parts = sa_entry.split(",")
            if len(parts) < 3:
                continue
            sa_chrom = parts[0]
            try:
                sa_pos = int(parts[1]) - 1  # 1-based to 0-based
            except ValueError:
                continue

            if sa_chrom not in chrom_starts:
                continue
            starts = chrom_starts[sa_chrom]
            sorted_regions = chrom_regions_sorted[sa_chrom]
            idx = bisect.bisect_right(starts, sa_pos) - 1
            if idx >= 0:
                t_chrom, t_start, t_end = sorted_regions[idx]
                if t_start <= sa_pos < t_end:
                    target_region = (t_chrom, t_start, t_end)
                    for p_region in primary_regions:
                        if p_region != target_region:
                            key = tuple(sorted(
                                [p_region, target_region],
                            ))
                            sa_bridges.setdefault(
                                key, set(),
                            ).add(qname)

    # Also link regions sharing query_names across primary/supplementary
    for qname, rset in read_to_regions.items():
        if len(rset) >= 2:
            rlist = sorted(rset)
            for i in range(len(rlist)):
                for j in range(i + 1, len(rlist)):
                    key = (rlist[i], rlist[j])
                    sa_bridges.setdefault(key, set()).add(qname)

    # Build links list
    links = []
    for region_a, region_b in sorted(sa_bridges):
        qnames = sa_bridges[(region_a, region_b)]
        sv_type = _infer_sv_type(region_a, region_b)
        links.append({
            "region_a": region_a,
            "region_b": region_b,
            "supporting_reads": qnames,
            "sv_type_hint": sv_type,
        })

    return annotations, links


def _write_bedpe(links, bedpe_path):
    """Write linked SV breakpoint pairs to a BEDPE file.

    Args:
        links: List of link dicts from ``_annotate_and_link_from_metadata()``.
        bedpe_path: Output BEDPE file path.
    """
    with open(bedpe_path, "w") as fh:
        fh.write(
            "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2"
            "\tsv_id\tsupporting_reads\tsv_type\n"
        )
        for idx, link in enumerate(links, 1):
            ra = link["region_a"]
            rb = link["region_b"]
            n_support = len(link["supporting_reads"])
            sv_type = link["sv_type_hint"]
            fh.write(
                f"{ra[0]}\t{ra[1]}\t{ra[2]}"
                f"\t{rb[0]}\t{rb[1]}\t{rb[2]}"
                f"\tSV_{idx}\t{n_support}\t{sv_type}\n"
            )
    logger.info("BEDPE file written: %s (%d links)", bedpe_path, len(links))


def _classify_regions(regions, region_annotations, sv_links):
    """Assign SV classification to each region.

    - ``SV``: split_reads >= 2 OR discordant_pairs >= 2 OR
      unmapped_mates >= 2 OR region is linked via sv_links
    - ``SMALL``: split_reads == 0 AND discordant_pairs == 0 AND
      unmapped_mates == 0 AND not linked
    - ``AMBIGUOUS``: otherwise

    Updates region_annotations in place with a ``class`` key.
    """
    linked_regions = set()
    for link in sv_links:
        linked_regions.add(link["region_a"])
        linked_regions.add(link["region_b"])

    for region_key in regions:
        ann = region_annotations.get(region_key, {})
        split_reads = ann.get("split_reads", 0)
        discordant_pairs = ann.get("discordant_pairs", 0)
        unmapped_mates = ann.get("unmapped_mates", 0)
        if (split_reads >= 2 or discordant_pairs >= 2
                or unmapped_mates >= 2
                or region_key in linked_regions):
            ann["class"] = "SV"
        elif split_reads == 0 and discordant_pairs == 0 and unmapped_mates == 0:
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
                # e.g. "chr11:55003995" "T>C" "24" "53" "24" "0.4528" ...
                variant = parts[0]  # chr:pos
                ref_alt = parts[1]  # R>A
                dku = int(parts[2])
                dkt = int(parts[3])
                dka = int(parts[4])
                dka_dkt = float(parts[5])
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


# ── Curated DNM region definitions (Sulovari et al. 2023) ──────────

#: Curated de novo mutation regions from Sulovari et al. 2023
#: (PMID: 36894594, PMC10006329).  Each tuple:
#: (chrom, position, size_bp_or_None, event_type)
SULOVARI_DNM_REGIONS = [
    ("chr17", 53340465, 107, "deletion"),
    ("chr14", 23280711, None, "microsatellite_expansion"),
    ("chr3", 85552367, 64, "sv_like"),
    ("chr5", 97089276, 43, "sv_like"),
    ("chr8", 125785998, 43, "sv_like"),
    ("chr18", 62805217, 34, "sv_like"),
    ("chr7", 142786222, 10607, "deletion"),
]


def _evaluate_dnm_regions(discovery_regions, region_detail,
                          dnm_regions=None):
    """Evaluate how well VCF-free discovery captures curated DNM regions.

    For each curated de novo mutation region from Sulovari et al. 2023,
    determines whether it was nominated by the discovery pipeline and
    collects quantitative k-mer and SV-signal evidence from the
    overlapping discovery region(s).

    The evaluation provides a simple genotype-like assessment per region:

    - **DETECTED**: ≥1 discovery region overlaps the curated locus with
      informative reads carrying proband-unique k-mers.
    - **NOT_DETECTED**: No overlapping discovery region found.

    For detected regions, a *k-mer signal score* summarises evidence
    strength as ``unique_kmers / region_size_bp``.  Higher density
    indicates more child-specific sequence variation concentrated in the
    region, consistent with a de novo SV.

    Args:
        discovery_regions: List of (chrom, start, end) tuples (0-based,
            half-open) from the discovery BED.
        region_detail: List of dicts with per-region metrics (the
            ``regions`` array from the discovery metrics JSON).
        dnm_regions: Optional list of (chrom, pos, size_or_None,
            event_type) tuples.  Defaults to ``SULOVARI_DNM_REGIONS``.

    Returns:
        List of dicts, one per curated region, with keys:

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
    if dnm_regions is None:
        dnm_regions = SULOVARI_DNM_REGIONS

    # Build an index from region tuple → detail dict
    detail_by_key = {}
    for rd in region_detail:
        key = (rd["chrom"], rd["start"], rd["end"])
        detail_by_key[key] = rd

    results = []
    for chrom, pos, size, event_type in dnm_regions:
        dnm_start = pos
        dnm_end = pos + (size if size else 1)  # point if no size

        # Find overlapping discovery regions
        matches = []
        for dr_key in discovery_regions:
            dr_chrom, dr_start, dr_end = dr_key
            if dr_chrom != chrom:
                continue
            # Overlap check (both 0-based half-open)
            if dr_start < dnm_end and dnm_start < dr_end:
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
        lines.append("-" * 120)
        lines.append(
            f"  {'Region':<35s} {'Size':>8s} {'Reads':>6s}"
            f" {'Unique K-mers':>14s}"
            f" {'Split':>6s} {'Disc':>5s} {'MaxClip':>8s}"
            f" {'UnmapMate':>10s} {'Class':>10s}"
        )
        lines.append(
            f"  {'------':<35s} {'----':>8s} {'-----':>6s}"
            f" {'-------------':>14s}"
            f" {'-----':>6s} {'----':>5s} {'-------':>8s}"
            f" {'---------':>10s} {'-----':>10s}"
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
                f" {ann.get('class', 'SMALL'):>10s}"
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
        lines.append(
            "Curated DNM Region Evaluation (Sulovari et al. 2023)"
        )
        lines.append("-" * 80)
        lines.append(f"  Curated DNM loci:            {n_total:>8}")
        lines.append(
            f"  Detected by discovery:       {n_detected:>8}"
            f" / {n_total} ({pct:.1f}%)"
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


def _write_informative_reads_discovery(
    child_bam, ref_fasta, proband_unique_kmers_or_path, kmer_size,
    output_bam,
):
    """Write child reads carrying proband-unique k-mers to a BAM file.

    Includes **all** primary, non-duplicate reads — mapped and unmapped,
    regardless of mapping quality — so that downstream re-alignment can
    rescue initially unmapped reads.

    Each output read is tagged with ``dk:i:1`` indicating it contains
    a proband-unique k-mer. Reads are sorted and indexed for IGV.

    Supports two scanning backends:
    - **Jellyfish** — when *proband_unique_kmers_or_path* is a ``.jf``
      file path; uses disk-backed queries.
    - **Aho-Corasick** — when a FASTA path or Python set is provided;
      builds an in-memory automaton.

    Args:
        child_bam: Path to the child BAM file.
        ref_fasta: Path to the reference FASTA.
        proband_unique_kmers_or_path: Set of canonical k-mer strings,
            a path to a FASTA file, or a path to a ``.jf`` index.
        kmer_size: K-mer size.
        output_bam: Path for the output BAM file.
    """
    _log_memory("before informative reads k-mer load")

    automaton = None
    jf_query = None

    if isinstance(proband_unique_kmers_or_path, str):
        if proband_unique_kmers_or_path.endswith(".jf"):
            jf_query = JellyfishKmerQuery(proband_unique_kmers_or_path)
            logger.info(
                "  Using jellyfish query for BAM writing: %s",
                proband_unique_kmers_or_path,
            )
        else:
            kmers = _load_kmers_from_fasta(proband_unique_kmers_or_path)
            logger.info(
                "  Loaded %d proband-unique k-mers from %s for BAM writing",
                len(kmers), proband_unique_kmers_or_path,
            )
            automaton = build_kmer_automaton(kmers)
            del kmers
    else:
        kmers = proband_unique_kmers_or_path or set()
        automaton = build_kmer_automaton(kmers)
        del kmers

    _log_memory("after informative reads scanner init")

    bam_in = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )

    unsorted_path = output_bam + ".unsorted.bam"
    bam_out = pysam.AlignmentFile(unsorted_path, "wb", header=bam_in.header)

    written = set()
    for read in bam_in.fetch():
        if read.is_secondary:
            continue
        if read.is_duplicate:
            continue

        seq = read.query_sequence
        if seq is None:
            continue

        has_unique = False
        if automaton is not None:
            for _end_idx, _canonical_kmer in automaton.iter(seq):
                has_unique = True
                break
        elif jf_query is not None:
            unique_in_read, _ = jf_query.scan_read(seq, kmer_size)
            has_unique = bool(unique_in_read)

        dedup_key = (read.query_name, read.is_supplementary)
        if has_unique and dedup_key not in written:
            read.set_tag("dk", 1, value_type="i")
            bam_out.write(read)
            written.add(dedup_key)

    bam_out.close()
    bam_in.close()
    if jf_query is not None:
        jf_query.close()


    pysam.sort("-o", output_bam, unsorted_path)
    pysam.index(output_bam)
    os.remove(unsorted_path)

    logger.info(
        "Informative reads BAM written: %s (%d reads)",
        output_bam, len(written),
    )


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
    logger.info("  Min base quality:  %d", args.min_baseq)
    logger.info("  Min distinct kmers/read: %d", min_dk_per_read)
    logger.info("  JF hash size:      %s", jf_hash_size or "(auto)")
    logger.info("  Threads:           %d", args.threads)
    logger.info(
        "  Memory limit:      %s",
        f"{memory_limit_gb:.1f} GB" if memory_limit_gb is not None
        else "(auto-detect)",
    )
    logger.info("  Tmp dir:           %s", getattr(args, 'tmp_dir', None) or "(auto)")
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
    tmp_root = _resolve_tmp_dir(args, out_dir)
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
            )
        )
        logger.info(
            "[Module 3] Complete (%s)",
            _format_elapsed(time.monotonic() - step_start),
        )
        _log_memory("after Module 3")

        # Write informative reads BAM using the jellyfish index.
        logger.info("[Module 4] Writing informative reads BAM: %s", info_bam_path)
        _write_informative_reads_discovery(
            args.child, args.ref_fasta, proband_jf,
            args.kmer_size, info_bam_path,
        )

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

    # Informative reads BAM was already written inside the tmpdir
    # context (before tmpdir cleanup) to avoid copying the FASTA.
    _log_memory("after informative reads BAM")

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

    # ── Curated DNM region evaluation ─────────────────────────────
    dnm_evaluation = _evaluate_dnm_regions(
        regions, metrics["regions"],
    )
    n_dnm_detected = sum(1 for e in dnm_evaluation if e["detected"])
    logger.info(
        "[Module 4] Curated DNM evaluation: %d / %d detected",
        n_dnm_detected, len(dnm_evaluation),
    )
    metrics["dnm_evaluation"] = {
        "total_loci": len(dnm_evaluation),
        "detected": n_dnm_detected,
        "detection_rate": (
            n_dnm_detected / len(dnm_evaluation)
        ) if dnm_evaluation else 0.0,
        "loci": dnm_evaluation,
    }

    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    logger.info("[Module 4] Metrics written to: %s", metrics_path)

    logger.info("[Module 4] Writing summary: %s", summary_path)
    summary_text = _write_discovery_summary(
        summary_path, regions, region_reads, region_kmers, metrics,
        candidate_comparison=candidate_comparison,
        region_annotations=region_annotations,
        dnm_evaluation=dnm_evaluation,
    )
    logger.info("\n%s", summary_text)

    logger.info(
        "[Module 4] Output complete (%s)",
        _format_elapsed(time.monotonic() - step_start),
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
    tmp_root = _resolve_tmp_dir(args, out_dir)
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

        dkt = len(read_kmers_list)
        running_reads += dkt
        dku = 0
        dka = 0
        informative_names = set()
        informative_alt_names = set()
        all_variant_kmers = set()
        alt_variant_kmers = set()
        for read_name, kmers, supports_alt in read_kmers_list:
            all_variant_kmers.update(kmers)
            if supports_alt:
                alt_variant_kmers.update(kmers)
            # A read is informative if it has at least one variant-spanning
            # k-mer that is absent from both parents.
            if not kmers.issubset(parent_kmer_set):
                dku += 1
                informative_names.add(read_name)
                if supports_alt:
                    dka += 1
                    informative_alt_names.add(read_name)

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

    # ── Kraken2 bacterial content flagging (VCF mode) ──────────────
    kraken2_result = None
    if kraken2_db is not None:
        step_start = time.monotonic()
        all_informative_names = set()
        for names in informative_reads_by_variant.values():
            all_informative_names.update(names)
        logger.info(
            "[Kraken2] Classifying %d informative reads for bacterial content",
            len(all_informative_names),
        )
        kraken2_result = _run_kraken2_on_reads(
            args.child, args.ref_fasta, all_informative_names,
            kraken2_db, confidence=kraken2_confidence,
            threads=args.threads,
            informative_reads_by_variant=informative_reads_by_variant,
        )
        logger.info(
            "[Kraken2] %s (%s)",
            kraken2_result.summary(),
            _format_elapsed(time.monotonic() - step_start),
        )
        bacterial_read_names = kraken2_result.bacterial_read_names
        for var_key, ann in annotations.items():
            dku_names = informative_reads_by_variant.get(var_key, set())
            dka_names = informative_alt_reads_by_variant.get(var_key, set())

            dku_bacterial = len(dku_names.intersection(bacterial_read_names))
            dka_bacterial = len(dka_names.intersection(bacterial_read_names))

            ann["dku_bacterial_fraction"] = (
                round(dku_bacterial / len(dku_names), _BACTERIAL_FRACTION_PRECISION)
                if dku_names else 0.0
            )
            ann["dka_bacterial_fraction"] = (
                round(dka_bacterial / len(dka_names), _BACTERIAL_FRACTION_PRECISION)
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
                "human_reads": kraken2_result.human_count,
                "root_reads": kraken2_result.root_count,
                "bacterial_fraction": kraken2_result.bacterial_fraction,
            }
        with open(args.metrics, "w") as fh:
            json.dump(metrics, fh, indent=2)
        logger.info("[Step 5/5] Metrics written to: %s", args.metrics)

    if args.summary:
        logger.info("[Step 5/5] Writing summary: %s", args.summary)
        summary_text = _write_summary(args.summary, variants, annotations)
        logger.info("\n%s", summary_text)

    logger.info(
        "[Step 5/5] Output complete (%s)",
        _format_elapsed(time.monotonic() - step_start),
    )

    # ── Done ───────────────────────────────────────────────────────
    logger.info(
        "Pipeline finished successfully in %s",
        _format_elapsed(time.monotonic() - pipeline_start),
    )
