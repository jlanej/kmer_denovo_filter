"""Main pipeline for de novo variant k-mer analysis."""

import collections
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
    _is_symbolic,
    canonicalize,
    extract_variant_spanning_kmers,
    read_supports_alt,
)

logger = logging.getLogger(__name__)


def _check_tool(name):
    """Check if an external tool is available on PATH."""
    return shutil.which(name) is not None


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


def _validate_inputs(args):
    """Validate pipeline inputs before starting computation.

    Raises SystemExit with a clear error message for any invalid input.
    """
    errors = []

    # Check required input files exist
    for label, path in [
        ("Child BAM/CRAM (--child)", args.child),
        ("Mother BAM/CRAM (--mother)", args.mother),
        ("Father BAM/CRAM (--father)", args.father),
        ("Input VCF (--vcf)", args.vcf),
    ]:
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
    if args.min_mapq < 0:
        errors.append(
            f"--min-mapq must be >= 0, got {args.min_mapq}"
        )
    if args.threads < 1:
        errors.append(
            f"--threads must be >= 1, got {args.threads}"
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
        var_key = f"{chrom}:{pos}"
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
                supports = read_supports_alt(read, pos, ref, alt,
                                             aligned_pairs=aligned_pairs, seq=seq)
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


def _write_kmer_fasta(kmers, filepath):
    """Write k-mers to a FASTA file for jellyfish ``--if``."""
    with open(filepath, "w") as fh:
        for i, kmer in enumerate(kmers):
            fh.write(f">{i}\n{kmer}\n")


def _scan_parent_jellyfish(
    parent_bam, ref_fasta, kmer_fasta, kmer_size, parent_dir, threads=4,
):
    """Scan a parent BAM and find which child k-mers are present.

    Uses ``jellyfish count`` with the ``--if`` filter so only child k-mers
    are tracked while the parent BAM is streamed.

    Returns:
        Set of canonical k-mers found in the parent.
    """
    os.makedirs(parent_dir, exist_ok=True)
    jf_output = os.path.join(parent_dir, "parent.jf")

    samtools_cmd = ["samtools", "fasta", "-F", "0x500", "-@", "2", parent_bam]
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
        result = subprocess.run(
            dump_cmd, capture_output=True, text=True, check=True,
        )
        for line in result.stdout.strip().split("\n"):
            if line:
                parts = line.split()
                found_kmers[parts[0]] = int(parts[1])

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

    When *proband_id* matches a sample in the VCF, DKU and DKT are written
    as FORMAT fields on that sample.  Otherwise they are written as INFO
    fields.

    Returns:
        The actual output path (with ``.gz`` suffix).
    """
    vcf_in = pysam.VariantFile(input_vcf)

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

    if not output_vcf.endswith(".gz"):
        output_vcf = output_vcf + ".gz"

    vcf_out = pysam.VariantFile(output_vcf, "wz", header=vcf_in.header)

    for rec in vcf_in:
        var_key = f"{rec.chrom}:{rec.start}"
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
        chrom, pos_str = var_key.rsplit(":", 1)
        pos = int(pos_str)
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
        var_key = f"{var['chrom']}:{var['pos']}"
        ann = annotations.get(var_key, {"dku": 0, "dkt": 0, "dka": 0, "dku_dkt": 0.0, "dka_dkt": 0.0, "max_pkc": 0, "avg_pkc": 0.0, "min_pkc": 0, "max_pkc_alt": 0, "avg_pkc_alt": 0.0, "min_pkc_alt": 0})
        ref = var["ref"]
        alts = var["alts"]
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
    logger.info("  Proband ID:        %s", args.proband_id or "(not set)")
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

    with tempfile.TemporaryDirectory(prefix="kmer_denovo_") as tmpdir:
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
        var_key = f"{var['chrom']}:{var['pos']}"
        read_kmers_list = variant_read_kmers.get(var_key, [])

        dkt = len(read_kmers_list)
        running_reads += dkt
        dku = 0
        dka = 0
        informative_names = set()
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
