"""Main pipeline for de novo variant k-mer analysis."""

import json
import logging
import os
import shutil
import statistics
import subprocess
import sys
import tempfile

import pysam

from kmer_denovo_filter.kmer_utils import (
    _is_symbolic,
    canonicalize,
    extract_variant_spanning_kmers,
)

logger = logging.getLogger(__name__)


def _check_tool(name):
    """Check if an external tool is available on PATH."""
    return shutil.which(name) is not None


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
            (read_name, kmer_set) tuples
    """
    bam = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )
    batch = set()
    total_written = 0
    fasta_fh = open(kmer_fasta, "w")
    variant_read_kmers = {}

    def _flush_batch():
        nonlocal total_written
        for kmer in batch:
            fasta_fh.write(f">{total_written}\n{kmer}\n")
            total_written += 1
        batch.clear()

    for var in variants:
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

            kmers = extract_variant_spanning_kmers(
                read, pos, kmer_size, min_baseq, ref=ref, alt=alt,
            )
            if kmers:
                read_kmers.append((read.query_name, kmers))
                batch.update(kmers)
                if len(batch) >= flush_threshold:
                    _flush_batch()

        variant_read_kmers[var_key] = read_kmers

        if debug_kmers:
            unique = (
                set().union(*(k for _, k in read_kmers)) if read_kmers
                else set()
            )
            logger.info(
                "Variant %s: %d reads, %d unique k-mers",
                var_key, len(read_kmers), len(unique),
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

    samtools_cmd = ["samtools", "fasta", "-F", "0x400", "-@", "2", parent_bam]
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

    logger.info("Scanning parent: %s", parent_bam)

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

    _, jf_stderr = p_jellyfish.communicate()
    p_samtools.communicate()

    if p_jellyfish.returncode != 0:
        raise RuntimeError(
            f"jellyfish count failed: {jf_stderr.decode()}"
        )

    found_kmers = set()
    if os.path.exists(jf_output):
        dump_cmd = ["jellyfish", "dump", "-c", "-L", "1", jf_output]
        result = subprocess.run(
            dump_cmd, capture_output=True, text=True, check=True,
        )
        for line in result.stdout.strip().split("\n"):
            if line:
                kmer = line.split()[0]
                found_kmers.add(kmer)

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
            else:
                rec.info["DKU"] = ann["dku"]
                rec.info["DKT"] = ann["dkt"]
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
        median_dku = statistics.median(dku_values)
        lines.append("Read Support Statistics")
        lines.append("-" * 40)
        lines.append(f"  DKU  mean:   {mean_dku:>6.1f}   median: {median_dku:>4}")
        lines.append(f"  DKT  mean:   {mean_dkt:>6.1f}")
        lines.append("")

    if dnm_dku:
        mean_dnm_dku = sum(dnm_dku) / len(dnm_dku)
        lines.append(f"  Avg DKU among likely DNMs:   {mean_dnm_dku:>6.1f}")
        lines.append("")

    lines.append("Per-Variant Results")
    lines.append("-" * 60)
    lines.append(f"  {'Variant':<30s} {'DKU':>5s} {'DKT':>5s}  Call")
    lines.append(f"  {'-------':<30s} {'---':>5s} {'---':>5s}  ----")

    for var in variants:
        var_key = f"{var['chrom']}:{var['pos']}"
        ann = annotations.get(var_key, {"dku": 0, "dkt": 0})
        ref = var["ref"]
        alts = var["alts"]
        alt = alts[0] if alts else "."
        label = f"{var['chrom']}:{var['pos'] + 1} {ref}>{alt}"
        call = "DE_NOVO" if ann["dku"] > 0 else "inherited"
        lines.append(f"  {label:<30s} {ann['dku']:>5d} {ann['dkt']:>5d}  {call}")

    lines.append("")
    lines.append("=" * 60)
    lines.append("")

    text = "\n".join(lines)

    with open(summary_path, "w") as fh:
        fh.write(text)

    return text


def run_pipeline(args):
    """Run the de novo k-mer analysis pipeline."""
    logging.basicConfig(
        level=logging.DEBUG if args.debug_kmers else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    for tool in ("samtools", "jellyfish"):
        if not _check_tool(tool):
            logger.error("%s not found in PATH", tool)
            sys.exit(1)

    # 1. Parse VCF
    logger.info("Parsing VCF: %s", args.vcf)
    variants = _parse_vcf_variants(args.vcf)
    logger.info("Found %d candidate variants", len(variants))

    if not variants:
        logger.warning("No variants found in VCF")
        _write_annotated_vcf(args.vcf, args.output, {}, args.proband_id)
        if args.metrics:
            with open(args.metrics, "w") as fh:
                json.dump({"total_variants": 0}, fh, indent=2)
        return

    # 2. Extract child k-mers
    logger.info("Extracting child k-mers (k=%d)", args.kmer_size)

    # 3. Scan parents using jellyfish
    parent_found_kmers = set()

    with tempfile.TemporaryDirectory(prefix="kmer_denovo_") as tmpdir:
        kmer_fasta = os.path.join(tmpdir, "child_kmers.fa")

        total_child_kmers, variant_read_kmers = _collect_child_kmers(
            args.child, args.ref_fasta, variants,
            args.kmer_size, args.min_baseq, args.min_mapq, args.debug_kmers,
            kmer_fasta,
        )
        logger.info(
            "Wrote %d child k-mers (partially deduplicated)",
            total_child_kmers,
        )

        # Scan mother
        logger.info("Scanning mother BAM: %s", args.mother)
        mother_kmers = _scan_parent_jellyfish(
            args.mother, args.ref_fasta, kmer_fasta, args.kmer_size,
            os.path.join(tmpdir, "mother"), args.threads,
        )
        parent_found_kmers.update(mother_kmers)
        logger.info("Found %d child k-mers in mother", len(mother_kmers))

        # Scan father
        logger.info("Scanning father BAM: %s", args.father)
        father_kmers = _scan_parent_jellyfish(
            args.father, args.ref_fasta, kmer_fasta, args.kmer_size,
            os.path.join(tmpdir, "father"), args.threads,
        )
        parent_found_kmers.update(father_kmers)
        logger.info("Found %d child k-mers in father", len(father_kmers))

    child_unique_kmers = max(0, total_child_kmers - len(parent_found_kmers))
    logger.info(
        "Child-unique k-mers (approx): %d / %d (batched dedup)",
        child_unique_kmers,
        total_child_kmers,
    )

    # 4. Count proband-unique reads per variant
    annotations = {}
    informative_reads_by_variant = {}
    for var in variants:
        var_key = f"{var['chrom']}:{var['pos']}"
        read_kmers_list = variant_read_kmers.get(var_key, [])

        dkt = len(read_kmers_list)
        dku = 0
        informative_names = set()
        for read_name, kmers in read_kmers_list:
            # A read is informative if it has at least one variant-spanning
            # k-mer that is absent from both parents.
            if kmers - parent_found_kmers:
                dku += 1
                informative_names.add(read_name)

        annotations[var_key] = {"dku": dku, "dkt": dkt}
        if informative_names:
            informative_reads_by_variant[var_key] = informative_names

        if args.debug_kmers:
            logger.info("Variant %s: DKU=%d DKT=%d", var_key, dku, dkt)

    # 5. Write annotated VCF
    logger.info("Writing annotated VCF: %s", args.output)
    actual_output = _write_annotated_vcf(
        args.vcf, args.output, annotations, args.proband_id,
    )

    # 6. Write informative reads BAM
    if args.informative_reads:
        logger.info(
            "Writing informative reads BAM: %s", args.informative_reads,
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
            "Wrote %d informative reads across %d variants",
            total_reads, len(informative_reads_by_variant),
        )

    # 7. Write metrics
    if args.metrics:
        metrics = {
            "total_variants": len(variants),
            "total_child_kmers": total_child_kmers,
            "parent_found_kmers": len(parent_found_kmers),
            "child_unique_kmers": child_unique_kmers,
            "variants_with_unique_reads": sum(
                1 for a in annotations.values() if a["dku"] > 0
            ),
        }
        with open(args.metrics, "w") as fh:
            json.dump(metrics, fh, indent=2)
        logger.info("Metrics written to: %s", args.metrics)

    # 8. Write summary
    if args.summary:
        logger.info("Writing summary: %s", args.summary)
        summary_text = _write_summary(args.summary, variants, annotations)
        logger.info("\n%s", summary_text)

    logger.info("Done")
