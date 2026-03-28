"""Shared test helper functions for creating synthetic BAM, VCF, and FASTA data."""

import pysam


def _create_ref_fasta(path, chrom="chr1", length=200):
    """Create a small reference FASTA with non-repetitive sequence.

    Uses MD5 solely to produce a deterministic pseudo-random base sequence
    for test data (not for any cryptographic purpose).
    """
    import hashlib
    seq = []
    bases = "ACGT"
    for i in range(length):
        h = hashlib.md5(str(i).encode()).hexdigest()
        seq.append(bases[int(h, 16) % 4])
    seq = "".join(seq)
    with open(path, "w") as fh:
        fh.write(f">{chrom}\n{seq}\n")
    pysam.faidx(path)
    return seq


def _create_bam(path, ref_fasta, chrom, reads):
    """Create a BAM file from a list of (name, pos, seq, quals) tuples.

    ``pos`` is 0-based.  An optional fifth element may supply a CIGAR
    tuple list; otherwise a simple all-M alignment is used.
    """
    header = pysam.AlignmentHeader.from_references(
        [chrom], [300],
    )
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for entry in reads:
            name, pos, seq, *rest = entry
            quals = rest[0] if rest else None
            cigar = rest[1] if len(rest) > 1 else [(0, len(seq))]
            seg = pysam.AlignedSegment()
            seg.query_name = name
            seg.query_sequence = seq
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = pos
            seg.mapping_quality = 60
            seg.cigar = cigar
            if quals is not None:
                seg.query_qualities = pysam.qualitystring_to_array(quals)
            else:
                seg.query_qualities = pysam.qualitystring_to_array(
                    "I" * len(seq)
                )
            bam.write(seg)
    pysam.sort("-o", path, path)
    pysam.index(path)


def _create_vcf(path, chrom, variants, sample="HG002"):
    """Create a VCF with given variants.

    ``variants`` is a list of (pos_1based, ref, alt) tuples.
    """
    header = pysam.VariantHeader()
    header.add_sample(sample)
    header.add_line(f"##contig=<ID={chrom},length=300>")
    with pysam.VariantFile(path, "w", header=header) as vcf:
        for pos, ref, alt in variants:
            rec = vcf.new_record(
                contig=chrom, start=pos - 1, stop=pos - 1 + len(ref),
                alleles=(ref, alt),
            )
            vcf.write(rec)


def _create_bam_with_supplementary(path, ref_fasta, chroms, chrom_lengths,
                                   reads):
    """Create a BAM with support for supplementary alignments and SA tags.

    ``reads`` is a list of dicts with keys:
        name: query name
        chrom_idx: reference index (0-based)
        pos: 0-based reference start
        seq: query sequence
        cigar: CIGAR tuples list (default all-M)
        flag: SAM flag (default 0)
        sa_tag: SA tag string (optional)
        mapq: mapping quality (default 60)
    """
    header = pysam.AlignmentHeader.from_references(chroms, chrom_lengths)
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for entry in reads:
            seg = pysam.AlignedSegment()
            seg.query_name = entry["name"]
            seg.query_sequence = entry["seq"]
            seg.flag = entry.get("flag", 0)
            seg.reference_id = entry.get("chrom_idx", 0)
            seg.reference_start = entry["pos"]
            seg.mapping_quality = entry.get("mapq", 60)
            seg.cigar = entry.get("cigar", [(0, len(entry["seq"]))])
            seg.query_qualities = pysam.qualitystring_to_array(
                "I" * len(entry["seq"])
            )
            if "sa_tag" in entry:
                seg.set_tag("SA", entry["sa_tag"])
            bam.write(seg)
    pysam.sort("-o", path, path)
    pysam.index(path)
