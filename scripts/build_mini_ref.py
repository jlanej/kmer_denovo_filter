#!/usr/bin/env python3
"""Build a minimal reference FASTA from BAM reads with no mismatches.

This script scans one or more BAM files and reconstructs reference
sequence from reads whose CIGAR string contains only matches ('M')
and no mismatches (NM:i:0 tag).  It produces a FASTA + .fai index
suitable for use as a mini reference in integration tests.

The resulting FASTA contains one sequence per chromosome with reads,
where each base is either reconstructed from perfect-match reads or
filled with 'N' where no coverage exists.

Usage
-----
    python scripts/build_mini_ref.py \
        -b tests/data/giab/HG002_child.bam \
        -b tests/data/giab/HG003_father.bam \
        -b tests/data/giab/HG004_mother.bam \
        -o tests/data/giab/mini_ref.fa
"""

import argparse
import collections
import os
import sys

import pysam


def _reads_with_no_mismatches(bam_path):
    """Yield (chrom, pos, seq) for reads with NM==0.

    Only primary, non-duplicate, mapped reads with the NM:i:0 tag
    and a CIGAR consisting entirely of M (BAM_CMATCH) operations are
    yielded.
    """
    bam = pysam.AlignmentFile(bam_path)
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_duplicate:
            continue
        if read.query_sequence is None:
            continue

        # Check NM tag — must be 0 (no mismatches)
        try:
            nm = read.get_tag("NM")
        except KeyError:
            continue
        if nm != 0:
            continue

        # CIGAR must be all M (match/mismatch) — skip if it has
        # insertions, deletions, soft clips, etc.
        cigar = read.cigartuples
        if cigar is None:
            continue
        if not all(op == 0 for op, _ in cigar):  # 0 = BAM_CMATCH
            continue

        yield (
            read.reference_name,
            read.reference_start,  # 0-based
            read.query_sequence,
        )
    bam.close()


def _cluster_intervals(positions, max_gap=1000):
    """Cluster sorted positions into contiguous intervals.

    Adjacent positions separated by more than *max_gap* bases start a
    new interval.

    Returns list of (start, end) tuples (0-based, half-open).
    """
    if not positions:
        return []

    sorted_pos = sorted(positions)
    intervals = []
    start = sorted_pos[0]
    prev = sorted_pos[0]

    for pos in sorted_pos[1:]:
        if pos - prev > max_gap:
            intervals.append((start, prev + 1))
            start = pos
        prev = pos

    intervals.append((start, prev + 1))
    return intervals


def build_mini_ref(bam_paths, output_fasta, padding=100):
    """Build a minimal reference FASTA from perfect-match BAM reads.

    For each covered region, creates a separate contig named
    ``{chrom}_{start}_{end}`` containing only that region's sequence
    (with small gaps filled by 'N').  This avoids multi-megabyte
    sequences for chromosomes with reads at distant loci.

    Args:
        bam_paths: List of BAM file paths to scan.
        output_fasta: Output FASTA path.
        padding: Extra bases of padding around each region (filled
            with 'N' if not covered).
    """
    # Collect per-chromosome base coverage
    chrom_bases = collections.defaultdict(dict)
    total_reads = 0

    for bam_path in bam_paths:
        print(f"Scanning {bam_path}...", file=sys.stderr)
        for chrom, start, seq in _reads_with_no_mismatches(bam_path):
            total_reads += 1
            bases = chrom_bases[chrom]
            for i, base in enumerate(seq):
                pos = start + i
                if pos not in bases:
                    bases[pos] = base

    print(f"Total perfect-match reads used: {total_reads}", file=sys.stderr)
    print(f"Chromosomes: {sorted(chrom_bases.keys())}", file=sys.stderr)

    # Write FASTA — one contig per clustered region
    chroms = sorted(chrom_bases.keys())
    with open(output_fasta, "w") as fh:
        for chrom in chroms:
            bases = chrom_bases[chrom]
            intervals = _cluster_intervals(list(bases.keys()))

            for iv_start, iv_end in intervals:
                # Add padding
                padded_start = max(0, iv_start - padding)
                padded_end = iv_end + padding

                # Build sequence
                seq_chars = []
                for pos in range(padded_start, padded_end):
                    seq_chars.append(bases.get(pos, "N"))
                seq = "".join(seq_chars)

                contig = f"{chrom}_{padded_start}_{padded_end}"
                fh.write(f">{contig}\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")

                covered = sum(1 for c in seq if c != "N")
                print(
                    f"  {contig}: {len(seq)} bp "
                    f"({covered} covered, {len(seq) - covered} gaps)",
                    file=sys.stderr,
                )

    # Index with samtools
    pysam.faidx(output_fasta)
    print(f"Written: {output_fasta} (+.fai)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Build a minimal reference FASTA from BAM reads "
                    "with no mismatches (NM:i:0, CIGAR all-M).",
    )
    parser.add_argument(
        "-b", "--bam", action="append", required=True,
        help="BAM file(s) to scan (can be specified multiple times)",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output FASTA path",
    )
    args = parser.parse_args()

    build_mini_ref(args.bam, args.output)


if __name__ == "__main__":
    main()
