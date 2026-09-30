#!/usr/bin/env python3
"""Simulate a small trio whose child carries de novo structural variants.

Writes to ``--outdir``:

* ``ref.fa`` (+ ``.fai``): a random reference with two chromosomes,
  chr1 (60 kb) and chr2 (30 kb).
* ``child.bam``, ``mother.bam``, ``father.bam`` (+ ``.bai``): 2x150 bp
  paired-end reads at 30x, coordinate-sorted and indexed.
* ``events.tsv``: the child's de novo events (``EVENTS``).

The parents carry only the reference.  One of the child's two haplotypes
carries every event: a deletion, a tandem duplication, an inversion, a
reciprocal translocation, two novel insertions and an SNV.

No aligner is run.  Each read's alignment is written from where the read
came from, following BWA-MEM's conventions for reads that cross a
breakpoint:

* A read part shorter than ``MIN_ALIGNED`` bp is soft-clipped.
* When two or more parts of a read align to places that are not
  collinear, the longest is the primary alignment (soft-clipped) and the
  others are supplementary alignments (hard-clipped), all with SA tags.
* An indel of at most ``MAX_CIGAR_INDEL`` bp, with at least
  ``MIN_INDEL_FLANK`` bp aligned on both sides, goes in the CIGAR.
* A read with no part of ``MIN_ALIGNED`` bp is unmapped and placed at its
  mate's position.
* A pair is proper when its reads face each other on one chromosome
  (forward-reverse) with an insert size within four standard deviations
  of the mean.

Usage::

    python examples/sv_demo/simulate_sv_trio.py --outdir sv_demo
    kmer-discovery --child sv_demo/child.bam --mother sv_demo/mother.bam \\
        --father sv_demo/father.bam --ref-fasta sv_demo/ref.fa \\
        --out-prefix sv_demo/sv_demo

Requires pysam, which kmer-denovo-filter already depends on.
"""

import argparse
import os
import random

import pysam

READ_LEN = 150
INSERT_MEAN = 350
INSERT_SD = 50
PROPER_MIN = INSERT_MEAN - 4 * INSERT_SD
PROPER_MAX = INSERT_MEAN + 4 * INSERT_SD
MIN_ALIGNED = 30      # bp: shorter read parts are clipped (BWA-MEM -T 30)
MIN_INDEL_FLANK = 20  # bp aligned on both sides of a CIGAR indel
MAX_CIGAR_INDEL = 100  # bp: longer indels split or clip the alignment
MAPQ = 60
BASE_QUALITY = 30

CHROM_SIZES = {"chr1": 60_000, "chr2": 30_000}

#: The child's de novo events: (chrom, 0-based start, end, type, note).
#: Insertions and translocation breakends have start == end.
EVENTS = [
    ("chr1", 10_000, 12_000, "DEL", "2 kb deletion"),
    ("chr1", 25_000, 26_500, "DUP", "1.5 kb tandem duplication"),
    ("chr1", 40_000, 42_000, "INV", "2 kb inversion"),
    ("chr1", 55_000, 55_000, "BND",
     "reciprocal translocation with chr2:15000"),
    ("chr2", 5_000, 5_001, "SNV", "single-base substitution"),
    ("chr2", 10_000, 10_000, "INS", "55 bp novel insertion"),
    ("chr2", 15_000, 15_000, "BND",
     "reciprocal translocation with chr1:55000"),
    ("chr2", 22_000, 22_000, "INS", "250 bp novel insertion"),
]

_CIGAR_CODE = {"M": 0, "I": 1, "D": 2, "S": 4, "H": 5}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def random_seq(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


class Haplotype:
    """A sequence assembled from reference segments and novel sequence.

    *blocks* are ``(chrom, start, end, strand)`` reference segments or
    novel sequence strings.  *ref* maps chromosome names to the sequence
    the reference segments are cut from.
    """

    def __init__(self, name, blocks, ref):
        self.name = name
        self.segments = []  # (hap_start, hap_end, chrom, start, end, strand)
        parts = []
        pos = 0
        for block in blocks:
            if isinstance(block, str):
                seq, segment = block, (None, 0, 0, "+")
            else:
                chrom, start, end, strand = block
                seq = ref[chrom][start:end]
                if strand == "-":
                    seq = revcomp(seq)
                segment = block
            self.segments.append((pos, pos + len(seq)) + segment)
            parts.append(seq)
            pos += len(seq)
        self.seq = "".join(parts)

    def read_parts(self, start, end, reverse):
        """Split the read from [start, end) into parts, in read order.

        Each part is ``(qstart, qend, chrom, rstart, rend, strand)``:
        read bases [qstart, qend) align to reference [rstart, rend) on
        *strand*, or are novel sequence when *chrom* is None.
        """
        parts = []
        for h0, h1, chrom, rstart, rend, strand in self.segments:
            s, e = max(start, h0), min(end, h1)
            if s >= e:
                continue
            if chrom is not None and strand == "+":
                rstart, rend = rstart + (s - h0), rstart + (e - h0)
            elif chrom is not None:
                rstart, rend = rend - (e - h0), rend - (s - h0)
            if reverse:
                qstart, qend = end - e, end - s
                strand = "-" if strand == "+" else "+"
            else:
                qstart, qend = s - start, e - start
            parts.append((qstart, qend, chrom, rstart, rend, strand))
        return sorted(parts)


def child_haplotypes(ref, rng):
    """Return the child's haplotypes: the reference and one with EVENTS."""
    alt = dict(ref)
    snv_pos = 5_000
    alt_base = rng.choice([b for b in "ACGT" if b != ref["chr2"][snv_pos]])
    alt["chr2"] = ref["chr2"][:snv_pos] + alt_base + ref["chr2"][snv_pos + 1:]
    ins55 = random_seq(rng, 55)
    ins250 = random_seq(rng, 250)
    # Derivative chromosomes of the chr1:55000 / chr2:15000 translocation
    der1 = Haplotype("der1", [
        ("chr1", 0, 10_000, "+"),         # deletion of 10,000-12,000
        ("chr1", 12_000, 26_500, "+"),
        ("chr1", 25_000, 40_000, "+"),    # second copy of 25,000-26,500
        ("chr1", 40_000, 42_000, "-"),    # inversion
        ("chr1", 42_000, 55_000, "+"),
        ("chr2", 15_000, 22_000, "+"),
        ins250,
        ("chr2", 22_000, 30_000, "+"),
    ], alt)
    der2 = Haplotype("der2", [
        ("chr2", 0, 10_000, "+"),         # SNV at 5,000
        ins55,
        ("chr2", 10_000, 15_000, "+"),
        ("chr1", 55_000, 60_000, "+"),
    ], alt)
    return [
        Haplotype(f"ref_{chrom}", [(chrom, 0, size, "+")], ref)
        for chrom, size in CHROM_SIZES.items()
    ] + [der1, der2]


def parental_haplotypes(ref):
    return [
        Haplotype(f"ref{copy}_{chrom}", [(chrom, 0, size, "+")], ref)
        for copy in (1, 2)
        for chrom, size in CHROM_SIZES.items()
    ]


def _merge_ops(ops):
    merged = []
    for op, length in ops:
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    return merged


def align(parts):
    """Group read parts into alignments, as an aligner would.

    Returns a list of dicts (primary first) with the read-order query
    interval, CIGAR operations in read order, reference interval and
    strand; empty when the read is unmapped.
    """
    groups = []
    prev = None
    novel = 0
    for part in parts:
        qstart, qend, chrom, rstart, rend, strand = part
        if chrom is None:
            novel += qend - qstart
            continue
        op = None
        if prev is not None and prev[2] == chrom and prev[5] == strand:
            gap = rstart - prev[4] if strand == "+" else prev[3] - rend
            size = gap or novel
            flanks = min(prev[1] - prev[0], qend - qstart)
            if gap >= 0 and not (gap and novel) and (
                    size == 0 or (size <= MAX_CIGAR_INDEL
                                  and flanks >= MIN_INDEL_FLANK)):
                op = ("D", gap) if gap else ("I", novel)
        if op is not None:
            group = groups[-1]
            if op[1]:
                group["ops"].append(op)
            group["ops"].append(("M", qend - qstart))
            group["parts"].append(part)
            group["qend"] = qend
        else:
            groups.append({
                "ops": [("M", qend - qstart)], "parts": [part],
                "qstart": qstart, "qend": qend,
            })
        prev = part
        novel = 0

    alignments = []
    for group in groups:
        aligned = sum(n for op, n in group["ops"] if op == "M")
        if aligned < MIN_ALIGNED:
            continue
        group["ops"] = _merge_ops(group["ops"])
        group["aligned"] = aligned
        group["chrom"] = group["parts"][0][2]
        group["strand"] = group["parts"][0][5]
        group["rstart"] = min(p[3] for p in group["parts"])
        group["rend"] = max(p[4] for p in group["parts"])
        alignments.append(group)
    # The longest alignment is the primary; ties go to the first in the read
    alignments.sort(key=lambda g: -g["aligned"])
    return alignments


def _ref_cigar(group, clip):
    """CIGAR string, in reference orientation, with the given clip op."""
    left, right = group["qstart"], READ_LEN - group["qend"]
    ops = group["ops"]
    if group["strand"] == "-":
        left, right, ops = right, left, ops[::-1]
    return ([(clip, left)] if left else []) + ops + (
        [(clip, right)] if right else [])


def _cigar_str(cigar):
    return "".join(f"{n}{op}" for op, n in cigar)


def _edit_distance(group, read_seq, ref):
    nm = sum(n for op, n in group["ops"] if op in "ID")
    for qstart, qend, chrom, rstart, rend, strand in group["parts"]:
        ref_seq = ref[chrom][rstart:rend]
        if strand == "-":
            ref_seq = revcomp(ref_seq)
        nm += sum(a != b for a, b in zip(read_seq[qstart:qend], ref_seq))
    return nm


def _pair_records(header, ref, name, reads):
    """Build BAM records for one read pair.

    *reads* is ``[(sequence, alignments), ...]`` for read 1 and read 2.
    """
    primaries = [alns[0] if alns else None for _, alns in reads]
    proper, tlen = False, [0, 0]
    p1, p2 = primaries
    if p1 and p2 and p1["chrom"] == p2["chrom"]:
        start = min(p1["rstart"], p2["rstart"])
        end = max(p1["rend"], p2["rend"])
        tlen = [end - start, start - end]
        if p2["rstart"] < p1["rstart"]:
            tlen.reverse()
        if p1["strand"] != p2["strand"]:
            fwd, rev = (p1, p2) if p1["strand"] == "+" else (p2, p1)
            proper = (fwd["rstart"] < rev["rend"]
                      and PROPER_MIN <= end - start <= PROPER_MAX)

    records = []
    for i, (seq, alns) in enumerate(reads):
        mate = primaries[1 - i]
        own = primaries[i]
        flag_base = 0x1 | (0x40 if i == 0 else 0x80) | (0x2 if proper else 0)
        if mate is None:
            flag_base |= 0x8
        elif mate["strand"] == "-":
            flag_base |= 0x20
        quals = [BASE_QUALITY] * READ_LEN

        if own is None:
            seg = pysam.AlignedSegment(header)
            seg.query_name = name
            seg.flag = flag_base | 0x4
            seg.query_sequence = seq
            seg.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(q + 33) for q in quals))
            if mate is not None:
                # Placed at the mapped mate's position
                seg.reference_name = mate["chrom"]
                seg.reference_start = mate["rstart"]
                seg.next_reference_name = mate["chrom"]
                seg.next_reference_start = mate["rstart"]
                seg.set_tag("MC", _cigar_str(_ref_cigar(mate, "S")))
                seg.set_tag("MQ", MAPQ)
            records.append(seg)
            continue

        sa_entries = [
            f"{g['chrom']},{g['rstart'] + 1},{g['strand']},"
            f"{_cigar_str(_ref_cigar(g, 'S'))},{MAPQ},"
            f"{_edit_distance(g, seq, ref)}"
            for g in alns
        ]
        for j, group in enumerate(alns):
            supplementary = j > 0
            cigar = _ref_cigar(group, "H" if supplementary else "S")
            ref_seq = seq if group["strand"] == "+" else revcomp(seq)
            ref_quals = quals
            if supplementary:
                left = cigar[0][1] if cigar[0][0] == "H" else 0
                right = cigar[-1][1] if cigar[-1][0] == "H" else 0
                ref_seq = ref_seq[left:READ_LEN - right]
                ref_quals = quals[left:READ_LEN - right]
            seg = pysam.AlignedSegment(header)
            seg.query_name = name
            seg.flag = (flag_base | (0x10 if group["strand"] == "-" else 0)
                        | (0x800 if supplementary else 0))
            seg.reference_name = group["chrom"]
            seg.reference_start = group["rstart"]
            seg.mapping_quality = MAPQ
            seg.cigartuples = [(_CIGAR_CODE[op], n) for op, n in cigar]
            seg.query_sequence = ref_seq
            seg.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(q + 33) for q in ref_quals))
            if mate is None:
                # An unmapped mate is placed here, so it points back here
                seg.next_reference_name = own["chrom"]
                seg.next_reference_start = own["rstart"]
            else:
                seg.next_reference_name = mate["chrom"]
                seg.next_reference_start = mate["rstart"]
                seg.set_tag("MC", _cigar_str(_ref_cigar(mate, "S")))
                seg.set_tag("MQ", MAPQ)
            if not supplementary:
                seg.template_length = tlen[i]
            seg.set_tag("NM", _edit_distance(group, seq, ref))
            if len(alns) > 1:
                seg.set_tag("SA", "".join(
                    entry + ";" for k, entry in enumerate(sa_entries)
                    if k != j
                ))
            records.append(seg)
    return records


def _add_errors(seq, rate, rng):
    if not rate:
        return seq
    bases = list(seq)
    for i, base in enumerate(bases):
        if rng.random() < rate:
            bases[i] = rng.choice([b for b in "ACGT" if b != base])
    return "".join(bases)


def simulate_sample(path, haplotypes, ref, header, depth, error_rate, rng,
                    sample):
    """Write a sorted, indexed BAM of read pairs from *haplotypes*.

    Each haplotype is sequenced to *depth* / 2, so a sample with two
    copies of each chromosome has *depth* overall.
    """
    records = []
    for hap in haplotypes:
        length = len(hap.seq)
        n_pairs = round(depth / 2 * length / (2 * READ_LEN))
        for n in range(n_pairs):
            insert = min(length, max(
                READ_LEN, round(rng.gauss(INSERT_MEAN, INSERT_SD))))
            start = rng.randrange(length - insert + 1)
            left = (start, start + READ_LEN, False)
            right = (start + insert - READ_LEN, start + insert, True)
            ends = [left, right] if rng.random() < 0.5 else [right, left]
            reads = []
            for s, e, reverse in ends:
                seq = hap.seq[s:e]
                if reverse:
                    seq = revcomp(seq)
                seq = _add_errors(seq, error_rate, rng)
                reads.append((seq, align(hap.read_parts(s, e, reverse))))
            records.extend(_pair_records(
                header, ref, f"{sample}:{hap.name}:{n:05d}", reads))

    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        for record in records:
            bam.write(record)
    pysam.sort("-o", path, unsorted)
    os.remove(unsorted)
    pysam.index(path)
    return len(records)


def write_reference(path, ref):
    with open(path, "w") as fh:
        for chrom, seq in ref.items():
            fh.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--outdir", required=True,
                        help="Output directory (created if missing)")
    parser.add_argument("--seed", type=int, default=1,
                        help="Random seed (default: 1)")
    parser.add_argument("--depth", type=float, default=30,
                        help="Sequencing depth per sample (default: 30)")
    parser.add_argument("--error-rate", type=float, default=0.001,
                        help="Per-base substitution error rate "
                             "(default: 0.001)")
    args = parser.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)
    ref = {chrom: random_seq(rng, size) for chrom, size in CHROM_SIZES.items()}
    write_reference(os.path.join(args.outdir, "ref.fa"), ref)

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": n} for c, n in CHROM_SIZES.items()],
        "PG": [{"ID": "simulate_sv_trio", "PN": "simulate_sv_trio.py"}],
    })
    samples = [
        ("child", child_haplotypes(ref, rng)),
        ("mother", parental_haplotypes(ref)),
        ("father", parental_haplotypes(ref)),
    ]
    for sample, haplotypes in samples:
        path = os.path.join(args.outdir, f"{sample}.bam")
        n = simulate_sample(path, haplotypes, ref, header, args.depth,
                            args.error_rate, rng, sample)
        print(f"{path}: {n} records")

    events_path = os.path.join(args.outdir, "events.tsv")
    with open(events_path, "w") as fh:
        fh.write("#chrom\tstart\tend\ttype\tnote\n")
        for event in EVENTS:
            fh.write("\t".join(map(str, event)) + "\n")
    print(f"{events_path}: {len(EVENTS)} de novo events")


if __name__ == "__main__":
    main()
