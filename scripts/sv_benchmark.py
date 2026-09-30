#!/usr/bin/env python3
"""Benchmark kmer-discovery's SV calls on a trio simulated from real sequence.

Unlike the demo in ``examples/sv_demo`` (random sequence, alignments written
by construction), this cuts real GRCh38 sequence, with its repeats, and
aligns simulated reads with BWA-MEM.

``prepare`` cuts chr20:30.5-34.5 Mb and chr21:30-32 Mb from a GRCh38 FASTA
and places events on the trio's haplotypes:

* de novo SVs in the child, heterozygous: deletions, tandem duplications,
  inversions, novel insertions, Alu insertions with target-site
  duplications, a reciprocal translocation, and a deletion between two Alu
  copies (as from non-allelic homologous recombination);
* de novo SNVs and small indels;
* mosaic de novo SVs, carried by half the cells of one haplotype (25% of
  reads);
* inherited SVs, carried by a parent and transmitted, which must not be
  called.

It simulates 2x150 bp reads at 30x per sample, with random substitution
errors, and aligns them with ``bwa mem``.  Run ``kmer-discovery`` on the
output, then ``score`` to compare its BED and BEDPE with the truth:

    python scripts/sv_benchmark.py prepare --reference GRCh38.fa --outdir bench
    kmer-discovery --child bench/child.bam --mother bench/mother.bam \\
        --father bench/father.bam --ref-fasta bench/ref.fa \\
        --out-prefix bench/disc
    python scripts/sv_benchmark.py score --outdir bench --prefix bench/disc

Requires samtools and bwa on PATH.
"""

import argparse
import collections
import json
import math
import os
import random
import re
import subprocess

SEGMENTS = [("chr20", 30_500_000, 34_500_000),
            ("chr21", 30_000_000, 32_000_000)]
READ_LEN = 150
INSERT_MEAN, INSERT_SD = 400, 80
ERROR_RATE = 0.002
TSD_LEN = 12
#: Alu 5' region, conserved across subfamilies; marks Alu copies
ALU_MOTIF = "GGCCGGGCGCGGTGGCTCA"
#: AluY consensus (without its poly(A) tail)
ALU_Y = (
    "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTC"
    "AGGAGATCGAGACCATCCCGGCTAAAACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGG"
    "GCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGA"
    "GGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCT"
    "CAAAAAAA"
)
#: Regions within this distance of a breakpoint count as calling it
MATCH_BP = 300

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


# -- prepare ------------------------------------------------------------------

def cut_reference(reference, outdir):
    ref = {}
    for chrom, start, end in SEGMENTS:
        out = subprocess.run(
            ["samtools", "faidx", reference, f"{chrom}:{start + 1}-{end}"],
            capture_output=True, text=True, check=True).stdout
        ref[chrom] = "".join(out.splitlines()[1:]).upper()
    path = os.path.join(outdir, "ref.fa")
    write_fasta(path, ref)
    subprocess.run(["samtools", "faidx", path], check=True)
    return ref


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def alu_copies(ref):
    """(chrom, start, strand) of Alu copies found by their 5' motif."""
    hits = []
    for chrom, seq in ref.items():
        for m in re.finditer(ALU_MOTIF, seq):
            hits.append((chrom, m.start(), "+"))
        for m in re.finditer(revcomp(ALU_MOTIF), seq):
            hits.append((chrom, m.end(), "-"))
    return hits


class Placer:
    """Picks event positions with room around them, away from Ns."""

    def __init__(self, ref, rng, spacing=25_000):
        self.ref, self.rng, self.spacing = ref, rng, spacing
        self.taken = collections.defaultdict(list)

    def free(self, chrom, start, end):
        seq = self.ref[chrom]
        if start < 5_000 or end > len(seq) - 5_000:
            return False
        if "N" in seq[start - 2_000:end + 2_000]:
            return False
        return all(end + self.spacing < s or start > e + self.spacing
                   for s, e in self.taken[chrom])

    def take(self, chrom, start, end):
        self.taken[chrom].append((start, end))

    def place(self, size, chrom=None):
        for _ in range(10_000):
            c = chrom or self.rng.choice(sorted(self.ref))
            if c == "chr21" and chrom is None and self.rng.random() < 0.6:
                continue  # keep most events on the larger segment
            start = self.rng.randrange(5_000, len(self.ref[c]) - size - 5_000)
            if self.free(c, start, start + size):
                self.take(c, start, start + size)
                return c, start
        raise RuntimeError("no room for another event")


def make_events(ref, rng):
    """The truth set: a list of event dicts, in reference coordinates."""
    placer = Placer(ref, rng)
    events = []

    def add(category, kind, size, haplotype, **extra):
        chrom, start = placer.place(max(size, 1), extra.pop("chrom", None))
        event = {"id": f"{kind}{len(events) + 1}", "category": category,
                 "type": kind, "chrom": chrom, "start": start,
                 "end": start + (0 if kind == "INS" else size),
                 "size": size, "haplotype": haplotype, **extra}
        events.append(event)
        return event

    for size in (50, 75, 100, 150, 300, 500, 1000, 2000, 5000, 10000):
        add("de_novo", "DEL", size, rng.choice("12"))
    for size in (50, 100, 300, 1000, 3000, 10000):
        add("de_novo", "DUP", size, rng.choice("12"))
    for size in (500, 1000, 3000, 10000):
        add("de_novo", "INV", size, rng.choice("12"))
    for size in (50, 100, 200, 400):
        add("de_novo", "INS", size, rng.choice("12"),
            seq="".join(rng.choice("ACGT") for _ in range(size)))
    # Alu insertions, as mobile-element insertions leave them: a young AluY
    # with a little divergence (so, like a real one, it matches many
    # reference copies well), a poly(A) tail and a target-site duplication
    for _ in range(4):
        alu = "".join(b if rng.random() > 0.02 else rng.choice("ACGT")
                      for b in ALU_Y) + "A" * rng.randrange(15, 30)
        if rng.random() < 0.5:
            alu = revcomp(alu)
        add("de_novo", "INS", len(alu), rng.choice("12"), seq=alu,
            tsd=TSD_LEN, note="AluY insertion")
    # A reciprocal translocation between the two segments
    c20, p = placer.place(1, "chr20")
    c21, q = placer.place(1, "chr21")
    events.append({"id": f"BND{len(events) + 1}", "category": "de_novo",
                   "type": "BND", "chrom": c20, "start": p, "end": p,
                   "chrom2": c21, "pos2": q, "size": 0, "haplotype": "2"})
    # Up to two deletions between Alu copies in the same orientation (the
    # junction lies inside their shared motif, as for homologous
    # recombination)
    by_chrom = collections.defaultdict(list)
    for chrom, pos, strand in alu_copies(ref):
        if strand == "+":
            by_chrom[chrom].append(pos)
    n_nahr = 0
    for chrom, starts in sorted(by_chrom.items()):
        starts.sort()
        for a, b in zip(starts, starts[1:]):
            if 2_000 <= b - a <= 8_000 and placer.free(chrom, a, b):
                placer.take(chrom, a, b)
                events.append({"id": f"DEL{len(events) + 1}",
                               "category": "de_novo_alu_mediated",
                               "type": "DEL", "chrom": chrom, "start": a + 10,
                               "end": b + 10, "size": b - a,
                               "haplotype": "1"})
                n_nahr += 1
                break
        if n_nahr == 2:
            break
    for _ in range(3):
        event = add("de_novo_small", "SNV", 1, rng.choice("12"))
        ref_base = ref[event["chrom"]][event["start"]]
        event["seq"] = rng.choice([b for b in "ACGT" if b != ref_base])
    for size in (2, 8):
        add("de_novo_small", "DEL", size, rng.choice("12"))
    add("de_novo_small", "INS", 6, rng.choice("12"),
        seq="".join(rng.choice("ACGT") for _ in range(6)))
    for kind, size in (("DEL", 1000), ("DUP", 1000)):
        add("mosaic", kind, size, "1m")
    for kind, size, parent in (("DEL", 1000, "m"), ("DEL", 5000, "f"),
                               ("DUP", 2000, "m"), ("INV", 3000, "f")):
        add(f"inherited_{parent}", kind, size, "1" if parent == "m" else "2")
    add("inherited_m", "INS", 200, "1",
        seq="".join(rng.choice("ACGT") for _ in range(200)))
    return events


def apply_edits(seq, events, offset=0):
    """Apply events (reference coordinates, minus *offset*) to *seq*."""
    for e in sorted(events, key=lambda e: e["start"], reverse=True):
        s, t = e["start"] - offset, e["end"] - offset
        kind = e["type"]
        if kind == "DEL":
            seq = seq[:s] + seq[t:]
        elif kind == "DUP":
            seq = seq[:t] + seq[s:t] + seq[t:]
        elif kind == "INV":
            seq = seq[:s] + revcomp(seq[s:t]) + seq[t:]
        elif kind == "INS":
            tsd = e.get("tsd", 0)
            seq = seq[:s + tsd] + e["seq"] + seq[s:]
        elif kind == "SNV":
            seq = seq[:s] + e["seq"] + seq[s + 1:]
    return seq


def haplotype(ref, events):
    """Contig sequences of a haplotype carrying *events*."""
    bnd = [e for e in events if e["type"] == "BND"]
    local = [e for e in events if e["type"] != "BND"]
    by_chrom = collections.defaultdict(list)
    for e in local:
        by_chrom[e["chrom"]].append(e)
    if not bnd:
        return {c: apply_edits(seq, by_chrom[c]) for c, seq in ref.items()}
    (b,) = bnd
    p, q = b["start"], b["pos2"]
    pieces = {}
    for chrom, cut in ((b["chrom"], p), (b["chrom2"], q)):
        left = [e for e in by_chrom[chrom] if e["end"] <= cut]
        right = [e for e in by_chrom[chrom] if e["start"] >= cut]
        pieces[chrom] = (apply_edits(ref[chrom][:cut], left),
                         apply_edits(ref[chrom][cut:], right, cut))
    (a, bb), (c, d) = pieces[b["chrom"]], pieces[b["chrom2"]]
    return {f"der_{b['chrom']}": a + d, f"der_{b['chrom2']}": c + bb}


def simulate_reads(seqs, depth, prefix, rng):
    """Write read pairs from *seqs* at *depth*; returns the FASTQ paths.

    Fragments are drawn uniformly, with normally distributed lengths, and
    read from both ends.  Sequencing errors are random substitutions.
    (wgsim is not used: its errors at a position always give the same
    base, which makes recurrent errors look like variants.)
    """
    fq1, fq2 = prefix + "_1.fq", prefix + "_2.fq"
    qual = "I" * READ_LEN
    log_keep = math.log(1 - ERROR_RATE)
    n = 0
    with open(fq1, "w") as out1, open(fq2, "w") as out2:
        for name, seq in seqs.items():
            n_pairs = round(depth * len(seq) / (2 * READ_LEN))
            for _ in range(n_pairs):
                size = max(READ_LEN, round(rng.gauss(INSERT_MEAN, INSERT_SD)))
                start = rng.randrange(len(seq) - size + 1)
                frag = seq[start:start + size]
                if "N" in frag[:READ_LEN] or "N" in frag[-READ_LEN:]:
                    continue
                reads = [frag[:READ_LEN], revcomp(frag[-READ_LEN:])]
                if rng.random() < 0.5:
                    reads.reverse()
                n += 1
                for out, read in zip((out1, out2), reads):
                    bases = list(read)
                    # Positions of errors, by skipping error-free runs
                    i = int(math.log(1 - rng.random()) / log_keep)
                    while i < READ_LEN:
                        bases[i] = rng.choice([b for b in "ACGT"
                                               if b != bases[i]])
                        i += 1 + int(math.log(1 - rng.random()) / log_keep)
                    out.write(f"@{name}_{n}\n{''.join(bases)}\n+\n{qual}\n")
    return fq1, fq2


def align(ref_fa, fastqs, sample, bam, threads):
    fq1 = bam + "_1.fq"
    fq2 = bam + "_2.fq"
    for out, idx in ((fq1, 0), (fq2, 1)):
        with open(out, "w") as fh:
            for h, pair in enumerate(fastqs):
                with open(pair[idx]) as src:
                    # Reads are numbered per haplotype: keep names unique
                    for n, line in enumerate(src):
                        fh.write(f"@h{h}_{line[1:]}" if n % 4 == 0 else line)
    rg = f"@RG\\tID:{sample}\\tSM:{sample}"
    bwa = subprocess.Popen(
        ["bwa", "mem", "-t", str(threads), "-R", rg, ref_fa, fq1, fq2],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "sort", "-@", "4", "-o", bam, "-"],
                   stdin=bwa.stdout, check=True, stderr=subprocess.DEVNULL)
    bwa.wait()
    subprocess.run(["samtools", "index", bam], check=True)
    for path in (fq1, fq2) + tuple(p for pair in fastqs for p in pair):
        os.remove(path)


def prepare(args):
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)
    ref = cut_reference(args.reference, args.outdir)
    ref_fa = os.path.join(args.outdir, "ref.fa")
    subprocess.run(["bwa", "index", ref_fa], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    events = make_events(ref, rng)
    with open(os.path.join(args.outdir, "truth.json"), "w") as fh:
        json.dump(events, fh, indent=1)

    # Maternal and paternal inherited SVs sit on the haplotype each parent
    # transmits; the child's de novo SVs on haplotype 1 or 2, and the
    # mosaic ones on a copy of haplotype 1 sequenced at half its depth.
    maternal = [e for e in events if e["category"] == "inherited_m"]
    paternal = [e for e in events if e["category"] == "inherited_f"]
    child1 = [e for e in events if e["haplotype"] == "1"]
    child2 = [e for e in events if e["haplotype"] == "2"]
    mosaic = [e for e in events if e["haplotype"] == "1m"]
    half = args.depth / 2
    samples = {
        "mother": [(haplotype(ref, maternal), half), (ref, half)],
        "father": [(haplotype(ref, paternal), half), (ref, half)],
        "child": [(haplotype(ref, child1), half / 2),
                  (haplotype(ref, child1 + mosaic), half / 2),
                  (haplotype(ref, child2), half)],
    }
    for sample, haps in samples.items():
        fastqs = []
        for i, (seqs, depth) in enumerate(haps):
            fastqs.append(simulate_reads(
                seqs, depth, os.path.join(args.outdir, f"{sample}_h{i}"),
                rng))
        align(ref_fa, fastqs, sample,
              os.path.join(args.outdir, f"{sample}.bam"), args.threads)
        print(f"{sample}.bam written")


# -- score --------------------------------------------------------------------

def breakpoints(e):
    if e["type"] == "BND":
        return [(e["chrom"], e["start"]), (e["chrom2"], e["pos2"])]
    if e["type"] in ("INS", "SNV") or e["end"] - e["start"] <= 1:
        return [(e["chrom"], e["start"])]
    return [(e["chrom"], e["start"]), (e["chrom"], e["end"])]


def _rows(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh
                if not line.startswith("#")]


def read_bed(path):
    return [{"chrom": f[0], "start": int(f[1]), "end": int(f[2]),
             "reads": int(f[3]), "class": f[9], "sv_type": f[12]}
            for f in _rows(path)]


def read_bedpe(path):
    return [{"a": (f[0], int(f[1]), int(f[2])),
             "b": (f[3], int(f[4]), int(f[5])),
             "reads": int(f[7]), "strands": f[8] + f[9], "sv_type": f[10]}
            for f in _rows(path)]


def near(region, chrom, pos, slack=MATCH_BP):
    return (region["chrom"] == chrom
            and region["start"] - slack <= pos <= region["end"] + slack)


def score(args):
    with open(os.path.join(args.outdir, "truth.json")) as fh:
        events = json.load(fh)
    regions = read_bed(args.prefix + ".bed")
    links = read_bedpe(args.prefix + ".sv.bedpe")
    rank = {"SV": 3, "AMBIGUOUS": 2, "SMALL": 1}
    used = set()
    rows = []
    for e in events:
        hits = []
        for chrom, pos in breakpoints(e):
            hits += [i for i, r in enumerate(regions) if near(r, chrom, pos)]
        hits = sorted(set(hits))
        used.update(hits)
        cls = max((regions[i]["class"] for i in hits), key=rank.get,
                  default="-")
        types = sorted({regions[i]["sv_type"] for i in hits})
        keys = {(regions[i]["chrom"], regions[i]["start"], regions[i]["end"])
                for i in hits}
        link = sorted({f"{l['sv_type']}{l['strands']}x{l['reads']}"
                       for l in links if l["a"] in keys and l["b"] in keys})
        rows.append((e, cls, types, link))

    print(f"{'id':<8}{'category':<22}{'type':<5}{'size':>7}  "
          f"{'class':<10}{'region types':<16}links")
    for e, cls, types, link in rows:
        print(f"{e['id']:<8}{e['category']:<22}{e['type']:<5}{e['size']:>7}  "
              f"{cls:<10}{','.join(types) or '-':<16}{' '.join(link) or '-'}")

    print("\nSummary (SV class, and region type matching the event):")
    tally = collections.defaultdict(collections.Counter)
    for e, cls, types, link in rows:
        key = (e["category"], e["type"])
        tally[key]["n"] += 1
        tally[key]["SV"] += cls == "SV"
        tally[key]["found"] += cls != "-"
        tally[key]["typed"] += e["type"] in types
    for (category, kind), t in sorted(tally.items()):
        print(f"  {category:<22}{kind:<5} found {t['found']}/{t['n']}  "
              f"SV {t['SV']}/{t['n']}  typed {t['typed']}/{t['n']}")
    extra = [r for i, r in enumerate(regions) if i not in used]
    classes = collections.Counter(r["class"] for r in extra)
    print(f"\nRegions away from every event: {len(extra)} "
          f"({', '.join(f'{n} {c}' for c, n in sorted(classes.items()))})")
    for r in extra:
        if r["class"] != "SMALL":
            print(f"  {r['chrom']}:{r['start']}-{r['end']} "
                  f"reads={r['reads']} {r['class']} {r['sv_type']}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare", help="Simulate and align the trio")
    p.add_argument("--reference", required=True,
                   help="GRCh38 FASTA with a .fai index")
    p.add_argument("--outdir", required=True)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--depth", type=float, default=30)
    p.add_argument("--threads", type=int, default=8)
    s = sub.add_parser("score", help="Compare discovery output to the truth")
    s.add_argument("--outdir", required=True)
    s.add_argument("--prefix", required=True,
                   help="kmer-discovery --out-prefix")
    args = parser.parse_args(argv)
    (prepare if args.command == "prepare" else score)(args)


if __name__ == "__main__":
    main()
