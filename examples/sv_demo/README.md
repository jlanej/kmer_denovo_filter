# Structural-variant demo

A small simulated trio whose child carries de novo structural variants. It
shows what `kmer-discovery` reports for each SV type, and runs in seconds.

The reference is random: chr1 is 60 kb and chr2 is 30 kb. The reads are
2×150 bp pairs at 30×. The parents carry only the reference. The child is
heterozygous for:

* a 2 kb deletion;
* a 1.5 kb tandem duplication;
* a 2 kb inversion;
* a reciprocal chr1/chr2 translocation;
* a 55 bp and a 250 bp novel insertion;
* an SNV.

See [`events.tsv`](expected/events.tsv) for the coordinates.

```bash
python examples/sv_demo/simulate_sv_trio.py --outdir sv_demo
kmer-discovery --child sv_demo/child.bam --mother sv_demo/mother.bam \
  --father sv_demo/father.bam --ref-fasta sv_demo/ref.fa \
  --out-prefix sv_demo/sv_demo
```

The simulator needs only pysam. It writes each read's alignment from the
read's known origin, following BWA-MEM's conventions, so no aligner is
needed. `kmer-discovery` needs samtools and Jellyfish.

The expected output is in [`expected/`](expected/):
[`sv_demo.bed`](expected/sv_demo.bed),
[`sv_demo.sv.bedpe`](expected/sv_demo.sv.bedpe) and
[`sv_demo.summary.txt`](expected/sv_demo.summary.txt).
[`tests/discovery/test_sv_demo.py`](../../tests/discovery/test_sv_demo.py)
checks that a fresh run matches it. The results by event:

| Event | Class | Type | BEDPE |
|---|---|---|---|
| Deletion | SV (both ends) | DEL | `+ -` DEL |
| Tandem duplication | SV (both ends) | DUP | `- +` DUP |
| Inversion | SV (both ends) | INV | `+ +` INV |
| Translocation | SV (both ends) | BND | `+ -` BND |
| 55 bp insertion | SV | INS (CIGAR) | – |
| 250 bp insertion | SV | `.` (clips and unmapped mates only) | – |
| SNV | SMALL | `.` | – |

For how each call is made, see
[Structural Variant Calling](../../docs/sv_calling.md#demo-a-simulated-trio).
