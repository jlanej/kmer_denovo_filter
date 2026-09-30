# Structural Variant Calling in Discovery Mode

`kmer-discovery` finds regions where the child's reads carry k-mers found in
neither parent nor the reference. At a structural variant (SV), those k-mers
span the SV's **junction**: the point where two pieces of the genome that are
apart in the reference are joined, or where new sequence is inserted. This
page explains how discovery mode turns the reads that carry them into SV
calls:

* evidence counts for each region;
* a class for each region (`SV`, `AMBIGUOUS` or `SMALL`);
* links between breakpoint regions (the BEDPE);
* an SV type from each breakpoint's orientation.

Two worked examples follow: a simulated trio that runs in seconds, and the
GIAB HG002 test data. The page ends with the limitations.

* [How it works](#how-it-works)
* [Evidence per region](#evidence-per-region)
* [Classes](#classes)
* [Linking breakpoints](#linking-breakpoints)
* [Breakpoint orientation and SV types](#breakpoint-orientation-and-sv-types)
* [Demo: a simulated trio](#demo-a-simulated-trio)
* [Real data: GIAB HG002](#real-data-giab-hg002)
* [Limitations](#limitations)
* [Working with the output](#working-with-the-output)

## How it works

```
child reads
  │  k-mers absent from the reference and from both parents ("proband-unique"):
  │  an SV junction makes k−1 of them, inserted sequence one more per base
  ▼
informative reads: at least --min-distinct-kmers-per-read proband-unique k-mers
  │  (default k/4); each one's alignment is recorded during the same scan:
  │  clips, SA tag, mate, CIGAR indels, strand, mapping quality
  ▼
regions: informative alignments within --cluster-distance (500 bp) of each other
  │  (--min-supporting-reads and --min-distinct-kmers are applied here)
  ▼
evidence per region ──► class: SV / AMBIGUOUS / SMALL
  ▼
links between regions (BEDPE) ──► SV type from breakpoint orientation
```

A k-mer that spans a junction holds bases from both sides of it, so it is
absent from the reference and from both parents. With k = 31, a read carries
min(*a*, *b*, 30) of these k-mers, where *a* and *b* are the read's bases on
each side of the junction. With the default of 7, a read that crosses a
junction with at least 7 bases on each side is informative. So is a read
that reaches at least 7 bases into inserted sequence.

The pipeline annotates, links, classifies and types the regions after
anchoring (Module 4). It reuses the alignment details recorded for each
informative read during anchoring, so it never reads the BAM a second time.

## Evidence per region

Each evidence column counts **molecules** (read names). A read pair or a
split read counts at most once per region, even when several of its
alignments are informative.

| Column | Molecules with… | Counted in |
|---|---|---|
| `split_reads` | a split alignment (an SA tag) | every region holding one of the molecule's informative alignments |
| `discordant_pairs` | a mapped mate, not properly paired | every region holding one of the molecule's informative alignments |
| `unmapped_mates` | an unmapped mate; or an unmapped informative read that the aligner placed at its mate's position | every region holding one of the molecule's informative alignments; for a placed read, the region at the placement |
| `breakpoint_reads` | a soft clip of at least 20 bp. This is the largest group clipped at positions within 5 bp of one another, and 0 when fewer than two, because a single long clip is often an adapter or a low-quality tail | the region of the clipped alignment |
| `large_indel_reads` | a CIGAR insertion or deletion of at least 50 bp | the region of the alignment |
| `max_clip_len` | (not a count) the longest soft clip | the region of the alignment |

Only informative reads count. A discordant pair whose two reads both lie
away from the junction carries no new k-mers, so it is never seen. The
counts therefore show how the junction-crossing reads align; they are not
everything an SV caller would find at the locus.

## Classes

| Class | Rule |
|---|---|
| `SV` | At least two molecules show one kind of evidence (`split_reads`, `discordant_pairs`, `unmapped_mates`, `breakpoint_reads` or `large_indel_reads` ≥ 2), or the region is linked to another in the BEDPE |
| `AMBIGUOUS` | Evidence from a single molecule only |
| `SMALL` | No evidence: the reads align end to end, as for SNVs and small indels |

`max_clip_len` is reported but does not affect the class. Clips mark a
breakpoint without measuring the event. An `SV` region can therefore hold an
indel shorter than 50 bp, like the 43 bp events in the
[GIAB example](#real-data-giab-hg002).

## Linking breakpoints

Two regions are linked when a molecule joins them. Each junction between
them is one line of `{prefix}.sv.bedpe` (see [Voting](#voting)). A molecule
can join two regions in three ways:

* **SA tag.** A supplementary alignment listed in the SA tag of an
  informative primary alignment lies in the other region.
* **Split molecule.** The molecule has informative alignments in both
  regions. Examples are both parts of a split read carrying junction k-mers,
  or both reads of a pair doing so.
* **Discordant pair.** The mate of an informative read that is not properly
  paired lies in the other region.

An SA or mate position within `--cluster-distance` of a region counts as in
it. Every alignment used for linking needs a mapping quality of at least 20:
the informative read, the SA entry, and the mate. The mate's MAPQ is checked
only when the BAM records it in the `MQ` tag. `supporting_reads` counts
the molecules that show the line's junction. A linked region is always
`SV`.

BWA-MEM hard-clips supplementary alignments by default, so they hold only
their aligned bases. With no junction k-mers they are not informative, and
the SA tag on the primary alignment makes the link. When supplementary
alignments are soft-clipped (`bwa mem -Y`), both parts of the split read are
informative.

## Breakpoint orientation and SV types

### Breakends

A junction joins two **breakends**. Each breakend is a position plus a side,
which says where the joined sequence lies:

```
                breakpoint
                    |
  "+"     ==========|             the joined sequence lies left of the breakpoint (ends there)
  "-"               |==========   the joined sequence lies right of it (starts there)
```

Each molecule that crosses the junction shows the side:

```
split read   an alignment clipped on its right     AAAAAAAAAAAA~~~~~~    -> "+" at its end
             an alignment clipped on its left            ~~~~~~CCCCCCCCCCCC  -> "-" at its start
read pair    a forward read: the breakpoint is ahead, to its right   -->     -> "+"
             a reverse read: the breakpoint is ahead, to its left    <--     -> "-"
```

`~` marks clipped bases, which can be soft or hard clips. An alignment
clipped at both ends takes the side of its longer clip. CIGARs are written
in reference orientation, so the clip rule holds on either strand. For the
other part of a split read, the SA tag gives its position, strand and CIGAR,
and hence its side as well. The pair rule assumes a standard forward-reverse
(FR) paired-end library.

The two breakends are ordered by chromosome name, then by position; for a
read pair, the position is where each read starts. On a tie, the `+` end
comes first. The ordered pair of sides gives the type:

| Sides (first, second) | Type | Meaning |
|---|---|---|
| `+ -` | `DEL` | The end of one piece joins the start of a piece further right, so the sequence between is lost |
| `- +` | `DUP` | The end of a piece joins back to an earlier start, so the sequence between is repeated (a tandem duplication) |
| `+ +` or `- -` | `INV` | A piece joins the reverse strand of another: one junction of an inversion |
| any, on two chromosomes | `BND` | A translocation, or another junction between chromosomes |
| unknown, on one chromosome | `INTRA` | Links only: the orientation is unknown, or the molecules disagree |

### What each type looks like

Capital letters are reference segments, and lowercase is a segment
reverse-complemented:

```
DEL  reference  AAAAAAAAAA BBBBBBBBBB CCCCCCCCCC
     child      AAAAAAAAAA CCCCCCCCCC
                         +|-                    end of A (+) joined to start of C (-)   -> (+, -)
     pairs across the junction face each other, farther apart than the library insert

DUP  reference  AAAAAAAAAA BBBBBBBBBB CCCCCCCCCC
     child      AAAAAAAAAA BBBBBBBBBB BBBBBBBBBB CCCCCCCCCC
                                    +|-         end of B (+) joined to start of B (-);
                                                on the reference the start comes first  -> (-, +)
     pairs across the junction face away from each other (reverse read first)

INV  reference  AAAAAAAAAA BBBBBBBBBB CCCCCCCCCC
     child      AAAAAAAAAA bbbbbbbbbb CCCCCCCCCC
                         +|+        -|-         junction 1: end of A (+) to end of B (+)      -> (+, +)
                                                junction 2: start of B (-) to start of C (-)  -> (-, -)
     both reads of a pair across a junction are on the same strand

BND  reference  chr1: AAAAAAAAAA BBBBBBBBBB     chr2: CCCCCCCCCC DDDDDDDDDD
     child      der(1): AAAAAAAAAA DDDDDDDDDD   end of A (+) to start of D (-)          -> (+, -)
                der(2): CCCCCCCCCC BBBBBBBBBB   end of C (+) to start of B (-);
                                                chr1 comes first                        -> (-, +)
```

**Insertions** have no second breakend. An aligner can place a short
insertion in the CIGAR, such as `45M55I50M`. A read with a CIGAR insertion
or deletion of at least 50 bp votes `INS` or `DEL` for its region. Reads
reaching into a longer insertion are clipped at the insertion point from
both sides, and reads entirely inside it are unmapped. That gives
`breakpoint_reads` and `unmapped_mates` but no type (see
[Limitations](#limitations)).

### Voting

Each molecule votes once for each pair of regions it joins, and once for
each region it touches:

* **Junctions (BEDPE lines).** The molecules joining two regions vote for an
  SV type, and each orientation of the majority type is one line. That line
  is a junction, and it counts the molecules that show it. Molecules of a
  minority type are outvoted. With no majority (no orientation known, or a
  tie), the two regions get one line with `.` strands: `INTRA` on one
  chromosome, `BND` across two.
* **Region type (BED column 13).** This is the majority over the region's
  molecules. Votes come from its links, from junctions that lie wholly
  inside the region (such as a deletion shorter than `--cluster-distance`),
  and from CIGAR indels of at least 50 bp (`DEL` or `INS`). It is `.` when
  there are no votes or they tie.
* **No-vote case.** A forward-reverse pair with both reads in one region
  does not vote. Such a pair is flagged discordant either because its
  insert is too long (a deletion) or too short (an insertion, such as one
  longer than the reads), and the pair alone cannot say which.
* **Balanced events.** Both junctions of a balanced inversion or a
  reciprocal translocation join the same two regions. Each gets its own
  line, so the event shows as two lines with the same coordinates and
  opposite orientations. A single line of type `INV` or `BND` means only one
  junction was seen, as for an unbalanced translocation.

## Demo: a simulated trio

[`examples/sv_demo/simulate_sv_trio.py`](../examples/sv_demo/simulate_sv_trio.py)
writes a trio on a random two-chromosome reference. chr1 is 60 kb and chr2
is 30 kb, and the reads are 2×150 bp pairs at 30×. The parents carry only
the reference. The child is heterozygous for seven events, six of them
structural variants:

| Event | Position (1-based) |
|---|---|
| 2 kb deletion | chr1:10,001–12,000 |
| 1.5 kb tandem duplication | chr1:25,001–26,500 |
| 2 kb inversion | chr1:40,001–42,000 |
| Reciprocal translocation | after chr1:55,000 ↔ after chr2:15,000 |
| SNV | chr2:5,001 |
| 55 bp novel insertion | after chr2:10,000 |
| 250 bp novel insertion | after chr2:22,000 |

No aligner is run. Each read's alignment is written from where it came
from, following BWA-MEM's conventions:

* Split reads get a soft-clipped primary and hard-clipped supplementary
  alignments, with SA tags.
* Parts shorter than 30 bp are soft-clipped.
* Indels of up to 100 bp go in the CIGAR when both flanks have at least
  20 bp.
* Unmapped reads are placed at their mate's position.
* Pairs are proper only when they face each other with an insert size
  within 4 SD of the mean.

It needs only pysam, plus samtools and Jellyfish for `kmer-discovery`:

```bash
python examples/sv_demo/simulate_sv_trio.py --outdir sv_demo
kmer-discovery --child sv_demo/child.bam --mother sv_demo/mother.bam \
  --father sv_demo/father.bam --ref-fasta sv_demo/ref.fa \
  --out-prefix sv_demo/sv_demo
```

The Docker image includes the script at `/app/examples/sv_demo/`:

```bash
docker run --rm -v "$PWD:/data" --entrypoint python \
  ghcr.io/jlanej/kmer_denovo_filter:latest \
  /app/examples/sv_demo/simulate_sv_trio.py --outdir /data/sv_demo
docker run --rm -v "$PWD:/data" --entrypoint kmer-discovery \
  ghcr.io/jlanej/kmer_denovo_filter:latest \
  --child /data/sv_demo/child.bam --mother /data/sv_demo/mother.bam \
  --father /data/sv_demo/father.bam --ref-fasta /data/sv_demo/ref.fa \
  --out-prefix /data/sv_demo/sv_demo
```

Both steps take a few seconds. The output is deterministic, and
[`tests/discovery/test_sv_demo.py`](../tests/discovery/test_sv_demo.py)
checks it against
[`examples/sv_demo/expected/`](../examples/sv_demo/expected/).

Eleven regions are found: one at each breakpoint position, one at the SNV
and one at each insertion. The table uses 1-based coordinates. BkptClip is
`breakpoint_reads` and Indel50 is `large_indel_reads`:

| Event | Region | Reads | Split | Disc | Unmapped mates | BkptClip | Indel50 | Class | Type |
|---|---|---|---|---|---|---|---|---|---|
| Deletion, left end | chr1:9,861–10,000 | 11 | 6 | 4 | 0 | 8 | 0 | SV | DEL |
| Deletion, right end | chr1:12,001–12,140 | 8 | 7 | 3 | 0 | 7 | 0 | SV | DEL |
| Duplication, start | chr1:25,001–25,131 | 10 | 7 | 4 | 0 | 9 | 0 | SV | DUP |
| Duplication, end | chr1:26,359–26,500 | 10 | 7 | 6 | 0 | 8 | 0 | SV | DUP |
| Inversion, left end | chr1:39,871–40,141 | 13 | 9 | 6 | 0 | 12 | 0 | SV | INV |
| Inversion, right end | chr1:41,891–42,133 | 11 | 9 | 5 | 0 | 10 | 0 | SV | INV |
| Translocation, chr1 | chr1:54,861–55,115 | 11 | 8 | 5 | 0 | 9 | 0 | SV | BND |
| SNV | chr2:4,873–5,124 | 10 | 0 | 0 | 0 | 0 | 0 | SMALL | . |
| 55 bp insertion | chr2:9,859–10,132 | 23 | 0 | 0 | 0 | 12 | 7 | SV | INS |
| Translocation, chr2 | chr2:14,862–15,141 | 12 | 6 | 6 | 0 | 10 | 0 | SV | BND |
| 250 bp insertion | chr2:21,867–22,078 | 10 | 0 | 1 | 7 | 9 | 0 | SV | . |

`sv_demo.sv.bedpe` (0-based, like BED):

```
#chrom1  start1  end1   chrom2  start2  end2   sv_id  supporting_reads  strand1  strand2  sv_type
chr1     9860    10000  chr1    12000   12140  SV_1   15                +        -        DEL
chr1     25000   25131  chr1    26358   26500  SV_2   18                -        +        DUP
chr1     39870   40141  chr1    41890   42133  SV_3   13                +        +        INV
chr1     39870   40141  chr1    41890   42133  SV_4   7                 -        -        INV
chr1     54860   55115  chr2    14861   15141  SV_5   10                +        -        BND
chr1     54860   55115  chr2    14861   15141  SV_6   8                 -        +        BND
```

How each event is called:

* **Deletion.** Reads across the junction split into a part that ends at
  10,000 and is clipped on its right (`+`), and a part that starts at
  12,001 and is clipped on its left (`-`). Pairs across it have a forward
  read on the left and a reverse read on the right. Both breakpoints are
  `SV` regions typed `DEL`, and SV_1 links them as `+ -`.
* **Tandem duplication.** Reads run off the end of the duplicated segment
  (26,500) and continue at its start (25,001). On the reference the start
  comes first, which gives `- +`, a `DUP`. Pairs across the junction face
  away from each other.
* **Inversion.** Junction 1 joins 40,000 to 42,000, giving `+ +`. There the
  second part of each split read aligns to the reverse strand and is
  clipped on its right. Junction 2 joins 40,001 to 42,001, giving `- -`.
  Both junctions join the same two regions: SV_3 is junction 1 (13
  molecules) and SV_4 is junction 2 (7 molecules).
* **Translocation.** der(1) joins chr1:55,000 to chr2:15,001: `+ -`, SV_5,
  10 molecules. der(2) joins chr2:15,000 to chr1:55,001, which reads `- +`
  with chr1 first: SV_6, 8 molecules.
* **SNV.** Its reads align end to end, so the region is `SMALL`.
* **55 bp insertion.** Seven molecules span it with a CIGAR insertion,
  which types the region `INS`. Twelve molecules are clipped at the
  insertion point.
* **250 bp insertion.** The insertion is too long for a read to align on
  both sides of it. Reads are clipped at the insertion point (9 molecules),
  and reads inside it are unmapped, so 7 molecules have an unmapped mate.
  The region is `SV` but untyped (`.`):
  * Clips alone cannot tell an insertion from one side of a translocation.
  * The one discordant pair is forward-reverse, with an insert 250 bp
    short. Its reads are in one region, so it does not vote.

## Real data: GIAB HG002

The integration test runs discovery on slices of the GIAB HG002 trio (see
[`tests/README.md`](../tests/README.md)), with 7 SV-like *de novo* events
curated by Sulovari et al. 2023. The committed output in
[`tests/example_output_discovery/`](../tests/example_output_discovery/)
shows three ways the SV evidence appears on real reads. These BAMs were
aligned with novoalign and carry no SA tags, so every link comes from read
pairs.

| Event | Region | Reads | Disc | Unmapped mates | BkptClip | Indel50 | Class | Type |
|---|---|---|---|---|---|---|---|---|
| 107 bp deletion at chr17:53,340,465 | chr17:53,340,222–53,341,164 | 18 | 0 | 1 | 6 | 3 | SV | DEL |
| 10.6 kb deletion at chr7:142,786,222 (TRB), left end | chr7:142,785,588–142,786,223 | 7 | 0 | 4 | 0 | 0 | SV | DEL |
| ″ inside the deleted interval | chr7:142,788,519–142,789,264 | 3 | 0 | 1 | 0 | 0 | AMBIGUOUS | . |
| ″ inside the deleted interval | chr7:142,792,552–142,792,982 | 2 | 0 | 1 | 0 | 0 | AMBIGUOUS | . |
| ″ right end | chr7:142,796,830–142,797,423 | 14 | 2 | 4 | 7 | 0 | SV | DEL |
| 43 bp SV-like event at chr5:97,089,276 | chr5:97,089,091–97,089,510 | 22 | 0 | 0 | 3 | 0 | SV | . |
| 43 bp insertion at chr8:125,785,998 | chr8:125,785,801–125,786,483 | 34 | 0 | 0 | 5 | 0 | SV | . |
| 34 bp insertion at chr18:62,805,217 | chr18:62,804,895–62,805,439 | 7 | 0 | 0 | 0 | 0 | SMALL | . |

The three patterns:

* **Deletion within one region.** Both ends of the 107 bp deletion on chr17
  fall in one region. Three reads carry it as `107D` in their CIGAR, which
  types the region `DEL`, and six are clipped at one breakpoint. There is no
  link, because there is no second region.
* **Deletion across two regions.** The 10.6 kb deletion on chr7 has its
  ends in separate regions, both typed `DEL`. The one BEDPE link joins them:

  ```
  chr7  142785587  142786223  chr7  142796829  142797423  SV_1  2  +  -  DEL
  ```

  Two reverse reads at the right end have forward mates at the left end
  (`+ -`, an 11 kb insert). Two more regions inside the deleted interval
  each have one unmapped mate, so they are `AMBIGUOUS`.
* **Insertions shorter than 50 bp.** The 43 bp events on chr5 and chr8 are
  `SV` from clips at one breakpoint. They are too short to count as
  `large_indel_reads`, and clips give no type. The 34 bp insertion on chr18
  has no two molecules clipped at one breakpoint, so it is `SMALL`.

## Limitations

* **Only novel junctions are visible.** An SV is found only if its junction
  creates k-mers absent from the reference and from both parents. This
  misses:
  * SVs whose breakpoints lie in identical repeat copies longer than k
    (for example, recombination between segmental duplications), so the
    junction sequence already exists in the reference;
  * SVs whose junction sequence occurs elsewhere in the genome.
* **Evidence comes only from informative reads.** Discordant pairs and
  clipped reads that do not carry junction k-mers are not counted (see
  [Evidence per region](#evidence-per-region)).
* **Pair orientation assumes forward-reverse (FR) libraries.** This is the
  standard Illumina paired-end layout. Mate-pair (RF) libraries would be
  typed wrongly.
* **Some insertions are untyped.** Insertions that the aligner leaves out of
  the CIGAR give `SV` with type `.`. This includes insertions longer than
  the reads, such as mobile-element insertions.
* **Nearby breakpoints share a region.** Breakpoints within
  `--cluster-distance` of each other fall in one region. The SV is typed but
  not linked.
* **Region filters run first.** The region filters (`--min-supporting-reads`,
  `--min-distinct-kmers`) run before annotation. A breakpoint region they
  remove cannot be linked.
* **Types are hints, not calls.** Discovery mode does not genotype SVs,
  refine breakpoints or assemble inserted sequence. Confirm candidates with
  a dedicated SV caller, or by viewing the informative reads.

## Working with the output

* **View the evidence in IGV.** Load `{prefix}.informative.bam` with the BED
  and BEDPE; IGV draws BEDPE files as arcs. Every read in the BAM carries
  proband-unique k-mers (`dk:i:1`).
* **Select SV regions.** Use the class column:
  `awk '$10 == "SV"' {prefix}.bed`. To also require a type, add
  `&& $13 != "."`.
* **Compare with other calls.** Use `bedtools pairtobed` or
  `bedtools pairtopair`. The BEDPE's first ten columns follow the standard
  layout, and each line is one junction.
* **Check known events.** `--dnm-regions` reports, for each listed event,
  the regions that overlap it or lie within `--cluster-distance` of it. A
  deletion's two breakpoint regions both count, although they flank the
  deleted interval rather than overlap it.
