# Kraken2 Non-Human Content Detection

This document describes how `kmer-denovo` uses [Kraken2](https://github.com/DerrickWood/kraken2)
to detect non-human content — bacteria, archaea, fungi, protists, viruses,
and synthetic sequencing vectors — in the child's sequencing reads, and why
Kraken2's k-mer–based classification approach is well suited to that goal.

> **Scope / where the engine lives.** The classification *engine* (Kraken2 LCA
> classification, lineage-aware domain assignment, the human-homology guard,
> UniVec-Core exclusion, and the non-human fraction) was extracted into the
> standalone [`nonhuman-screen`](../packages/nonhuman-screen) package and is
> **authoritatively documented there** —
> [methodology](../packages/nonhuman-screen/docs/methodology.md) and
> [database setup](../packages/nonhuman-screen/docs/database.md). This document
> covers only how `kmer-denovo` *integrates* that engine into the trio/VCF de
> novo workflow: informative-read selection, the `DKU_*`/`DKA_*` VCF
> annotations, the Kraken2 BED outputs, and the `--kraken2*` flags. Engine
> internals below are summarized for context; the package docs are the source
> of truth.

---

## Why Non-Human Content Detection Matters

`kmer-denovo` identifies *de novo* variants by looking for k-mers present in a
child but absent from both parents. When the child's sample contains non-human
contamination (bacterial, archaeal, fungal, protist, viral, or other),
sequencing reads from these organisms can carry k-mers that are truly absent
from parental genomes — not because a *de novo* mutation occurred, but simply
because the non-human sequence has no counterpart in the human reference or in
the parents. Without a contamination check, these reads would be
indistinguishable from genuine *de novo* signals.

Reads flagged as non-human can be used to compute per-domain fraction
annotations (e.g. **DKU_BF**, **DKU_AF**, **DKU_FF**, **DKU_PF**,
**DKU_VF**) and a consolidated **DKU_NHF** (non-human fraction), which
indicate what proportion of the informative reads supporting a candidate
variant appear to derive from a non-human organism rather than from the human
child genome.  **DKU_UCF** separately tracks the fraction of reads classified
as UniVec Core (synthetic vectors/adapters); these reads are excluded from
DKU_NHF because they are artificial constructs, not biological contamination.
A high non-human fraction is a strong indicator of a false-positive *de novo*
call.

---

## How Kraken2 Classification Works

Kraken2 assigns each read a taxon via k-mer–based LCA classification, gated by a
`--confidence` threshold (exposed here as `--kraken2-confidence`, default
`0.0`). See the package
[methodology §1](../packages/nonhuman-screen/docs/methodology.md) for the full
algorithm and confidence-threshold semantics.

The one detail this integration depends on directly is the **per-read output
format** the host parses into BED columns — Kraken2 emits one line per read:

```
C/U  read_name  taxid  length  kmer_detail_string
```

- `C` = classified, `U` = unclassified
- `taxid` = NCBI taxonomy ID of the LCA classification
- `kmer_detail_string` = space-separated `taxid:count` tokens showing how many
  k-mers voted for each taxid (paired-end reads use `|:|` as a mate delimiter)

---

## The PrackenDB Reference Database

Database acquisition, layout, required files, and the k-mer length are
documented in the package
[database setup guide](../packages/nonhuman-screen/docs/database.md); the
bundled `download_kraken2_db.sh` fetches and validates **PrackenDB** (a curated,
pre-built Kraken2 database with `taxonomy/nodes.dmp` and `taxonomy/names.dmp`).

**Why PrackenDB for this workflow**: it uses a single reference genome per
species (with a couple of exceptions such as normal and pathogenic *E. coli*),
so a k-mer shared across species is LCA-elevated to a genus/family node rather
than an unrelated lineage. This keeps k-mer counting per species unambiguous —
which matters because `kmer-denovo` reasons about per-read, per-species
evidence.

---

## How `kmer-denovo` Uses Kraken2 Output

### Step 1 — Identify informative reads

The pipeline first identifies *informative* child reads: reads that carry at
least one variant-spanning k-mer absent from both parents (DKU reads). These
are the reads that support candidate *de novo* variants.

### Step 2 — Classify informative reads with Kraken2

Only informative reads are passed to Kraken2 (not the entire BAM). The pipeline
extracts their sequences from the child BAM/CRAM and writes a temporary FASTQ.
This is substantially more efficient than classifying all reads, and ensures
the fraction annotations are computed on exactly the reads contributing to each
variant's evidence.

The host does not call the engine per read; `_run_kraken2_on_reads`
(`src/kmer_denovo_filter/vcf/pipeline.py`) translates the informative
`chrom:pos → read-names` map into loci and delegates extraction +
classification to the package:

```python
from nonhuman_screen.bam import classify_reads_from_bam
result = classify_reads_from_bam(child_bam, kraken2_db, read_names=..., loci=...)
```

### Step 3 — Lineage-aware multi-domain classification

The engine assigns each read a single taxid and maps it to one or more domains
by walking the NCBI taxonomy (`nodes.dmp`). The domains and their root taxids
(Bacteria 2, Archaea 2157, Fungi 4751, Viruses 10239, UniVec-Core 81077, and
Protist = Eukaryota − Metazoa − Fungi − Viridiplantae) — plus the exact-taxid
fallback when `nodes.dmp` is missing — are documented in the package
[methodology §2 and §7](../packages/nonhuman-screen/docs/methodology.md).

What matters for the integration: each informative read contributes to the
per-domain **DKU_\*/DKA_\*** fractions below according to its assigned domain,
and lineage-aware matching means a read assigned to *E. coli* (562) or
*S. aureus* (1280) counts toward the bacterial fraction — not only reads whose
LCA is exactly taxid 2. UniVec-Core (81077) reads are tracked separately and
excluded from the non-human fraction.

### Step 4 — Human homology guard

The engine's **human-homology guard** drops any read with human (taxid 9606)
k-mer evidence from *every* non-human numerator — a conservative measure that
avoids over-flagging human reads carrying non-human-like k-mers. The mechanism,
and why it matters for integrating viruses (ERVs, HBV, HPV), is documented in
the package
[methodology §3](../packages/nonhuman-screen/docs/methodology.md).

For this integration the consequences are:

- the guard status is surfaced per read in the detail BED (`guard_status`
  column, value `HHG` for a guard-excluded read); and
- **DKU_VF reflects only exogenous, non-integrated viral contamination** — a
  read spanning an integrating virus's integration junction carries human
  k-mers and is conservatively excluded from the viral count.

### Step 5 — Conservative non-human fraction (NHF)

The consolidated **non-human fraction** (DKU_NHF / DKA_NHF) counts a read as
non-human only when it is classified, off the human→root lineage, outside the
human clade, outside UniVec-Core, and clears the human-homology guard. The full
read-inclusion definition, worked taxid examples, and the four-way partition
(`nonhuman + univec_core + human_lineage + unclassified = 1`) are documented in
the package
[methodology §5](../packages/nonhuman-screen/docs/methodology.md).

The host maps those per-domain and consolidated fractions onto the `DKU_*`
(all informative reads) and `DKA_*` (alt-supporting reads) VCF tags described in
[Output Annotations](#output-annotations) below.

---

## Output Annotations

The classification results are added to the output VCF as per-variant fields:

### Domain-specific fractions

| Field | Description |
|-------|-------------|
| **DKU_BF** | Fraction of DKU fragments classified as **bacterial** by Kraken2. Denominator = DKU (both are fragment-based, i.e. unique read names). Always in [0.0, 1.0]. |
| **DKA_BF** | Fraction of DKA fragments (DKU fragments that also directly support the alternate allele) classified as **bacterial**. DKA fragments are always a subset of DKU. Always in [0.0, 1.0]. |
| **DKU_AF** | Fraction of DKU fragments classified as **archaeal** by Kraken2. |
| **DKA_AF** | Fraction of DKA fragments classified as **archaeal**. |
| **DKU_FF** | Fraction of DKU fragments classified as **fungal** by Kraken2. |
| **DKA_FF** | Fraction of DKA fragments classified as **fungal**. |
| **DKU_PF** | Fraction of DKU fragments classified as **protist** by Kraken2. |
| **DKA_PF** | Fraction of DKA fragments classified as **protist**. |
| **DKU_VF** | Fraction of DKU fragments classified as **viral** by Kraken2 (RefSeq viral genomes in PrackenDB). Reads with any human k-mer evidence are conservatively excluded, which handles viruses that can integrate into human DNA. |
| **DKA_VF** | Fraction of DKA fragments classified as **viral**. |
| **DKU_UCF** | Fraction of DKU fragments classified as **UniVec Core** (synthetic sequencing-vector and adapter sequences, taxid 81077) by Kraken2. Reads with any human k-mer evidence are excluded. UniVec Core reads are **not** included in DKU_NHF because they are artificial constructs, not biological contamination. |
| **DKA_UCF** | Fraction of DKA fragments classified as **UniVec Core**. |

### Consolidated non-human fraction

| Field | Description |
|-------|-------------|
| **DKU_NHF** | Fraction of DKU fragments classified as **non-human** by Kraken2 (consolidated across all non-human domains). UniVec Core reads are excluded. |
| **DKA_NHF** | Fraction of DKA fragments classified as **non-human**. UniVec Core reads are excluded. |

### Unclassified fraction

| Field | Description |
|-------|-------------|
| **DKU_UF** | Fraction of DKU fragments that received **no taxonomic assignment** from Kraken2 (status "U"). |
| **DKA_UF** | Fraction of DKA fragments that received **no taxonomic assignment**. |

### Human lineage fraction

| Field | Description |
|-------|-------------|
| **DKU_HLF** | Fraction of DKU fragments in the **human lineage**: classified reads that are neither definitively non-human (NHF) nor UniVec Core (UCF). Includes reads directly classified as human, reads cleared by the human homology guard (HHG), and reads assigned to broad taxonomic ranks on the human-to-root path (e.g. Eukaryota, Root). |
| **DKA_HLF** | Fraction of DKA fragments in the **human lineage**. |

### Sum-to-one property

For every variant, the four partition fractions sum to 1.0:

```
DKU_NHF + DKU_UCF + DKU_HLF + DKU_UF = 1.0
DKA_NHF + DKA_UCF + DKA_HLF + DKA_UF = 1.0
```

This provides a complete accounting of all variant-supporting reads: every read
is in exactly one of non-human (NHF), UniVec Core (UCF), human lineage (HLF),
or unclassified (UF).

These are per-variant fractions: each fraction is computed from the intersection
of that variant's informative reads with the global set of domain-specific or
non-human read names returned by Kraken2.

**Interpretation**:

- `DKU_NHF` close to `1.0` — essentially all evidence for this variant comes
  from reads classified as non-human; strong indicator of contamination artifact
- `DKU_BF` close to `1.0` — specifically bacterial contamination
- `DKU_VF` close to `1.0` — specifically exogenous viral contamination (reads
  with integrated-virus k-mer signatures are excluded via the human homology guard)
- `DKU_UCF` close to `1.0` — reads are classified as synthetic sequencing
  vectors/adapters; likely library-preparation artifacts, not biological contamination
- `DKA_NHF` close to `1.0` — reads that directly support the alternate allele
  sequence are predominantly non-human; high-confidence contamination flag
- All fractions near `0.0` — no detectable non-human content among the
  supporting reads; the candidate variant is more likely genuine
- `DKU_NHF` ≥ `DKU_BF` — the non-human fraction is always at least as large as
  any individual domain fraction, since it consolidates all non-human categories
- `DKU_UCF` is separate from `DKU_NHF` — UniVec Core reads are tracked
  independently and never inflate the non-human contamination signal
- `DKU_HLF` close to `1.0` — essentially all variant-supporting reads are in
  the human lineage; the variant is likely a genuine human variant
- `DKU_UF` close to `1.0` — nearly all variant-supporting reads could not be
  classified; may warrant investigation of read quality or database coverage
- `DKU_NHF + DKU_UCF + DKU_HLF + DKU_UF = 1.0` — a quick sanity check that
  all reads are accounted for in exactly one category

---

## Per-Read Classification Detail (BED)

In addition to the VCF summary fractions, the pipeline writes a companion
**BED file** with one row per (variant, read) pair.  This exposes the full
per-read Kraken2 classification detail so users can audit any individual
variant.

### Why BED and not VCF INFO fields

Per-read detail is **not appropriate for VCF INFO fields** because VCF
parsers expect well-typed, fixed-schema fields — not multi-kilobyte
free-text blobs.  The BED file is directly queryable with standard
genomics tools (`tabix`, `bedtools`), loadable in pandas/R, and joinable
to the VCF on the variant key.

### File naming

When `--kraken2-db` is provided in VCF mode, the BED file is written
alongside the output VCF:

- If `--output my_trio.annotated.vcf.gz`, the BED is
  `my_trio.annotated.kraken2_reads.bed.gz` (with a `.tbi` tabix index).
- Override with `--kraken2-read-detail <path>`.

### Columns

The file uses a `#`-prefixed header so that downstream tools can parse the
schema.  The first three columns are standard BED coordinates; subsequent
columns carry classification detail.

| Column | Type | Description |
|--------|------|-------------|
| `#chrom` | String | Chromosome name. |
| `chromStart` | Integer | 0-based start position (same as the internal pipeline `pos`). |
| `chromEnd` | Integer | Exclusive end position (`chromStart + len(ref)`). |
| `variant` | String | Variant key `chr:pos:ref:alt` (0-based `pos`). Join key to the VCF — add 1 to convert to VCF POS. |
| `read_name` | String | Fragment name (SAM QNAME). |
| `read_set` | Enum | `DKU` (informative-only) or `DKA` (also supports alt allele). |
| `kraken2_status` | Enum | `C` (classified) or `U` (unclassified). |
| `assigned_taxid` | Integer | NCBI taxonomy ID assigned by Kraken2. `0` for unclassified. |
| `assigned_taxon` | String | `.` if the read is unclassified; otherwise the scientific name from `names.dmp` (spaces → underscores), falling back to the numeric taxid when `names.dmp` is unavailable. |
| `domain` | String | `Bacteria`, `Archaea`, `Fungi`, `Protist`, `Viruses`, `UniVec_Core`, `Human`, `Root`, `Unclassified`, or `Ambiguous_Ancestor`. |
| `guard_status` | String | `PASS`, `HHG` (human homology guard), `UVC` (UniVec Core), `HUMAN`, or `UNCLASSIFIED`. |
| `is_nonhuman` | Boolean | `true` if counted in NHF numerator after all guards. |
| `kmer_votes` | String | Top k-mer vote summary: `taxid1:count1;taxid2:count2;...` (top 10, descending). |
| `kmer_votes_named` | String | Same with scientific names: `Escherichia_coli:25;Bacteria:8;unclassified:3`. |
| `total_kmers` | Integer | Total classified + unclassified k-mers in the read. |
| `human_kmer_count` | Integer | K-mers that voted for taxid 9606 (human). |

### Sort order

Rows are sorted by chromosome (lexicographic), then position (numeric),
then read name (lexicographic).

### Tabix indexing

The output is bgzipped and tabix-indexed (`-p bed`), so regions can be
queried efficiently:

```bash
tabix my_trio.annotated.kraken2_reads.bed.gz chr1:1000-2000
```

### Example

```
#chrom	chromStart	chromEnd	variant	read_name	read_set	kraken2_status	assigned_taxid	assigned_taxon	domain	guard_status	is_nonhuman	kmer_votes	kmer_votes_named	total_kmers	human_kmer_count
chr1	1000	1001	chr1:1000:A:T	read001	DKA	C	562	Escherichia_coli	Bacteria	PASS	true	562:25;2:8;0:3	Escherichia_coli:25;Bacteria:8;unclassified:3	36	0
chr1	1000	1001	chr1:1000:A:T	read002	DKU	C	9606	Homo_sapiens	Human	HUMAN	false	9606:40;0:2	Homo_sapiens:40;unclassified:2	42	40
```

### `names.dmp` requirement

The `assigned_taxon` and `kmer_votes_named` columns use `taxonomy/names.dmp`
(or `names.dmp` at the DB root for PrackenDB) to look up scientific names.
If `names.dmp` is not present, the BED file falls back to numeric taxids.
The download script (`download_kraken2_db.sh`) validates the presence of
`names.dmp` and warns if missing.  PrackenDB includes this file.

---

## Genomic Span BED File

In addition to the per-read classification detail BED described above, the
pipeline writes a second **species-annotated genomic span BED file** that maps
each informative read's **aligned reference span** to its Kraken2-assigned
species.  This enables:

1. **Visual audit in IGV**: Load the BED alongside the informative reads BAM
   (`--informative-reads`) to see species labels on each read's footprint.
2. **Spatial clustering detection**: Non-human reads clustered at a single locus
   suggest a contamination pile-up (likely false positive).  Scattered non-human
   reads suggest low-level background contamination.
3. **Clipping-aware interpretation**: A read with large soft clips classified as
   *Bacteria* may represent a chimeric molecule (human + bacterial junction from
   library-prep contamination), distinct from a fully-bacterial read that happened
   to map due to low-complexity sequence.

### Important Scientific Scope Limitation

**Kraken2 classifies the full read, not sub-regions.**  The BED coordinates
represent the read's **aligned reference span** (from `reference_start` to
`reference_end` in pysam), and the species label applies to the **entire read**.
You cannot conclude "bases chr1:1000–1050 are bacterial and 1050–1100 are human"
from this file.  The spatial information here is about **where the read aligns**,
not where the non-human sequence within the read begins or ends.

For split reads (SA tag), both alignment segments are emitted as separate BED
intervals with the same classification, linked by `read_name`.

### File naming

Written alongside the output VCF when `--kraken2-db` is provided:

- Default: derived from `--output` (e.g. `my_trio.annotated.kraken2_spans.bed.gz`)
- Override with `--kraken2-span-bed <path>`

### Columns

Standard BED4+ format (tab-delimited, 0-based half-open coordinates), bgzipped
and tabix-indexed.

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `chrom` | String | Reference contig name. |
| 2 | `start` | Integer | 0-based aligned start (`reference_start`). |
| 3 | `end` | Integer | 0-based exclusive aligned end (`reference_end`). |
| 4 | `taxon_name` | String | Scientific name (underscored) from Kraken2 LCA assignment. `Unclassified` for unclassified reads. `Unknown_taxid_NNN` if `names.dmp` is unavailable. |
| 5 | `domain` | String | `Bacteria`, `Archaea`, `Fungi`, `Protist`, `Viruses`, `UniVec_Core`, `Human`, `Root`, `Unclassified`, or `Ambiguous_Ancestor`. |
| 6 | `guard_status` | String | `PASS`, `HHG`, `UVC`, `HUMAN`, or `UNCLASSIFIED`. |
| 7 | `is_nonhuman` | String | `true` or `false`. Final NHF determination after all guards. |
| 8 | `read_name` | String | SAM QNAME. |
| 9 | `variant` | String | Comma-separated `chr:pos:ref:alt` variant key(s) this read is informative for. |
| 10 | `read_set` | String | `DKA` or `DKU`. |
| 11 | `mapq` | Integer | Mapping quality of this alignment. |
| 12 | `softclip_left` | Integer | Number of soft-clipped bases at the left (5′ aligned) end. |
| 13 | `softclip_right` | Integer | Number of soft-clipped bases at the right (3′ aligned) end. |
| 14 | `is_split` | String | `true` if the read has an SA (supplementary alignment) tag; `false` otherwise. |
| 15 | `is_supplementary` | String | `true` if this specific alignment record is the supplementary; `false` if it is the primary. |

### Sort order & indexing

Sorted by `chrom` (reference order), then `start`.  Bgzipped and tabix-indexed
(`tabix -p bed`) for efficient region lookup:

```bash
tabix my_trio.annotated.kraken2_spans.bed.gz chr1:100000-200000
```

### Example

```bed
chr1	100000	100150	Escherichia_coli	Bacteria	PASS	true	read001	chr1:100050:A:T	DKA	60	0	5	false	false
chr1	100020	100170	Homo_sapiens	Human	HUMAN	false	read002	chr1:100050:A:T	DKU	60	0	0	false	false
chr1	100030	100100	Bacteria	Bacteria	HHG	false	read003	chr1:100050:A:T	DKA	25	42	0	true	false
chr7	50000	50070	Bacteria	Bacteria	HHG	false	read003	chr1:100050:A:T	DKA	0	0	42	true	true
```

In this example, `read003` is a split read: the primary alignment (chr1) has
42 bp left-clipped and is classified as *Bacteria* but was intercepted by the
human homology guard.  Its supplementary alignment lands on chr7.  Both records
carry the same Kraken2 classification (because Kraken2 classified the full read
sequence, not each alignment segment independently).

### Interpretation Guidance

**Clipping patterns and what they suggest:**

| Pattern | Likely Interpretation |
|---------|----------------------|
| Non-human read, no clips, high MAPQ | Genuinely non-human sequence that happens to align to the human reference (low-complexity or conserved region). Check the k-mer votes — if overwhelming bacterial/viral, likely real contamination. |
| Non-human read, large soft clips (>30 bp), moderate MAPQ | Possible chimeric molecule: the aligned portion maps to human, the clipped portion is non-human. Common in library-prep contamination where bacterial and human DNA ligate. |
| Non-human read, split alignment (`is_split=true`) | The aligner split the read across two locations. If one segment is in a known bacterial integration site, this may be a genuine insertion. If the segments are on different chromosomes with no biological rationale, likely an artifact. |
| HHG-intercepted read, low human k-mer count (1–3) | Possibly a genuine non-human read with a single conserved k-mer that happened to match human. Consider the ratio: 2 human k-mers out of 80 total is very different from 30/80. The companion per-read detail BED's `human_kmer_count` column provides this detail. |
| Cluster of non-human reads at one locus | Strong signal of contamination pile-up or a genuine non-human insertion. Cross-reference with the variant's `DKA_NHF` — if near 1.0, the variant is likely a contamination artifact. |
| Scattered non-human reads across many loci | Low-level sample contamination. Individual variants may still be genuine if their specific `DKA_NHF` is low. |

---

## Expanded Genomic Span BED File

In addition to the standard span BED described above, the pipeline writes an
**expanded span BED file** that naively extends each read's BED coordinates by
the observed soft-clip lengths.  The expanded coordinates hypothesize the
genomic span as if all soft-clipped bases were aligned contiguously to the
reference at the mapped location:

```
expanded_start = max(0, reference_start - softclip_left)
expanded_end   = reference_end + softclip_right
```

**The expanded spans are for visualization only** — they do not represent
verified reference alignments.  The actual mapped region is preserved in
the `aligned_start` and `aligned_end` columns.

### Scientific Rationale

Contamination and library-prep chimeras often manifest as reads with partial
human alignments and substantial soft-clipped ends representing non-human
(e.g. bacterial) sequence.  Visually, clusters of soft-clipped non-human
reads with congruent expanded spans strongly suggest local contamination,
integration, or library chimera events.  The expanded span provides a
"maximum hypothetical coverage" window to focus curation or pileup analysis.

### File naming

Written alongside the standard span BED when `--kraken2-db` is provided
(unless `--no-expanded-bed` is specified):

- Default: derived from `--output` (e.g. `my_trio.annotated.kraken2_spans_expanded.bed.gz`)
- Disable with `--no-expanded-bed`

### Columns

The expanded BED contains all 15 columns from the standard span BED, plus
two additional columns referencing the original mapped coordinates.

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `chrom` | String | Reference contig name. |
| 2 | `start` | Integer | 0-based **expanded** start: `max(0, reference_start - softclip_left)`. |
| 3 | `end` | Integer | 0-based exclusive **expanded** end: `reference_end + softclip_right`. May exceed chromosome length. |
| 4 | `taxon_name` | String | Scientific name (underscored) from Kraken2 LCA assignment. |
| 5 | `domain` | String | Domain classification. |
| 6 | `guard_status` | String | Human homology guard status. |
| 7 | `is_nonhuman` | String | `true` or `false`. |
| 8 | `read_name` | String | SAM QNAME. |
| 9 | `variant` | String | Comma-separated variant key(s). |
| 10 | `read_set` | String | `DKA` or `DKU`. |
| 11 | `mapq` | Integer | Mapping quality. |
| 12 | `softclip_left` | Integer | Left soft-clip length. |
| 13 | `softclip_right` | Integer | Right soft-clip length. |
| 14 | `is_split` | String | `true` if the read has an SA tag. |
| 15 | `is_supplementary` | String | `true` if this record is supplementary. |
| 16 | `aligned_start` | Integer | Original 0-based aligned start (`reference_start`). |
| 17 | `aligned_end` | Integer | Original 0-based exclusive aligned end (`reference_end`). |

### Example

```bed
#chrom	start	end	taxon_name	domain	guard_status	is_nonhuman	read_name	variant	read_set	mapq	softclip_left	softclip_right	is_split	is_supplementary	aligned_start	aligned_end
chr1	99980	100155	Escherichia_coli	Bacteria	PASS	true	read001	chr1:100050:A:T	DKA	60	20	5	false	false	100000	100150
```

Here the original aligned span is `chr1:100000–100150` (columns 16–17), and
the expanded span extends 20 bp left (soft-clip) and 5 bp right.  Compare
with the standard span BED row for the same read:

```bed
chr1	100000	100150	Escherichia_coli	Bacteria	PASS	true	read001	chr1:100050:A:T	DKA	60	20	5	false	false
```

### Comparing Standard and Expanded BED Tracks

Loading both BED tracks in IGV (or a similar genome browser) enables
powerful visual contamination auditing:

| Pattern | Standard BED | Expanded BED | Interpretation |
|---------|-------------|-------------|----------------|
| Clustered non-human soft-clipped reads | Small, congruent aligned regions | Large, consistently expanded intervals spanning a locus | Possible library/prep contamination or integration breakpoint |
| Fully mapped non-human read | Standard and expanded spans nearly identical | Standard and expanded spans nearly identical | Likely genuine non-human sequence mapping to conserved/low-complexity region |
| Mixed/ambiguous read | Small aligned span | Expanded covers ambiguous region | Carefully review; may indicate integration or clipped artifact |
| Split read with large clips | Two small intervals on different chroms | Each interval extended by its clips | Cross-chromosome chimera; if congruent across reads, likely systematic contamination |

---

## Why Kraken2 Is Well Suited to This Task

| Property | Benefit |
|----------|---------|
| **K-mer–based, alignment-free** | No need to align reads to a non-human reference; classification runs in seconds even for thousands of reads |
| **Taxonomic LCA over the full tree** | Correctly identifies reads from any bacterium, archaeon, fungus, or protist — not just species explicitly in the database — a read from an unknown strain will be classified at the correct genus or family |
| **Per-read k-mer detail output** | Enables the human homology guard: per-read k-mer votes reveal when a classified read also matches human sequence |
| **PrackenDB coverage** | One genome per species across all NCBI reference bacteria, archaea, protists, fungi, human, and viruses — captures a broad contamination landscape |
| **Confidence threshold** | `--kraken2-confidence` allows tuning sensitivity vs. specificity without rerunning the database build |
| **Scalable** | Multi-threaded with `--threads`; only informative reads are classified, so runtimes are proportional to variant evidence, not total sequencing depth |

---

## Configuration Reference

| CLI argument | Default | Effect |
|---|---|---|
| `--kraken2-db` | *(disabled)* | Path to the Kraken2 database directory; enables non-human fraction annotations (DKU_BF/DKA_BF, DKU_AF/DKA_AF, DKU_FF/DKA_FF, DKU_PF/DKA_PF, DKU_VF/DKA_VF, DKU_UCF/DKA_UCF, DKU_NHF/DKA_NHF) in VCF mode |
| `--kraken2-confidence` | `0.0` | LCA confidence threshold (0.0–1.0); higher values reduce sensitivity, increase specificity |
| `--kraken2-memory-mapping` | `false` | Pass `--memory-mapping` to Kraken2 so the database index is memory-mapped rather than loaded into RAM (much lower resident memory, slower classification) |
| `--kraken2-read-detail` | *(auto-derived)* | Output path for the per-read classification detail BED file. Auto-derived from `--output` when `--kraken2-db` is provided (e.g. `my_trio.annotated.kraken2_reads.bed.gz`). |
| `--kraken2-span-bed` | *(auto-derived)* | Output path for the species-annotated genomic span BED file. Auto-derived from `--output` when `--kraken2-db` is provided (e.g. `my_trio.annotated.kraken2_spans.bed.gz`). |
| `--no-expanded-bed` | `false` | When set, disables generation of the expanded span BED file. By default both standard and expanded span BEDs are produced. |

See [Kraken2 Database Setup Helper](../README.md#kraken2-database-setup-helper)
for instructions on downloading PrackenDB.

See the [Kraken2 manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
for full documentation on the confidence parameter, database building, and
output formats.
