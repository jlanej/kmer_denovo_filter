# Kraken2 Non-Human Content Detection

This document describes how `kmer-denovo` uses [Kraken2](https://github.com/DerrickWood/kraken2)
to detect non-human content — bacteria, archaea, fungi, protists, viruses,
and synthetic sequencing vectors — in the child's sequencing reads, and why
Kraken2's k-mer–based classification approach is well suited to that goal.

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

Kraken2 uses a **k-mer–based taxonomic classification** algorithm
([Wood & Salzberg, 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46);
[Wood et al., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)):

1. **Database construction** — Every k-mer (default k=35 for Kraken2) from all
   reference sequences in the database is mapped to its *Lowest Common Ancestor*
   (LCA) in the NCBI taxonomy tree. If a k-mer appears in genomes from multiple
   species, its LCA is promoted upward in the tree to the most specific node
   that covers all contributing species.

2. **Per-read classification** — For each input read, Kraken2 extracts every
   k-mer in a sliding window and looks up each k-mer in the database. Each
   k-mer lookup returns the taxid at whose LCA that k-mer was stored. All
   retrieved taxids are fed into an LCA vote: the read's final classification
   is the deepest taxon in the NCBI tree that is consistent with a sufficient
   fraction of the k-mer votes.

3. **Confidence threshold** — The `--confidence` parameter (exposed as
   `--kraken2-confidence`, default `0.0`) sets a minimum fraction of k-mers
   that must vote at or below the assigned clade. Higher values give more
   specific but potentially fewer classifications; lower values are more
   sensitive. A value of `0.0` assigns the LCA of all k-mer votes regardless
   of consistency.

4. **Per-read output** — Kraken2 emits one output line per read:

   ```
   C/U  read_name  taxid  length  kmer_detail_string
   ```

   - `C` = classified, `U` = unclassified
   - `taxid` = NCBI taxonomy ID of the LCA classification
   - `kmer_detail_string` = space-separated `taxid:count` tokens showing
     how many k-mers voted for each taxid (paired-end reads use `|:|` as a
     mate delimiter)

---

## The PrackenDB Reference Database

`kmer-denovo` downloads **PrackenDB** — a curated, pre-built Kraken2 database
published by the Kraken2 project
([CCB JHU downloads](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)).

As of January 2026, PrackenDB contains all NCBI reference assemblies (GenBank
and RefSeq) of bacteria, archaea, protists, and fungi as of October 7, 2025.
It also includes the human genome, RefSeq viral genomes, and UniVec Core. A
key difference from other Kraken2 databases is that PrackenDB has only a single
reference genome per species (with a couple of exceptions such as normal and
pathogenic *E. coli*), which is useful for methods that count k-mers per
species.

**Why one genome per species matters**: Because each species contributes exactly
one reference, a k-mer that appears in multiple species is LCA-elevated to a
genus or family node — not to an unrelated lineage. This preserves the
specificity of classification while avoiding inflation from redundant genomes.
It also makes k-mer counting per species unambiguous.

**Taxonomy files**: PrackenDB includes `taxonomy/nodes.dmp` (used for
lineage-aware classification) and `taxonomy/names.dmp` (used to map taxonomy
IDs to scientific names in the per-read detail BED file).  Both are validated
by the download script.

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

```python
kr = Kraken2Runner(kraken2_db, confidence=confidence, threads=threads)
result = kr.classify_sequences(sequences, tmpdir=tmpdir)
```

### Step 3 — Lineage-aware multi-domain classification

Kraken2 assigns each read a single taxid. The pipeline classifies each read
into one or more biological domains by traversing the NCBI taxonomy tree loaded
from `taxonomy/nodes.dmp` (or `nodes.dmp` at the database root for PrackenDB)
in the database directory.

The following domain-specific taxid sets are computed:

| Domain | Root taxid | Description |
|--------|-----------|-------------|
| **Bacteria** | 2 | All descendants of the Bacteria domain |
| **Archaea** | 2157 | All descendants of the Archaea domain |
| **Fungi** | 4751 | All descendants of the Fungi kingdom |
| **Protist** | (computed) | Eukaryota (2759) descendants **minus** Metazoa (33208), Fungi (4751), and Viridiplantae (33090) descendants |
| **Viruses** | 10239 | All descendants of the Viruses superkingdom (PrackenDB includes RefSeq viral genomes) |
| **UniVec Core** | 81077 | Synthetic sequencing-vector and adapter sequences — **excluded** from non-human counts (see below) |

```python
taxid_sets = Kraken2Runner._load_all_taxid_sets(db_path)
# taxid_sets["bacterial"]   = set of ALL taxids descending from taxid 2
# taxid_sets["archaeal"]    = set of ALL taxids descending from taxid 2157
# taxid_sets["fungal"]      = set of ALL taxids descending from taxid 4751
# taxid_sets["protist"]     = eukaryota - metazoa - fungi - viridiplantae
# taxid_sets["viral"]       = set of ALL taxids descending from taxid 10239
# taxid_sets["univec_core"] = set of ALL taxids descending from taxid 81077
```

This lineage-aware check correctly classifies reads assigned to a specific
genus, family, or order — not just reads whose LCA happens to be the domain
root. Without this tree traversal, reads assigned to *Escherichia coli*
(taxid 562) or *Staphylococcus aureus* (taxid 1280) would be missed because
those taxids are not equal to 2.

If `taxonomy/nodes.dmp` is missing or unreadable, `Kraken2Runner` logs a
warning and falls back to exact taxid matching only (taxid == 2 for bacteria,
taxid == 2157 for archaea, etc.). This is a less sensitive fallback: reads
assigned to specific species below the domain root will be missed. The download
script warns when `nodes.dmp` is absent after extraction; PrackenDB does
include this file.

### Step 4 — Human homology guard

Some non-human k-mers share sequence with the human genome (e.g. highly
conserved ribosomal sequences, mobile genetic elements, or horizontal gene
transfer events). A read that contains such shared k-mers could be assigned a
non-human LCA by Kraken2 even though it actually originated from human DNA.

To reduce false flagging, `kmer-denovo` applies a **human homology guard** to
**all** non-human categories (bacterial, archaeal, fungal, protist, viral,
UniVec Core, and consolidated non-human):

```python
# kmer_taxids = taxids voting in the per-read kmer_detail_string
has_human_kmer = _HUMAN_TAXID in kmer_taxids

# Applied to EVERY non-human category:
if has_human_kmer:
    is_bacterial = False
    is_archaeal = False
    is_fungal = False
    is_protist = False
    is_viral = False
    is_univec_core = False
    is_nonhuman = False
```

If Kraken2's per-read k-mer detail string includes any k-mer that voted for
human (taxid 9606), the read is conservatively excluded from **every**
non-human numerator. This means:

- A read assigned to *Bacteria* LCA but with some human k-mer evidence →
  **not counted as bacterial or non-human**
- A read assigned to *Archaea* LCA with some human k-mer evidence →
  **not counted as archaeal or non-human**
- A read assigned to a virus with some human k-mer evidence →
  **not counted as viral or non-human**
- A read assigned to *Bacteria* LCA with no human k-mer evidence →
  **counted as bacterial and non-human**

This is deliberately conservative: it may slightly undercount non-human reads
that happen to contain a human-matching k-mer, but it avoids over-flagging
human reads with non-human-like k-mers as contamination.

#### Viral reads and human DNA integration

Viruses receive the same human homology guard as all other domains, but the
guard is **especially important** for viral reads because some viruses can
integrate into or co-evolve with the human genome:

- **Endogenous retroviruses (ERVs)** — ERV sequences make up ~8% of the human
  genome. Reads from known ERV loci are already covered by the human reference
  and will be classified as human (not viral) by Kraken2 without any special
  handling. Exogenous retroviruses or novel ERV insertions may share k-mers
  with both viral references and the human reference, making the human homology
  guard essential.
- **HBV and HPV** — Hepatitis B virus and human papillomavirus can integrate
  into host chromosomes. A read spanning an integration junction would contain
  both viral and human k-mers, and the human homology guard conservatively
  excludes it from the viral count.
- **UniVec Core** — PrackenDB includes UniVec Core (sequencing vector and
  adapter sequences, taxid 81077). These synthetic constructs are handled in
  two layers: (1) the human homology guard excludes any UniVec-classified read
  that also has human k-mer evidence, and (2) UniVec Core reads are
  *unconditionally* excluded from the consolidated non-human fraction (NHF),
  because they are artificial sequences, not biological organisms, and their
  k-mers can overlap with real human genomic sequence.  See
  [Step 5](#step-5--conservative-non-human-fraction-nhf) for details.

In practice, reads from stably integrated viral sequences are expected to
produce human k-mer evidence and be excluded from the viral count, meaning
**DKU_VF reflects only reads from exogenous, non-integrated viral contamination**.
This is the conservative behavior intended by the design.

### Step 5 — Conservative non-human fraction (NHF)

In addition to domain-specific fractions, the pipeline computes a consolidated
**non-human fraction** (DKU_NHF / DKA_NHF). A read is counted as "non-human"
only if:

1. It is classified (not unclassified)
2. Its assigned taxid is **not** on the human lineage (the path from
   taxid 9606 up to root) — this excludes ambiguous ranks like Eukaryota
   (2759), Metazoa (33208), or root (1) where the read could plausibly
   be human
3. Its assigned taxid is **not** a descendant of human (9606) — this
   excludes human subspecies and populations
4. Its assigned taxid is **not** under UniVec Core (taxid 81077) — synthetic
   vector/adapter sequences are unconditionally excluded from NHF because
   they are artificial constructs and may share k-mers with human DNA,
   meaning a human read misclassified as UniVec Core would otherwise produce
   a false positive
5. It has **no human k-mer evidence** in the k-mer detail string (the
   human homology guard)

This conservative definition means:

- A read classified as *E. coli* (562) with no human k-mers → **counted as non-human** ✓
- A read classified as *Bacteria* (2) with no human k-mers → **counted as non-human** ✓
- A read classified as *Eukaryota* (2759) → **not counted as non-human** (ambiguous ancestor of human)
- A read classified as *Metazoa* (33208) → **not counted as non-human** (ancestor of human)
- A read classified as *Drosophila melanogaster* (7227) → **counted as non-human** ✓ (not in human lineage)
- A read classified as *Homo sapiens* (9606) → **not counted as non-human**
- A read classified as UniVec Core (81077) → **not counted as non-human** (synthetic construct)

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
| `assigned_taxon` | String | Scientific name from `names.dmp` (spaces → underscores). `.` if unclassified or name unavailable. |
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
| `--kraken2-read-detail` | *(auto-derived)* | Output path for the per-read classification detail BED file. Auto-derived from `--output` when `--kraken2-db` is provided (e.g. `my_trio.annotated.kraken2_reads.bed.gz`). |
| `--kraken2-span-bed` | *(auto-derived)* | Output path for the species-annotated genomic span BED file. Auto-derived from `--output` when `--kraken2-db` is provided (e.g. `my_trio.annotated.kraken2_spans.bed.gz`). |
| `--no-expanded-bed` | `false` | When set, disables generation of the expanded span BED file. By default both standard and expanded span BEDs are produced. |

See [Kraken2 Database Setup Helper](../README.md#kraken2-database-setup-helper)
for instructions on downloading PrackenDB.

See the [Kraken2 manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
for full documentation on the confidence parameter, database building, and
output formats.
