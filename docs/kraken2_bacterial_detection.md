# Kraken2 Bacterial and Non-Human Content Detection

This document describes how `kmer-denovo` uses [Kraken2](https://github.com/DerrickWood/kraken2)
to detect bacterial and non-human content in the child's sequencing reads, and
why Kraken2's k-mer–based classification approach is well suited to that goal.

---

## Why Bacterial Content Detection Matters

`kmer-denovo` identifies *de novo* variants by looking for k-mers present in a
child but absent from both parents. When the child's sample contains bacterial
contamination, sequencing reads from bacteria can carry k-mers that are truly
absent from parental genomes — not because a *de novo* mutation occurred, but
simply because the bacterial sequence has no counterpart in the human reference
or in the parents. Without a contamination check, these reads would be
indistinguishable from genuine *de novo* signals.

Reads flagged as bacterial can be used to compute **DKU_BF** and **DKA_BF**
(bacterial-fraction annotations), which indicate what proportion of the
informative reads supporting a candidate variant appear to derive from a
bacterial organism rather than from the human child genome. A high bacterial
fraction is a strong indicator of a false-positive *de novo* call.

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
PrackenDB contains exactly one NCBI reference genome per species, covering:

- All NCBI reference assemblies: bacteria, archaea, protists, fungi
- The human genome (GRCh38)
- RefSeq viral genomes
- UniVec Core (sequencing adapter / vector contamination sequences)

**Why one genome per species matters**: Because each species contributes exactly
one reference, a k-mer that appears in multiple bacterial species is LCA-elevated
to a genus or family node — not to an unrelated non-bacterial lineage. This
preserves the specificity of bacterial classification while avoiding inflation
from redundant genomes. It also makes k-mer counting per species unambiguous.

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
the bacterial-fraction annotations are computed on exactly the reads contributing
to each variant's evidence.

```python
kr = Kraken2Runner(kraken2_db, confidence=confidence, threads=threads)
result = kr.classify_sequences(sequences, tmpdir=tmpdir)
```

### Step 3 — Lineage-aware bacterial classification

Kraken2 assigns each read a single taxid. A read is counted as **bacterial** if
that taxid falls anywhere within the Bacteria domain (NCBI taxid 2). This is
determined by traversing the NCBI taxonomy tree loaded from
`taxonomy/nodes.dmp` in the database directory:

```python
bacterial_taxids = Kraken2Runner._load_bacterial_taxids(db_path)
# bacterial_taxids is the set of ALL taxids that descend from taxid 2
is_bacterial = taxid in bacterial_taxids
```

This lineage-aware check correctly classifies reads assigned to a specific
bacterial genus, family, or order — not just reads whose LCA happens to be the
root Bacteria node (taxid 2). Without this tree traversal, reads assigned to
*Escherichia coli* (taxid 562) or *Staphylococcus aureus* (taxid 1280) would
be missed, because those taxids are not equal to 2.

If `taxonomy/nodes.dmp` is missing or unreadable, `Kraken2Runner` logs a
warning and falls back to exact taxid==2 matching only. This is a less sensitive
fallback: reads assigned to specific bacterial species below the root Bacteria
node will be missed. The download script warns when `nodes.dmp` is absent after
extraction; PrackenDB does include this file.

### Step 4 — Human homology guard

Some bacterial k-mers share sequence with the human genome (e.g. highly
conserved ribosomal sequences, mobile genetic elements, or horizontal gene
transfer events). A read that contains such shared k-mers could be assigned a
bacterial LCA by Kraken2 even though it actually originated from human DNA.

To reduce false flagging, `kmer-denovo` applies a **human homology guard**:

```python
# kmer_taxids = taxids voting in the per-read kmer_detail_string
if is_bacterial and _HUMAN_TAXID in kmer_taxids:
    is_bacterial = False
```

If Kraken2's per-read k-mer detail string includes any k-mer that voted for
human (taxid 9606), the read is conservatively excluded from the bacterial
numerator. This means:

- A read assigned to *Bacteria* LCA but with some human k-mer evidence →
  **not counted as bacterial**
- A read assigned to *Bacteria* LCA with no human k-mer evidence →
  **counted as bacterial**

This is deliberately conservative: it may slightly undercount bacterial reads
that happen to contain a human-matching k-mer, but it avoids over-flagging
human reads with bacterial-like k-mers as contamination.

---

## Output Annotations

The bacterial classification results are added to the output VCF as two
per-variant fields:

| Field | Description |
|-------|-------------|
| **DKU_BF** | Fraction of DKU fragments classified as bacterial by Kraken2. Denominator = DKU (both are fragment-based, i.e. unique read names). Always in [0.0, 1.0]. |
| **DKA_BF** | Fraction of DKA fragments (DKU fragments that also directly support the alternate allele) classified as bacterial. DKA fragments are always a subset of DKU. Always in [0.0, 1.0]. |

These are per-variant fractions: each fraction is computed from the intersection
of that variant's informative reads with the global set of bacterial read names
returned by Kraken2.

**Interpretation**:

- `DKU_BF` close to `1.0` — essentially all evidence for this variant comes
  from reads classified as bacterial; strong indicator of contamination artifact
- `DKA_BF` close to `1.0` — reads that directly support the alternate allele
  sequence are predominantly bacterial; high-confidence contamination flag
- `DKU_BF` and `DKA_BF` both near `0.0` — no detectable bacterial content
  among the supporting reads; the candidate variant is more likely genuine

---

## Why Kraken2 Is Well Suited to This Task

| Property | Benefit |
|----------|---------|
| **K-mer–based, alignment-free** | No need to align reads to a bacterial reference; classification runs in seconds even for thousands of reads |
| **Taxonomic LCA over the full tree** | Correctly identifies reads from any bacterium, not just species explicitly in the database — a read from an unknown strain will be classified at the correct genus or family |
| **Per-read k-mer detail output** | Enables the human homology guard: per-read k-mer votes reveal when a classified read also matches human sequence |
| **PrackenDB coverage** | One genome per species across all NCBI reference bacteria, archaea, protists, fungi, human, and viruses — captures a broad contamination landscape |
| **Confidence threshold** | `--kraken2-confidence` allows tuning sensitivity vs. specificity without rerunning the database build |
| **Scalable** | Multi-threaded with `--threads`; only informative reads are classified, so runtimes are proportional to variant evidence, not total sequencing depth |

---

## Configuration Reference

| CLI argument | Default | Effect |
|---|---|---|
| `--kraken2-db` | *(disabled)* | Path to the Kraken2 database directory; enables DKU_BF/DKA_BF annotations in VCF mode |
| `--kraken2-confidence` | `0.0` | LCA confidence threshold (0.0–1.0); higher values reduce sensitivity, increase specificity |

See [Kraken2 Database Setup Helper](../README.md#kraken2-database-setup-helper)
for instructions on downloading PrackenDB.

See the [Kraken2 manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
for full documentation on the confidence parameter, database building, and
output formats.
