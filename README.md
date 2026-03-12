# kmer_denovo_filter

De novo variant curation using k-mer analysis.

## Overview

`kmer-denovo` identifies candidate *de novo* variants in a child by
comparing k-mers from the child's sequencing reads against both parents.
K-mers present in the child but absent from both parents signal potential
*de novo* mutations.

The tool supports two modes:

* **VCF mode** – Given a VCF of candidate variants, annotate each variant
  with k-mer–based evidence of *de novo* origin.
* **VCF-free discovery mode** – Without any input VCF, scan the child's
  entire genome for regions harboring proband-unique k-mers and report
  candidate *de novo* regions.

## Algorithm

### VCF Mode

1. **Extract child k-mers** – For each candidate variant, child reads
   overlapping the position are fetched and only k-mers whose genomic span
   includes the variant position are kept. K-mers are canonicalized
   (lexicographically smaller of k-mer and its reverse complement).
   For multiallelic variants, only the first ALT allele is evaluated.

2. **Scan parents** – All child k-mers are collected into a single set.
   Each parent's entire BAM/CRAM is streamed through
   [Jellyfish](https://github.com/gmarcais/Jellyfish) `count` with the
   `--if` filter so only the child k-mers are tracked. No mapping-quality
   filter is applied to parent reads. K-mer counts from both parents are
   summed into a single combined count.

3. **Count proband-unique fragments** – For each variant, a child fragment
   is counted as *unique* when at least one of its variant-spanning k-mers
   is absent from both parents. DKU, DKT, and DKA all count unique
   fragment names: when both mates of a paired-end read span the variant,
   they share a read name and are counted as one fragment.

### VCF-Free Discovery Mode

1. **Index reference** – Build (or reuse) a Jellyfish k-mer index of the
   reference genome so that reference k-mers can be subtracted later.

2. **Extract & filter child k-mers** – Count all k-mers in the child
   BAM/CRAM with Jellyfish. Keep only k-mers with at least
   `--min-child-count` occurrences (default 3), then subtract any k-mers
   present in the reference. The remaining k-mers are non-reference child
   k-mers. Intermediate data (the child Jellyfish index and candidates
   FASTA) is removed as soon as each step completes to keep disk and
   memory usage low.

3. **Filter against parents** – Stream each parent's BAM/CRAM through
   Jellyfish to count the non-reference child k-mers. Remove any k-mer
   found in either parent. The survivors are *proband-unique* k-mers.
   Each parent's Jellyfish index is removed immediately after it is
   queried.

4. **Build proband-unique index** – Build a Jellyfish hash index from the
   proband-unique k-mers FASTA. This index is memory-mapped and shared
   across workers via the OS page cache, so N workers ≈ 1× the index
   size in RAM (typically 2–10 GB for WGS data).

5. **Anchor & cluster** – Scan each child read with a sliding k-mer
   window. Each k-mer is canonicalized and queried against the
   proband-unique Jellyfish index via a long-lived ``jellyfish query``
   subprocess. Reads with at least `--min-distinct-kmers-per-read`
   (default k/4) distinct proband-unique k-mers are retained. Cluster
   nearby reads (within `--cluster-distance` bp) into candidate genomic
   regions. Regions are then filtered by `--min-supporting-reads` and
   `--min-distinct-kmers`. The number of parallel workers is dynamically
   capped based on available system memory.

6. **Output** – Write a BED file of candidate regions, k-mer coverage
   bedGraph, read coverage BED, informative-reads BAM, SV breakpoints
   BEDPE, a metrics JSON file, and a human-readable summary.

## Prerequisites

* Python ≥ 3.9
* [samtools](https://www.htslib.org/) on `PATH` [1]
* [Jellyfish ≥ 2](https://github.com/gmarcais/Jellyfish) on `PATH` [2]
* Optional for VCF-mode non-human fraction annotations:
  * [Kraken2](https://github.com/DerrickWood/kraken2) on `PATH` [3]
  * A Kraken2 database — the helper script downloads the pre-built
    [PrackenDB](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)
    database (see [Kraken2 Database Setup Helper](#kraken2-database-setup-helper))

## Installation

```bash
pip install .
```

## Usage

### VCF Mode

Annotate a VCF of candidate *de novo* variants:

```bash
kmer-denovo \
  --child   child.bam \
  --mother  mother.bam \
  --father  father.bam \
  --vcf     candidates.vcf \
  --output  annotated.vcf \
  --proband-id HG002 \
  --metrics summary.json \
  --summary summary.txt \
  --kmer-size 31
```

For CRAM files, supply the reference FASTA:

```bash
kmer-denovo \
  --child   child.cram \
  --mother  mother.cram \
  --father  father.cram \
  --ref-fasta GRCh38.fa \
  --vcf     candidates.vcf \
  --output  annotated.vcf \
  --proband-id HG002
```

### VCF-Free Discovery Mode

Discover candidate *de novo* regions without an input VCF. A reference
FASTA (or a precomputed Jellyfish index) is required so that reference
k-mers can be subtracted:

```bash
kmer-denovo \
  --child   child.bam \
  --mother  mother.bam \
  --father  father.bam \
  --ref-fasta GRCh38.fa \
  --out-prefix discovery_output \
  --min-child-count 3 \
  --kmer-size 31 \
  --threads 8
```

This produces seven output files (see [Discovery Mode Output](#discovery-mode-output)):

* `discovery_output.bed`
* `discovery_output.informative.bam`
* `discovery_output.sv.bedpe`
* `discovery_output.kmer_coverage.bedgraph`
* `discovery_output.read_coverage.bed`
* `discovery_output.metrics.json`
* `discovery_output.summary.txt`

To skip the reference indexing step on subsequent runs, pass a precomputed
Jellyfish index:

```bash
kmer-denovo \
  --child   child.bam \
  --mother  mother.bam \
  --father  father.bam \
  --ref-jf  GRCh38.fa.k31.jf \
  --out-prefix discovery_output
```

Optionally, compare discovered regions against high-quality candidates from
a previous VCF-mode run:

```bash
kmer-denovo \
  --child   child.bam \
  --mother  mother.bam \
  --father  father.bam \
  --ref-fasta GRCh38.fa \
  --out-prefix discovery_output \
  --candidate-summary vcf_mode_summary.txt
```

### Arguments

| Argument | Default | Description |
|---|---|---|
| **Common** | | |
| `--child` | *required* | Child BAM/CRAM file (indexed) |
| `--mother` | *required* | Mother BAM/CRAM file (indexed) |
| `--father` | *required* | Father BAM/CRAM file (indexed) |
| `--ref-fasta` / `-r` | – | Reference FASTA with `.fai` index (required for CRAM; required for discovery mode unless `--ref-jf` is provided) |
| `--kmer-size` / `-k` | 31 | K-mer size (must be odd, 3–201) |
| `--min-baseq` | 20 | Minimum base quality for read k-mers |
| `--threads` / `-t` | 4 | Number of threads for Jellyfish and parallel anchoring workers |
| `--memory` | auto | Available memory in GB. On HPC (e.g. SLURM), set this to the allocated memory so worker counts and hash sizes are tuned correctly. When omitted, auto-detected from the system |
| `--debug-kmers` | false | Enable per-variant debug output |
| `--kraken2-db` | – | Optional Kraken2 database path. In VCF mode, enables non-human fraction annotations (DKU_BF/DKA_BF, DKU_AF/DKA_AF, DKU_FF/DKA_FF, DKU_PF/DKA_PF, DKU_VF/DKA_VF, DKU_UCF/DKA_UCF, DKU_NHF/DKA_NHF). Ignored in discovery mode. **Memory note:** the standard Kraken2 DB typically needs ~50–100 GB RAM to load/classify; provision memory accordingly to avoid OOM |
| `--kraken2-confidence` | 0.0 | Kraken2 LCA confidence threshold (0.0–1.0) |
| `--kraken2-memory-mapping` | false | Passes Kraken2 `--memory-mapping` so DB files are memory-mapped from disk to reduce RAM footprint (usually slower, but helpful on RAM-constrained nodes) |
| **VCF mode** | | |
| `--vcf` | – | Input VCF with candidate variants (activates VCF mode) |
| `--output` / `-o` | – | Output annotated VCF (required with `--vcf`) |
| `--min-mapq` | 20 | Minimum mapping quality for child reads (VCF mode only; discovery mode scans all primary reads regardless of mapping quality) |
| `--proband-id` | – | Sample ID of the proband in the VCF. When provided and matching a VCF sample, annotations are written as FORMAT fields; otherwise INFO fields |
| `--informative-reads` | – | Output BAM of reads carrying child-unique k-mers (tagged with `DV`) |
| `--metrics` | – | Output summary metrics JSON file |
| `--summary` | – | Output human-readable summary text file |
| **Discovery mode** | | |
| `--out-prefix` | – | Output prefix for discovery mode files (activates discovery mode) |
| `--ref-jf` | – | Precomputed Jellyfish reference index; defaults to `[ref-fasta].k[kmer-size].jf` |
| `--min-child-count` | 3 | Minimum k-mer occurrences in the child to be considered a candidate |
| `--cluster-distance` | 500 | Maximum gap (bp) for merging adjacent regions |
| `--min-distinct-kmers-per-read` | k/4 | Minimum distinct proband-unique k-mers a read must carry to be retained. Applied before region-level and bedGraph filters (see [Filtering Flow](#discovery-mode-filtering-flow)) |
| `--min-supporting-reads` | 1 | Minimum number of supporting reads per region |
| `--min-distinct-kmers` | 1 | Minimum number of distinct proband-unique k-mers per region |
| `--min-bedgraph-reads` | 3 | Minimum distinct reads at a position for inclusion in the bedGraph and read coverage BED |
| `--parent-max-count` | 0 | Maximum k-mer count in a parent before the k-mer is considered parental; k-mers with count > this value in either parent are removed |
| `--candidate-summary` | – | Path to a VCF-mode `summary.txt` for candidate comparison. High-quality *de novos* (DKA\_DKT > 0.25, DKA > 10) are checked against discovered regions |
| `--sv-bedpe` | – | Output BEDPE file for linked SV breakpoint pairs (default: `[out-prefix].sv.bedpe`) |
| `--jf-hash-size` | auto | Initial hash size for `jellyfish count` (e.g. `2G`, `500M`). Estimated from the child BAM file size by default. A larger value avoids hash overflow (which creates multi-file indexes requiring more memory to merge/dump) |
| `--tmp-dir` | auto | Directory for temporary files (jellyfish indexes, intermediate FASTA files). Defaults to a subdirectory next to the output files. Avoid RAM-backed filesystems like tmpfs (`/tmp` on many HPC systems), as intermediate files can exceed 100 GB for WGS data |

### VCF Mode Output

The output VCF is bgzipped (`.vcf.gz`) with a tabix index (`.vcf.gz.tbi`),
and annotated with the following fields:

* **DKU** – Number of child fragments (unique read names) with at least one
  variant-spanning k-mer unique to the child (absent from both parents).
  When both mates of a paired-end read span the variant, they share a
  read name and are counted as one fragment.
* **DKT** – Total child fragments (unique read names) with variant-spanning
  k-mers.
* **DKA** – Number of child fragments with at least one unique k-mer that
  also exactly supports the first ALT allele. DKA ≤ DKU by construction.
  For multiallelic variants only the first ALT allele is evaluated.
* **DKU_DKT** – Proportion of child fragments with unique k-mers (DKU / DKT).
  A value of 1.0 means all fragments spanning the variant carry child-unique
  k-mers. When DKT is 0 the value is 0.0.
* **DKA_DKT** – Proportion of child fragments with unique allele-supporting
  k-mers (DKA / DKT). When DKT is 0 the value is 0.0.
* **MAX_PKC** / **AVG_PKC** / **MIN_PKC** – Maximum, average (2 dp), and
  minimum combined k-mer count (summed across both parents) among
  variant-spanning k-mers that were found in at least one parent.
  Returns 0 when no variant-spanning k-mers appear in the parents.
* **MAX_PKC_ALT** / **AVG_PKC_ALT** / **MIN_PKC_ALT** – Same as above but
  restricted to k-mers from reads that directly support the alternate allele.
* **DKU_BF** / **DKA_BF** *(optional; when `--kraken2-db` is provided in VCF mode)* –
  Fraction of DKU/DKA fragments classified as bacterial by Kraken2.
* **DKU_AF** / **DKA_AF** – Fraction classified as archaeal.
* **DKU_FF** / **DKA_FF** – Fraction classified as fungal.
* **DKU_PF** / **DKA_PF** – Fraction classified as protist
  (eukaryotic but not metazoan, fungal, or plant).
* **DKU_VF** / **DKA_VF** – Fraction classified as viral (RefSeq viral
  genomes are included in PrackenDB). Reads with any human k-mer evidence
  are conservatively excluded — this specifically handles viruses that can
  integrate into the human genome (e.g. endogenous retroviruses, HBV, HPV),
  ensuring integrated viral sequences are not counted as exogenous contamination.
* **DKU_UCF** / **DKA_UCF** – Fraction classified as UniVec Core (synthetic
  sequencing-vector and adapter sequences, taxid 81077). Reads with any
  human k-mer evidence are conservatively excluded. UniVec Core reads are
  **not** included in the non-human fraction (DKU_NHF/DKA_NHF) because they
  are artificial constructs, not biological contamination.
* **DKU_NHF** / **DKA_NHF** – Consolidated non-human fraction: any read
  definitively classified outside the human lineage. UniVec Core reads are
  excluded from this fraction.

  Since DKU, DKA, and the fraction denominators all count unique fragment names,
  every fraction is always in [0.0, 1.0]. To reduce over-flagging from shared
  human homology, reads with explicit human taxid k-mer evidence in Kraken2's
  per-read output are conservatively excluded from **all** non-human numerators.
  Reads classified at ambiguous taxonomic ranks that include human (e.g.
  Eukaryota, Metazoa, root) are also excluded from the non-human fraction.
  See [Kraken2 Non-Human Content Detection](docs/kraken2_bacterial_detection.md)
  for a detailed description of the classification method.

When `--proband-id` is provided and matches a sample in the input VCF,
annotations are written as **FORMAT** (per-sample) fields on that sample.
Otherwise they are written as **INFO** fields.

The optional `--metrics` JSON file provides a summary including total
variants, child-unique k-mer counts, and the number of variants with unique
reads.

The optional `--summary` text file provides a human-readable overview
including variant counts, read-support statistics, and a per-variant table
showing all annotation values and *de novo* calls.

The optional `--informative-reads` BAM file contains child reads that carry
at least one variant-spanning k-mer absent from both parents. Each read is
tagged with `DV` indicating which variant(s) it supports. The BAM is sorted
and indexed for direct visualization in IGV.

### Kraken2 Database Setup Helper

To download the pre-built [PrackenDB](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)
Kraken2 database (NCBI reference assemblies — one genome per species):

```bash
./scripts/download_kraken2_db.sh --db /path/to/kraken2_db
```

PrackenDB contains all NCBI reference assemblies (GenBank and RefSeq)
of bacteria, archaea, protists, and fungi as of October 7, 2025, plus
the human genome, RefSeq viral genomes, and UniVec Core.  A key
difference from other Kraken2 databases is that PrackenDB has only a
single reference genome per species (with a couple of exceptions such
as normal and pathogenic *E. coli*), which is useful for methods that
count k-mers per species.

> **Kraken2 k-mer size:** Kraken2 uses a k-mer length that is baked into
> the database at build time and stored in the database's `opts.k2d` file.
> The PrackenDB pre-built database uses Kraken2's default **k = 35**
> (35-mers).  This is distinct from the Jellyfish k-mer size (`--kmer-size`,
> default 31) used for de-novo variant k-mer counting.  The pipeline reads
> the database k-mer length from `opts.k2d` at startup and logs it as
> `[Kraken2] database k-mer length: 35`.

This helper validates that required Kraken2 DB files are present, including
`taxonomy/nodes.dmp` used for lineage-aware non-human classification.

The script requires only `wget` — it downloads and extracts a pre-built
database archive and does not need `k2` or `kraken2-build`.

See the [Kraken2 manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
for full database build options.

You can also run the helper inside the published container:

```bash
docker run --rm \
  -v "$PWD:/work" \
  ghcr.io/jlanej/kmer_denovo_filter:latest \
  bash /work/scripts/download_kraken2_db.sh --db /work/kraken2_db
```

Or via Apptainer on HPC (see [Running on HPC with Apptainer](#running-on-hpc-with-apptainer)):

```bash
apptainer exec kmer_denovo.sif \
  bash /app/scripts/download_kraken2_db.sh --db /scratch/kraken2_db
```

### Discovery Mode Output

Discovery mode always produces seven files based on `--out-prefix`:

#### BED file (`{prefix}.bed`)

Tab-delimited file with candidate genomic regions. Each row represents a
cluster of reads carrying proband-unique k-mers. The file begins with a
`#filters:` comment line that records the filter parameters used to
generate the file:

| Column | Description |
|---|---|
| chrom | Chromosome |
| start | 0-based start coordinate |
| end | End coordinate (exclusive) |
| read_count | Number of unique reads with proband-unique k-mers in this region |
| kmer_count | Number of distinct proband-unique k-mers in this region |
| split_reads | Number of split-read alignments (SA-tag evidence) in this region |
| discordant_pairs | Number of discordant read pairs in this region |
| max_clip_len | Maximum soft-clip length among reads in this region |
| unmapped_mates | Number of reads whose mate is unmapped |
| class | SV classification: `SV`, `AMBIGUOUS`, or `SMALL` |

#### K-mer coverage bedGraph (`{prefix}.kmer_coverage.bedgraph`)

4-column bedGraph of de novo k-mer reference coverage. For each reference
position covered by at least `--min-bedgraph-reads` distinct informative
reads, reports the total number of unique k-mer base overlaps across all
reads. Adjacent positions with identical coverage are merged into
intervals.

#### Read coverage BED (`{prefix}.read_coverage.bed`)

Per-position read support BED file. For each reference position meeting
the `--min-bedgraph-reads` threshold, reports:

| Column | Description |
|---|---|
| chrom | Chromosome |
| start | 0-based start coordinate |
| end | End coordinate (exclusive) |
| read_count | Number of distinct reads touching the position with at least one de novo k-mer |
| avg_kmers_per_read | Average number of unique k-mer overlaps per read at this position |

#### Informative BAM (`{prefix}.informative.bam`)

Sorted, indexed BAM containing all child reads with at least one
proband-unique k-mer. Each read carries a `dk:i:1` tag. Both mapped and
unmapped informative reads are included to support downstream
re-alignment analysis.

#### Metrics JSON (`{prefix}.metrics.json`)

Machine-readable pipeline statistics:

```json
{
  "mode": "discovery",
  "child_candidate_kmers": 150000,
  "non_ref_kmers": 45000,
  "proband_unique_kmers": 1200,
  "informative_reads": 350,
  "unmapped_informative_reads": 5,
  "candidate_regions": 42,
  "filters": {
    "min_distinct_kmers_per_read": 7,
    "min_supporting_reads": 1,
    "min_distinct_kmers": 1,
    "min_bedgraph_reads": 3
  },
  "regions": [
    {
      "chrom": "chr1",
      "start": 1000,
      "end": 5500,
      "size": 4500,
      "reads": 12,
      "unique_kmers": 35,
      "split_reads": 0,
      "discordant_pairs": 1,
      "max_clip_len": 50,
      "unmapped_mates": 0,
      "class": "AMBIGUOUS"
    }
  ]
}
```

When `--candidate-summary` is used, a `candidate_comparison` object is
included with capture rate and per-candidate details.

#### Summary text (`{prefix}.summary.txt`)

Human-readable overview including:

* K-mer filtering statistics (child candidates → non-reference → proband-unique)
* Region counts and informative read totals
* Region size statistics (mean, median, max)
* Per-region results table with coordinates, size, read count, k-mer count, SV annotations, and classification
* Candidate comparison results (when `--candidate-summary` is provided)

#### SV breakpoints BEDPE (`{prefix}.sv.bedpe`)

Tab-delimited [BEDPE](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
file [4] listing linked SV breakpoint pairs identified from split-read and
discordant-pair evidence across discovery regions:

| Column | Description |
|---|---|
| chrom1 | Chromosome of the first breakpoint region |
| start1 | 0-based start of the first breakpoint region |
| end1 | End of the first breakpoint region |
| chrom2 | Chromosome of the second breakpoint region |
| start2 | 0-based start of the second breakpoint region |
| end2 | End of the second breakpoint region |
| sv_id | Identifier for the SV link (e.g. `SV_1`) |
| supporting_reads | Number of reads supporting the link |
| sv_type | SV type hint: `INTRA` (intra-chromosomal) or `BND` (inter-chromosomal) |

When no linked breakpoints are found the file contains only the header line.

### Discovery Mode Filtering Flow

Discovery mode applies filters at four levels, in this order:

1. **K-mer–level**:
   * `--min-child-count` — Minimum k-mer occurrences in the child
     BAM/CRAM (default 3). K-mers below this count are discarded.
   * Reference subtraction — K-mers present in the reference genome are
     removed.
   * `--parent-max-count` — K-mers with count > this value in either
     parent are removed (default 0, meaning any parental occurrence
     removes the k-mer).

    After filtering, a Jellyfish hash index is built from the surviving
    proband-unique k-mers and shared across workers via the OS page cache.

2. **Read-level** (anchoring):
   * `--min-distinct-kmers-per-read` — A read must carry at least this
     many distinct proband-unique k-mers to be considered informative.
     The default is `k/4` (e.g. 7 for k=31).  Reads below the threshold
     are excluded from **all** downstream outputs: regions, bedGraph,
     read coverage BED, and the informative BAM's coverage signal.

3. **Region-level** (post-clustering):
   * `--min-supporting-reads` — Minimum number of reads in a region
     (default 1).
   * `--min-distinct-kmers` — Minimum number of distinct proband-unique
     k-mers in a region (default 1).

   These two filters control which regions appear in the BED file and
   metrics JSON.

4. **Position-level** (output):
   * `--min-bedgraph-reads` — Minimum number of distinct reads at a
     single reference position for inclusion in the bedGraph and read
     coverage BED (default 3).

The BED file header records the applied filter values so that output
provenance is self-documenting.

## Notes & Limitations

* **Multiallelic variants** – When a VCF record contains multiple ALT
  alleles, only the **first ALT allele** is evaluated. DKA, DKA_DKT,
  and the `_ALT` PKC metrics all reflect only the first ALT. A warning
  is logged for each multiallelic record. To evaluate all alleles,
  decompose the VCF beforehand (e.g. `bcftools norm -m-`).

* **Fragment-based counting** – DKU, DKT, and DKA all count unique
  fragment names (read names). When both mates of a paired-end read span
  the variant, they share a read name and are counted as one fragment.
  The non-human fraction metrics (DKU_BF, DKA_BF, DKU_AF, DKA_AF,
  DKU_FF, DKA_FF, DKU_PF, DKA_PF, DKU_VF, DKA_VF, DKU_UCF, DKA_UCF,
  DKU_NHF, DKA_NHF) use the same fragment-based denominator, so all
  counts and fractions are consistent and all fraction values are always
  in [0.0, 1.0].

* **PKC is a combined parental count** – MAX_PKC, AVG_PKC, and MIN_PKC
  represent the **sum** of mother and father k-mer counts for each
  variant-spanning k-mer. They do not distinguish which parent
  contributes the count. Only k-mers that were found in at least one
  parent contribute to the statistics; k-mers absent from both parents
  are excluded from the average.

* **Symbolic alleles are skipped** – Variants whose first ALT is a
  symbolic allele (`<DEL>`, `<INS>`, breakend notation, `*`) are
  skipped during k-mer extraction and receive zero-valued annotations.

## Docker

A Docker image is published to GitHub Container Registry on every push to
`main`. The image includes `samtools`, `jellyfish`, `kraken2`, and `wget`
for downloading the pre-built PrackenDB database.  The helper
script is available inside the container at `/app/scripts/download_kraken2_db.sh`:

```bash
# VCF mode
docker run --rm -v $PWD:/data ghcr.io/jlanej/kmer_denovo_filter:latest \
  --child /data/child.bam \
  --mother /data/mother.bam \
  --father /data/father.bam \
  --vcf /data/candidates.vcf \
  --output /data/annotated.vcf \
  --kraken2-db /data/kraken2_db \
  --proband-id HG002

# Discovery mode
docker run --rm -v $PWD:/data ghcr.io/jlanej/kmer_denovo_filter:latest \
  --child /data/child.bam \
  --mother /data/mother.bam \
  --father /data/father.bam \
  --ref-fasta /data/GRCh38.fa \
  --out-prefix /data/discovery_output
```

## Running on HPC with Apptainer

[Apptainer](https://apptainer.org/) (formerly Singularity) can pull the
Docker image directly and run it on HPC clusters without root access [5].

### Pull the container image

```bash
apptainer pull kmer_denovo.sif docker://ghcr.io/jlanej/kmer_denovo_filter:latest
```

### Interactive usage

```bash
# VCF mode
apptainer exec kmer_denovo.sif kmer-denovo \
  --child   child.bam \
  --mother  mother.bam \
  --father  father.bam \
  --vcf     candidates.vcf \
  --output  annotated.vcf \
  --proband-id HG002 \
  --threads 8

# Discovery mode
apptainer exec kmer_denovo.sif kmer-denovo \
  --child   child.bam \
  --mother  mother.bam \
  --father  father.bam \
  --ref-fasta GRCh38.fa \
  --out-prefix discovery_output \
  --threads 8
```

### Downloading the Kraken2 database via Apptainer

The helper script downloads the pre-built PrackenDB database:

```bash
apptainer exec --bind /scratch kmer_denovo.sif \
  bash /app/scripts/download_kraken2_db.sh \
    --db /scratch/kraken2_db
```

> **Note:** The PrackenDB download is approximately 50 GB.  Ensure your
> scratch directory and job allocation are sized accordingly.  At
> runtime, Kraken2 loads the database into memory so provision
> sufficient RAM as well.

### Example SLURM batch script

```bash
#!/bin/bash
#SBATCH --job-name=kmer_denovo
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=kmer_denovo_%j.log

module load apptainer   # or: module load singularity

SIF=/path/to/kmer_denovo.sif

# VCF mode (16 GB is typically sufficient)
apptainer exec --bind /data,/scratch "$SIF" kmer-denovo \
  --child   /data/trio/child.bam \
  --mother  /data/trio/mother.bam \
  --father  /data/trio/father.bam \
  --vcf     /data/trio/candidates.vcf \
  --output  /scratch/annotated.vcf \
  --proband-id HG002 \
  --metrics /scratch/metrics.json \
  --summary /scratch/summary.txt \
  --kmer-size 31 \
  --threads ${SLURM_CPUS_PER_TASK} \
  --memory $(( ${SLURM_MEM_PER_NODE:-32768} / 1024 ))

# Discovery mode (uncomment to use instead)
# For WGS discovery, request at least 64 GB; 128 GB recommended.
# The child k-mer counting step (Module 1) may need 80–120 GB
# for large WGS BAMs. Module 3 (anchoring) uses a disk-backed
# Jellyfish index, keeping per-worker memory low (~2–10 GB shared
# via page cache). The --tmp-dir flag defaults to a subdirectory
# next to --out-prefix; on HPC systems, point it to a fast scratch
# filesystem (avoid RAM-backed /tmp or tmpfs).
# Use --memory to tell the tool how much RAM your SLURM job has
# so that worker counts are tuned correctly (system-reported memory
# may reflect the full node, not your allocation).
# apptainer exec --bind /data,/scratch "$SIF" kmer-denovo \
#   --child   /data/trio/child.bam \
#   --mother  /data/trio/mother.bam \
#   --father  /data/trio/father.bam \
#   --ref-fasta /data/reference/GRCh38.fa \
#   --out-prefix /scratch/discovery_output \
#   --min-child-count 3 \
#   --kmer-size 31 \
#   --threads ${SLURM_CPUS_PER_TASK} \
#   --memory 128
```

> **Tip:** Use `--bind` to make host directories visible inside the
> container. Bind your input data directory and a writable scratch directory
> for output.

### Using CRAM files on HPC

When working with CRAM files, the reference FASTA must also be accessible
inside the container:

```bash
apptainer exec --bind /data,/references "$SIF" kmer-denovo \
  --child   /data/trio/child.cram \
  --mother  /data/trio/mother.cram \
  --father  /data/trio/father.cram \
  --ref-fasta /references/GRCh38.fa \
  --vcf     /data/trio/candidates.vcf \
  --output  /data/trio/annotated.vcf \
  --proband-id HG002
```

## Testing

```bash
pip install pytest
pytest
```

## References

1. Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience*. 2021;10(2):giab008. https://doi.org/10.1093/gigascience/giab008
2. Marçais G, Kingsford C. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. *Bioinformatics*. 2011;27(6):764-770. https://doi.org/10.1093/bioinformatics/btr011
3. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biology*. 2019;20:257. https://doi.org/10.1186/s13059-019-1891-0
4. Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26(6):841-842. https://doi.org/10.1093/bioinformatics/btq033
5. Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. *PLOS ONE*. 2017;12(5):e0177459. https://doi.org/10.1371/journal.pone.0177459
