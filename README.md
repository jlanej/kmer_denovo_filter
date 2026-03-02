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

2. **Scan parents** – All child k-mers are collected into a single set.
   Each parent's entire BAM/CRAM is streamed through
   [Jellyfish](https://github.com/gmarcais/Jellyfish) `count` with the
   `--if` filter so only the child k-mers are tracked. No mapping-quality
   filter is applied to parent reads.

3. **Count proband-unique reads** – For each variant, a child read is counted
   as *unique* when at least one of its variant-spanning k-mers is absent
   from both parents.

### VCF-Free Discovery Mode

1. **Index reference** – Build (or reuse) a Jellyfish k-mer index of the
   reference genome so that reference k-mers can be subtracted later.

2. **Extract & filter child k-mers** – Count all k-mers in the child
   BAM/CRAM with Jellyfish. Keep only k-mers with at least
   `--min-child-count` occurrences (default 3), then subtract any k-mers
   present in the reference. The remaining k-mers are non-reference child
   k-mers.

3. **Filter against parents** – Stream each parent's BAM/CRAM through
   Jellyfish to count the non-reference child k-mers. Remove any k-mer
   found in either parent. The survivors are *proband-unique* k-mers.

4. **Anchor & cluster** – Scan the child BAM/CRAM for reads carrying
   proband-unique k-mers. Only reads with at least
   `--min-distinct-kmers-per-read` (default k/4) distinct proband-unique
   k-mers are retained. Cluster nearby reads (within
   `--cluster-distance` bp) into candidate genomic regions. Regions are
   then filtered by `--min-supporting-reads` and `--min-distinct-kmers`.

5. **Output** – Write a BED file of candidate regions, k-mer coverage
   bedGraph, read coverage BED, informative-reads BAM, SV breakpoints
   BEDPE, a metrics JSON file, and a human-readable summary.

## Prerequisites

* Python ≥ 3.9
* [samtools](https://www.htslib.org/) on `PATH`
* [Jellyfish ≥ 2](https://github.com/gmarcais/Jellyfish) on `PATH`

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
| `--threads` / `-t` | 4 | Number of threads for Jellyfish |
| `--debug-kmers` | false | Enable per-variant debug output |
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

### VCF Mode Output

The output VCF is bgzipped (`.vcf.gz`) with a tabix index (`.vcf.gz.tbi`),
and annotated with the following fields:

* **DKU** – Number of child reads with at least one variant-spanning k-mer
  unique to the child (absent from both parents).
* **DKT** – Total child reads with variant-spanning k-mers.
* **DKA** – Number of child reads with at least one unique k-mer that also
  exactly supports the candidate allele.
* **DKU_DKT** – Proportion of child reads with unique k-mers (DKU / DKT).
  A value of 1.0 means all reads spanning the variant carry child-unique
  k-mers. When DKT is 0 the value is 0.0.
* **DKA_DKT** – Proportion of child reads with unique allele-supporting
  k-mers (DKA / DKT). When DKT is 0 the value is 0.0.
* **MAX_PKC** / **AVG_PKC** / **MIN_PKC** – Maximum, average (2 dp), and
  minimum k-mer count among variant-spanning k-mers found in the parents.
* **MAX_PKC_ALT** / **AVG_PKC_ALT** / **MIN_PKC_ALT** – Same as above but
  restricted to k-mers from reads that directly support the alternate allele.

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
file listing linked SV breakpoint pairs identified from split-read and
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

Discovery mode applies filters at three levels, in this order:

1. **K-mer–level** (Modules 1–2):
   * `--min-child-count` — Minimum k-mer occurrences in the child
     BAM/CRAM (default 3). K-mers below this count are discarded.
   * Reference subtraction — K-mers present in the reference genome are
     removed.
   * `--parent-max-count` — K-mers with count > this value in either
     parent are removed (default 0, meaning any parental occurrence
     removes the k-mer).

2. **Read-level** (Module 3, anchoring):
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

## Docker

A Docker image is published to GitHub Container Registry on every push to
`main`:

```bash
# VCF mode
docker run --rm -v $PWD:/data ghcr.io/jlanej/kmer_denovo_filter:latest \
  --child /data/child.bam \
  --mother /data/mother.bam \
  --father /data/father.bam \
  --vcf /data/candidates.vcf \
  --output /data/annotated.vcf \
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
Docker image directly and run it on HPC clusters without root access.

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

### Example SLURM batch script

```bash
#!/bin/bash
#SBATCH --job-name=kmer_denovo
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=kmer_denovo_%j.log

module load apptainer   # or: module load singularity

SIF=/path/to/kmer_denovo.sif

# VCF mode
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
  --threads ${SLURM_CPUS_PER_TASK}

# Discovery mode (uncomment to use instead)
# apptainer exec --bind /data,/scratch "$SIF" kmer-denovo \
#   --child   /data/trio/child.bam \
#   --mother  /data/trio/mother.bam \
#   --father  /data/trio/father.bam \
#   --ref-fasta /data/reference/GRCh38.fa \
#   --out-prefix /scratch/discovery_output \
#   --min-child-count 3 \
#   --kmer-size 31 \
#   --threads ${SLURM_CPUS_PER_TASK}
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
