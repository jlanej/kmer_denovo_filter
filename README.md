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
   proband-unique k-mers. Cluster nearby reads (within 500 bp) into
   candidate genomic regions.

5. **Output** – Write a BED file of candidate regions, an informative-reads
   BAM, a metrics JSON file, and a human-readable summary.

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

This produces four output files (see [Discovery Mode Output](#discovery-mode-output)):

* `discovery_output.bed`
* `discovery_output.informative.bam`
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
| `--min-mapq` | 20 | Minimum mapping quality for child reads |
| `--threads` / `-t` | 4 | Number of threads for Jellyfish |
| `--debug-kmers` | false | Enable per-variant debug output |
| **VCF mode** | | |
| `--vcf` | – | Input VCF with candidate variants (activates VCF mode) |
| `--output` / `-o` | – | Output annotated VCF (required with `--vcf`) |
| `--proband-id` | – | Sample ID of the proband in the VCF. When provided and matching a VCF sample, annotations are written as FORMAT fields; otherwise INFO fields |
| `--informative-reads` | – | Output BAM of reads carrying child-unique k-mers (tagged with `DV`) |
| `--metrics` | – | Output summary metrics JSON file |
| `--summary` | – | Output human-readable summary text file |
| **Discovery mode** | | |
| `--out-prefix` | – | Output prefix for discovery mode files (activates discovery mode) |
| `--ref-jf` | – | Precomputed Jellyfish reference index; defaults to `[ref-fasta].k[kmer-size].jf` |
| `--min-child-count` | 3 | Minimum k-mer occurrences in the child to be considered a candidate |
| `--cluster-distance` | 500 | Maximum gap (bp) for merging adjacent regions |
| `--min-supporting-reads` | 1 | Minimum number of supporting reads per region |
| `--min-distinct-kmers` | 1 | Minimum number of distinct proband-unique k-mers per region |
| `--parent-max-count` | 0 | Maximum k-mer count in a parent before the k-mer is considered parental; k-mers with count > this value in either parent are removed |
| `--candidate-summary` | – | Path to a VCF-mode `summary.txt` for candidate comparison. High-quality *de novos* (DKA\_DKT > 0.25, DKA > 10) are checked against discovered regions |

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

Discovery mode always produces four files based on `--out-prefix`:

#### BED file (`{prefix}.bed`)

Tab-delimited file with candidate genomic regions. Each row represents a
cluster of reads carrying proband-unique k-mers:

| Column | Description |
|---|---|
| chrom | Chromosome |
| start | 0-based start coordinate |
| end | End coordinate (exclusive) |
| read_count | Number of unique reads with proband-unique k-mers in this region |
| kmer_count | Number of distinct proband-unique k-mers in this region |

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
  "regions": [
    {
      "chrom": "chr1",
      "start": 1000,
      "end": 5500,
      "size": 4500,
      "reads": 12,
      "unique_kmers": 35
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
* Per-region results table with coordinates, size, read count, and k-mer count
* Candidate comparison results (when `--candidate-summary` is provided)

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
