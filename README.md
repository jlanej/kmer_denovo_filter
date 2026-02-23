# kmer_denovo_filter

Methods to filter for high quality de novo variants using k-mers.

## Overview

`kmer-denovo` curates candidate de novo variants by comparing k-mers from
child reads against both parents. For each variant the tool reports how many
child reads carry a variant-spanning k-mer that is absent from both parents.

### Algorithm

1. **Extract child k-mers** – For each candidate variant, child reads
   overlapping the position are fetched and only k-mers whose genomic span
   includes the variant position are kept. K-mers are canonicalized
   (lexicographically smaller of k-mer and its reverse complement).

2. **Scan parents (whole-file)** – All child k-mers are collected into a
   single set. Each parent's entire BAM/CRAM is streamed through
   [Jellyfish](https://github.com/gmarcais/Jellyfish) `count` with the
   `--if` filter so only the child k-mers are tracked. No mapping-quality
   filter is applied to parent reads.

3. **Count proband-unique reads** – For each variant, a child read is counted
   as *unique* when at least one of its variant-spanning k-mers is absent
   from both parents.

## Prerequisites

* Python ≥ 3.9
* [samtools](https://www.htslib.org/) on `PATH`
* [Jellyfish ≥ 2](https://github.com/gmarcais/Jellyfish) on `PATH`

## Installation

```bash
pip install .
```

## Usage

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

### Arguments

| Argument | Default | Description |
|---|---|---|
| `--child` | *required* | Child BAM/CRAM file (indexed) |
| `--mother` | *required* | Mother BAM/CRAM file (indexed) |
| `--father` | *required* | Father BAM/CRAM file (indexed) |
| `--ref-fasta` / `-r` | – | Reference FASTA with `.fai` index (required for CRAM) |
| `--vcf` | *required* | Input VCF with candidate variants |
| `--output` / `-o` | *required* | Output annotated VCF |
| `--proband-id` | – | Sample ID of the proband in the VCF (see [Output](#output)) |
| `--informative-reads` | – | Output BAM with reads carrying informative k-mers (for IGV) |
| `--metrics` | – | Output summary metrics JSON file |
| `--summary` | – | Output human-readable summary of variant stats and likely DNMs |
| `--kmer-size` / `-k` | 31 | K-mer size |
| `--min-baseq` | 20 | Minimum base quality for read k-mers |
| `--min-mapq` | 20 | Minimum mapping quality for child reads |
| `--threads` / `-t` | 4 | Number of threads for jellyfish |
| `--debug-kmers` | false | Enable per-variant debug output |

### Output

The output VCF is annotated with five fields, **DKU**, **DKT**, **DKA**,
**MAX_PKC**, and **AVG_PKC**:

* **DKU** – Number of child reads with at least one variant-spanning k-mer
  unique to the child (absent from both parents).
* **DKT** – Total child reads with variant-spanning k-mers.
* **DKA** – Number of child reads with at least one unique k-mer that also
  exactly support the candidate allele. This helps distinguish real de novo
  signal from spurious noise.
* **MAX_PKC** – Maximum k-mer count among variant-spanning k-mers found in
  the parents. A value of 0 means none of the variant k-mers were observed
  in either parent.
* **AVG_PKC** – Average k-mer count (rounded to 2 decimal places) across
  variant-spanning k-mers found in the parents.

When `--proband-id` is provided and the given ID matches a sample in the
input VCF, DKU, DKT, DKA, MAX_PKC, and AVG_PKC are written as **FORMAT**
(per-sample) fields on that sample. If `--proband-id` is omitted or does
not match any VCF sample, they are written as **INFO** fields instead.

The optional `--metrics` JSON file provides a summary including total
variants, child-unique k-mer counts, and the number of variants with unique
reads.

The optional `--summary` text file provides a human-readable overview
including variant counts, read-support statistics (DKU, DKT, DKA, MAX_PKC,
AVG_PKC), and a per-variant table showing all annotation values and de novo
calls.

The optional `--informative-reads` BAM file contains child reads that carry
at least one variant-spanning k-mer absent from both parents. Each read is
tagged with `DV` indicating which variant(s) it supports. The BAM is sorted
and indexed for direct visualization in IGV.

## Docker

A Docker image is published to GitHub Container Registry on every push to
`main`:

```bash
docker run --rm -v $PWD:/data ghcr.io/jlanej/kmer_denovo_filter:latest \
  --child /data/child.bam \
  --mother /data/mother.bam \
  --father /data/father.bam \
  --vcf /data/candidates.vcf \
  --output /data/annotated.vcf \
  --proband-id HG002
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
apptainer exec kmer_denovo.sif kmer-denovo \
  --child   child.bam \
  --mother  mother.bam \
  --father  father.bam \
  --vcf     candidates.vcf \
  --output  annotated.vcf \
  --proband-id HG002 \
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
