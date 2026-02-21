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
  --output  annotated.vcf
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
| `--informative-reads` | – | Output BAM with reads carrying informative k-mers (for IGV) |
| `--metrics` | – | Output summary metrics JSON file |
| `--summary` | – | Output human-readable summary of variant stats and likely DNMs |
| `--kmer-size` / `-k` | 31 | K-mer size |
| `--min-baseq` | 20 | Minimum base quality for read k-mers |
| `--min-mapq` | 20 | Minimum mapping quality for child reads |
| `--threads` / `-t` | 4 | Number of threads for jellyfish |
| `--debug-kmers` | false | Enable per-variant debug output |

### Output

The output VCF contains two additional INFO fields:

* **DKU** – Number of child reads with at least one variant-spanning k-mer
  unique to the child (absent from both parents).
* **DKT** – Total child reads with variant-spanning k-mers.

The optional `--metrics` JSON file provides a summary including total
variants, child-unique k-mer counts, and the number of variants with unique
reads.

The optional `--summary` text file provides a human-readable overview
including variant counts, read-support statistics, and a per-variant table
showing DKU/DKT values and de novo calls.

The optional `--informative-reads` BAM file contains child reads that carry
at least one variant-spanning k-mer absent from both parents. Each read is
tagged with `DV` indicating which variant(s) it supports. The BAM is sorted
and indexed for direct visualization in IGV.

## Testing

```bash
pip install pytest
pytest
```
