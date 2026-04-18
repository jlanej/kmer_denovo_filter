# HG002 Trio – End-to-End De Novo Variant Filtering Example

This directory contains a complete, batteries-included example for running
`kmer-denovo` on the [GIAB](https://www.nist.gov/programs-projects/genome-bottle)
HG002 Ashkenazi trio using Apptainer (Singularity) on a SLURM-managed HPC
cluster.

## Trio

| Sample | Role   | GIAB ID   | Coriell ID |
|--------|--------|-----------|------------|
| HG002  | Child  | NA24385   | GM24385    |
| HG003  | Father | NA24149   | GM24149    |
| HG004  | Mother | NA24143   | GM24143    |

## Data Sources

| Data           | Description                                               |
|----------------|-----------------------------------------------------------|
| Alignment files | NIST Illumina 2×250 bp WGS, novoalign, GRCh38            |
| Benchmark VCFs | GIAB v4.2.1 (GRCh38, chr1-22), high-confidence calls      |
| Download method | Aspera (`ascp`) with HTTPS/`wget` fallback                |

All files are hosted on the NCBI FTP:
`ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/`

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────┐
│  Step 1: Download GIAB trio BAMs + VCFs via Aspera/HTTPS    │
│          (idempotent – skips files already downloaded)       │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│  Step 2: Pull/cache Apptainer container                     │
│          (docker://ghcr.io/jlanej/kmer_denovo_filter:latest)│
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│  Step 3: Identify putative de novo variants                 │
│          (bcftools isec – child-private sites from trio)     │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│  Step 4: Run kmer-denovo via Apptainer                      │
│          (k-mer evidence annotation of candidates)          │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│  Step 5: Extract mini CRAMs/BAMs for IGV review             │
│          (±1 kb around each candidate site, per trio member) │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│  Step 6: Create IGV variant review TSV                      │
│          (variant table + mini CRAM + annotated VCF paths)  │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│  Step 7: Report results                                     │
│          (annotated VCF, metrics, summary, reads BAM)       │
└─────────────────────────────────────────────────────────────┘
```

## Quick Start

### SLURM submission (recommended)

```bash
sbatch --partition=compute --account=mylab \
    examples/HG002_trio/run_hg002_trio.sh \
    --data-dir    /scratch/$USER/hg002_data \
    --results-dir /scratch/$USER/hg002_results
```

### Interactive run

```bash
bash examples/HG002_trio/run_hg002_trio.sh \
    --data-dir    /scratch/$USER/hg002_data \
    --results-dir /scratch/$USER/hg002_results \
    --threads 8 \
    --memory 32
```

### Using a pre-built container

```bash
# Pull the container first
apptainer pull kmer_denovo.sif docker://ghcr.io/jlanej/kmer_denovo_filter:latest

# Then run with the local .sif
bash examples/HG002_trio/run_hg002_trio.sh \
    --data-dir /scratch/$USER/hg002_data \
    --results-dir /scratch/$USER/hg002_results \
    --container kmer_denovo.sif
```

### Resuming after interruption

The script is idempotent. Simply re-run the same command; downloads that
completed successfully are skipped and the pipeline resumes from the
next incomplete step.

```bash
# Re-run – completed downloads are skipped
sbatch --partition=compute \
    examples/HG002_trio/run_hg002_trio.sh \
    --data-dir /scratch/$USER/hg002_data \
    --results-dir /scratch/$USER/hg002_results
```

### Skip downloads (data already staged)

```bash
bash examples/HG002_trio/run_hg002_trio.sh \
    --skip-download \
    --data-dir /data/giab/hg002 \
    --results-dir /scratch/$USER/hg002_results
```

## Prerequisites

| Tool       | Version | Notes                                         |
|------------|---------|-----------------------------------------------|
| Apptainer  | ≥ 1.1   | Or Singularity ≥ 3.8; `module load apptainer` |
| bcftools   | ≥ 1.10  | For de novo candidate identification           |
| tabix      | –       | Usually bundled with htslib/bcftools           |
| bgzip      | –       | Usually bundled with htslib/bcftools           |
| ascp       | –       | Aspera CLI; optional (falls back to wget)      |
| wget       | –       | Used as fallback if Aspera is unavailable      |

### Installing Aspera CLI

Aspera provides the fastest transfers from NCBI. Install via conda:

```bash
conda install -c conda-forge aspera-cli
```

Or install the [IBM Aspera Connect](https://www.ibm.com/products/aspera/downloads)
client manually. The script auto-discovers the SSH key from common locations.

If Aspera is not available, the script falls back to `wget` over HTTPS.

## Parameters

### Data & Output

| Parameter         | Default            | Description                                  |
|-------------------|--------------------|----------------------------------------------|
| `--data-dir`      | `./hg002_data`     | Download directory for BAMs and VCFs         |
| `--results-dir`   | `./hg002_results`  | Results directory                            |
| `--tmp-dir`       | `RESULTS_DIR/tmp`  | Temp dir for jellyfish hashes (avoid tmpfs!) |

### Compute

| Parameter      | Default | Description                              |
|----------------|---------|------------------------------------------|
| `--threads`    | 16      | Thread count for jellyfish / kmer-denovo |
| `--memory`     | 64      | Available memory in GB                   |
| `--kmer-size`  | 31      | K-mer size                               |

### Container

| Parameter       | Default                                           | Description                    |
|-----------------|---------------------------------------------------|--------------------------------|
| `--container`   | `docker://ghcr.io/jlanej/kmer_denovo_filter:latest` | Apptainer image URI or .sif  |
| `--sif-cache`   | `./containers`                                    | Directory for cached .sif files|

### Download

| Parameter           | Default | Description                                |
|---------------------|---------|--------------------------------------------|
| `--aspera-key`      | auto    | Path to Aspera SSH key                     |
| `--aspera-max-rate` | `500m`  | Aspera max transfer rate                   |
| `--skip-download`   | off     | Skip downloads; use pre-existing files     |
| `--force-download`  | off     | Force re-download even if files exist      |

### Analysis

| Parameter         | Default | Description                                      |
|-------------------|---------|--------------------------------------------------|
| `--ref-fasta`     | –       | Reference FASTA (required for CRAM input)        |
| `--variant-types` | all     | Variant types for de novo scan (e.g. `snps`)     |
| `--proband-id`    | `HG002` | Proband sample ID in the VCF                     |
| `--extra-args`    | –       | Additional arguments passed to `kmer-denovo`     |
| `--mini-cram-padding` | 1000 | Padding in bp for mini CRAM extraction (±bp)  |

Parameters can also be set via environment variables (e.g. `DATA_DIR`,
`THREADS`, `MEMORY_GB`). Command-line arguments take precedence.

## Resource Requirements

| Resource  | Estimate               | Notes                                       |
|-----------|------------------------|---------------------------------------------|
| Disk      | ~700 GB                | BAMs (~480 GB) + VCFs + indices + working   |
| Memory    | 64 GB recommended      | For jellyfish hash tables                   |
| CPUs      | 16 cores recommended   | Parallelizes jellyfish counting             |
| Wall time | 6–24 hours             | Depends on network, cluster load, I/O speed |

## Output Files

| File                             | Description                                    |
|----------------------------------|------------------------------------------------|
| `putative_denovos.vcf.gz`       | Child-private variants (input to kmer-denovo)  |
| `putative_denovos.vcf.gz.tbi`   | Tabix index                                    |
| `HG002_denovo_annotated.vcf`    | Annotated VCF with DKU/DKT/DKA scores         |
| `HG002_metrics.json`            | Per-variant metrics in JSON format             |
| `HG002_summary.txt`             | Human-readable summary of results              |
| `HG002_informative_reads.bam`   | BAM with reads carrying child-unique k-mers    |
| `mini_crams/`                   | Directory with mini alignment files for IGV    |
| `mini_crams/HG002_trio_child.*` | Child reads ±1 kb around each candidate        |
| `mini_crams/HG002_trio_father.*`| Father reads ±1 kb around each candidate       |
| `mini_crams/HG002_trio_mother.*`| Mother reads ±1 kb around each candidate       |
| `mini_crams/HG002_trio_regions.bed` | Extraction regions BED file                |
| `mini_crams/HG002_trio_regions_merged.bed` | Merged extraction regions         |
| `HG002_igv_review.tsv`                  | IGV variant review server input (variant table with mini CRAM + annotated VCF paths) |

## Scripts

### `run_hg002_trio.sh`

Main driver script. Orchestrates the full pipeline from download through
analysis. Can be submitted directly as a SLURM batch job or run
interactively.

### `identify_putative_denovos.sh`

Reusable helper script that identifies child-private variants from trio
VCFs using `bcftools isec`. Performs a two-pass complement:

1. Remove variants shared with the father
2. Remove variants shared with the mother

The remaining variants are child-private (putative de novos). This script
can be used independently on any trio:

```bash
bash examples/HG002_trio/identify_putative_denovos.sh \
    --child-vcf  child.vcf.gz  \
    --father-vcf father.vcf.gz \
    --mother-vcf mother.vcf.gz \
    --output     putative_denovos.vcf.gz
```

### `create_igv_review_tsv.sh`

Generates a variant TSV file ready for the
[IGV de novo variant review server](https://github.com/jlanej/igv.js/tree/master/server).

The TSV contains all required IGV server columns:

- `chrom`, `pos`, `ref`, `alt` — required variant coordinates
- `quality`, `filter`, `child_gt` — from the annotated VCF
- `dku`, `dkt`, `dka`, `dku_dkt`, `dka_dkt` — kmer-denovo evidence scores
  (plus any Kraken2 contamination fractions if the pipeline was run with a
  Kraken2 database)
- `inheritance` — set to `de_novo` for every row
- `child_file` / `father_file` / `mother_file` + index — absolute paths to
  the mini CRAM/BAM files from `extract_mini_crams.sh`
- `child_vcf` / `child_vcf_index` / `child_vcf_id` — bgzipped annotated VCF
  (bgzipped + tabix-indexed automatically if the input is a plain `.vcf`)

Population frequency, variant impact, and gene annotations are not available
from the kmer-denovo pipeline and are therefore omitted.

```bash
bash examples/HG002_trio/create_igv_review_tsv.sh \
    --vcf         HG002_denovo_annotated.vcf.gz \
    --mini-dir    mini_crams/                   \
    --output      HG002_igv_review.tsv          \
    --prefix      HG002_trio                    \
    --proband-id  HG002
```

Start the review server (requires [igv.js](https://github.com/jlanej/igv.js)):

```bash
node /path/to/igv.js/server/server.js \
    --variants  HG002_igv_review.tsv \
    --data-dir  /results             \
    --genome    hg38                 \
    --port      3000
# Open: http://127.0.0.1:3000
```

Reusable helper script that extracts small alignment files (CRAM or BAM)
containing only reads within ±padding of candidate variant sites.  These
"mini" files are ideal for IGV review without requiring hundreds of
gigabytes of full-genome alignment data.

When `--ref-fasta` is provided the output is CRAM (highly compressed);
otherwise BAM is produced. Overlapping regions are automatically merged
to reduce redundant extraction.

```bash
bash examples/HG002_trio/extract_mini_crams.sh \
    --vcf         candidates.vcf.gz \
    --child-bam   child.bam         \
    --father-bam  father.bam        \
    --mother-bam  mother.bam        \
    --output-dir  mini_crams/       \
    --ref-fasta   GRCh38.fa         \
    --padding     1000
```

## Troubleshooting

### Aspera transfer fails

Set `--aspera-max-rate` to a lower value (e.g. `100m`) or let the script
fall back to HTTPS automatically. Ensure your Aspera SSH key is readable.

### Container pull fails

Pre-pull the container on a node with internet access:

```bash
apptainer pull kmer_denovo.sif docker://ghcr.io/jlanej/kmer_denovo_filter:latest
```

Then pass `--container kmer_denovo.sif` to the run script.

### "No putative de novo variants found"

This typically indicates the VCF files were not downloaded correctly.
Verify with:

```bash
bcftools view -H data/vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | head
```

### Out of memory

Increase `--memory` and the SLURM `--mem` allocation. For the full HG002
trio, 64 GB is recommended. You can also set `--extra-args "--jf-hash-size 30G"`
to control the jellyfish hash size for the **discovery mode** child k-mer
counting step. Note that `--jf-hash-size` has no effect in VCF mode; VCF mode
parent scan hash sizes are automatically sized from the number of child k-mers.

### Disk space errors

The pipeline needs ~700 GB total. Ensure your scratch filesystem has
sufficient space. Use `--tmp-dir` to point to a directory on a large
filesystem. **Avoid RAM-backed filesystems** (e.g. `/tmp` on many HPC
systems) for `--tmp-dir`.

## References

- **GIAB**: Zook JM et al. "A robust benchmark for detection of germline
  large deletions and insertions." *Nature Biotechnology*, 2020.
  https://doi.org/10.1038/s41587-020-0538-8

- **GIAB Benchmark VCFs v4.2.1**:
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/

- **kmer_denovo_filter**: https://github.com/jlanej/kmer_denovo_filter
