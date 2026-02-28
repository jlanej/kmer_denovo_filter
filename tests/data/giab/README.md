# GIAB HG002 Trio – Child-Private Variant Test Data

Small BAM slices around child-private SNVs from the GIAB Ashkenazi trio,
for integration testing of `kmer_denovo_filter`.

Child-private variants are SNVs present in the HG002 (child) GIAB v4.2.1
benchmark VCF but absent from both HG003 (father) and HG004 (mother)
benchmark VCFs.  These serve as realistic candidates for testing de novo
mutation filtering workflows.

## Samples

| File              | Sample | Role   | GIAB ID   |
|-------------------|--------|--------|-----------|
| HG002_child.bam   | HG002  | Child  | NA24385   |
| HG003_father.bam  | HG003  | Father | NA24149   |
| HG004_mother.bam  | HG004  | Mother | NA24143   |

## Data sources

- **BAMs**: NIST Illumina 2×250 bp WGS, aligned with novoalign to GRCh38
  (queried via HTTPS random access – never downloaded in full)
- **Variant discovery**: SNVs streamed from HG002 GIAB v4.2.1 benchmark VCF
  across small chromosomal windows, filtered for absence in both parents
- **Verification**: Each position verified against HG003 and HG004 GIAB
  v4.2.1 benchmark VCFs over HTTPS

## Files

- `HG002_child.bam` / `.bai` – Child reads around variant sites
- `HG003_father.bam` / `.bai` – Father reads around variant sites
- `HG004_mother.bam` / `.bai` – Mother reads around variant sites
- `candidates.vcf.gz` / `.tbi` – VCF with child-private SNVs
- `mini_ref.fa` / `.fai` – Mini reference FASTA built from BAM reads with
  no mismatches (NM:i:0, CIGAR all-M), for discovery-mode testing
- `mini_ref.fa.k31.jf` – Precomputed Jellyfish k-mer index of mini reference

## Usage with kmer-denovo

### VCF mode (candidate VCF required)

```bash
kmer-denovo \
    --child   tests/data/giab/HG002_child.bam   \
    --father  tests/data/giab/HG003_father.bam  \
    --mother  tests/data/giab/HG004_mother.bam  \
    --vcf     tests/data/giab/candidates.vcf.gz  \
    --output  output.vcf
```

### Discovery mode (no VCF needed)

```bash
kmer-denovo \
    --child   tests/data/giab/HG002_child.bam   \
    --father  tests/data/giab/HG003_father.bam  \
    --mother  tests/data/giab/HG004_mother.bam  \
    --ref-fasta tests/data/giab/mini_ref.fa      \
    --ref-jf  tests/data/giab/mini_ref.fa.k31.jf \
    --out-prefix discovery_output \
    --min-child-count 3 --kmer-size 31
```

## Building the mini reference

The mini reference can be regenerated from the BAMs:

```bash
python scripts/build_mini_ref.py \
    -b tests/data/giab/HG002_child.bam \
    -b tests/data/giab/HG003_father.bam \
    -b tests/data/giab/HG004_mother.bam \
    -o tests/data/giab/mini_ref.fa

jellyfish count -m 31 -s 10M -t 4 -C \
    tests/data/giab/mini_ref.fa \
    -o tests/data/giab/mini_ref.fa.k31.jf
```
