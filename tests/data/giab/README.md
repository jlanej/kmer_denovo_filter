# GIAB HG002 Trio – Child-Private Variant Test Data

Small BAM slices around child-private SNVs and curated SV-like de novo
mutation candidates from the GIAB Ashkenazi trio, for integration testing
of `kmer_denovo_filter`.

## Variant sources

### 1. Discovered child-private SNVs

Child-private variants are SNVs present in the HG002 (child) GIAB v4.2.1
benchmark VCF but absent from both HG003 (father) and HG004 (mother)
benchmark VCFs.  These serve as realistic candidates for testing de novo
mutation filtering workflows.

### 2. Curated SV-like de novo mutation candidates

BAM regions are always extracted around SV-like de novo mutations reported
in **Sulovari et al. 2023** (PMID: 36894594, PMC10006329), regardless of
whether they appear in the benchmark VCF.  This ensures read-level data is
available for tuning discovery methods.  Any HG002 benchmark VCF variants
overlapping these regions are also checked against the parental VCFs; only
confirmed child-private variants are included in the candidates VCF.

| Locus              | Event type                     | Size (bp) | Padding           |
|--------------------|--------------------------------|-----------|--------------------|
| chr17:53340465     | Deletion                       | 107       | ±500 bp            |
| chr14:23280711     | Microsatellite repeat expansion| –         | ±500 bp            |
| chr3:85552367      | SV-like event                  | 64        | ±500 bp            |
| chr5:97089276      | SV-like event                  | 43        | ±500 bp            |
| chr8:125785998     | SV-like event                  | 43        | ±500 bp            |
| chr18:62805217     | SV-like event                  | 34        | ±500 bp            |
| chr7:142786222     | Deletion (TRB locus, lower-confidence) | 10,607 | −1 kb / +11 kb |

**Reference:**
> Sulovari A, Li R, Audano PA, et al. "Human-specific tandem repeat
> expansion and deletion variants show widespread signatures of
> positive selection." *Science*, 2023; 380(6645):eabn6358.
> https://pmc.ncbi.nlm.nih.gov/articles/PMC10006329/

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
- **SV-like DNMs**: Curated regions from Sulovari et al. 2023 (PMC10006329);
  BAM slices extracted regardless of VCF content, but VCF candidates verified
  for child-private status against parental benchmark VCFs
- **Verification**: Each SNV position verified against HG003 and HG004 GIAB
  v4.2.1 benchmark VCFs over HTTPS

## Files

- `HG002_child.bam` / `.bai` – Child reads around variant sites
- `HG003_father.bam` / `.bai` – Father reads around variant sites
- `HG004_mother.bam` / `.bai` – Mother reads around variant sites
- `candidates.vcf.gz` / `.tbi` – VCF with child-private SNVs and verified
  de novo variants from curated SV-like DNM regions

## Usage with kmer-denovo

```bash
kmer-denovo \
    --child   tests/data/giab/HG002_child.bam   \
    --father  tests/data/giab/HG003_father.bam  \
    --mother  tests/data/giab/HG004_mother.bam  \
    --vcf     tests/data/giab/candidates.vcf.gz  \
    --output  output.vcf
```
