# nonhuman-screen

Classify sequencing reads by **non-human taxonomic content** using
[kraken2](https://github.com/DerrickWood/kraken2), and reduce the per-read
verdicts to a **non-human fraction (NHF)** plus per-domain breakdowns
(bacterial / archaeal / fungal / protist / viral / UniVec-Core).

The headline capability is the **allele-based NHF**: given a BAM/CRAM and a
variant, compute the non-human fraction of the reads that support that
variant's ALT allele — in one call, no boilerplate:

```python
from nonhuman_screen import classify_variant_alt_reads

v = classify_variant_alt_reads(
    "sample.bam", "kraken2_db", "chr1", 12345, "A", "T", ref_fasta="ref.fa",
)
v.nonhuman_fraction   # 0.0–1.0 over ALT-supporting reads
v.fractions.bacterial # per-domain breakdown
v.supporting_reads    # denominator
```

`nonhuman-screen` is **sample-agnostic**: it classifies whatever reads you hand
it and never asks whose they are. The same call works on a proband, a parent, a
tumour, or a plain FASTQ — so it can back higher-level analyses (e.g. flagging
inherited variant calls that are actually parental contamination) without
knowing anything about pedigrees.

> Extracted from [`kmer-denovo-filter`](https://github.com/jlanej/kmer_denovo_filter),
> where this methodology was originally developed, so it can be reused
> independently.

## How it works (in one paragraph)

Reads are written to a temporary FASTQ and classified by kraken2's LCA
algorithm. Each read's assigned taxid is mapped to a domain using the
database's `nodes.dmp` taxonomy. A read counts as **non-human** only if it is
classified and lies outside the human lineage, outside the human clade, and
outside UniVec-Core. Two conservative guards reduce false positives: a
**human-homology guard** (any k-mer voting for human, taxid 9606, drops the read
from every non-human numerator — important for integrating viruses such as
HBV/HPV/ERVs) and **UniVec-Core exclusion** (synthetic vector/adapter sequences,
taxid 81077, are tracked separately and never counted as contamination). See
[docs/methodology.md](docs/methodology.md).

## Install

```bash
pip install nonhuman-screen           # core engine (standard library only)
pip install 'nonhuman-screen[bam]'    # + BAM/CRAM & allele-based NHF (pysam)
pip install 'nonhuman-screen[plots]'  # + optional plotting helpers (plotly)
```

You also need the `kraken2` binary on `PATH` and a kraken2 database — see
[docs/database.md](docs/database.md).

## Usage

### Allele-based NHF for many variants (batched — one kraken2 call)

```python
from nonhuman_screen import classify_variants_alt_reads

variants = [("chr1", 12345, "A", "T"), ("chr2", 999, "G", "GAC")]  # pos 0-based
for v in classify_variants_alt_reads("sample.bam", "kraken2_db", variants,
                                      ref_fasta="ref.fa"):
    print(v.variant_key, v.nonhuman_fraction, v.supporting_reads)
```

The reads supporting all variants' ALT alleles are gathered and classified in a
**single** kraken2 invocation (its database load dominates runtime), then split
back per variant. You never re-implement fetch/allele-filter/classify yourself.

### Classify an arbitrary set of reads

```python
from nonhuman_screen import Kraken2Runner

result = Kraken2Runner("kraken2_db").classify_sequences(
    {"read1": "ACGT...", "read2": "TTGCA..."}
)
result.nonhuman_fraction        # 0.0–1.0
result.fractions()              # TaxonomicFractions(nonhuman=…, bacterial=…, …)
result.taxonomy_available       # False if nodes.dmp was missing (fractions unreliable)
result.nonhuman_read_names      # set of read names
```

### Command line

```bash
# Per-variant allele NHF table (headline)
nonhuman-screen classify \
    --bam sample.bam --kraken2-db kraken2_db --ref-fasta ref.fa \
    --variants candidates.vcf.gz --out-prefix sample_contam
# -> sample_contam.variant_nhf.tsv, sample_contam.summary.json

# Whole-BAM summary
nonhuman-screen classify --bam sample.bam --kraken2-db kraken2_db
```

See [docs/cli.md](docs/cli.md).

## Coordinate convention

All `pos` arguments in the Python API are **0-based** (pysam reference
coordinates); variant keys are `"{chrom}:{pos}:{ref}:{alt}"`. The CLI reads a
VCF and handles the conversion for you.

## API summary

| Symbol | Needs `[bam]` | Purpose |
|---|:---:|---|
| `Kraken2Runner.classify_sequences(seqs)` | no | Classify `{name: seq}` → `ClassificationResult` |
| `ClassificationResult` | no | Per-domain read-name sets, counts, `fractions()`, `nonhuman_fraction`, `taxonomy_available` |
| `TaxonomicFractions` | no | Per-domain fractions; `from_result` / `over_reads` |
| `read_supports_alt(read, pos, ref, alt)` | no | Does an aligned read carry the ALT allele? |
| `parse_kmer_votes(kmer_string)` | no | Parse kraken2 per-read k-mer detail into taxid votes |
| `classify_variant_alt_reads(...)` | yes | Allele-based NHF for one variant → `VariantNHF` |
| `classify_variants_alt_reads(...)` | yes | Batched allele-based NHF for many variants |
| `classify_reads_from_bam(...)` | yes | Classify reads by name and/or locus |
| `reads_supporting_alt(...)` | yes | Read names supporting an ALT allele |

## License

MIT — see [LICENSE](LICENSE).
