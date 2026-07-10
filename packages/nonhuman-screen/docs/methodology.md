# Methodology

`nonhuman-screen` classifies reads by non-human taxonomic content using
kraken2's lowest-common-ancestor (LCA) algorithm and reduces the per-read
verdicts to a conservative **non-human fraction (NHF)** plus per-domain
breakdowns.

## 1. Classification

Reads are written to a temporary FASTQ (with dummy qualities) and classified by
the `kraken2` binary. kraken2 breaks each read into k-mers, looks up the LCA of
each k-mer against the database, and assigns the read a single taxid.

The `--confidence` threshold (0.0–1.0, default 0.0) controls how strict the LCA
must be: a value of 0.2 requires ≥20% of a read's k-mers to agree on the
assigned clade. Higher values reduce spurious assignments at the cost of
sensitivity.

## 2. Lineage-aware domain assignment

The database's `nodes.dmp` taxonomy is parsed into a child→parent map, and
descendant sets are computed for each domain from its NCBI root taxid:

| Domain | Root taxid | Definition |
|---|---|---|
| Bacteria | 2 | all descendants |
| Archaea | 2157 | all descendants |
| Fungi | 4751 | all descendants |
| Protist | — | `eukaryota(2759) − metazoa(33208) − fungi − viridiplantae(33090)` |
| Viruses | 10239 | all descendants |
| UniVec-Core | 81077 | synthetic vector/adapter sequences (see §4) |
| Human clade | 9606 | *Homo sapiens* and any sub-taxa |
| Human lineage | — | ancestors of 9606 up to root (Eukaryota, Root, …) |

A read is **non-human** iff it is classified **and** its taxid is *not* on the
human→root lineage, *not* in the human clade, and *not* under UniVec-Core.

## 3. Human-homology guard (HHG)

Some non-human references share k-mers with the human genome (notably
integrating viruses — endogenous retroviruses, HBV, HPV). To avoid over-calling,
the per-read k-mer detail string is inspected: **if any k-mer voted for human
(taxid 9606), the read is removed from every non-human numerator**, regardless
of its LCA assignment. Such reads are marked with guard status `HHG`.

## 4. UniVec-Core exclusion

UniVec-Core (taxid 81077) is a curated set of synthetic sequencing-vector and
adapter sequences. These are laboratory artefacts, not biological contamination,
and can share k-mers with human DNA. Reads assigned to UniVec-Core are tracked
separately (`univec_core_read_names`) and **unconditionally excluded** from the
non-human fraction. This is a second safety net applied even when the HHG does
not fire.

## 5. The non-human fraction and its partition

Over any set of reads, the classified+unclassified reads are partitioned into
four disjoint sets whose fractions sum to 1.0 (modulo rounding):

```
nonhuman + univec_core + human_lineage + unclassified = 1.0
```

- **nonhuman** — definitively outside the human lineage (after HHG + UniVec).
- **univec_core** — synthetic vectors.
- **human_lineage** — human clade, HHG-cleared reads, and reads assigned to
  ranks too broad to call (Root, Eukaryota, …).
- **unclassified** — kraken2 returned no assignment.

The per-domain breakdowns (`bacterial`, `archaeal`, `fungal`, `protist`,
`viral`) are sub-partitions of `nonhuman`.

## 6. Allele-based NHF

For variant-level screening, the reads supporting a variant's **ALT allele** are
identified by `read_supports_alt`, which extracts the read sequence aligned to
the variant's reference span (via pysam's aligned pairs) and compares it
strictly to the ALT — handling SNPs, MNPs, insertions, deletions, and complex
indels. Only ALT-supporting reads are classified, and the NHF is computed over
that set. When many variants are screened, the union of their ALT-supporting
reads is classified in a single kraken2 invocation and split back per variant.

**Caveat — whole-read scope:** kraken2 classifies the *entire* read, not the
sub-region overlapping the variant. A read is judged non-human as a whole; the
allele-based NHF therefore measures "how many ALT-supporting reads are non-human
reads," which is the intended signal for contamination screening but should not
be read as locus-level taxonomy.

## 7. When taxonomy is unavailable

If `nodes.dmp` cannot be read, classification falls back to **exact-taxid
matching only** (no lineage walk), and the results become unreliable in *both*
directions:

- The **per-domain breakdowns** (`bacterial`, `viral`, …) **under-count**: only
  reads assigned to the exact domain-root taxid match, and `protist` cannot be
  computed at all, so all descendant taxa are missed.
- The **consolidated non-human fraction over-counts**: the fallback treats a
  read as non-human unless its taxid is exactly human (9606), root (1), or
  UniVec-Core, so human-lineage ancestors (e.g. genus *Homo*, Eukaryota) and
  non-9606 human subspecies fall through and are mis-counted as non-human.

In this mode the engine logs a warning and sets
`ClassificationResult.taxonomy_available = False`. **Consumers that gate
decisions on non-human content should treat `taxonomy_available is False` as
"unknown" — neither "clean" nor "contaminated"** — since a database missing its
taxonomy dumps otherwise silently corrupts the signal. Always verify your
database includes `nodes.dmp` (and `names.dmp` for taxon names); see
[database.md](database.md).
