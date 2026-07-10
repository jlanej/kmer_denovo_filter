# Changelog

All notable changes to `nonhuman-screen` are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/), and this project
adheres to [Semantic Versioning](https://semver.org/).

## [0.1.0] — unreleased

Initial release. Extracted from
[`kmer-denovo-filter`](https://github.com/jlanej/kmer_denovo_filter), where the
methodology was originally developed.

### Added
- `Kraken2Runner` engine: kraken2 LCA classification of reads into a non-human
  fraction plus per-domain breakdowns, with a human-homology guard and
  UniVec-Core exclusion (stdlib-only core).
- `ClassificationResult`, `TaxonomicFractions`, `VariantNHF` result types;
  `ClassificationResult.taxonomy_available` flag for fail-closed consumers.
- Allele-based non-human fraction (`[bam]` extra): `classify_variant_alt_reads`,
  `classify_variants_alt_reads` (batched, single kraken2 invocation),
  `classify_reads_from_bam`, `reads_supporting_alt`, `read_supports_alt`.
- `parse_kmer_votes` helper.
- `nonhuman-screen classify` CLI (requires the `[bam]` extra).
- Docs: methodology, database setup, CLI reference; `download_kraken2_db.sh`;
  Dockerfile pinning kraken2 v2.17.1; CI + PyPI-publish workflows.
