# Tests

## Unit Tests

- **test_cli.py** – CLI argument parsing and validation.
- **test_kmer_utils.py** – K-mer utility functions (canonicalization, extraction).
- **vcf/test_pipeline.py** – VCF-mode pipeline integration tests using synthetic BAM/VCF data.
- **discovery/test_pipeline.py** – Discovery-mode pipeline integration tests using synthetic BAM data.
- **helpers.py** – Shared test helper functions for creating synthetic BAM, VCF, and FASTA data.
- **test_example_output.py** – Regression tests that fail when committed example
  output changes (metrics, summary, VCF annotations). Shows a unified diff on
  failure.
- **test_example_output_discovery.py** – Regression tests for discovery-mode
  output (BED, metrics, summary). Shows a unified diff on failure.
- **test_integration_comparison.py** – Integration tests comparing VCF-mode
  candidates to discovery-mode regions, verifying that high-quality de novo
  candidates are captured within discovered genomic regions.  Also evaluates
  the 7 curated Sulovari et al. 2023 DNM SV regions
  (`examples/HG002_trio/sulovari2023_dnm_regions.tsv`) against discovery output.
- **conftest.py** – Shared fixtures, including session-scoped
  `generated_example_output` and `generated_discovery_output` fixtures that
  run the GIAB pipeline once and return paths to all output files for reuse
  by multiple tests.

```bash
pytest tests/test_cli.py tests/test_kmer_utils.py tests/vcf/test_pipeline.py tests/discovery/test_pipeline.py -v
```

## GIAB Integration Test

The CI workflow ([integration-test.yml](../.github/workflows/integration-test.yml))
runs a full end-to-end pipeline on real data from the
[Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
HG002 trio. The test data includes:

1. **Child-private SNVs** – SNVs present in HG002 but absent from both
   parents' GIAB v4.2.1 benchmark VCFs, discovered across multiple chromosomes.
2. **Curated SV-like de novo mutation candidates** – 7 structural variant-like
   events from Sulovari et al. 2023 (PMC10006329), including deletions,
   microsatellite expansions, and SV-like events ranging from 34 bp to ~10.6 kb.
   BAM regions are always extracted for these loci; any HG002 benchmark VCF
   variants overlapping these regions are verified as child-private against
   parental VCFs before inclusion in the candidates VCF.

Together these yield 22 candidate variants in the VCF with BAM slices for
the child (HG002), father (HG003), and mother (HG004).

### Example Output

Up-to-date example output from the latest successful integration test on
`main` is committed automatically to
[`tests/example_output/`](example_output/) (VCF mode) and
[`tests/example_output_discovery/`](example_output_discovery/) (discovery
mode). These directories are refreshed by CI on every push to `main`, so
they always reflect the current state of the tool.

#### VCF-mode output files

| File | Description |
|---|---|
| `annotated.vcf.gz` | Input VCF annotated with k-mer–based de novo metrics |
| `annotated.vcf.gz.tbi` | Tabix index for the annotated VCF |
| `metrics.json` | Summary counts (total variants, unique k-mers, etc.) |
| `summary.txt` | Human-readable report with per-variant results |

#### Discovery-mode output files

| File | Description |
|---|---|
| `giab_discovery.bed` | Candidate regions with read/k-mer counts and SV annotations |
| `giab_discovery.metrics.json` | Per-region detail with SV classification |
| `giab_discovery.summary.txt` | Human-readable discovery summary with candidate comparison |
| `giab_discovery.sv.bedpe` | Linked breakpoint pairs (BEDPE format) |

### Result Highlights

The GIAB test set contains 22 candidate variants (child-private SNVs plus
verified de novo variants from curated SV-like DNM regions). Key results
from the example output:

- **12 of 22** candidates are classified as likely de novo (`DKU > 0`).
- **10 of 22** candidates show no child-unique k-mers and are marked as
  inherited.
- The average DKU among likely de novo variants is **6.7**, indicating
  strong child-unique read support.
- The SV-like insertion at chr8:125785997 (43 bp) shows strong de novo
  signal with DKA=24 and DKA_DKT=0.30.
- Metrics show **1,484 total child k-mers** extracted, of which **190**
  (≈13%) were absent from both parents.

Discovery mode identifies **27 candidate regions** from the same data,
with **3 high-quality candidates** (DKA_DKT > 0.25, DKA > 10) captured
at 100% rate.

### Curated DNM Region Evaluation

When run with `--dnm-regions examples/HG002_trio/sulovari2023_dnm_regions.tsv`
(as the integration test and the `generated_discovery_output` fixture do), the
discovery pipeline evaluates its regions against the 7 curated de novo SV loci
from Sulovari et al. 2023 (PMC10006329).  5 of the 7 are detected:

| Locus | Event | Size | Reads | K-mers | Signal | MaxClip | Class | Status |
|---|---|---|---|---|---|---|---|---|
| chr17:53340465 | Deletion | 107 bp | 18 | 39 | 0.0414 | 108 | SV | DETECTED |
| chr14:23280711 | Microsatellite expansion | – | – | – | – | – | – | NOT_DETECTED |
| chr3:85552367 | SV-like event | 64 bp | – | – | – | – | – | NOT_DETECTED |
| chr5:97089276 | SV-like event | 43 bp | 22 | 30 | 0.0714 | 50 | SV | DETECTED |
| chr8:125785998 | SV-like event | 43 bp | 34 | 53 | 0.0776 | 54 | SV | DETECTED |
| chr18:62805217 | SV-like event | 34 bp | 7 | 39 | 0.0716 | 27 | SMALL | DETECTED |
| chr7:142786222 | Deletion (TRB) | 10,607 bp | 12 | 112 | 0.0100 | 94 | SV | DETECTED |

**Observations:**
- The chr17 107 bp deletion is `SV` from both kinds of breakpoint evidence:
  three reads carry `107D` in their CIGAR, which also types the region
  `DEL`, and six are soft-clipped at one breakpoint.
- The chr7 TRB locus 10.6 kb deletion is captured by 3 separate discovery
  regions, the expected breakpoint pattern for a large deletion. Two
  discordant pairs join its breakpoints, which is the one link in
  `giab_discovery.sv.bedpe`: a forward read at the left breakpoint and a
  reverse read at the right one give orientation `+ -`, a `DEL` (the test
  BAMs are aligned with novoalign and carry no SA tags, so every link must
  come from mates).
- The chr18 34 bp event (DKU=0, inherited in VCF mode) still shows 39
  proband-unique k-mers in discovery mode, illustrating that k-mer-based
  discovery can surface variants missed by VCF-guided annotation.

### Keeping Output Up to Date

The integration test workflow automatically commits updated output to
`tests/example_output/` and `tests/example_output_discovery/` after every
successful run on the `main` branch. This means the example output in this
repository always matches the latest version of the tool. Workflow artifacts
for individual runs are also available in the
[Actions tab](../../actions/workflows/integration-test.yml).
