# Tests

## Unit Tests

- **test_cli.py** – CLI argument parsing and validation.
- **test_kmer_utils.py** – K-mer utility functions (canonicalization, extraction).
- **test_pipeline.py** – Pipeline integration tests using synthetic BAM/VCF data.
- **test_example_output.py** – Regression tests that fail when committed example
  output changes (metrics, summary, VCF annotations). Shows a unified diff on
  failure.
- **test_example_output_discovery.py** – Regression tests for discovery-mode
  output (BED, metrics, summary). Shows a unified diff on failure.
- **test_integration_comparison.py** – Integration tests comparing VCF-mode
  candidates to discovery-mode regions, verifying that high-quality de novo
  candidates are captured within discovered genomic regions.  Also evaluates
  all 7 curated Sulovari et al. 2023 DNM SV regions against discovery output.
- **conftest.py** – Shared fixtures, including session-scoped
  `generated_example_output` and `generated_discovery_output` fixtures that
  run the GIAB pipeline once and return paths to all output files for reuse
  by multiple tests.

```bash
pytest tests/test_cli.py tests/test_kmer_utils.py tests/test_pipeline.py -v
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

The discovery pipeline automatically evaluates its regions against the 7
curated de novo SV loci from Sulovari et al. 2023 (PMC10006329).  All 7
are detected (100% sensitivity) with the following evidence profile:

| Locus | Event | Size | Reads | K-mers | Signal | MaxClip | Class |
|---|---|---|---|---|---|---|---|
| chr17:53340465 | Deletion | 107 bp | 21 | 40 | 0.0330 | 108 | AMBIGUOUS |
| chr14:23280711 | Microsatellite expansion | – | 10 | 16 | 0.0182 | 29 | SV |
| chr3:85552367 | SV-like event | 64 bp | 4 | 7 | 0.0205 | 9 | SMALL |
| chr5:97089276 | SV-like event | 43 bp | 22 | 30 | 0.0690 | 50 | SMALL |
| chr8:125785998 | SV-like event | 43 bp | 37 | 56 | 0.0538 | 78 | SMALL |
| chr18:62805217 | SV-like event | 34 bp | 14 | 42 | 0.0515 | 34 | SMALL |
| chr7:142786222 | Deletion (TRB) | 10,607 bp | 14 | 112 | 0.0100 | 94 | SV |

**Observations:**
- The chr5 and chr8 loci show the strongest k-mer signal (≥0.05 k-mers/bp),
  consistent with insertion events that create novel sequence not present in
  the reference or parents.
- The chr17 107 bp deletion is classified as AMBIGUOUS with a 108 bp max clip
  length matching the expected event size, confirming breakpoint-spanning reads.
- The chr14 microsatellite expansion is correctly classified as SV based on
  5 unmapped mates and 1 discordant pair, consistent with reads failing to
  align across an expanded repeat.
- The chr7 TRB locus 10.6 kb deletion is captured by 3 separate discovery
  regions with SV classification, reflecting the expected breakpoint pattern
  for a large deletion.
- The chr18 34 bp event (DKU=0, inherited in VCF mode) still shows 42
  proband-unique k-mers in discovery mode, illustrating that k-mer-based
  discovery can surface variants missed by VCF-guided annotation.

### Keeping Output Up to Date

The integration test workflow automatically commits updated output to
`tests/example_output/` and `tests/example_output_discovery/` after every
successful run on the `main` branch. This means the example output in this
repository always matches the latest version of the tool. Workflow artifacts
for individual runs are also available in the
[Actions tab](../../actions/workflows/integration-test.yml).
