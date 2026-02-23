# Tests

## Unit Tests

- **test_cli.py** – CLI argument parsing and validation.
- **test_kmer_utils.py** – K-mer utility functions (canonicalization, extraction).
- **test_pipeline.py** – Pipeline integration tests using synthetic BAM/VCF data.

```bash
pytest tests/test_cli.py tests/test_kmer_utils.py tests/test_pipeline.py -v
```

## GIAB Integration Test

The CI workflow ([integration-test.yml](../.github/workflows/integration-test.yml))
runs a full end-to-end pipeline on real data from the
[Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
HG002 trio. It uses 20 curated candidate de novo variants across multiple
chromosomes with BAM slices for the child (HG002), father (HG003), and
mother (HG004).

### Example Output

Up-to-date example output from the latest successful integration test on
`main` is committed automatically to
[`tests/example_output/`](example_output/). This directory is refreshed by
CI on every push to `main`, so it always reflects the current state of the
tool.

The output files are:

| File | Description |
|---|---|
| `annotated.vcf.gz` | Input VCF annotated with k-mer–based de novo metrics |
| `annotated.vcf.gz.tbi` | Tabix index for the annotated VCF |
| `metrics.json` | Summary counts (total variants, unique k-mers, etc.) |
| `summary.txt` | Human-readable report with per-variant results |

### Result Highlights

The GIAB test set contains 20 candidate variants. Key results from the
example output:

- **11 of 20** candidates are classified as likely de novo (`DKU > 0`).
- **9 of 20** candidates show no child-unique k-mers and are marked as
  inherited.
- The average DKU among likely de novo variants is **4.5**, indicating
  strong child-unique read support.
- Parent k-mer counts (`MAX_PKC`, `AVG_PKC`) vary widely, helping
  distinguish genuine de novo events from regions with high parental
  background.
- Metrics show **1,291 total child k-mers** extracted, of which **157**
  (≈12%) were absent from both parents.

### Keeping Output Up to Date

The integration test workflow automatically commits updated output to
`tests/example_output/` after every successful run on the `main` branch.
This means the example output in this repository always matches the latest
version of the tool. Workflow artifacts for individual runs are also
available in the
[Actions tab](../../actions/workflows/integration-test.yml).
