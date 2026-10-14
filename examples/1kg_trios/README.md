# 1000 Genomes trios: an at-scale de novo SV screen on SLURM

This directory runs `kmer-discovery` on real 30x whole-genome trios from the
1000 Genomes Project, as many of them as you like, on a SLURM cluster, and
summarises the structural-variant (SV) calls across the cohort. It is the
scale test for discovery mode: the GIAB example in
[`examples/HG002_trio`](../HG002_trio/) is one trio on one node, this is a
cohort.

Everything it reads is public, and nothing in it is specific to one cluster.
The scheduler specifics (partition, account, QoS) go in one variable, and
every path hangs off `WORK_DIR`.

> **Cite the data:** Byrska-Bishop, M. et al. (2022). High-coverage
> whole-genome sequencing of the expanded 1000 Genomes Project cohort
> including 602 trios. *Cell* 185(18), 3426–3440.
> https://doi.org/10.1016/j.cell.2022.08.004. The data-reuse terms are in
> the collection's
> [README](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20200526_1000G_2504plus698_high_cov_data_reuse_README.txt).

* [What it does](#what-it-does)
* [Data sources](#data-sources)
* [Prerequisites](#prerequisites)
* [Quick start](#quick-start)
* [Step by step](#step-by-step)
* [Outputs](#outputs)
* [Reading the results](#reading-the-results)
* [Configuration reference](#configuration-reference)
* [Monitoring and troubleshooting](#monitoring-and-troubleshooting)
* [Testing](#testing)

## What it does

```
00_setup.sh            container, reference, NYGC indexes + pedigree, trio manifest,
                       SV call set subset; submits 00b (reference Jellyfish index)
        │
01_stage_and_discover.sh   one manager job: downloads each trio's three CRAMs
        │                  (S3 mirror / ENA, MD5-checked), then submits ...
        ├──► 01b_discover_trio.sh   one job per trio: child-only SVs from the
        ├──► 01b_discover_trio.sh   call set -> --dnm-regions; kmer-discovery;
        └──► 01b_discover_trio.sh   timing; cleanup of the CRAMs
        │
02_summarize.sh        cohort tables: per-trio counts, recurrent loci across
                       children, the de novo SV candidate list
```

The 3,202 NYGC samples include 602 complete trios (603 in the pedigree; one
GWD child has no CRAM). `select_trios.py` picks the ones to run: a seeded
random subset, a population, or an explicit list. For each trio the
discovery run reports regions of the child's genome carrying k-mers absent
from the reference and from both parents, classes each as `SV`, `AMBIGUOUS`
or `SMALL`, links breakpoints and types SVs (see
[Structural Variant Calling](../../docs/sv_calling.md)). Across trios, the
summary drops loci that recur in unrelated children, which are reference
and alignment artefacts rather than de novo events, and leaves a candidate
list to review.

Two things make this a *screen* rather than just a batch run:

* **Known SVs.** The published 3,202-sample SV call set genotypes every
  sample at every SV. For each trio, the SVs genotyped non-reference in the
  child and homozygous reference in both parents are passed to discovery as
  `--dnm-regions`, so the metrics say, per locus, whether discovery found
  reads there. Most child-only calls in a population call set are parents
  genotyped wrongly rather than de novo SVs, which is why the filters are
  strict by default and both counts are reported.
* **Recurrence.** A region at the same locus in two or more children of
  different families cannot be de novo in all of them. With a few hundred
  trios the recurrent set is a map of discovery's systematic false
  positives on this reference, and what remains per child is the real
  candidate list.

## Data sources

| Data | Where | Size |
|---|---|---|
| CRAMs (30x NovaSeq, GRCh38) | ENA FTP (`ftp.sra.ebi.ac.uk`), listed in the two NYGC `sequence.index` files; mirrored byte-identically on the AWS Open Data bucket `1000genomes` | ~15 GB per sample |
| Sequence indexes | `1000G_2504_high_coverage.sequence.index`, `1000G_698_related_high_coverage.sequence.index` (URL, MD5, sample, population, library) | 1.2 MB |
| Pedigree | `20130606_g1k_3202_samples_ped_population.txt` | 95 KB |
| Reference | `GRCh38_full_analysis_set_plus_decoy_hla.fa` (+ `.fai`), what the CRAMs were aligned to | 3.1 GB |
| SV call set | `1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz` (+ `.tbi`): GATK-SV, svtools and Absinthe calls on all 3,202 samples, with allele frequencies | 1.7 GB |

All under
`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/`
except the CRAMs and the reference (`.../technical/reference/GRCh38_reference_genome/`).

Downloads follow what worked for the
[NGS-PCA 1000G example](https://github.com/jlanej/NGS-PCA/tree/main/example/1000G_highcov):
the AWS mirror first, then ENA over HTTPS, with `aria2c` (16 range streams),
then parallel `curl` ranges, then `wget`, and every CRAM checked against the
index MD5 before a job sees it. Aspera is deliberately absent: the current
ENA server needs a newer `ascp` and a new public key, and on one cluster
every Aspera transfer arrived corrupt while HTTPS ran clean. NGS-PCA keeps a
working Aspera path if you need one. For bulk staging,
`globus_batch_from_manifest.sh` writes Globus batch files for the cohort
and the manager picks the staged files up.

## Prerequisites

| Requirement | Notes |
|---|---|
| Apptainer ≥ 1.1 (or Singularity) | `module load apptainer`; pulls `ghcr.io/jlanej/kmer_denovo_filter:latest`, which has `kmer-discovery`, samtools, bcftools, jellyfish, aria2c and these helper scripts |
| SLURM | `sbatch`/`squeue`; site flags via `SBATCH_EXTRA_ARGS` |
| Internet | On the node that runs `00_setup.sh` and on the node that runs the manager (a DTN or login node with `MANAGER_LOCAL=1`, or a compute partition with outbound access) |
| python3 ≥ 3.9 | On the setup/summary node (standard library only; the container's python is used otherwise) |
| Disk | Reference + index ~25 GB; SV call set 1.7 GB; **~45 GB per staged trio × `MAX_STAGED_TRIOS`**; **100+ GB of scratch per running discovery job** (`TMP_ROOT`, not a RAM-backed `/tmp`); results ~1 GB per trio (informative BAM and report included) |
| Memory | 64 GB for the reference index job (23 GB used); `DISCOVERY_MEM` (default 128 GB) per trio job |

## Quick start

```bash
cd examples/1kg_trios
export WORK_DIR=/scratch/$USER/1kg_trios
export SBATCH_EXTRA_ARGS="--partition=compute --account=mylab"   # your site

# 1. One-time setup: ~5 GB of downloads, then submits the reference index job
N_TRIOS=10 bash 00_setup.sh

# 2. Once $WORK_DIR/reference/*.k31.jf exists: stage and run the trios
bash 01_stage_and_discover.sh            # TRIO_LIMIT=1 for a smoke test

# 3. Any time, as trios finish
bash 02_summarize.sh
cat $WORK_DIR/results/cohort/cohort_summary.md
```

`N_TRIOS=0` selects all 602. Re-running any script continues where it
stopped.

## Step by step

### 0. Setup

`00_setup.sh` is idempotent and runs on a login node:

1. creates the tree under `WORK_DIR`;
2. pulls the container to `SIF_IMAGE` (or symlinks a local `.sif` given as
   `CONTAINER_URI`);
3. downloads the reference and its `.fai`;
4. downloads the two sequence indexes and the pedigree, and runs
   `select_trios.py`, which writes `trios.tsv` (one row per trio) and
   `trio_files.tsv` (three rows per trio: role, sample, CRAM and CRAI
   URLs, MD5, release batch). Delete both to select a different cohort;
   the results directory is keyed by trio, so a new selection only adds
   or drops trios from the sweep;
5. downloads the SV call set and subsets it to the cohort's samples with
   bcftools (`cohort_trios.sv.vcf.gz`), redone whenever the manifest is
   newer. `USE_SV_CALLSET=0` skips this;
6. submits `00b_build_ref_index.sh`, which builds the canonical k=31
   Jellyfish index of the reference once (`SETUP_LOCAL=1` builds it in
   place instead). Every trio job reuses it through `--ref-jf`.

Selection knobs (all environment variables, see `config.sh`):

| Variable | Effect |
|---|---|
| `N_TRIOS` | how many (0 = all complete trios) |
| `TRIO_SEED` | seeded shuffle before taking `N_TRIOS` (export an empty value for sorted order) |
| `TRIO_POPULATIONS`, `TRIO_SUPERPOPULATIONS` | e.g. `GBR,CHS` or `EUR,EAS` |
| `TRIO_CHILDREN` | explicit child IDs, comma-separated or a file, one per line |
| `TRIO_ALL_SIBLINGS` | 1 keeps every child of a family (default: one per parent pair) |

### 1. Staging and discovery

`01_stage_and_discover.sh` submits itself as one long-running manager job
(or runs in place with `MANAGER_LOCAL=1`, e.g. on a data-transfer node
under `tmux`). The manager:

* skips trios that are done or whose discovery job is already queued or
  running (one `squeue` snapshot at start);
* stages `TRIO_SLOTS` trios at a time, each downloading its child, father
  and mother CRAMs in turn with `DOWNLOAD_CONNECTIONS` streams per file;
* verifies every CRAM against the index MD5, deleting and re-downloading
  once on a mismatch;
* reuses CRAMs already on disk (a dead run, a Globus transfer, a sibling's
  parents), still verifying them;
* pauses downloads while `MAX_STAGED_TRIOS` trios sit on disk awaiting
  discovery, so a stalled queue cannot fill scratch;
* submits `01b_discover_trio.sh` for each verified trio with
  `DISCOVERY_CPUS`, `DISCOVERY_MEM` and `DISCOVERY_TIME`;
* exits non-zero listing the trios that failed. Re-run it to retry them.

A parent shared by two selected trios (siblings) is downloaded once and
kept until both are done.

### 1b. One trio

`01b_discover_trio.sh TRIO` runs on a compute node:

1. `child_only_svs.py` lists the trio's child-only SVs from the cohort
   call set, keeping calls that are `PASS`, have a population allele
   frequency of at most `SV_MAX_AF` (0.01) and genotype quality of at
   least `SV_MIN_GQ` (20) in all three members. Multiallelic `<CNV>`
   records are skipped. The result is `<trio>.child_only_svs.tsv` in
   `--dnm-regions` format and a detail table of every child-only call,
   filtered or not. With at least one event the run gets `--dnm-regions`.
2. `kmer-discovery` runs in the container with `--ref-fasta`, `--ref-jf`,
   `--threads` from the allocation and `--memory` from
   `SLURM_MEM_PER_NODE` (or `DISCOVERY_MEM`), `--tmp-dir` under
   `TMP_ROOT/<trio>` (or `DISCOVERY_TMP_DIR`), and `--report`. The log is
   kept as `<trio>.log`.
3. On success it writes `<trio>.timing.tsv` (wall seconds, node, threads,
   memory, job id, image) and `<trio>.done`, then deletes the trio's CRAMs
   unless `KEEP_CRAMS=1` or another selected trio still needs one. On
   failure the CRAMs stay for a retry.

Pass `verify` as a second argument to check the MD5s first, for a trio
staged by hand.

Sizing: Module 1 counts the child's k-mers with a Bloom filter (about 1.2
bytes per distinct k-mer, 10–15 GB for a 30x child) and a hash sized to
the memory budget, counted with `--disk` so it cannot outgrow it. 128 GB is
comfortable; the README's
[SLURM section](../../README.md#example-slurm-batch-script) has the
details. Wall time is dominated by streaming the three CRAMs through
samtools and jellyfish, so it depends on the node's I/O; the per-module
times in `trio_summary.tsv` tell you where it goes on your system, and
`DISCOVERY_TIME` should be set from the first few runs.

### 2. Summary

`02_summarize.sh` runs `summarize_trios.py` over every finished trio
(standard library only, seconds). Regions from all trios are clustered
within `RECURRENT_SLACK_BP` (500 bp); a cluster holding regions from at
least `RECURRENT_MIN_TRIOS` (2) children is recurrent.

## Outputs

Per trio, `$WORK_DIR/results/<trio>/<trio>.*`:

| File | Contents |
|---|---|
| `.bed`, `.sv.bedpe`, `.metrics.json`, `.summary.txt`, `.kmer_coverage.bedgraph`, `.read_coverage.bed`, `.informative.bam` | the `kmer-discovery` outputs ([Discovery Mode Output](../../README.md#discovery-mode-output)) |
| `.report.html` | the interactive report |
| `.child_only_svs.tsv` | the `--dnm-regions` file: child-only SVs passing the filters |
| `.child_only_svs.detail.tsv` | every child-only call with ID, type, size, AF, filters, the three genotypes and GQs, and `PASS_FILTERS`/`FILTERED` |
| `.log` | the full discovery log, with per-module times |
| `.timing.tsv`, `.done` | wall time and resources; the completion marker the manager keys on |

Cohort, `$WORK_DIR/results/cohort/`:

| File | Contents |
|---|---|
| `cohort_summary.md` | headline numbers and a per-trio table |
| `trio_summary.tsv` | one row per trio: regions by class, `SV`-class regions by type, non-recurrent SVs, junctions, k-mer and read totals, child-only SVs (all, filtered, tested, supported), wall and per-module hours |
| `recurrent_regions.bed` | the recurrent loci: interval, number of children, their IDs, classes and types |
| `sv_candidates.tsv` | every non-recurrent `SV`-class region with its evidence columns and junctions (`partner|strands|type|reads`) |
| `sv_candidates.bed`, `sv_candidates.bedpe` | the same for IGV and bedtools, with the trio in the name / last column |

## Reading the results

* **Candidate de novo SVs** are `sv_candidates.tsv`. Expect few per child:
  de novo SVs are rare, and discovery sees only junctions that create new
  k-mers (see [Limitations](../../docs/sv_calling.md#limitations)). A
  candidate with a `DEL`/`DUP`/`INV`/`BND` junction supported by several
  molecules in `junctions` is the strongest kind; an `INS` with no
  junction is reads clipped from both sides at one point.
* **Review** a candidate in IGV: load the trio's `<trio>.informative.bam`
  (every read carries proband-unique k-mers) with `sv_candidates.bed` and
  `sv_candidates.bedpe`. The CRAMs are gone by then unless `KEEP_CRAMS=1`;
  the S3 URLs in `trio_files.tsv` work directly in IGV and samtools for
  a region.
* **Recurrent loci** are the pipeline's systematic false positives on
  GRCh38 with this call set of parents: compare `recurrent_regions.bed`
  with known problematic regions, and consider a cohort-derived blacklist
  for other trios.
* **Known SVs**: `known_svs_tested` and `known_svs_detected` in
  `trio_summary.tsv` count the published child-only SVs passing the
  filters and how many discovery supported. A supported one is either a
  true de novo SV or a parental genotyping miss where discovery sees the
  child's junction reads; check the parents at the locus. An unsupported
  one is most often a parent genotyping miss too, and the detail table's
  `EVIDENCE` and `ALGORITHMS` columns say what called it.
* **Scale**: per-trio rows let you see depth, population and node effects
  on region counts and runtime. Siblings (`TRIO_ALL_SIBLINGS=1`) share no
  de novo SVs, which is a cheap specificity check.

## Configuration reference

All in [`config.sh`](config.sh); export to override.

| Variable | Default | Description |
|---|---|---|
| `WORK_DIR` | `/scratch/$USER/1kg_trios` | root of the working tree |
| `CRAM_DIR`, `RESULTS_DIR`, `LOG_DIR`, `REF_DIR`, `STATE_DIR`, `TMP_ROOT` | under `WORK_DIR` | staged CRAMs, outputs, logs, reference, manager state, discovery scratch |
| `CONTAINER_URI` | `docker://ghcr.io/jlanej/kmer_denovo_filter:latest` | image to pull, or a local `.sif` |
| `SBATCH_EXTRA_ARGS` | empty | site flags for every `sbatch` (partition, account, QoS) |
| `N_TRIOS`, `TRIO_SEED`, `TRIO_POPULATIONS`, `TRIO_SUPERPOPULATIONS`, `TRIO_CHILDREN`, `TRIO_ALL_SIBLINGS` | 10, 1, –, –, –, 0 | trio selection |
| `USE_SV_CALLSET`, `SV_MAX_AF`, `SV_MIN_GQ`, `SV_PASS_ONLY` | 1, 0.01, 20, 1 | the child-only SV screen |
| `S3_HTTPS_BASE`, `ENA_HTTPS_BASE`, `DOWNLOAD_CONNECTIONS` | AWS mirror, ENA, 16 | download sources and streams per file |
| `TRIO_SLOTS`, `MAX_STAGED_TRIOS`, `TRIO_LIMIT` | 2, 4, 0 | trios staged at once; disk backpressure; sweep limit |
| `KEEP_CRAMS` | 0 | 1 keeps CRAMs after discovery |
| `DISCOVERY_CPUS`, `DISCOVERY_MEM`, `DISCOVERY_TIME` | 16, 128G, 36:00:00 | per-trio job resources |
| `MIN_CHILD_COUNT`, `KMER_SIZE`, `DISCOVERY_EXTRA_ARGS`, `DISCOVERY_TMP_DIR` | 3, 31, –, – | discovery parameters |
| `REF_INDEX_CPUS`, `REF_INDEX_MEM`, `REF_INDEX_TIME`, `REF_JF_HASH_SIZE` | 8, 64G, 06:00:00, 4G | reference index job |
| `RECURRENT_MIN_TRIOS`, `RECURRENT_SLACK_BP` | 2, 500 | recurrence filter |
| `JOB_PREFIX` | `1kgtrio` | job names: `1kgtrio_manager`, `1kgtrio_refidx`, `1kgtrio_disc_<trio>` |

## Monitoring and troubleshooting

```bash
squeue -u $USER -n 1kgtrio_manager,1kgtrio_refidx        # the manager and the index build
squeue -u $USER | grep 1kgtrio_disc                        # discovery jobs
tail -f $WORK_DIR/logs/manager_<JOBID>.out                 # one line per dispatch
cat $WORK_DIR/logs/stage_<TRIO>.log                        # a trio's download story
cat $WORK_DIR/logs/discover_<TRIO>_<JOBID>.out             # a trio's job
ls $WORK_DIR/results/*/*.done | wc -l                      # finished trios
grep -h "download complete" $WORK_DIR/logs/stage_*.log | sort | uniq -c   # transport census
```

* **"reference Jellyfish index missing"**: the index job has not finished;
  `squeue -n 1kgtrio_refidx`, or build it in place with
  `bash 00b_build_ref_index.sh` on a node with 32 GB.
* **MD5 mismatches on every file**: a transport corrupting payloads. The
  census above says which; set `S3_HTTPS_BASE=""` to use ENA only, or stage
  with Globus.
* **No aria2c completions**: the host lacks `aria2c` and the container
  could not be run on the manager node; parallel `curl` carries the load,
  more slowly. `module load aria2` on the manager node, or run the manager
  where Apptainer works.
* **Out of memory in a trio job**: raise `DISCOVERY_MEM`; the child count
  sizes itself to `--memory`, so the log's "hash ... entries" line says
  what it asked for.
* **Scratch full**: lower `MAX_STAGED_TRIOS`, or point `TMP_ROOT` at a
  larger filesystem. Never a RAM-backed `/tmp`.
* **Changing the cohort**: delete `trios.tsv` and `trio_files.tsv`, re-run
  setup; finished trios not in the new selection keep their results.
* **Two managers at once** fight over the same files; run one, and let it
  resweep.

## Testing

[`tests/test_1kg_trios_helpers.py`](../../tests/test_1kg_trios_helpers.py)
tests the three Python helpers on synthetic inputs: trio selection, the
child-only SV extraction (whose `--dnm-regions` file is loaded by the
pipeline's own parser), and the cohort summary.
[`tests/test_1kg_trios_pipeline.py`](../../tests/test_1kg_trios_pipeline.py)
runs the SLURM scripts end to end with `sbatch`, `squeue`, `curl`, `wget`
and `apptainer` shimmed and no network: three trios, two of them siblings
and one with a corrupt CRAM, through setup, staging, the shared-parent
bookkeeping, a failed trio healed on a resweep, cleanup and the summary.
Both run in CI and need only pysam and bcftools.

```bash
pytest tests/test_1kg_trios_helpers.py tests/test_1kg_trios_pipeline.py -v
```
