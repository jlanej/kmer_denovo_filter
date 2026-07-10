# kraken2 binary and database setup

`nonhuman-screen` shells out to the `kraken2` binary and reads a kraken2
database directory. Neither is bundled.

## The kraken2 binary

Install kraken2 and make sure `kraken2` is on `PATH`. This package is developed
and tested against **kraken2 v2.17.1** (pinned in the package `Dockerfile`).

```bash
git clone --depth 1 --branch v2.17.1 https://github.com/DerrickWood/kraken2.git
cd kraken2 && ./install_kraken2.sh /usr/local/bin
kraken2 --version
```

Or use the package's Docker image, which installs the pinned version.

## The database

A kraken2 database directory must contain the index files `hash.k2d`,
`opts.k2d`, and `taxo.k2d`. It should **also** contain the taxonomy dumps
`nodes.dmp` and `names.dmp` (either under `taxonomy/` or at the database root):
these are strongly recommended but technically optional — without them the
engine still runs but falls back to exact-taxid matching, which corrupts the
non-human signal in both directions (see [methodology.md §7](methodology.md)).
`nodes.dmp` drives lineage-aware classification and `names.dmp` supplies taxon
names.

### PrackenDB (recommended)

PrackenDB is a curated, pre-built kraken2 database of NCBI reference assemblies
across all domains, built with kraken2's default k-mer length (**35**) and
shipping the taxonomy dumps. Download and validate it with the bundled script:

```bash
scripts/download_kraken2_db.sh --db /path/to/kraken_db
# override the source with --url <tarball>
```

The script downloads the tarball
(`k2_NCBI_reference_20251007.tar.gz` by default), extracts it (handling the
versioned subdirectory layout), and validates that `hash.k2d`, `opts.k2d`, and
`taxo.k2d` are present (failing if any is missing), warning separately if
`nodes.dmp` or `names.dmp` is absent.

### Confirming the database k-mer length

```python
from nonhuman_screen import Kraken2Runner
print(Kraken2Runner.read_kmer_length("/path/to/kraken_db"))  # e.g. 35
```

This reads the `k` field from the database's `opts.k2d` (an `IndexOptions`
struct). It is a quick way to confirm the database is intact and which k-mer
size kraken2 will use.

## Memory

Large databases are memory-hungry because kraken2 loads the hash table into
RAM. Pass `memory_mapping=True` (Python) or `--memory-mapping` (CLI) to memory-map
the index files instead, trading speed for a much smaller resident footprint.
