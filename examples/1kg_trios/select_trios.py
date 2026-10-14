#!/usr/bin/env python3
"""Select 1000 Genomes trios and write the manifests the trio screen reads.

Inputs are the public NYGC 30x files (00_setup.sh downloads them):

* the pedigree of the 3,202 samples (``20130606_g1k_3202_samples_ped_population.txt``):
  whitespace-separated ``FamilyID SampleID FatherID MotherID Sex Population
  Superpopulation``, one header line;
* one or more ``sequence.index`` files (the 2,504 unrelated and the 698
  related samples): tab-separated, ``#`` header lines, with the ENA CRAM
  path in column 1, its MD5 in column 2 and the sample name in column 10.

A trio is *complete* when the child and both parents have a CRAM.  Of the
603 trios in the pedigree, 602 are complete.  Siblings share parents, so by
default one child per parent pair is kept (``--all-siblings`` keeps them
all; the manager then stages a shared parent once and the finished trio
keeps it until the other is done).

Two manifests are written:

``trios.tsv``
    ``TRIO_ID CHILD FATHER MOTHER FAMILY_ID CHILD_SEX POPULATION SUPERPOPULATION``,
    one row per trio.  The trio ID is the child's sample ID.
``trio_files.tsv``
    ``TRIO_ID ROLE SAMPLE_ID CRAM_FTP_URL CRAI_FTP_URL CRAM_MD5 RELEASE_BATCH``,
    three rows per trio (child, father, mother), which the download manager
    and the per-trio job read.

Selection is deterministic: trios are sorted by child ID, shuffled with
``--seed`` when one is given, and the first ``--n`` are taken.  Pass
``--children`` to name the children instead.
"""

import argparse
import os
import random
import re
import sys

PED_COLUMNS = ("FamilyID", "SampleID", "FatherID", "MotherID", "Sex",
               "Population", "Superpopulation")
TRIOS_HEADER = ("TRIO_ID", "CHILD", "FATHER", "MOTHER", "FAMILY_ID",
                "CHILD_SEX", "POPULATION", "SUPERPOPULATION")
FILES_HEADER = ("TRIO_ID", "ROLE", "SAMPLE_ID", "CRAM_FTP_URL",
                "CRAI_FTP_URL", "CRAM_MD5", "RELEASE_BATCH")
ROLES = ("child", "father", "mother")
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
MD5_RE = re.compile(r"^[0-9a-f]{32}$")
ENA_PREFIX = "ftp://ftp.sra.ebi.ac.uk/"
#: ENA study accession of each release batch (column 4 of the index)
STUDY_BATCHES = {"ERP114329": "2504", "ERP120144": "698"}


def read_pedigree(path):
    """Rows of the pedigree as dicts keyed by PED_COLUMNS."""
    rows = []
    with open(path) as fh:
        header = fh.readline().split()
        if header[:4] != list(PED_COLUMNS[:4]):
            raise ValueError(
                f"{path}: unexpected pedigree header {header[:4]} "
                f"(expected {list(PED_COLUMNS[:4])})")
        for lineno, line in enumerate(fh, 2):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < len(PED_COLUMNS):
                raise ValueError(
                    f"{path}:{lineno}: expected {len(PED_COLUMNS)} columns, "
                    f"found {len(fields)}")
            rows.append(dict(zip(PED_COLUMNS, fields)))
    return rows


def read_sequence_index(path, batch=None):
    """{sample: (cram_url, md5, batch)} for the CRAM rows of an index file.

    The batch labels the row's release (the manager maps a URL to its S3
    mirror path by batch).  Unless *batch* is given it comes from each row's
    study accession (STUDY_BATCHES), falling back to ``698`` when the file
    name says so and ``2504`` otherwise, which matches the NYGC file names.
    """
    file_batch = batch or ("698" if "698" in os.path.basename(path) else "2504")
    samples = {}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith("#"):
                continue
            fields = [f.strip() for f in line.rstrip("\r\n").split("\t")]
            if len(fields) < 10 or not fields[0].endswith(".cram"):
                continue
            cram, md5, sample = fields[0], fields[1].lower(), fields[9]
            if not sample or not cram:
                continue
            row_batch = batch or STUDY_BATCHES.get(fields[3], file_batch)
            # The EBI index has shipped a sample name with a trailing space
            # ("HG03366 "); the strip above removes it, this catches worse.
            if not SAMPLE_ID_RE.match(sample):
                raise ValueError(f"{path}:{lineno}: bad sample name {sample!r}")
            if not cram.startswith(ENA_PREFIX):
                raise ValueError(f"{path}:{lineno}: CRAM is not on the ENA "
                                 f"FTP: {cram}")
            if not MD5_RE.match(md5):
                raise ValueError(f"{path}:{lineno}: bad MD5 {md5!r}")
            if sample in samples and samples[sample][0] != cram:
                raise ValueError(f"{path}:{lineno}: {sample} listed twice "
                                 f"with different CRAMs")
            samples[sample] = (cram, md5, row_batch)
    return samples


def complete_trios(pedigree, samples):
    """(trio rows with CRAMs for all three members, rows without)."""
    complete, incomplete = [], []
    for row in pedigree:
        if row["FatherID"] == "0" or row["MotherID"] == "0":
            continue
        members = (row["SampleID"], row["FatherID"], row["MotherID"])
        if all(m in samples for m in members):
            complete.append(row)
        else:
            incomplete.append(row)
    return complete, incomplete


def _split_list(value):
    if not value:
        return []
    if os.path.isfile(value):
        with open(value) as fh:
            items = [line.split("#")[0].strip() for line in fh]
    else:
        items = value.replace(",", " ").split()
    return [item for item in items if item]


def select(trios, n=0, seed=None, populations=(), superpopulations=(),
           children=(), exclude=(), all_siblings=False):
    """The trios to run, in manifest order."""
    chosen = list(trios)
    if populations:
        chosen = [t for t in chosen if t["Population"] in set(populations)]
    if superpopulations:
        chosen = [t for t in chosen
                  if t["Superpopulation"] in set(superpopulations)]
    if exclude:
        drop = set(exclude)
        chosen = [t for t in chosen if t["SampleID"] not in drop]
    if children:
        by_child = {t["SampleID"]: t for t in chosen}
        missing = [c for c in children if c not in by_child]
        if missing:
            raise ValueError(
                "not a complete trio child (or filtered out): "
                + ", ".join(missing))
        return [by_child[c] for c in children]
    chosen.sort(key=lambda t: t["SampleID"])
    if not all_siblings:
        seen, unique = set(), []
        for t in chosen:
            parents = (t["FatherID"], t["MotherID"])
            if parents in seen:
                continue
            seen.add(parents)
            unique.append(t)
        chosen = unique
    if seed is not None:
        random.Random(seed).shuffle(chosen)
    if n and n > 0:
        chosen = chosen[:n]
    return chosen


def write_manifests(trios, samples, trios_path, files_path):
    with open(trios_path, "w") as out:
        out.write("\t".join(TRIOS_HEADER) + "\n")
        for t in trios:
            out.write("\t".join((
                t["SampleID"], t["SampleID"], t["FatherID"], t["MotherID"],
                t["FamilyID"], t["Sex"], t["Population"],
                t["Superpopulation"])) + "\n")
    with open(files_path, "w") as out:
        out.write("\t".join(FILES_HEADER) + "\n")
        for t in trios:
            members = (t["SampleID"], t["FatherID"], t["MotherID"])
            for role, sample in zip(ROLES, members):
                cram, md5, batch = samples[sample]
                out.write("\t".join((
                    t["SampleID"], role, sample, cram, cram + ".crai", md5,
                    batch)) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ped", required=True, help="3,202-sample pedigree")
    parser.add_argument("--index", required=True, action="append",
                        help="sequence.index file (repeat for both batches)")
    parser.add_argument("--out-trios", required=True)
    parser.add_argument("--out-files", required=True)
    parser.add_argument("--n", type=int, default=0,
                        help="trios to select (0 = all)")
    parser.add_argument("--seed", type=lambda v: int(v) if v != "" else None,
                        default=None, metavar="INT",
                        help="shuffle with this seed before taking --n "
                             "(empty: keep the sorted order)")
    parser.add_argument("--populations", default="",
                        help="comma-separated population codes, or a file")
    parser.add_argument("--superpopulations", default="",
                        help="comma-separated superpopulation codes, or a file")
    parser.add_argument("--children", default="",
                        help="child sample IDs to select (comma-separated, "
                             "or a file, one per line); overrides --n/--seed")
    parser.add_argument("--exclude", default="",
                        help="child sample IDs to leave out")
    parser.add_argument("--all-siblings", action="store_true",
                        help="keep every child of a family, not one per "
                             "parent pair")
    parser.add_argument("--expected-complete", type=int, default=0,
                        help="warn when the number of complete trios differs")
    args = parser.parse_args(argv)

    pedigree = read_pedigree(args.ped)
    samples = {}
    for path in args.index:
        for sample, entry in read_sequence_index(path).items():
            if sample in samples and samples[sample][0] != entry[0]:
                raise SystemExit(f"{sample} has different CRAMs in two indexes")
            samples[sample] = entry
    complete, incomplete = complete_trios(pedigree, samples)
    print(f"pedigree: {len(pedigree)} samples, {len(complete) + len(incomplete)}"
          f" trios, {len(complete)} with CRAMs for all three members",
          file=sys.stderr)
    for row in incomplete:
        print(f"  incomplete trio skipped: child {row['SampleID']} "
              f"(family {row['FamilyID']})", file=sys.stderr)
    if args.expected_complete and len(complete) != args.expected_complete:
        print(f"WARNING: expected {args.expected_complete} complete trios, "
              f"found {len(complete)}; check the pedigree and index files",
              file=sys.stderr)

    try:
        chosen = select(
            complete, n=args.n, seed=args.seed,
            populations=_split_list(args.populations),
            superpopulations=_split_list(args.superpopulations),
            children=_split_list(args.children),
            exclude=_split_list(args.exclude),
            all_siblings=args.all_siblings)
    except ValueError as exc:
        raise SystemExit(f"error: {exc}")
    if not chosen:
        raise SystemExit("error: no trio matches the selection")
    write_manifests(chosen, samples, args.out_trios, args.out_files)
    print(f"selected {len(chosen)} trios -> {args.out_trios}, "
          f"{args.out_files}", file=sys.stderr)


if __name__ == "__main__":
    main()
