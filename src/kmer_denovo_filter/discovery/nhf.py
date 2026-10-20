"""The non-human fraction (NHF) screen of a discovery run.

A read of a bacterium, a virus, a fungus or a sequencing adapter carries
k-mers that are in neither parent nor the reference, so discovery keeps it
as informative by construction; aligned with clips to a human locus, or
placed at its mate's position, it looks like structural-variant evidence.
This module classifies a run's informative reads with the engine VCF mode
already uses for ``DKU_NHF`` / ``DKA_NHF`` (the ``nonhuman-screen`` package
around kraken2) and writes what it finds beside each region.  Nothing is
filtered on it.

**The reads** are the read ends of ``{prefix}.informative.bam``: every read
carrying proband-unique k-mers, mapped or not.  The two mates of a pair are
two sequences and are classified as two reads.

**The classes** are kraken2's taxon per read, reduced by the engine to a
four-way partition: non-human (classified, off the human lineage, outside
UniVec Core and with no k-mer voting human: the human-homology guard),
UniVec Core (vectors and adapters), human lineage, unclassified; the four
sum to 1.  NHF is the non-human share, and the non-human part is broken
down by domain (bacterial, archaeal, fungal, protist, viral).

**Two levels** are measured, with one kraken2 run:

* the run: a deterministic sample of up to ``--nhf-sample-reads``
  informative read ends - those with the smallest hash of their name - the
  child's informative reads as a whole;
* each ``SV`` and ``AMBIGUOUS`` region, and each region at a BEDPE junction
  end: the informative alignments overlapping it, the other read end of each
  of those molecules, and the unmapped reads placed within
  ``--cluster-distance`` of it (an insertion's interior read is unmapped
  and placed at its mate's position, up to an insert size away: the
  ``unmapped_mates`` evidence).  ``SMALL`` regions are not screened.

**The status** is ``screened``, ``disabled`` (no ``--kraken2-db``) or
``failed`` (kraken2, the database or its taxonomy), recorded in the metrics
and never a crash of the run.  A failed screen is never a clean one, so an
unscreened region shows ``.`` in the BED, not 0.

The same functions serve ``kmer-discovery`` (during the run) and
``kmer-nhf`` (afterwards, on the run's files), so the two give the same
columns.
"""

import bisect
import collections
import heapq
import importlib.metadata
import json
import logging
import os
import shutil
import subprocess
import zlib

import pysam

from nonhuman_screen.engine import Kraken2Runner
from nonhuman_screen.result import DOMAINS, TaxonomicFractions

try:  # nonhuman-screen >= 0.1.2 raises on a failed kraken2 run
    from nonhuman_screen.engine import Kraken2Error
except ImportError:  # pragma: no cover - an older engine returns an empty result instead
    class Kraken2Error(RuntimeError):
        """Stand-in for engines without it; the completeness check catches those."""

logger = logging.getLogger(__name__)

#: The two BED columns, appended after ``mapq0_reads``.
NHF_COLUMNS = ("nhf_reads", "nhf")
#: What an unscreened region shows in both columns: not 0.
UNSCREENED = "."
DEFAULT_SAMPLE_READS = 200_000
#: Region classes screened (every region at a BEDPE junction end is too).
SCREENED_CLASSES = frozenset({"SV", "AMBIGUOUS"})
STATUS_SCREENED = "screened"
STATUS_DISABLED = "disabled"
STATUS_FAILED = "failed"
#: The metrics key the block lives under.
METRICS_KEY = "nonhuman_fraction"
SUMMARY_TITLE = "Non-Human Fraction Screen"
_DB_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")
_NONHUMAN_DOMAINS = ("bacterial", "archaeal", "fungal", "protist", "viral")


# ── the database and the engine ────────────────────────────────────


def _nodes_dmp(db_path):
    for candidate in (os.path.join(db_path, "taxonomy", "nodes.dmp"),
                      os.path.join(db_path, "nodes.dmp")):
        if os.path.isfile(candidate):
            return candidate
    return None


def kraken2_db_problems(db_path):
    """What stops the screen before kraken2 runs: a list of reasons, empty when none.

    kraken2 must be on PATH; *db_path* must be a directory holding
    ``hash.k2d``, ``opts.k2d`` and ``taxo.k2d``; and its taxonomy
    (``taxonomy/nodes.dmp``, or ``nodes.dmp`` at the root as PrackenDB has
    it) must be readable, because without it the engine falls back to
    exact-taxid matching and under-counts non-human reads.
    """
    problems = []
    if shutil.which("kraken2") is None:
        problems.append("kraken2 is not on PATH")
    if not db_path:
        problems.append("no Kraken2 database given")
        return problems
    if not os.path.isdir(db_path):
        problems.append(f"Kraken2 database directory not found: {db_path}")
        return problems
    missing = [f for f in _DB_FILES
               if not os.path.isfile(os.path.join(db_path, f))]
    if missing:
        problems.append(
            f"{db_path} lacks {', '.join(missing)} (a Kraken2 database "
            f"directory holds hash.k2d, opts.k2d and taxo.k2d)")
    if _nodes_dmp(db_path) is None:
        problems.append(
            f"no taxonomy/nodes.dmp (nor nodes.dmp) under {db_path}: without "
            f"the taxonomy the engine falls back to exact-taxid matching and "
            f"under-counts non-human reads")
    return problems


def engine_info(db_path=None):
    """Versions of the engine and of kraken2, and the database's k-mer length."""
    info = {}
    try:
        info["nonhuman_screen"] = importlib.metadata.version("nonhuman-screen")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        info["nonhuman_screen"] = "unknown"
    try:
        proc = subprocess.run(["kraken2", "--version"], capture_output=True,
                              text=True, timeout=60)
        first = (proc.stdout or proc.stderr).strip().splitlines()
        info["kraken2"] = first[0] if first else "unknown"
    except (OSError, subprocess.SubprocessError):
        info["kraken2"] = "not found"
    if db_path and os.path.isdir(db_path):
        info["db_kmer_len"] = Kraken2Runner.read_kmer_length(db_path)
    return info


# ── the result block ───────────────────────────────────────────────


def _block(status, reason="", db=None, params=None, engine=None,
           informative_reads=None, sampled=None, regions=()):
    regions = list(regions)
    return {
        "status": status,
        "reason": reason,
        "db": db,
        "engine": dict(engine or {}),
        "params": dict(params or {}),
        "informative_reads": informative_reads,
        "sampled": sampled,
        "regions_screened": len(regions),
        "regions": regions,
    }


def disabled_block(params=None):
    """No ``--kraken2-db``: the reads were not screened, and that is recorded."""
    return _block(STATUS_DISABLED,
                  "no --kraken2-db: the informative reads were not screened",
                  params=params)


def failed_block(reason, db=None, params=None, engine=None):
    """The screen could not be completed; *reason* says why."""
    return _block(STATUS_FAILED, reason, db=db, params=params, engine=engine)


def empty_run_block(db, params=None, engine=None):
    """A run that stopped before anchoring: nothing to classify."""
    return _block(STATUS_SCREENED, "no informative read to classify", db=db,
                  params=params, engine=engine, informative_reads=0,
                  sampled={"reads": 0, "fractions": None})


def screen_params(cluster_distance, sample_reads, confidence, memory_mapping,
                  threads):
    return {
        "cluster_distance": int(cluster_distance),
        "sample_reads": int(sample_reads),
        "confidence": float(confidence),
        "memory_mapping": bool(memory_mapping),
        "threads": int(threads),
    }


def region_values(block):
    """``{(chrom, start, end): (nhf_reads, nhf)}`` for the regions screened."""
    if block.get("status") != STATUS_SCREENED:
        return {}
    return {
        (r["chrom"], r["start"], r["end"]): (r["reads"], r["fractions"]["nonhuman"])
        for r in block.get("regions", [])
    }


def _fmt_fraction(value):
    return f"{value:.4f}"


# ── the reads ──────────────────────────────────────────────────────


def _read_end(read):
    """One id per read end: the two mates of a pair are two sequences."""
    return f"{read.query_name}/{2 if read.is_read2 else 1}"


def _hash_key(text):
    """The sampling key: a hash of the read end's name, so the sample is
    the same whoever draws it and however the BAM is ordered."""
    return zlib.crc32(text.encode()) & 0xFFFFFFFF


def _region_key(region):
    return (region["chrom"], region["start"], region["end"])


def regions_to_screen(regions, junction_ends=()):
    """The regions the screen measures: ``SV``, ``AMBIGUOUS``, or at a junction end."""
    ends = set(junction_ends)
    return [r for r in regions
            if r.get("class") in SCREENED_CLASSES or _region_key(r) in ends]


def collect_reads(informative_bam, regions, cluster_distance, sample_reads):
    """Two passes over the informative BAM: the ids, then the sequences.

    Returns ``(n_ends, sample_ids, region_ids, region_placed, sequences)``:
    the number of read ends in the BAM; the sampled ids (the
    *sample_reads* ends with the smallest hash of their name, or all of
    them); per region (by key) the ids of its reads and of the unmapped
    reads placed in or near it; and ``{id: sequence}`` for every id in the
    sample or a region.  A read end's sequence comes from its primary
    record, or from a supplementary (hard-clipped) record only when that
    is all there is.
    """
    ids = set()
    by_name = {}
    placed_pos = collections.defaultdict(list)
    placed_ids = collections.defaultdict(list)
    with pysam.AlignmentFile(informative_bam, "rb") as bam:
        contigs = set(bam.references)
        bam.reset()
        for read in bam.fetch(until_eof=True):
            rid = _read_end(read)
            ids.add(rid)
            by_name.setdefault(read.query_name, set()).add(rid)
            if read.is_unmapped and read.reference_id >= 0:
                placed_pos[read.reference_name].append(read.reference_start)
                placed_ids[read.reference_name].append(rid)
        for chrom in placed_pos:
            order = sorted(range(len(placed_pos[chrom])),
                           key=placed_pos[chrom].__getitem__)
            placed_pos[chrom] = [placed_pos[chrom][i] for i in order]
            placed_ids[chrom] = [placed_ids[chrom][i] for i in order]

        n_ends = len(ids)
        if sample_reads >= n_ends:
            sample_ids = set(ids)
        else:
            sample_ids = set(heapq.nsmallest(sample_reads, ids, key=_hash_key))

        region_ids = {}
        region_placed = {}
        for region in regions:
            key = _region_key(region)
            got = set()
            placed = set()
            if region["chrom"] in contigs:
                for read in bam.fetch(region["chrom"], region["start"],
                                      region["end"]):
                    rid = _read_end(read)
                    got.add(rid)
                    if read.is_unmapped:
                        placed.add(rid)
                positions = placed_pos.get(region["chrom"], ())
                lo = bisect.bisect_left(positions,
                                        region["start"] - cluster_distance)
                hi = bisect.bisect_left(positions,
                                        region["end"] + cluster_distance)
                for rid in placed_ids[region["chrom"]][lo:hi]:
                    got.add(rid)
                    placed.add(rid)
                for rid in list(got):
                    got |= by_name[rid.rsplit("/", 1)[0]]
            region_ids[key] = got
            region_placed[key] = placed

        wanted = set(sample_ids)
        for got in region_ids.values():
            wanted |= got
        sequences = {}
        supplementary_only = set()
        bam.reset()
        for read in bam.fetch(until_eof=True):
            rid = _read_end(read)
            if rid not in wanted:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            if read.is_supplementary:
                if rid not in sequences:
                    sequences[rid] = seq
                    supplementary_only.add(rid)
            elif rid not in sequences or rid in supplementary_only:
                sequences[rid] = seq
                supplementary_only.discard(rid)
    return n_ends, sample_ids, region_ids, region_placed, sequences


def _has_index(bam_path):
    stem = bam_path[:-4] if bam_path.endswith(".bam") else bam_path
    return any(os.path.isfile(p) for p in (bam_path + ".bai", bam_path + ".csi",
                                            stem + ".bai", stem + ".csi"))


# ── the screen ─────────────────────────────────────────────────────


def _fractions(result, ids):
    return dict(TaxonomicFractions.over_reads(result, ids).to_dict()) if ids else None


def _top_taxon(result, ids, name_map):
    counts = collections.Counter()
    for rid in ids & result.nonhuman_read_names:
        detail = result.per_read_detail.get(rid) or {}
        counts[detail.get("taxid", 0)] += 1
    if not counts:
        return "."
    taxid, n = counts.most_common(1)[0]
    return f"{name_map.get(taxid, str(taxid))}({n})"


def screen(informative_bam, regions, junction_ends, kraken2_db, *,
           cluster_distance=500, sample_reads=DEFAULT_SAMPLE_READS,
           confidence=0.0, threads=1, memory_mapping=False, tmpdir=None):
    """Classify a run's informative reads; return the metrics block.

    Args:
        informative_bam: ``{prefix}.informative.bam`` (indexed).
        regions: The run's regions as dicts with ``chrom``, ``start``,
            ``end``, ``class`` and ``sv_type``, in BED order.
        junction_ends: ``(chrom, start, end)`` of every BEDPE end.
        kraken2_db: The Kraken2 database directory.
        cluster_distance: How far (bp) from a region an unmapped read may
            be placed and still belong to it (the run's
            ``--cluster-distance``).
        sample_reads: At most this many read ends for the run-level
            fractions.
        confidence, threads, memory_mapping: Passed to kraken2 as VCF
            mode passes them.
        tmpdir: Where kraken2's temporary FASTQ goes (disk-backed).

    Returns:
        The block ``kmer-discovery`` writes under ``nonhuman_fraction`` in
        the metrics: ``status`` ``screened`` or ``failed`` with the
        ``reason``, the ``db``, ``engine`` versions, ``params``,
        ``informative_reads`` (read ends), ``sampled`` (``reads`` and
        ``fractions``), ``regions_screened`` and ``regions`` (one row per
        screened region: its coordinates, class and type, ``reads``
        classified, ``placed_unmapped`` among them, the ``fractions`` and
        the ``top_taxon`` among its non-human reads).
    """
    params = screen_params(cluster_distance, sample_reads, confidence,
                           memory_mapping, threads)
    engine = engine_info(kraken2_db)
    problems = kraken2_db_problems(kraken2_db)
    if problems:
        return failed_block("; ".join(problems), db=kraken2_db, params=params,
                            engine=engine)
    if not os.path.isfile(informative_bam):
        return failed_block(f"informative BAM not found: {informative_bam}",
                            db=kraken2_db, params=params, engine=engine)
    if not _has_index(informative_bam):
        return failed_block(f"informative BAM has no index: {informative_bam}",
                            db=kraken2_db, params=params, engine=engine)

    targets = regions_to_screen(regions, junction_ends)
    try:
        n_ends, sample_ids, region_ids, region_placed, sequences = collect_reads(
            informative_bam, targets, cluster_distance, sample_reads)
    except (OSError, ValueError) as exc:
        return failed_block(f"could not read {informative_bam}: {exc}",
                            db=kraken2_db, params=params, engine=engine)
    wanted = set(sample_ids)
    for got in region_ids.values():
        wanted |= got
    missing = len(wanted) - len(sequences)
    if missing:  # a short tally would understate the fraction: never silent
        return failed_block(
            f"{missing} of {len(wanted)} reads to classify have no sequence in "
            f"{informative_bam}", db=kraken2_db, params=params, engine=engine)
    logger.info(
        "[NHF] %d informative read ends; sampling %d; %d regions to screen "
        "with %d reads; %d sequences to classify",
        n_ends, len(sample_ids), len(targets),
        sum(len(v) for v in region_ids.values()), len(sequences))

    if sequences:
        runner = Kraken2Runner(kraken2_db, confidence=confidence,
                               threads=threads, memory_mapping=memory_mapping)
        try:
            result = runner.classify_sequences(sequences, tmpdir=tmpdir)
        except (Kraken2Error, OSError, subprocess.SubprocessError) as exc:
            return failed_block(f"kraken2 failed: {str(exc)[:1500]}",
                                db=kraken2_db, params=params, engine=engine)
        not_classified = sum(1 for rid in sequences
                             if rid not in result.per_read_detail)
        if not_classified:
            return failed_block(
                f"kraken2 returned no classification for {not_classified} of "
                f"{len(sequences)} reads", db=kraken2_db, params=params,
                engine=engine)
        if not result.taxonomy_available:
            return failed_block(
                "the engine could not read the database's taxonomy "
                "(nodes.dmp): the fractions would under-count",
                db=kraken2_db, params=params, engine=engine)
        engine["taxonomy_available"] = True
    else:
        result = Kraken2Runner.Result()
    name_map = Kraken2Runner._load_name_map(kraken2_db) or {}

    rows = []
    for region in targets:
        key = _region_key(region)
        ids = region_ids[key]
        rows.append({
            "chrom": region["chrom"],
            "start": region["start"],
            "end": region["end"],
            "class": region.get("class", "."),
            "sv_type": region.get("sv_type", "."),
            "reads": len(ids),
            "placed_unmapped": len(region_placed[key]),
            "fractions": _fractions(result, ids) or _fractions_of_nothing(),
            "top_taxon": _top_taxon(result, ids, name_map),
        })
    block = _block(
        STATUS_SCREENED, db=kraken2_db, params=params, engine=engine,
        informative_reads=n_ends,
        sampled={"reads": len(sample_ids),
                 "fractions": _fractions(result, sample_ids)},
        regions=rows)
    logger.info("[NHF] %s", describe(block))
    return block


def _fractions_of_nothing():
    """The fractions of an empty read set: every share 0 over 0 reads."""
    return {"total": 0, **{d: 0.0 for d in DOMAINS}}


def run_screen(informative_bam, regions, junction_ends, kraken2_db, *,
               preflight_problems=None, **options):
    """:func:`screen` behind its safety net.

    No database gives the ``disabled`` block; *preflight_problems* found
    at start-up give ``failed`` without running kraken2; and any exception
    the screen raises becomes ``failed`` with its message, because the
    optional screen must never take a finished run down with it.
    """
    params = screen_params(options.get("cluster_distance", 500),
                           options.get("sample_reads", DEFAULT_SAMPLE_READS),
                           options.get("confidence", 0.0),
                           options.get("memory_mapping", False),
                           options.get("threads", 1))
    if kraken2_db is None:
        return disabled_block(params)
    if preflight_problems:
        return failed_block("; ".join(preflight_problems), db=kraken2_db,
                            params=params, engine=engine_info(kraken2_db))
    try:
        return screen(informative_bam, regions, junction_ends, kraken2_db,
                      **options)
    except Exception as exc:  # noqa: BLE001 - recorded, never fatal
        logger.exception("[NHF] the screen failed unexpectedly")
        return failed_block(f"{type(exc).__name__}: {exc}", db=kraken2_db,
                            params=params, engine=engine_info(kraken2_db))


def describe(block):
    """One line on the block, for the log."""
    status = block["status"]
    if status != STATUS_SCREENED:
        return f"non-human fraction screen {status}: {block['reason']}"
    sampled = block["sampled"] or {}
    fractions = sampled.get("fractions") or {}
    positive = sum(1 for r in block["regions"]
                   if r["fractions"]["nonhuman"] > 0)
    nhf = (_fmt_fraction(fractions["nonhuman"]) if fractions else ".")
    return (f"non-human fraction screen: {block['informative_reads']} "
            f"informative read ends; sampled {sampled.get('reads', 0)}, "
            f"non-human {nhf}; {block['regions_screened']} regions screened, "
            f"{positive} with non-human reads")


# ── the run's files ────────────────────────────────────────────────


def read_bed_regions(bed_path):
    """The regions of a discovery BED as dicts keyed by its header names."""
    header = None
    rows = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#chrom"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                raise ValueError(f"{bed_path}: no '#chrom' header line before "
                                 f"the first region")
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            rows.append(row)
    return rows


def read_bedpe_ends(bedpe_path):
    """``(chrom, start, end)`` of both ends of every junction; empty when missing."""
    ends = set()
    if not bedpe_path or not os.path.isfile(bedpe_path):
        return ends
    with open(bedpe_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            ends.add((f[0], int(f[1]), int(f[2])))
            ends.add((f[3], int(f[4]), int(f[5])))
    return ends


def _write_atomic(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def annotate_bed(bed_path, block):
    """Rewrite the ``nhf_reads`` and ``nhf`` columns of a discovery BED in place.

    A BED written before the columns existed gets them appended; every
    region not in *block* (unscreened, or a disabled or failed screen)
    shows ``.`` in both.
    """
    values = region_values(block)
    with open(bed_path) as fh:
        lines = fh.read().splitlines()
    n_base = None
    out = []
    for line in lines:
        if line.startswith("#chrom"):
            columns = [c for c in line[1:].split("\t") if c not in NHF_COLUMNS]
            n_base = len(columns)
            out.append("#" + "\t".join(columns + list(NHF_COLUMNS)))
            continue
        if line.startswith("#") or not line.strip():
            out.append(line)
            continue
        if n_base is None:
            raise ValueError(f"{bed_path}: no '#chrom' header line before "
                             f"the first region")
        fields = line.split("\t")[:n_base]
        found = values.get((fields[0], int(fields[1]), int(fields[2])))
        cells = ((str(found[0]), _fmt_fraction(found[1])) if found
                 else (UNSCREENED, UNSCREENED))
        out.append("\t".join(fields + list(cells)))
    _write_atomic(bed_path, "\n".join(out) + "\n")


def annotate_metrics(metrics, block):
    """Put *block* under ``nonhuman_fraction`` and the two values on each region."""
    metrics[METRICS_KEY] = block
    values = region_values(block)
    for region in metrics.get("regions", []):
        found = values.get((region["chrom"], region["start"], region["end"]))
        region["nhf_reads"] = found[0] if found else None
        region["nhf"] = found[1] if found else None
    return metrics


def summary_lines(block):
    """The screen's section of ``{prefix}.summary.txt``."""
    lines = [SUMMARY_TITLE, "-" * 40]
    status = block["status"]
    if status == STATUS_DISABLED:
        lines.append(f"  Status:                      disabled ({block['reason']})")
    elif status == STATUS_FAILED:
        lines.append("  Status:                      FAILED - the reads were not "
                     "screened; a blank in the BED is not 0")
        lines.append(f"  Reason:                      {block['reason']}")
        if block.get("db"):
            lines.append(f"  Database:                    {block['db']}")
    else:
        sampled = block.get("sampled") or {}
        fractions = sampled.get("fractions") or {}
        lines.append("  Status:                      screened")
        lines.append(f"  Database:                    {block.get('db')}")
        lines.append(f"  Informative read ends:       "
                     f"{block.get('informative_reads') or 0:>8}")
        lines.append(f"  Sampled read ends:           "
                     f"{sampled.get('reads', 0):>8}")
        if fractions:
            lines.append(
                f"  Sampled non-human fraction:  "
                f"{fractions['nonhuman']:>8.4f}   ("
                + ", ".join(f"{d} {fractions[d]:.4f}" for d in _NONHUMAN_DOMAINS)
                + ")")
            lines.append(
                f"  Sampled UniVec Core / human lineage / unclassified: "
                f"{fractions['univec_core']:.4f} / {fractions['human_lineage']:.4f}"
                f" / {fractions['unclassified']:.4f}")
        positive = [r for r in block["regions"] if r["fractions"]["nonhuman"] > 0]
        lines.append(f"  Regions screened:            "
                     f"{block['regions_screened']:>8}   (SV and AMBIGUOUS "
                     f"regions, and junction ends; one read end = one read)")
        lines.append(f"  Regions with non-human reads:{len(positive):>8}")
        if block["regions"]:
            lines.append("")
            lines.append(
                f"  {'Region':<30s} {'Class':>9s} {'Type':>5s} {'Reads':>6s}"
                f" {'Placed':>6s} {'NHF':>7s} {'Bact':>7s} {'Arch':>7s}"
                f" {'Fungi':>7s} {'Prot':>7s} {'Viral':>7s} {'UniVec':>7s}"
                f" {'Human':>7s} {'Unclass':>7s}  Top non-human taxon")
            lines.append(
                f"  {'------':<30s} {'-----':>9s} {'----':>5s} {'-----':>6s}"
                f" {'------':>6s} {'---':>7s} {'----':>7s} {'----':>7s}"
                f" {'-----':>7s} {'----':>7s} {'-----':>7s} {'------':>7s}"
                f" {'-----':>7s} {'-------':>7s}  -------------------")
            for r in block["regions"]:
                fr = r["fractions"]
                label = f"{r['chrom']}:{r['start'] + 1}-{r['end']}"
                lines.append(
                    f"  {label:<30s} {r['class']:>9s} {r['sv_type']:>5s}"
                    f" {r['reads']:>6d} {r['placed_unmapped']:>6d}"
                    f" {fr['nonhuman']:>7.4f} {fr['bacterial']:>7.4f}"
                    f" {fr['archaeal']:>7.4f} {fr['fungal']:>7.4f}"
                    f" {fr['protist']:>7.4f} {fr['viral']:>7.4f}"
                    f" {fr['univec_core']:>7.4f} {fr['human_lineage']:>7.4f}"
                    f" {fr['unclassified']:>7.4f}  {r['top_taxon']}")
    lines.append("")
    return lines


def update_summary(summary_path, block):
    """Rewrite the screen's section of an existing summary in place.

    The section sits last, before the closing rule; an older summary
    without one gets it there.  The result is what the run itself would
    have written.
    """
    with open(summary_path) as fh:
        lines = fh.read().splitlines()
    rule = "=" * 60
    closing = max((i for i, line in enumerate(lines) if line == rule),
                  default=None)
    if closing is None or closing == 0:
        raise ValueError(f"{summary_path}: not a discovery summary "
                         f"(no closing rule)")
    start = next((i for i, line in enumerate(lines[:closing])
                  if line == SUMMARY_TITLE), None)
    if start is None:
        start = closing
    head = lines[:start]
    if head and head[-1] != "":
        head.append("")
    new_lines = head + summary_lines(block) + lines[closing:]
    _write_atomic(summary_path, "\n".join(new_lines) + "\n")


def rescreen_run(out_prefix, kraken2_db, *, bedpe_path=None, **options):
    """``kmer-nhf``: the screen on an existing run's files, rewritten in place.

    Reads ``{prefix}.bed``, ``{prefix}.sv.bedpe`` (or *bedpe_path*),
    ``{prefix}.metrics.json`` and ``{prefix}.informative.bam``; rewrites
    the BED's two columns, the metrics block (and the two values on each
    region) and the summary's section.  Returns the block.
    """
    bed_path = f"{out_prefix}.bed"
    metrics_path = f"{out_prefix}.metrics.json"
    summary_path = f"{out_prefix}.summary.txt"
    informative_bam = f"{out_prefix}.informative.bam"
    bedpe_path = bedpe_path or f"{out_prefix}.sv.bedpe"
    for required in (bed_path, metrics_path):
        if not os.path.isfile(required):
            raise FileNotFoundError(required)
    regions = read_bed_regions(bed_path)
    junction_ends = read_bedpe_ends(bedpe_path)
    with open(metrics_path) as fh:
        metrics = json.load(fh)

    if not os.path.isfile(informative_bam) and not regions \
            and not metrics.get("informative_reads"):
        # A run that stopped before anchoring wrote no BAM: nothing to classify.
        params = screen_params(options.get("cluster_distance", 500),
                               options.get("sample_reads", DEFAULT_SAMPLE_READS),
                               options.get("confidence", 0.0),
                               options.get("memory_mapping", False),
                               options.get("threads", 1))
        problems = kraken2_db_problems(kraken2_db)
        block = (failed_block("; ".join(problems), db=kraken2_db, params=params,
                              engine=engine_info(kraken2_db)) if problems
                 else empty_run_block(kraken2_db, params, engine_info(kraken2_db)))
    else:
        block = run_screen(informative_bam, regions, junction_ends, kraken2_db,
                           **options)
    annotate_bed(bed_path, block)
    annotate_metrics(metrics, block)
    _write_atomic(metrics_path, json.dumps(metrics, indent=2) + "\n")
    if os.path.isfile(summary_path):
        update_summary(summary_path, block)
    return block
