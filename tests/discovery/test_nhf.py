"""The non-human fraction (NHF) screen of discovery mode.

The screen classifies a run's informative reads with kraken2 through the
``nonhuman-screen`` engine (what VCF mode uses for ``DKU_NHF`` /
``DKA_NHF``) and annotates every ``SV`` and ``AMBIGUOUS`` region with the
share of its reads that are non-human.  These tests run it on the
simulated trio of ``examples/sv_demo`` against a tiny Kraken2 database
built here from the demo's own sequences: the two reference chromosomes
as *Homo sapiens* (taxid 9606) and the child's 250 bp novel insertion as
*Escherichia coli* (taxid 562).  A read of a bacterium carries k-mers in
neither parent nor the reference, so it is informative by construction and
looks like SV evidence; here the reads inside the insertion are exactly
that, and the insertion's region must be the one region with non-human
reads.

``kraken2`` and ``kraken2-build`` must be on PATH for the tests that
classify; without them those are skipped with the reason.  The tests of
the read-set definition, the BED and summary rewriting and the failure
paths that never reach kraken2 run everywhere.
"""

import importlib.util
import json
import os
import random
import shutil
import subprocess

import pysam
import pytest

from kmer_denovo_filter.cli import (
    nhf_main,
    parse_discovery_args,
)
from kmer_denovo_filter.discovery import nhf
from kmer_denovo_filter.discovery.pipeline import run_discovery_pipeline

DEMO_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, "examples", "sv_demo",
)
HUMAN_TAXID = 9606
ECOLI_TAXID = 562
#: A minimal NCBI-style taxonomy: root, Bacteria with E. coli, Eukaryota
#: with Metazoa and Homo sapiens.  The engine's human-lineage rules need
#: the path from 9606 to the root.
TAXONOMY_NODES = [
    (1, 1, "no rank"), (131567, 1, "no rank"),
    (2, 131567, "superkingdom"), (ECOLI_TAXID, 2, "species"),
    (2759, 131567, "superkingdom"), (33208, 2759, "kingdom"),
    (HUMAN_TAXID, 33208, "species"),
]
TAXONOMY_NAMES = {
    1: "root", 131567: "cellular organisms", 2: "Bacteria",
    ECOLI_TAXID: "Escherichia coli", 2759: "Eukaryota", 33208: "Metazoa",
    HUMAN_TAXID: "Homo sapiens",
}
#: The 250 bp insertion's region in the demo's expected BED (0-based).
INSERTION_REGION = ("chr2", 21866, 22078)

needs_kraken2 = pytest.mark.skipif(
    shutil.which("kraken2") is None or shutil.which("kraken2-build") is None,
    reason="kraken2 and kraken2-build are needed to build and query the "
           "test database",
)


# ── fixtures ───────────────────────────────────────────────────────


def _simulator():
    spec = importlib.util.spec_from_file_location(
        "simulate_sv_trio", os.path.join(DEMO_DIR, "simulate_sv_trio.py"),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def demo_trio(tmp_path_factory):
    """The simulated trio of examples/sv_demo."""
    outdir = tmp_path_factory.mktemp("sv_demo")
    _simulator().main(["--outdir", str(outdir)])
    return outdir


@pytest.fixture(scope="module")
def demo_sequences(demo_trio):
    """The demo's reference and its 250 bp novel insertion, from the simulator's seed."""
    sim = _simulator()
    rng = random.Random(1)
    ref = {c: sim.random_seq(rng, n) for c, n in sim.CHROM_SIZES.items()}
    with pysam.FastaFile(str(demo_trio / "ref.fa")) as fa:
        assert all(fa.fetch(c) == ref[c] for c in ref)
    novel = [h.seq[s:e] for h in sim.child_haplotypes(ref, rng)
             for (s, e, chrom, *_rest) in h.segments if chrom is None]
    (ins250,) = [seq for seq in novel if len(seq) == 250]
    return ref, ins250


def _build_kraken2_db(db_dir, sequences):
    """A Kraken2 database of *sequences* (``{header: seq}``) with the toy taxonomy."""
    taxonomy = db_dir / "taxonomy"
    taxonomy.mkdir(parents=True)
    with open(taxonomy / "nodes.dmp", "w") as fh:
        for taxid, parent, rank in TAXONOMY_NODES:
            fh.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
    with open(taxonomy / "names.dmp", "w") as fh:
        for taxid, name in TAXONOMY_NAMES.items():
            fh.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    library = db_dir / "library.fa"
    with open(library, "w") as fh:
        for header, seq in sequences.items():
            fh.write(f">{header}\n{seq}\n")
    for cmd in (
        ["kraken2-build", "--add-to-library", str(library), "--db", str(db_dir),
         "--no-masking"],
        ["kraken2-build", "--build", "--db", str(db_dir), "--threads", "1"],
    ):
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            pytest.fail(f"{' '.join(cmd)} failed:\n{proc.stdout}\n{proc.stderr}")
    assert nhf.kraken2_db_problems(str(db_dir)) == []
    return str(db_dir)


@pytest.fixture(scope="module")
def kraken2_db(tmp_path_factory, demo_sequences):
    """The demo reference as human and its 250 bp insertion as E. coli."""
    if shutil.which("kraken2") is None or shutil.which("kraken2-build") is None:
        pytest.skip("kraken2 and kraken2-build are needed to build the test "
                    "database")
    ref, ins250 = demo_sequences
    sequences = {f"{c}|kraken:taxid|{HUMAN_TAXID} human": ref[c] for c in ref}
    sequences[f"ins250|kraken:taxid|{ECOLI_TAXID} insertion"] = ins250
    return _build_kraken2_db(tmp_path_factory.mktemp("kraken2_db"), sequences)


@pytest.fixture(scope="module")
def human_only_db(tmp_path_factory, demo_sequences):
    """The demo reference alone: nothing in it is non-human."""
    if shutil.which("kraken2") is None or shutil.which("kraken2-build") is None:
        pytest.skip("kraken2 and kraken2-build are needed to build the test "
                    "database")
    ref, _ins250 = demo_sequences
    sequences = {f"{c}|kraken:taxid|{HUMAN_TAXID} human": ref[c] for c in ref}
    return _build_kraken2_db(tmp_path_factory.mktemp("human_db"), sequences)


def _run_discovery(trio, prefix, *extra):
    run_discovery_pipeline(parse_discovery_args([
        "--child", str(trio / "child.bam"),
        "--mother", str(trio / "mother.bam"),
        "--father", str(trio / "father.bam"),
        "--ref-fasta", str(trio / "ref.fa"),
        "--out-prefix", str(prefix),
        "--threads", "2",
        *extra,
    ]))
    return str(prefix)


@pytest.fixture(scope="module")
def screened_run(demo_trio, kraken2_db, tmp_path_factory):
    """kmer-discovery with --kraken2-db and --report."""
    outdir = tmp_path_factory.mktemp("screened")
    return _run_discovery(demo_trio, outdir / "sv_demo",
                          "--kraken2-db", kraken2_db,
                          "--report", str(outdir / "sv_demo.report.html"))


@pytest.fixture(scope="module")
def disabled_run(demo_trio, tmp_path_factory):
    """kmer-discovery without --kraken2-db."""
    outdir = tmp_path_factory.mktemp("disabled")
    return _run_discovery(demo_trio, outdir / "sv_demo")


# ── helpers ────────────────────────────────────────────────────────


def _bed(prefix):
    """Rows of the BED as dicts keyed by the header names."""
    return nhf.read_bed_regions(prefix + ".bed")


def _metrics(prefix):
    with open(prefix + ".metrics.json") as fh:
        return json.load(fh)


def _read(path):
    with open(path) as fh:
        return fh.read()


def _copy_run(prefix, dest_dir):
    """Copy a run's output files to *dest_dir*; return the new prefix."""
    src_dir, name = os.path.split(prefix)
    for entry in os.listdir(src_dir):
        if entry.startswith(name + "."):
            shutil.copy(os.path.join(src_dir, entry), dest_dir / entry)
    return str(dest_dir / name)


# ── the screened run ───────────────────────────────────────────────


@needs_kraken2
class TestScreenedRun:
    def test_bed_columns(self, screened_run):
        rows = _bed(screened_run)
        header = _read(screened_run + ".bed").splitlines()[1]
        assert header.endswith("\tmax_mapq\tmapq0_reads\tnhf_reads\tnhf")
        assert len(rows) == 11
        for row in rows:
            if row["class"] == "SMALL":
                assert (row["nhf_reads"], row["nhf"]) == (".", ".")
                continue
            # a screened region's reads: at least the BED's own count, as
            # mates and placed unmapped reads come on top of it
            assert int(row["nhf_reads"]) >= int(row["reads"])
            assert 0.0 <= float(row["nhf"]) <= 1.0
        by_key = {(r["chrom"], r["start"], r["end"]): r for r in rows}
        insertion = by_key[INSERTION_REGION]
        assert float(insertion["nhf"]) > 0.5, insertion
        # 7 reads placed in the region plus 7 within --cluster-distance of it
        assert int(insertion["nhf_reads"]) > int(insertion["reads"]) + 7
        for key, row in by_key.items():
            if key != INSERTION_REGION and row["class"] != "SMALL":
                assert row["nhf"] == "0.0000", row

    def test_metrics_block(self, screened_run, kraken2_db):
        metrics = _metrics(screened_run)
        block = metrics["nonhuman_fraction"]
        assert block["status"] == "screened" and block["reason"] == ""
        assert block["db"] == kraken2_db
        assert block["engine"]["taxonomy_available"] is True
        assert block["engine"]["db_kmer_len"] == 35
        assert "kraken2" in block["engine"] and "nonhuman_screen" in block["engine"]
        assert block["params"] == {
            "cluster_distance": 500, "sample_reads": 200000, "confidence": 0.0,
            "memory_mapping": False, "threads": 2,
        }
        assert block["informative_reads"] > 0
        # fewer informative reads than the sample size: all are sampled
        assert block["sampled"]["reads"] == block["informative_reads"]
        fractions = block["sampled"]["fractions"]
        assert fractions["total"] == block["sampled"]["reads"]
        assert abs(fractions["nonhuman"] + fractions["univec_core"]
                   + fractions["human_lineage"] + fractions["unclassified"]
                   - 1.0) < 1e-3
        assert fractions["nonhuman"] == fractions["bacterial"] > 0

        screened = [r for r in metrics["regions"] if r["class"] != "SMALL"]
        assert block["regions_screened"] == len(block["regions"]) == len(screened)
        # the block's rows and the BED's two columns agree, region by region
        bed = {(r["chrom"], r["start"], r["end"]): r for r in _bed(screened_run)}
        for row in block["regions"]:
            key = (row["chrom"], row["start"], row["end"])
            assert bed[key]["nhf_reads"] == str(row["reads"])
            assert bed[key]["nhf"] == f"{row['fractions']['nonhuman']:.4f}"
            assert row["fractions"]["total"] == row["reads"]
        positive = [r for r in block["regions"] if r["fractions"]["nonhuman"] > 0]
        assert [(r["chrom"], r["start"], r["end"]) for r in positive] == [INSERTION_REGION]
        assert positive[0]["top_taxon"].startswith("Escherichia_coli(")
        assert positive[0]["placed_unmapped"] == 14
        # the two values ride on the metrics' own region rows too
        for region in metrics["regions"]:
            key = (region["chrom"], region["start"], region["end"])
            if region["class"] == "SMALL":
                assert region["nhf_reads"] is None and region["nhf"] is None
            else:
                assert region["nhf_reads"] == int(bed[key]["nhf_reads"])
                assert region["nhf"] == float(bed[key]["nhf"])

    def test_summary_and_report(self, screened_run):
        summary = _read(screened_run + ".summary.txt")
        section = summary[summary.index("Non-Human Fraction Screen"):]
        assert "Status:                      screened" in section
        assert "Regions with non-human reads:       1" in section
        assert "Escherichia_coli(" in section
        assert section.rstrip().endswith("=" * 60)
        report = _read(screened_run + ".report.html")
        assert "Non-human fraction screen" in report
        assert "Escherichia_coli(" in report
        assert "badge badge-orange" in report


# ── disabled and failed ────────────────────────────────────────────


class TestDisabledAndFailed:
    def test_disabled(self, disabled_run):
        rows = _bed(disabled_run)
        assert rows and all((r["nhf_reads"], r["nhf"]) == (".", ".") for r in rows)
        metrics = _metrics(disabled_run)
        block = metrics["nonhuman_fraction"]
        assert block["status"] == "disabled"
        assert "no --kraken2-db" in block["reason"]
        assert block["db"] is None and block["regions"] == []
        assert all(r["nhf"] is None and r["nhf_reads"] is None
                   for r in metrics["regions"])
        summary = _read(disabled_run + ".summary.txt")
        assert "Status:                      disabled (no --kraken2-db" in summary

    def test_unusable_db_is_a_failed_screen_not_a_crash(self, demo_trio, tmp_path,
                                                       monkeypatch):
        """An empty directory as the database: recorded, and the run exits 0."""
        empty_db = tmp_path / "not_a_db"
        empty_db.mkdir()

        def _never(*_args, **_kwargs):
            raise AssertionError("kraken2 must not run on an unusable database")
        monkeypatch.setattr(nhf.Kraken2Runner, "classify_sequences", _never)
        prefix = _run_discovery(demo_trio, tmp_path / "sv_demo",
                                "--kraken2-db", str(empty_db))
        block = _metrics(prefix)["nonhuman_fraction"]
        assert block["status"] == "failed"
        assert "hash.k2d" in block["reason"] and "nodes.dmp" in block["reason"]
        assert block["db"] == str(empty_db)
        assert all((r["nhf_reads"], r["nhf"]) == (".", ".") for r in _bed(prefix))
        summary = _read(prefix + ".summary.txt")
        assert "FAILED - the reads were not screened" in summary
        assert "hash.k2d" in summary

    @needs_kraken2
    @pytest.mark.parametrize("failure", ["raises", "returns nothing"])
    def test_kraken2_failing_is_a_failed_screen(self, demo_trio, kraken2_db,
                                                tmp_path, monkeypatch, failure):
        """kraken2 raising, or an engine that returns no classification."""
        if failure == "raises":
            def _classify(self, sequences, tmpdir=None):
                raise nhf.Kraken2Error("kraken2 exited with code 2")
        else:
            def _classify(self, sequences, tmpdir=None):
                result = nhf.Kraken2Runner.Result()
                result.total = len(sequences)
                return result
        monkeypatch.setattr(nhf.Kraken2Runner, "classify_sequences", _classify)
        prefix = _run_discovery(demo_trio, tmp_path / "sv_demo",
                                "--kraken2-db", kraken2_db)
        block = _metrics(prefix)["nonhuman_fraction"]
        assert block["status"] == "failed"
        if failure == "raises":
            assert block["reason"] == "kraken2 failed: kraken2 exited with code 2"
        else:
            assert block["reason"].startswith(
                "kraken2 returned no classification for")
        assert all((r["nhf_reads"], r["nhf"]) == (".", ".") for r in _bed(prefix))

    def test_early_exit_records_the_screen(self, demo_trio, tmp_path):
        """A child with no proband-unique k-mers: nothing to classify, recorded."""
        prefix = str(tmp_path / "none")
        run_discovery_pipeline(parse_discovery_args([
            "--child", str(demo_trio / "mother.bam"),
            "--mother", str(demo_trio / "mother.bam"),
            "--father", str(demo_trio / "father.bam"),
            "--ref-fasta", str(demo_trio / "ref.fa"),
            "--out-prefix", prefix, "--threads", "2",
        ]))
        metrics = _metrics(prefix)
        assert metrics["candidate_regions"] == 0
        assert metrics["nonhuman_fraction"]["status"] == "disabled"
        header = _read(prefix + ".bed").splitlines()[-1]
        assert header.endswith("\tnhf_reads\tnhf")
        assert "Non-Human Fraction Screen" in _read(prefix + ".summary.txt")


# ── kmer-nhf: the screen after the fact ────────────────────────────


@needs_kraken2
class TestRescreen:
    def test_reproduces_the_run(self, disabled_run, screened_run, kraken2_db,
                                tmp_path):
        """kmer-nhf on a run made without --kraken2-db gives the run's own output."""
        prefix = _copy_run(disabled_run, tmp_path)
        nhf_main(["--out-prefix", prefix, "--kraken2-db", kraken2_db,
                  "--threads", "2", "--report", prefix + ".report.html"])
        assert _read(prefix + ".bed") == _read(screened_run + ".bed")
        assert _read(prefix + ".summary.txt") == _read(screened_run + ".summary.txt")
        assert _metrics(prefix) == _metrics(screened_run)
        assert "Escherichia_coli(" in _read(prefix + ".report.html")

    def test_a_new_database_rewrites_in_place(self, screened_run, human_only_db,
                                              tmp_path):
        """Against a database without the insertion, its reads are unclassified."""
        prefix = _copy_run(screened_run, tmp_path)
        nhf_main(["--out-prefix", prefix, "--kraken2-db", human_only_db,
                  "--threads", "2"])
        rows = {(r["chrom"], r["start"], r["end"]): r for r in _bed(prefix)}
        assert rows[INSERTION_REGION]["nhf"] == "0.0000"
        block = _metrics(prefix)["nonhuman_fraction"]
        assert block["db"] == human_only_db
        (insertion,) = [r for r in block["regions"]
                        if (r["chrom"], r["start"], r["end"]) == INSERTION_REGION]
        assert insertion["fractions"]["unclassified"] > 0.5
        assert insertion["top_taxon"] == "."
        # the header is not duplicated, and the summary section is replaced
        bed_lines = _read(prefix + ".bed").splitlines()
        assert bed_lines[1].count("nhf_reads") == 1
        assert _read(prefix + ".summary.txt").count("Non-Human Fraction Screen") == 1

    def test_an_old_run_gets_the_columns(self, disabled_run, screened_run,
                                         kraken2_db, tmp_path):
        """A BED and summary from before the screen existed: columns appended."""
        prefix = _copy_run(disabled_run, tmp_path)
        bed_lines = _read(prefix + ".bed").splitlines()
        with open(prefix + ".bed", "w") as fh:
            for line in bed_lines:
                fields = line.split("\t")
                fh.write("\t".join(fields[:-2] if line.startswith("#chrom")
                                   or not line.startswith("#") else fields) + "\n")
        summary = _read(prefix + ".summary.txt")
        start = summary.index("\nNon-Human Fraction Screen")
        end = summary.rindex("=" * 60)
        with open(prefix + ".summary.txt", "w") as fh:
            fh.write(summary[:start] + "\n" + summary[end:])
        metrics = _metrics(prefix)
        del metrics["nonhuman_fraction"]
        for region in metrics["regions"]:
            del region["nhf"], region["nhf_reads"]
        with open(prefix + ".metrics.json", "w") as fh:
            json.dump(metrics, fh, indent=2)
        assert "nhf" not in _read(prefix + ".bed").splitlines()[1]
        nhf_main(["--out-prefix", prefix, "--kraken2-db", kraken2_db,
                  "--threads", "2"])
        assert _read(prefix + ".bed") == _read(screened_run + ".bed")
        assert _read(prefix + ".summary.txt") == _read(screened_run + ".summary.txt")
        assert _metrics(prefix) == _metrics(screened_run)

    def test_failed_screen_exits_1_and_is_recorded(self, disabled_run, tmp_path):
        prefix = _copy_run(disabled_run, tmp_path)
        empty_db = tmp_path / "not_a_db"
        empty_db.mkdir()
        with pytest.raises(SystemExit) as exc:
            nhf_main(["--out-prefix", prefix, "--kraken2-db", str(empty_db)])
        assert exc.value.code == 1
        block = _metrics(prefix)["nonhuman_fraction"]
        assert block["status"] == "failed" and "hash.k2d" in block["reason"]
        assert "FAILED" in _read(prefix + ".summary.txt")

    def test_not_a_run_exits_2(self, tmp_path, kraken2_db):
        with pytest.raises(SystemExit) as exc:
            nhf_main(["--out-prefix", str(tmp_path / "nothing"),
                      "--kraken2-db", kraken2_db])
        assert exc.value.code == 2


# ── the definitions, without kraken2 ───────────────────────────────


def _placed_unmapped(name, chrom_id, pos, read2=False):
    seg = pysam.AlignedSegment()
    seg.query_name = name
    seg.query_sequence = "ACGT" * 20
    seg.query_qualities = pysam.qualitystring_to_array("I" * 80)
    seg.flag = 0x1 | 0x4 | (0x80 if read2 else 0x40)
    seg.reference_id = chrom_id
    seg.reference_start = pos
    seg.next_reference_id = chrom_id
    seg.next_reference_start = pos
    seg.mapping_quality = 0
    return seg


def _mapped(name, chrom_id, pos, read2=False, mate_pos=None, supplementary=False):
    seg = pysam.AlignedSegment()
    seg.query_name = name
    seg.query_sequence = "ACGT" * 20
    seg.query_qualities = pysam.qualitystring_to_array("I" * 80)
    seg.flag = 0x1 | (0x80 if read2 else 0x40) | (0x800 if supplementary else 0)
    seg.reference_id = chrom_id
    seg.reference_start = pos
    seg.cigar = [(0, 80)]
    seg.mapping_quality = 60
    seg.next_reference_id = chrom_id
    seg.next_reference_start = mate_pos if mate_pos is not None else pos
    return seg


def _write_bam(path, segments):
    header = pysam.AlignmentHeader.from_references(["chr1"], [20000])
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for seg in segments:
            bam.write(seg)
    pysam.sort("-o", path, path)
    pysam.index(path)


class _FakeRunner:
    """Classifies by read name: a name holding 'bug' is non-human."""

    calls = []

    def __init__(self, db_path, **kwargs):
        self.db_path = db_path

    def classify_sequences(self, sequences, tmpdir=None):
        self.calls.append(dict(sequences))
        result = nhf.Kraken2Runner.Result()
        result.total = len(sequences)
        for rid in sequences:
            if "bug" in rid:
                result.nonhuman_read_names.add(rid)
                result.bacterial_read_names.add(rid)
                taxid, domain, guard, nonhuman = ECOLI_TAXID, "Bacteria", "PASS", True
            else:
                result.human_lineage_read_names.add(rid)
                taxid, domain, guard, nonhuman = HUMAN_TAXID, "Human", "HUMAN", False
            result.classified += 1
            result.per_read_detail[rid] = {
                "status": "C", "taxid": taxid, "domain": domain,
                "guard_status": guard, "is_nonhuman": nonhuman, "kmer_string": "",
            }
        result.nonhuman_count = len(result.nonhuman_read_names)
        result.bacterial_count = result.nonhuman_count
        result.human_lineage_count = len(result.human_lineage_read_names)
        result.taxonomy_available = True
        return result

    Result = nhf.Kraken2Runner.Result

    @staticmethod
    def _load_name_map(db_path):
        return {ECOLI_TAXID: "Escherichia_coli"}


@pytest.fixture()
def fake_engine(monkeypatch):
    """The screen's definitions without kraken2 or a database."""
    _FakeRunner.calls = []
    monkeypatch.setattr(nhf, "Kraken2Runner", _FakeRunner)
    monkeypatch.setattr(nhf, "kraken2_db_problems", lambda db: [])
    monkeypatch.setattr(nhf, "engine_info", lambda db=None: {"fake": True})
    return _FakeRunner


class TestReadSets:
    def test_a_region_s_reads(self, tmp_path, fake_engine):
        """Overlapping alignments, their mates, and placed reads within the slack."""
        bam = str(tmp_path / "informative.bam")
        _write_bam(bam, [
            _mapped("inside", 0, 1050),                       # overlaps the region
            _mapped("pair", 0, 1100, mate_pos=5000),           # overlaps; its mate is far
            _mapped("pair", 0, 5000, read2=True, mate_pos=1100),  # the other end: comes along
            _placed_unmapped("bug_near", 0, 1300),             # placed within 500 bp
            _placed_unmapped("bug_far", 0, 2000),              # placed beyond 500 bp
            _mapped("elsewhere", 0, 8000),                     # not this region
            _mapped("bug_other", 0, 8100, supplementary=True),  # a lone supplementary record
        ])
        regions = [
            {"chrom": "chr1", "start": 1000, "end": 1200, "class": "SV", "sv_type": "INS"},
            {"chrom": "chr1", "start": 7900, "end": 8200, "class": "SMALL", "sv_type": "."},
            {"chrom": "chr1", "start": 7900, "end": 8200, "class": "AMBIGUOUS", "sv_type": "."},
        ]
        block = nhf.screen(bam, regions[:2], set(), "fake-db", cluster_distance=500,
                           sample_reads=200000)
        assert block["status"] == "screened"
        assert block["informative_reads"] == 7
        assert block["sampled"]["reads"] == 7
        assert block["sampled"]["fractions"]["nonhuman"] == pytest.approx(3 / 7, abs=1e-4)
        (row,) = block["regions"]           # the SMALL region is not screened
        assert row["reads"] == 4 and row["placed_unmapped"] == 1
        assert row["fractions"]["nonhuman"] == 0.25
        assert row["top_taxon"] == "Escherichia_coli(1)"
        classified = fake_engine.calls[-1]
        assert set(classified) == {"inside/1", "pair/1", "pair/2", "bug_near/1",
                                   "bug_far/1", "elsewhere/1", "bug_other/1"}

        # a junction end is screened whatever its class; AMBIGUOUS always is
        block = nhf.screen(bam, regions[:2], {("chr1", 7900, 8200)}, "fake-db")
        assert [r["class"] for r in block["regions"]] == ["SV", "SMALL"]
        block = nhf.screen(bam, [regions[0], regions[2]], set(), "fake-db")
        assert [r["class"] for r in block["regions"]] == ["SV", "AMBIGUOUS"]
        assert block["regions"][1]["reads"] == 2            # elsewhere/1, bug_other/1

    def test_the_sample_is_deterministic_and_exact(self, tmp_path, fake_engine):
        bam = str(tmp_path / "informative.bam")
        _write_bam(bam, [_mapped(f"read{i}", 0, 100 * i) for i in range(1, 41)])
        blocks = [nhf.screen(bam, [], set(), "fake-db", sample_reads=10)
                  for _ in range(2)]
        assert blocks[0] == blocks[1]
        assert blocks[0]["sampled"]["reads"] == 10
        assert blocks[0]["informative_reads"] == 40
        assert set(fake_engine.calls[-1]) == set(fake_engine.calls[-2])
        assert len(fake_engine.calls[-1]) == 10
        # the same reads whatever the BAM's order: a hash of the name picks them
        shuffled = str(tmp_path / "shuffled.bam")
        _write_bam(shuffled, [_mapped(f"read{i}", 0, 100 * (41 - i))
                              for i in range(1, 41)])
        nhf.screen(shuffled, [], set(), "fake-db", sample_reads=10)
        assert set(fake_engine.calls[-1]) == set(fake_engine.calls[-2])
        for n, expected in ((40, 40), (100, 40)):
            full = nhf.screen(bam, [], set(), "fake-db", sample_reads=n)
            assert full["sampled"]["reads"] == expected

    def test_missing_bam_or_index_is_a_failed_screen(self, tmp_path, fake_engine):
        block = nhf.screen(str(tmp_path / "none.bam"), [], set(), "fake-db")
        assert block["status"] == "failed" and "not found" in block["reason"]
        bam = str(tmp_path / "informative.bam")
        _write_bam(bam, [_mapped("read", 0, 100)])
        os.remove(bam + ".bai")
        block = nhf.screen(bam, [], set(), "fake-db")
        assert block["status"] == "failed" and "no index" in block["reason"]

    def test_run_screen_never_raises(self, tmp_path, fake_engine, monkeypatch):
        def _boom(*_args, **_kwargs):
            raise RuntimeError("disk on fire")
        monkeypatch.setattr(nhf, "collect_reads", _boom)
        bam = str(tmp_path / "informative.bam")
        _write_bam(bam, [_mapped("read", 0, 100)])
        block = nhf.run_screen(bam, [], set(), "fake-db")
        assert block["status"] == "failed"
        assert block["reason"] == "RuntimeError: disk on fire"
        assert nhf.run_screen(bam, [], set(), None)["status"] == "disabled"
        block = nhf.run_screen(bam, [], set(), "fake-db",
                               preflight_problems=["kraken2 is not on PATH"])
        assert block["status"] == "failed" and block["reason"] == "kraken2 is not on PATH"


class TestDatabaseCheck:
    def test_problems_are_named(self, tmp_path, monkeypatch):
        monkeypatch.setattr(nhf.shutil, "which", lambda name: "/usr/bin/kraken2")
        assert nhf.kraken2_db_problems(None) == ["no Kraken2 database given"]
        assert nhf.kraken2_db_problems(str(tmp_path / "nope")) == [
            f"Kraken2 database directory not found: {tmp_path / 'nope'}"]
        db = tmp_path / "db"
        (db / "taxonomy").mkdir(parents=True)
        for name in ("hash.k2d", "opts.k2d"):
            (db / name).write_bytes(b"x")
        problems = nhf.kraken2_db_problems(str(db))
        assert len(problems) == 2
        assert problems[0].startswith(f"{db} lacks taxo.k2d")
        assert problems[1].startswith("no taxonomy/nodes.dmp")
        (db / "taxo.k2d").write_bytes(b"x")
        (db / "taxonomy" / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n")
        assert nhf.kraken2_db_problems(str(db)) == []
        monkeypatch.setattr(nhf.shutil, "which", lambda name: None)
        assert nhf.kraken2_db_problems(str(db)) == ["kraken2 is not on PATH"]


class TestFileRewrites:
    BED = ("#filters: min_distinct_kmers=1\n"
           "#chrom\tstart\tend\treads\tclass\tsv_type\n"
           "chr1\t100\t200\t5\tSV\tDEL\n"
           "chr1\t900\t950\t2\tSMALL\t.\n")

    def _block(self):
        return {
            "status": "screened", "reason": "", "db": "db", "engine": {},
            "params": {}, "informative_reads": 7,
            "sampled": {"reads": 7, "fractions": None}, "regions_screened": 1,
            "regions": [{
                "chrom": "chr1", "start": 100, "end": 200, "class": "SV",
                "sv_type": "DEL", "reads": 6, "placed_unmapped": 0,
                "fractions": {"total": 6, **{d: 0.0 for d in nhf.DOMAINS},
                              "nonhuman": 0.5, "viral": 0.5,
                              "human_lineage": 0.5},
                "top_taxon": "Virus(3)",
            }],
        }

    def test_annotate_bed_appends_then_replaces(self, tmp_path):
        bed = tmp_path / "run.bed"
        bed.write_text(self.BED)
        nhf.annotate_bed(str(bed), self._block())
        lines = bed.read_text().splitlines()
        assert lines[0] == "#filters: min_distinct_kmers=1"
        assert lines[1] == "#chrom\tstart\tend\treads\tclass\tsv_type\tnhf_reads\tnhf"
        assert lines[2] == "chr1\t100\t200\t5\tSV\tDEL\t6\t0.5000"
        assert lines[3] == "chr1\t900\t950\t2\tSMALL\t.\t.\t."
        nhf.annotate_bed(str(bed), nhf.disabled_block())
        lines = bed.read_text().splitlines()
        assert lines[1].count("nhf") == 2
        assert lines[2] == "chr1\t100\t200\t5\tSV\tDEL\t.\t."

    def test_annotate_metrics(self):
        metrics = {"regions": [
            {"chrom": "chr1", "start": 100, "end": 200, "class": "SV"},
            {"chrom": "chr1", "start": 900, "end": 950, "class": "SMALL"},
        ]}
        nhf.annotate_metrics(metrics, self._block())
        assert metrics["nonhuman_fraction"]["status"] == "screened"
        assert (metrics["regions"][0]["nhf_reads"], metrics["regions"][0]["nhf"]) == (6, 0.5)
        assert (metrics["regions"][1]["nhf_reads"], metrics["regions"][1]["nhf"]) == (None, None)

    def test_update_summary_inserts_or_replaces(self, tmp_path):
        rule = "=" * 60
        old = "\n".join([rule, "  kmer-denovo  —  Discovery Mode Summary", rule, "",
                         "Per-Region Results", "  chr1:101-200 ...", rule, ""])
        summary = tmp_path / "run.summary.txt"
        summary.write_text(old)
        nhf.update_summary(str(summary), nhf.disabled_block())
        text = summary.read_text()
        assert text.count("Non-Human Fraction Screen") == 1
        assert "  chr1:101-200 ...\n\nNon-Human Fraction Screen\n" in text
        assert text.endswith("disabled (no --kraken2-db: the informative reads "
                             "were not screened)\n\n" + rule + "\n")
        nhf.update_summary(str(summary), self._block())
        text = summary.read_text()
        assert text.count("Non-Human Fraction Screen") == 1
        assert "disabled" not in text and "Virus(3)" in text
        assert text.endswith(rule + "\n")
