"""SV calling end to end on the simulated trio in ``examples/sv_demo``.

The child is heterozygous for a deletion, a tandem duplication, an
inversion, a reciprocal translocation, two novel insertions and an SNV,
with reads crossing the breakpoints as split reads (SA tags), clipped
reads and discordant pairs.  Each event must be classed, typed and linked
as expected both with the simulator's own alignments and, when bwa is
installed, with the same reads aligned by BWA-MEM.

With the simulator's alignments the output must also match the files
committed in ``examples/sv_demo/expected/``, which ``docs/sv_calling.md``
walks through.  To refresh them after an intended change::

    python examples/sv_demo/simulate_sv_trio.py --outdir sv_demo
    kmer-discovery --child sv_demo/child.bam --mother sv_demo/mother.bam \\
        --father sv_demo/father.bam --ref-fasta sv_demo/ref.fa \\
        --out-prefix sv_demo/sv_demo
    cp sv_demo/events.tsv sv_demo/sv_demo.bed sv_demo/sv_demo.sv.bedpe \\
        sv_demo/sv_demo.summary.txt examples/sv_demo/expected/
"""

import difflib
import importlib.util
import os
import shutil
import subprocess

import pytest

from kmer_denovo_filter.cli import parse_discovery_args
from kmer_denovo_filter.discovery.pipeline import run_discovery_pipeline

DEMO_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, "examples", "sv_demo",
)
EXPECTED_DIR = os.path.join(DEMO_DIR, "expected")

#: A position in each expected region: (chrom, pos, class, sv_type)
EXPECTED_REGIONS = [
    ("chr1", 9_999, "SV", "DEL"),     # deletion, left breakpoint
    ("chr1", 12_000, "SV", "DEL"),    # deletion, right breakpoint
    ("chr1", 25_000, "SV", "DUP"),    # duplication, start
    ("chr1", 26_499, "SV", "DUP"),    # duplication, end
    ("chr1", 40_000, "SV", "INV"),    # inversion, both junctions
    ("chr1", 42_000, "SV", "INV"),
    ("chr1", 55_000, "SV", "BND"),    # translocation, both junctions
    ("chr2", 5_000, "SMALL", "."),    # SNV
    ("chr2", 10_000, "SV", "INS"),    # 55 bp insertion: CIGAR or split reads
    ("chr2", 15_000, "SV", "BND"),
    ("chr2", 22_000, "SV", "INS"),    # 250 bp insertion: clipped both sides
]

#: Junctions: (position in region 1, in region 2, sv_type, strands)
EXPECTED_LINKS = [
    (("chr1", 9_999), ("chr1", 12_000), "DEL", ("+", "-")),
    (("chr1", 25_000), ("chr1", 26_499), "DUP", ("-", "+")),
    # Both junctions of the balanced inversion and translocation
    (("chr1", 40_000), ("chr1", 42_000), "INV", ("+", "+")),
    (("chr1", 40_000), ("chr1", 42_000), "INV", ("-", "-")),
    (("chr1", 55_000), ("chr2", 15_000), "BND", ("+", "-")),
    (("chr1", 55_000), ("chr2", 15_000), "BND", ("-", "+")),
]


def _load_simulator():
    spec = importlib.util.spec_from_file_location(
        "simulate_sv_trio", os.path.join(DEMO_DIR, "simulate_sv_trio.py"),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_discovery(outdir):
    run_discovery_pipeline(parse_discovery_args([
        "--child", str(outdir / "child.bam"),
        "--mother", str(outdir / "mother.bam"),
        "--father", str(outdir / "father.bam"),
        "--ref-fasta", str(outdir / "ref.fa"),
        "--out-prefix", str(outdir / "sv_demo"),
        "--threads", "2",
    ]))


def _align_with_bwa(src, dest):
    """Align the reads of the trio BAMs in *src* again, with BWA-MEM."""
    ref = str(dest / "ref.fa")
    shutil.copy(src / "ref.fa", ref)
    subprocess.run(["bwa", "index", ref], check=True, capture_output=True)
    subprocess.run(["samtools", "faidx", ref], check=True)
    for sample in ("child", "mother", "father"):
        fastq = [str(dest / f"{sample}_{i}.fq") for i in (1, 2)]
        collated = subprocess.run(
            ["samtools", "collate", "-u", "-O", str(src / f"{sample}.bam")],
            check=True, capture_output=True).stdout
        subprocess.run(
            ["samtools", "fastq", "-n", "-1", fastq[0], "-2", fastq[1],
             "-0", "/dev/null", "-s", "/dev/null", "-"],
            input=collated, check=True, capture_output=True)
        sam = subprocess.run(["bwa", "mem", "-t", "2", ref, *fastq],
                             check=True, capture_output=True).stdout
        bam = str(dest / f"{sample}.bam")
        subprocess.run(["samtools", "sort", "-o", bam, "-"], input=sam,
                       check=True, capture_output=True)
        subprocess.run(["samtools", "index", bam], check=True)


@pytest.fixture(scope="module")
def sv_demo(tmp_path_factory):
    """Simulate the demo trio, run discovery on it; return the out dir."""
    outdir = tmp_path_factory.mktemp("sv_demo")
    _load_simulator().main(["--outdir", str(outdir)])
    _run_discovery(outdir)
    return outdir


@pytest.fixture(scope="module")
def sv_demo_bwa(sv_demo, tmp_path_factory):
    """The same reads aligned with BWA-MEM, then run through discovery."""
    if shutil.which("bwa") is None:
        pytest.skip("bwa is not installed")
    outdir = tmp_path_factory.mktemp("sv_demo_bwa")
    _align_with_bwa(sv_demo, outdir)
    _run_discovery(outdir)
    return outdir


@pytest.fixture(params=["simulated alignments", "bwa-mem"])
def demo_run(request):
    """Discovery output for both ways of aligning the demo reads."""
    return request.getfixturevalue(
        "sv_demo" if request.param == "simulated alignments"
        else "sv_demo_bwa")


def _rows(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh
                if not line.startswith("#")]


def _containing(regions, chrom, pos):
    hits = [r for r in regions if r[0] == chrom and r[1] <= pos < r[2]]
    assert len(hits) == 1, f"{len(hits)} regions contain {chrom}:{pos}"
    return hits[0]


@pytest.mark.parametrize("name", [
    "events.tsv", "sv_demo.bed", "sv_demo.sv.bedpe", "sv_demo.summary.txt",
])
def test_output_matches_expected(sv_demo, name):
    with open(os.path.join(EXPECTED_DIR, name)) as fh:
        expected = fh.read().splitlines()
    with open(sv_demo / name) as fh:
        generated = fh.read().splitlines()
    if generated != expected:
        diff = "\n".join(difflib.unified_diff(
            expected, generated, f"expected/{name}", f"generated/{name}",
            lineterm="",
        ))
        pytest.fail(f"{name} differs from the committed example (see this "
                    f"module's docstring to refresh it):\n{diff}")


def test_each_event_is_classified_and_typed(demo_run):
    rows = _rows(demo_run / "sv_demo.bed")
    regions = [(r[0], int(r[1]), int(r[2]), r[9], r[12]) for r in rows]
    called = {_containing(regions, chrom, pos)[:3]: (chrom, pos, cls, sv_type)
              for chrom, pos, cls, sv_type in EXPECTED_REGIONS}
    assert len(called) == len(regions) == len(EXPECTED_REGIONS)
    for region in regions:
        chrom, pos, cls, sv_type = called[region[:3]]
        assert region[3:] == (cls, sv_type), f"{chrom}:{pos}"


def test_each_junction_is_linked_with_its_type(demo_run):
    regions = [(r[0], int(r[1]), int(r[2]))
               for r in _rows(demo_run / "sv_demo.bed")]
    links = sorted(
        ((r[0], int(r[1]), int(r[2])), (r[3], int(r[4]), int(r[5])),
         r[10], (r[8], r[9]))
        for r in _rows(demo_run / "sv_demo.sv.bedpe")
    )
    expected = sorted(
        (_containing(regions, *end1), _containing(regions, *end2),
         sv_type, strands)
        for end1, end2, sv_type, strands in EXPECTED_LINKS
    )
    assert links == expected
