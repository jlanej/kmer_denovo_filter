"""Tests for the Python helpers of ``examples/1kg_trios`` (no network).

* ``select_trios.py`` on a synthetic pedigree and sequence indexes;
* ``child_only_svs.py`` on a synthetic multi-sample SV VCF, with the
  ``--dnm-regions`` file it writes checked by the discovery pipeline's own
  loader;
* ``summarize_trios.py`` on fabricated per-trio outputs.
"""

import importlib.util
import json
import os
import textwrap

import pytest

from kmer_denovo_filter.discovery.pipeline import _load_dnm_regions

EXAMPLE_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir, "examples", "1kg_trios",
)

INDEX_HEADER = (
    "#ENA_FILE_PATH\tMD5SUM\tRUN_ID\tSTUDY_ID\tSTUDY_NAME\tCENTER_NAME\t"
    "SUBMISSION_ID\tSUBMISSION_DATE\tSAMPLE_ID\tSAMPLE_NAME\tPOPULATION\t"
    "EXPERIMENT_ID\tINSTRUMENT_PLATFORM\tINSTRUMENT_MODEL\tLIBRARY_NAME\t"
    "RUN_NAME\tINSERT_SIZE\tLIBRARY_LAYOUT\tPAIRED_FASTQ\tREAD_COUNT\t"
    "BASE_COUNT\tANALYSIS_GROUP"
)

PED = textwrap.dedent("""\
    FamilyID SampleID FatherID MotherID Sex Population Superpopulation
    HG00096 HG00096 0 0 1 GBR EUR
    SH001 HG00403 0 0 1 CHS EAS
    SH001 HG00404 0 0 2 CHS EAS
    SH001 HG00405 HG00403 HG00404 2 CHS EAS
    SH001 HG00406 HG00403 HG00404 1 CHS EAS
    Y001 NA19238 0 0 2 YRI AFR
    Y001 NA19239 0 0 1 YRI AFR
    Y001 NA19240 NA19239 NA19238 2 YRI AFR
    GB20 HG02567 0 0 1 GWD AFR
    GB20 HG02568 0 0 2 GWD AFR
    GB20 HG02569 HG02567 HG02568 1 GWD AFR
    """)

#: (sample as written in the index, population, batch); HG02569 has no CRAM
INDEX_SAMPLES = [
    ("HG00096", "GBR", "2504"), ("HG00403", "CHS", "2504"),
    (" HG00404 ", "CHS", "2504"), ("HG00405", "CHS", "698"),
    ("HG00406", "CHS", "698"), ("NA19238", "YRI", "2504"),
    ("NA19239", "YRI", "2504"), ("NA19240", "YRI", "698"),
    ("HG02567", "GWD", "2504"), ("HG02568", "GWD", "2504"),
]


def load_module(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(EXAMPLE_DIR, name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def md5_of(sample):
    import hashlib
    return hashlib.md5(sample.strip().encode()).hexdigest()


def index_row(sample, pop, batch, md5=None, study=None):
    s = sample.strip()
    run = "ERR%07d" % (hash(s) % 10_000_000)
    url = f"ftp://ftp.sra.ebi.ac.uk/vol1/run/{run[:6]}/{run}/{s}.final.cram"
    if study is None:
        study = "ERP120144" if batch == "698" else "ERP114329"
    fields = [url, md5 or md5_of(s), run, study, "30X WGS", "NYGC",
              "ERA1", "2020-03-13 00:00", "ERS1", sample, pop, "ERX1",
              "ILLUMINA", "Illumina NovaSeq 6000", s, "run", "450", "PAIRED",
              "", "1", "2", "high_coverage"]
    return "\t".join(fields)


def write_inputs(tmp_path, samples=INDEX_SAMPLES, ped=PED):
    ped_path = tmp_path / "ped.txt"
    ped_path.write_text(ped)
    idx2504 = tmp_path / "1000G_2504_high_coverage.sequence.index"
    idx698 = tmp_path / "1000G_698_related_high_coverage.sequence.index"
    for path, batch in ((idx2504, "2504"), (idx698, "698")):
        rows = [index_row(s, p, b) for s, p, b in samples if b == batch]
        path.write_text("##FileDate=20200820\n" + INDEX_HEADER + "\n"
                        + "\n".join(rows) + "\n")
    return str(ped_path), [str(idx2504), str(idx698)]


def run_select(tmp_path, *extra):
    ped, indexes = write_inputs(tmp_path)
    out_trios = tmp_path / "trios.tsv"
    out_files = tmp_path / "trio_files.tsv"
    argv = ["--ped", ped, "--index", indexes[0], "--index", indexes[1],
            "--out-trios", str(out_trios), "--out-files", str(out_files)]
    load_module("select_trios").main(argv + list(extra))
    trios = [l.split("\t") for l in out_trios.read_text().splitlines()]
    files = [l.split("\t") for l in out_files.read_text().splitlines()]
    return trios, files


# -- select_trios -------------------------------------------------------------

class TestSelectTrios:
    def test_complete_trios_one_per_family(self, tmp_path):
        trios, files = run_select(tmp_path, "--seed", "")
        assert trios[0] == ["TRIO_ID", "CHILD", "FATHER", "MOTHER", "FAMILY_ID",
                            "CHILD_SEX", "POPULATION", "SUPERPOPULATION"]
        # HG02569 has no CRAM; HG00406 is HG00405's sibling (same parents)
        assert [t[0] for t in trios[1:]] == ["HG00405", "NA19240"]
        assert trios[1][1:] == ["HG00405", "HG00403", "HG00404", "SH001", "2",
                                "CHS", "EAS"]
        assert files[0] == ["TRIO_ID", "ROLE", "SAMPLE_ID", "CRAM_FTP_URL",
                            "CRAI_FTP_URL", "CRAM_MD5", "RELEASE_BATCH"]
        assert len(files) == 1 + 3 * 2
        roles = [(f[0], f[1], f[2]) for f in files[1:4]]
        assert roles == [("HG00405", "child", "HG00405"),
                         ("HG00405", "father", "HG00403"),
                         ("HG00405", "mother", "HG00404")]
        mother = files[3]
        assert mother[3].endswith("/HG00404.final.cram")   # trimmed name
        assert mother[4] == mother[3] + ".crai"
        assert mother[5] == md5_of("HG00404")
        assert mother[6] == "2504"
        assert files[1][6] == "698"                        # the child

    def test_all_siblings(self, tmp_path):
        trios, files = run_select(tmp_path, "--all-siblings", "--seed", "")
        assert [t[0] for t in trios[1:]] == ["HG00405", "HG00406", "NA19240"]
        assert len(files) == 1 + 9

    def test_n_and_seed_are_deterministic(self, tmp_path):
        first, _ = run_select(tmp_path, "--all-siblings", "--n", "2",
                              "--seed", "7")
        second, _ = run_select(tmp_path, "--all-siblings", "--n", "2",
                               "--seed", "7")
        assert first == second
        assert len(first) == 3
        sorted_order, _ = run_select(tmp_path, "--all-siblings", "--n", "2",
                                     "--seed", "")
        assert [t[0] for t in sorted_order[1:]] == ["HG00405", "HG00406"]

    def test_population_filters(self, tmp_path):
        trios, _ = run_select(tmp_path, "--populations", "YRI")
        assert [t[0] for t in trios[1:]] == ["NA19240"]
        trios, _ = run_select(tmp_path, "--superpopulations", "EAS")
        assert [t[0] for t in trios[1:]] == ["HG00405"]
        # excluding a child lets its sibling stand in for the family
        trios, _ = run_select(tmp_path, "--exclude", "HG00405")
        assert [t[0] for t in trios[1:]] == ["HG00406", "NA19240"]
        trios, _ = run_select(tmp_path, "--exclude", "HG00405,HG00406")
        assert [t[0] for t in trios[1:]] == ["NA19240"]

    def test_explicit_children(self, tmp_path):
        children = tmp_path / "children.txt"
        children.write_text("HG00406  # a sibling\nNA19240\n")
        trios, _ = run_select(tmp_path, "--children", str(children))
        assert [t[0] for t in trios[1:]] == ["HG00406", "NA19240"]
        with pytest.raises(SystemExit, match="HG02569"):
            run_select(tmp_path, "--children", "HG02569")

    def test_bad_index_rows_are_rejected(self, tmp_path):
        mod = load_module("select_trios")
        bad = tmp_path / "bad.index"
        bad.write_text(INDEX_HEADER + "\n"
                       + index_row("HG00405", "CHS", "698", md5="nothex")
                       + "\n")
        with pytest.raises(ValueError, match="bad MD5"):
            mod.read_sequence_index(str(bad))
        bad.write_text(INDEX_HEADER + "\n"
                       + index_row("HG00405", "CHS", "698").replace(
                           "ftp://ftp.sra.ebi.ac.uk", "ftp://elsewhere")
                       + "\n")
        with pytest.raises(ValueError, match="ENA"):
            mod.read_sequence_index(str(bad))

    def test_batch_from_study_then_file_name(self, tmp_path):
        mod = load_module("select_trios")
        idx = tmp_path / "1000G_698_related_high_coverage.sequence.index"
        idx.write_text(INDEX_HEADER + "\n"
                       + index_row("HG00405", "CHS", "698") + "\n"
                       + index_row("NA12718", "CEU", "2504") + "\n"
                       + index_row("HG00096", "GBR", "698", study="ERP999999") + "\n")
        samples = mod.read_sequence_index(str(idx))
        assert samples["HG00405"][2] == "698"
        assert samples["NA12718"][2] == "2504"        # by its study accession
        assert samples["HG00096"][2] == "698"         # unknown study: file name
        assert mod.read_sequence_index(str(idx), batch="x")["NA12718"][2] == "x"

    def test_no_match_exits(self, tmp_path):
        with pytest.raises(SystemExit, match="no trio"):
            run_select(tmp_path, "--populations", "XXX")


# -- child_only_svs -----------------------------------------------------------

SV_VCF_HEADER = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##contig=<ID=chr1,length=248956422>
    ##contig=<ID=chr2,length=242193529>
    ##contig=<ID=chr3,length=198295559>
    ##ALT=<ID=DEL,Description="Deletion">
    ##ALT=<ID=DUP,Description="Duplication">
    ##ALT=<ID=INS:ME:ALU,Description="Alu element insertion">
    ##ALT=<ID=INV,Description="Inversion">
    ##ALT=<ID=CNV,Description="Copy Number Polymorphism">
    ##FILTER=<ID=LowQual,Description="Variant failed lightGBM filter">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV Type">
    ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
    ##INFO=<ID=END,Number=1,Type=Integer,Description="End">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">
    ##INFO=<ID=EVIDENCE,Number=.,Type=String,Description="Evidence">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    """)


def sv_record(chrom, pos, svid, alt, filt, info, gts):
    return "\t".join([chrom, str(pos), svid, "N", alt, "100", filt, info,
                      "GT:GQ"] + gts) + "\n"


def write_sv_vcf(path, samples=("C1", "F1", "M1", "X1")):
    """A small call set: C1 is the child, F1/M1 the parents, X1 unrelated."""
    def gts(child="0/1:60", father="0/0:60", mother="0/0:60", other="0/0:60"):
        return [child, father, mother, other]

    body = SV_VCF_HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" \
        + "\t".join(samples) + "\n"
    rows = [
        # kept: a child-only PASS deletion, rare
        sv_record("chr1", 1000, "DEL1", "<DEL>", "PASS",
                  "SVTYPE=DEL;SVLEN=-500;END=1500;AF=0.001;ALGORITHMS=manta;EVIDENCE=PE,SR",
                  gts()),
        # kept: a child-only Alu insertion; size unknown, labelled from ALT
        sv_record("chr1", 5000, "INS1", "<INS:ME:ALU>", "PASS",
                  "SVTYPE=INS;SVLEN=300;END=5001;AF=0.002", gts()),
        # filtered: common in the population (a parent genotyping miss)
        sv_record("chr1", 9000, "DEL2", "<DEL>", "PASS",
                  "SVTYPE=DEL;SVLEN=-200;END=9200;AF=0.3", gts()),
        # not child-only: inherited from the father
        sv_record("chr1", 12000, "DUP1", "<DUP>", "PASS",
                  "SVTYPE=DUP;SVLEN=3000;END=15000;AF=0.001",
                  gts(father="0/1:60")),
        # not child-only: the father's genotype is missing
        sv_record("chr1", 15000, "DEL3", "<DEL>", "PASS",
                  "SVTYPE=DEL;SVLEN=-100;END=15100;AF=0.001",
                  gts(father="./.:1")),
        # filtered: failed the call set's filters
        sv_record("chr1", 18000, "DEL4", "<DEL>", "LowQual",
                  "SVTYPE=DEL;SVLEN=-100;END=18100;AF=0.001", gts()),
        # filtered: child GQ below the threshold
        sv_record("chr1", 21000, "INV1", "<INV>", "PASS",
                  "SVTYPE=INV;SVLEN=2000;END=23000;AF=0.001",
                  gts(child="0/1:5")),
        # skipped: multiallelic CNV
        sv_record("chr1", 25000, "CNV1", "<CNV>", "PASS",
                  "SVTYPE=CNV;SVLEN=5000;END=30000;AF=0.001", gts()),
        # kept: a breakend; no size
        sv_record("chr2", 100, "BND1", "N[chr3:500[", "PASS",
                  "SVTYPE=BND;AF=0.0", gts()),
        # kept: size from END when SVLEN is absent
        sv_record("chr2", 3000, "DEL5", "<DEL>", "PASS",
                  "SVTYPE=DEL;END=3400;AF=0.001", gts()),
    ]
    with open(path, "w") as fh:
        fh.write(body + "".join(rows))
    return path


def run_child_only(tmp_path, *extra):
    vcf = write_sv_vcf(str(tmp_path / "sv.vcf"))
    regions = tmp_path / "regions.tsv"
    table = tmp_path / "detail.tsv"
    load_module("child_only_svs").main([
        "--vcf", vcf, "--child", "C1", "--father", "F1", "--mother", "M1",
        "--out-regions", str(regions), "--out-table", str(table),
    ] + list(extra))
    return regions, table


class TestChildOnlySvs:
    def test_regions_file_loads_in_the_pipeline(self, tmp_path):
        regions, table = run_child_only(tmp_path)
        events = _load_dnm_regions(str(regions))
        assert events == [
            ("chr1", 1000, 500, "DEL"),
            ("chr1", 5000, None, "INS:ME:ALU"),
            ("chr2", 100, None, "BND"),
            ("chr2", 3000, 400, "DEL"),
        ]
        rows = [l.split("\t") for l in table.read_text().splitlines()]
        assert rows[0][:3] == ["ID", "CHROM", "POS"]
        by_id = {r[0]: r for r in rows[1:]}
        # every child-only call is in the table, filtered or not
        assert set(by_id) == {"DEL1", "INS1", "DEL2", "DEL4", "INV1", "BND1",
                              "DEL5"}
        assert by_id["DEL2"][-1] == "FILTERED"
        assert by_id["DEL1"][-1] == "PASS_FILTERS"
        assert by_id["DEL1"][9:11] == ["0/1", "60"]      # child GT, GQ
        assert by_id["DEL1"][15] == "manta"              # ALGORITHMS

    def test_filters_can_be_relaxed(self, tmp_path):
        regions, _ = run_child_only(tmp_path, "--no-pass-only", "--max-af",
                                    "-1", "--min-gq", "0")
        events = _load_dnm_regions(str(regions))
        assert [e[0] + ":" + str(e[1]) for e in events] == [
            "chr1:1000", "chr1:5000", "chr1:9000", "chr1:18000", "chr1:21000",
            "chr2:100", "chr2:3000"]

    def test_sv_types_and_cnv(self, tmp_path):
        regions, _ = run_child_only(tmp_path, "--sv-types", "DEL",
                                    "--include-cnv")
        events = _load_dnm_regions(str(regions))
        assert [e[3] for e in events] == ["DEL", "DEL"]
        regions, _ = run_child_only(tmp_path, "--sv-types", "CNV,DEL",
                                    "--include-cnv")
        assert [e[3] for e in _load_dnm_regions(str(regions))] == [
            "DEL", "CNV", "DEL"]

    def test_unknown_sample(self, tmp_path):
        vcf = write_sv_vcf(str(tmp_path / "sv.vcf"))
        with pytest.raises(SystemExit, match="not in the VCF"):
            load_module("child_only_svs").main([
                "--vcf", vcf, "--child", "NOPE", "--father", "F1",
                "--mother", "M1", "--out-regions", str(tmp_path / "r.tsv")])

    def test_no_events_leaves_an_empty_regions_file(self, tmp_path):
        regions, table = run_child_only(tmp_path, "--sv-types", "INV")
        data = [l for l in regions.read_text().splitlines()
                if not l.startswith("#")]
        assert data == []
        # which the pipeline's loader rejects, so the job must not pass it
        with pytest.raises(ValueError, match="no de novo events"):
            _load_dnm_regions(str(regions))


# -- summarize_trios ----------------------------------------------------------

REGION_KEYS = ("chrom", "start", "end", "size", "reads", "unique_kmers",
               "split_reads", "discordant_pairs", "max_clip_len",
               "unmapped_mates", "breakpoint_reads", "large_indel_reads",
               "sv_type", "class")


def region(chrom, start, end, cls, sv_type, reads=10):
    values = (chrom, start, end, end - start, reads, 30, 4, 2, 60, 0, 5, 0,
              sv_type, cls)
    return dict(zip(REGION_KEYS, values))


def write_trio_results(results_dir, trio, regions, links=(), dnm=None,
                       detail_rows=(), wall_s=7200, done=True):
    d = results_dir / trio
    d.mkdir(parents=True)
    p = d / trio
    metrics = {"mode": "discovery", "child_candidate_kmers": 1000,
               "non_ref_kmers": 200, "proband_unique_kmers": 50,
               "informative_reads": 40, "unmapped_informative_reads": 1,
               "candidate_regions": len(regions), "regions": regions}
    if dnm:
        metrics["dnm_evaluation"] = dnm
    (d / f"{trio}.metrics.json").write_text(json.dumps(metrics))
    (d / f"{trio}.bed").write_text("#chrom\tstart\tend\n" + "".join(
        f"{r['chrom']}\t{r['start']}\t{r['end']}\n" for r in regions))
    (d / f"{trio}.sv.bedpe").write_text(
        "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tsupporting_reads"
        "\tstrand1\tstrand2\tsv_type\n"
        + "".join("\t".join(str(v) for v in l) + "\n" for l in links))
    (d / f"{trio}.summary.txt").write_text("summary\n")
    (d / f"{trio}.timing.tsv").write_text(
        "trio\tchild\twall_s\tstarted\tfinished\tnode\tthreads\tmemory_gb\t"
        "job_id\timage\tregions\n"
        f"{trio}\t{trio}\t{wall_s}\t1\t2\tnode1\t16\t128\t42\timg\t{len(regions)}\n")
    (d / f"{trio}.log").write_text(
        "[Module 1] Complete (1h 0m 0s)\n[Module 2] Complete (20m 0.0s)\n"
        "[Module 2b] Complete (30.0s, index: 1 MB)\n"
        "[Module 3] Complete (39m 30s)\nPipeline finished in 2h 0m 0s\n")
    if detail_rows:
        (d / f"{trio}.child_only_svs.detail.tsv").write_text(
            "ID\tCHROM\tSTATUS\n" + "".join(
                f"id{i}\tchr1\t{s}\n" for i, s in enumerate(detail_rows)))
    if done:
        (d / f"{trio}.done").write_text("2026-10-14\n")


class TestSummarizeTrios:
    def test_cohort_summary(self, tmp_path):
        results = tmp_path / "results"
        trios = tmp_path / "trios.tsv"
        trios.write_text(
            "TRIO_ID\tCHILD\tFATHER\tMOTHER\tFAMILY_ID\tCHILD_SEX\tPOPULATION\tSUPERPOPULATION\n"
            "A\tA\tFA\tMA\tfam1\t1\tCHS\tEAS\n"
            "B\tB\tFB\tMB\tfam2\t2\tYRI\tAFR\n"
            "C\tC\tFC\tMC\tfam3\t1\tGBR\tEUR\n")
        write_trio_results(
            results, "A",
            [region("chr1", 1000, 1300, "SV", "DEL"),
             region("chr1", 5000, 5200, "SV", "DEL"),
             region("chr2", 100, 400, "SMALL", "."),
             region("chr3", 777, 999, "SV", "INS")],
            links=[("chr1", 1000, 1300, "chr1", 5000, 5200, "SV_1", 6, "+", "-",
                    "DEL")],
            dnm={"source": "x", "slack_bp": 500, "total_loci": 2,
                 "detected": 1, "detection_rate": 0.5, "loci": []},
            detail_rows=("PASS_FILTERS", "PASS_FILTERS", "FILTERED"))
        write_trio_results(
            results, "B",
            [region("chr3", 800, 1100, "SV", "INS"),
             region("chr4", 1, 200, "AMBIGUOUS", ".")],
            wall_s=3600)
        write_trio_results(results, "C", [region("chr5", 1, 2, "SV", "DEL")],
                           done=False)
        out = tmp_path / "cohort"
        load_module("summarize_trios").main([
            "--trios", str(trios), "--results-dir", str(results),
            "--out-dir", str(out), "--recurrent-min-trios", "2",
            "--slack", "500"])

        rows = [l.split("\t") for l in
                (out / "trio_summary.tsv").read_text().splitlines()]
        header = rows[0]
        by_trio = {r[0]: dict(zip(header, r)) for r in rows[1:]}
        assert by_trio["A"]["status"] == "done"
        assert by_trio["A"]["regions"] == "4"
        assert by_trio["A"]["SV"] == "3"
        assert by_trio["A"]["SMALL"] == "1"
        assert by_trio["A"]["SV_DEL"] == "2"
        assert by_trio["A"]["SV_INS"] == "1"
        assert by_trio["A"]["SV_nonrecurrent"] == "2"   # chr3 is shared with B
        assert by_trio["A"]["junctions"] == "1"
        assert by_trio["A"]["child_only_svs"] == "3"
        assert by_trio["A"]["child_only_svs_filtered"] == "2"
        assert by_trio["A"]["known_svs_tested"] == "2"
        assert by_trio["A"]["known_svs_detected"] == "1"
        assert by_trio["A"]["wall_hours"] == "2.00"
        assert by_trio["A"]["module1_hours"] == "1.00"
        assert by_trio["A"]["module3_hours"] == "0.66"
        assert by_trio["B"]["SV_nonrecurrent"] == "0"
        assert by_trio["B"]["AMBIGUOUS"] == "1"
        assert by_trio["C"]["status"] == "missing"

        recurrent = [l.split("\t") for l in
                     (out / "recurrent_regions.bed").read_text().splitlines()
                     if not l.startswith("#")]
        assert recurrent == [["chr3", "777", "1100", "2", "2", "SV:2", "INS:2",
                              "A,B"]]

        cands = [l.split("\t") for l in
                 (out / "sv_candidates.tsv").read_text().splitlines()]
        cand_header = cands[0]
        cands = [dict(zip(cand_header, c)) for c in cands[1:]]
        assert [(c["chrom"], c["start"], c["trio"]) for c in cands] == [
            ("chr1", "1000", "A"), ("chr1", "5000", "A")]
        assert cands[0]["junctions"] == "chr1:5000-5200|+-|DEL|6"
        assert cands[1]["junctions"] == "chr1:1000-1300|+-|DEL|6"
        assert cands[0]["n_trios_at_locus"] == "1"
        bedpe = (out / "sv_candidates.bedpe").read_text().splitlines()
        assert len(bedpe) == 2 and bedpe[1].endswith("\tDEL\tA")
        bed = (out / "sv_candidates.bed").read_text().splitlines()
        assert bed[1] == "chr1\t1000\t1300\tA|DEL|SV\t10\t."

        md = (out / "cohort_summary.md").read_text()
        assert "Trios selected: 3; finished: 2; missing: 1" in md
        assert "De novo SV candidates after dropping recurrent loci: 2" in md
        assert "supported by discovery: 1 (50%)" in md
        assert "| C | GBR | (missing) |" in md

    def test_no_finished_trio(self, tmp_path):
        trios = tmp_path / "trios.tsv"
        trios.write_text("TRIO_ID\tCHILD\tFATHER\tMOTHER\tFAMILY_ID\tCHILD_SEX\t"
                         "POPULATION\tSUPERPOPULATION\nA\tA\tF\tM\tf\t1\tCHS\tEAS\n")
        out = tmp_path / "cohort"
        load_module("summarize_trios").main([
            "--trios", str(trios), "--results-dir", str(tmp_path / "none"),
            "--out-dir", str(out)])
        assert "No finished trio yet." in (out / "cohort_summary.md").read_text()
        assert (out / "trio_summary.tsv").read_text() == "trio\tstatus\nA\tmissing\n"

    def test_parse_elapsed(self):
        mod = load_module("summarize_trios")
        assert mod.parse_elapsed("1h 2m 3s") == 3723
        assert mod.parse_elapsed("5m 1.5s") == 301.5
        assert mod.parse_elapsed("12.5s") == 12.5
