"""Module 3's scan: batches, windows and the shared k-mer set; and
Module 1's streamed reference subtraction.

The scan of a genome in windows, by several workers, in batches of reads,
must give exactly what one worker scanning each contig whole, read by
read, gives: the same informative reads (de-duplicated the same way when
both mates of a pair are informative), the same hit positions, the same
coverage and the same informative BAM.
"""

import os
import random
import subprocess

import pysam
import pytest

import kmer_denovo_filter.core.bam_scanner as bam_scanner_mod
from kmer_denovo_filter.core.kmer_set import KmerSet
from kmer_denovo_filter.discovery.pipeline import (
    _anchor_and_cluster,
    _scan_windows,
    _stream_records_to_pipes,
    _subtract_reference_kmers,
    _window_label,
)
from kmer_denovo_filter.kmer_utils import canonicalize

from tests.helpers import _create_bam_with_supplementary

K = 11


def _genome(rng, chroms):
    return {chrom: "".join(rng.choice("ACGT") for _ in range(length))
            for chrom, length in chroms.items()}


def _trio_like_bam(tmp_path, rng, genome, n_reads=600, read_len=60,
                   novel_every=10):
    """Paired reads over *genome*; every *novel_every*-th pair carries a
    substitution in both mates, so its mates are informative; a few are
    placed-unmapped or unplaced.  Returns (bam, proband-unique k-mers)."""
    chroms = list(genome)
    reads = []
    novel = set()
    for i in range(n_reads // 2):
        chrom = rng.choice(chroms)
        chrom_idx = chroms.index(chrom)
        length = len(genome[chrom])
        pos1 = rng.randrange(0, length - 2 * read_len - 100)
        pos2 = pos1 + read_len + rng.randrange(0, 100)
        seqs = [genome[chrom][pos1:pos1 + read_len],
                genome[chrom][pos2:pos2 + read_len]]
        if i % novel_every == 0:
            for m in range(2):
                j = read_len // 2
                base = rng.choice([b for b in "ACGT" if b != seqs[m][j]])
                seqs[m] = seqs[m][:j] + base + seqs[m][j + 1:]
                for w in range(j - K + 1, j + 1):
                    novel.add(canonicalize(seqs[m][w:w + K]))
        if i % 97 == 3:
            # the second mate unmapped, placed at its mate's position
            reads.append({"name": f"p{i}", "chrom_idx": chrom_idx,
                          "pos": pos1, "seq": seqs[0],
                          "flag": 0x1 | 0x8 | 0x40})
            reads.append({"name": f"p{i}", "chrom_idx": chrom_idx,
                          "pos": pos1, "seq": seqs[1],
                          "flag": 0x1 | 0x4 | 0x80, "cigar": [], "mapq": 0})
            continue
        reads.append({"name": f"p{i}", "chrom_idx": chrom_idx, "pos": pos1,
                      "seq": seqs[0], "flag": 0x1 | 0x40})
        reads.append({"name": f"p{i}", "chrom_idx": chrom_idx, "pos": pos2,
                      "seq": seqs[1], "flag": 0x1 | 0x80 | 0x10})
    for i in range(3):
        seq = "".join(rng.choice("ACGT") for _ in range(read_len))
        novel.update(canonicalize(seq[w:w + K]) for w in range(20))
        reads.append({"name": f"u{i}", "chrom_idx": -1, "pos": -1,
                      "seq": seq, "flag": 0x4, "cigar": [], "mapq": 0})
    bam = str(tmp_path / "child.bam")
    _create_bam_with_supplementary(
        bam, None, chroms, [len(genome[c]) for c in chroms], reads)
    return bam, novel


def _comparable(result):
    """The parts of _anchor_and_cluster's result that must not change."""
    (regions, region_reads, total, region_kmers, unmapped, sv_meta,
     kmer_cov, read_cov) = result
    return (regions, region_reads, total,
            {r: len(ks) for r, ks in region_kmers.items()}, unmapped,
            sv_meta, {c: dict(v) for c, v in kmer_cov.items()},
            {c: dict(v) for c, v in read_cov.items()})


def _bam_reads(path):
    with pysam.AlignmentFile(path) as fh:
        return [(r.query_name, r.flag, r.reference_name, r.reference_start,
                 r.get_tag("dk")) for r in fh.fetch(until_eof=True)]


class TestWindows:
    def test_windows_cover_each_contig_once(self):
        windows = _scan_windows(["a", "b", "c"], [25, 10, 0], 10)
        assert windows == [
            ("a", 0, 9), ("a", 9, 18), ("a", 18, 25),
            ("b", None, None), ("c", None, None), (None, None, None),
        ]
        assert _scan_windows(["a"], [25], None) == [
            ("a", None, None), (None, None, None)]
        assert _window_label("a", 9, 18) == "a:10-18"
        assert _window_label("a", None, None) == "a"
        assert _window_label(None, None, None) == "(unmapped)"

    @pytest.mark.parametrize("threads", [1, 3])
    def test_windowed_scan_equals_whole_contig_scan(self, tmp_path, threads,
                                                    monkeypatch):
        """Windows of 1 kb (reads and pairs straddle their boundaries) and
        batches of 7 reads give exactly the whole-contig, one-batch scan."""
        rng = random.Random(21)
        genome = _genome(rng, {"chrA": 6000, "chrB": 2500})
        bam, kmers = _trio_like_bam(tmp_path, rng, genome)
        whole = _anchor_and_cluster(
            bam, None, kmers, K, threads=1, min_distinct_kmers_per_read=3,
            tmpdir=str(tmp_path), informative_bam=str(tmp_path / "whole.bam"),
            window_size=None,
        )
        monkeypatch.setattr(bam_scanner_mod, "_SCAN_BATCH_READS", 7)
        windowed = _anchor_and_cluster(
            bam, None, kmers, K, threads=threads,
            min_distinct_kmers_per_read=3, tmpdir=str(tmp_path),
            informative_bam=str(tmp_path / "windowed.bam"), window_size=1000,
        )
        assert whole[2] > 20 and whole[4] >= 1, "the fixture is informative"
        assert _comparable(windowed) == _comparable(whole)
        assert _bam_reads(tmp_path / "windowed.bam") == _bam_reads(
            tmp_path / "whole.bam")

    def test_a_read_belongs_to_the_window_its_alignment_starts_in(
            self, tmp_path):
        """A read overlapping a window boundary is scanned once, by the
        window of its start; a pair straddling it keeps the first mate."""
        rng = random.Random(4)
        genome = _genome(rng, {"chrA": 400})
        seq_a = genome["chrA"][190:250]   # starts before the boundary at 200
        seq_b = genome["chrA"][205:265]
        kmers = {canonicalize(s[w:w + K]) for s in (seq_a, seq_b)
                 for w in range(20)}
        bam = str(tmp_path / "child.bam")
        _create_bam_with_supplementary(bam, None, ["chrA"], [400], [
            {"name": "pair", "chrom_idx": 0, "pos": 190, "seq": seq_a,
             "flag": 0x1 | 0x40},
            {"name": "pair", "chrom_idx": 0, "pos": 205, "seq": seq_b,
             "flag": 0x1 | 0x80},
        ])
        bam_scanner_mod._init_scan_worker(kmers, K, 1)
        try:
            first = bam_scanner_mod._scan_contig_for_hits(
                bam, None, "chrA", start=0, end=200)
            second = bam_scanner_mod._scan_contig_for_hits(
                bam, None, "chrA", start=200, end=400)
            whole = bam_scanner_mod._scan_contig_for_hits(bam, None, "chrA")
        finally:
            bam_scanner_mod._reset_scan_worker()
        assert first[3] == 1 and second[3] == 1 and whole[3] == 2
        assert [h[1] for h in first[0]] == [190]
        assert [h[1] for h in second[0]] == [205]
        # both windows report the pair's key; the whole scan keeps the first
        assert [h[1] for h in whole[0]] == [190]
        result = _anchor_and_cluster(bam, None, kmers, K, threads=1,
                                     window_size=200)
        assert result[2] == 1 and result[1] == {("chrA", 190, 250): {"pair"}}


class TestBatches:
    def test_batches_end_at_the_read_or_base_limit(self, monkeypatch):
        """Reads are scanned in batches; the batch size changes nothing."""

        class _Read:
            def __init__(self, name, seq):
                self.query_name = name
                self.query_sequence = seq
                self.is_secondary = self.is_duplicate = False
                self.is_supplementary = False
                self.is_unmapped = True
                self.reference_id = -1
                self.reference_start = -1

        class _Bam:
            def __init__(self, *args, **kwargs):
                self.header = None
                self._reads = [_Read(f"r{i}", "ACGTT" * 4 + "AAAAA" * (i % 2))
                               for i in range(23)]

            def fetch(self, *args, **kwargs):
                return iter(self._reads)

            def close(self):
                pass

        batches = []
        real_scan = KmerSet.scan

        def counting_scan(self, seqs):
            batches.append(len(seqs))
            return real_scan(self, seqs)

        monkeypatch.setattr(bam_scanner_mod.pysam, "AlignmentFile", _Bam)
        monkeypatch.setattr(KmerSet, "scan", counting_scan)
        bam_scanner_mod._init_scan_worker(["AAAAA"], 5, 1)
        try:
            monkeypatch.setattr(bam_scanner_mod, "_SCAN_BATCH_READS", 5)
            by_reads = bam_scanner_mod._scan_contig_for_hits("x", None, None)
            assert batches == [5, 5, 5, 5, 3]
            batches.clear()
            monkeypatch.setattr(bam_scanner_mod, "_SCAN_BATCH_READS", 1000)
            monkeypatch.setattr(bam_scanner_mod, "_SCAN_BATCH_BASES", 100)
            by_bases = bam_scanner_mod._scan_contig_for_hits("x", None, None)
            assert batches == [5, 5, 5, 5, 3]
            batches.clear()
            monkeypatch.setattr(bam_scanner_mod, "_SCAN_BATCH_BASES", 10**9)
            at_once = bam_scanner_mod._scan_contig_for_hits("x", None, None)
            assert batches == [23]
        finally:
            bam_scanner_mod._reset_scan_worker()
        # 11 reads carry AAAAA (the odd ones), all unmapped
        assert by_reads[2] == by_bases[2] == at_once[2] == 11
        assert by_reads[1] == by_bases[1] == at_once[1] == {
            (f"r{i}", False) for i in range(1, 23, 2)}
        assert by_reads[3] == 23


class TestStreamedSubtraction:
    """``jellyfish dump`` streamed through parallel ``jellyfish query``."""

    def _records(self, n):
        rng = random.Random(n)
        return [(f"{i}", "".join(rng.choice("ACGT") for _ in range(31)))
                for i in range(n)]

    @pytest.mark.parametrize("n_sinks,block", [(1, 64), (3, 64), (4, 1 << 20)])
    def test_records_arrive_whole_and_all_once(self, tmp_path, n_sinks, block):
        records = self._records(2000)
        fasta = "".join(f">{h}\n{s}\n" for h, s in records).encode()
        src_r, src_w = os.pipe()
        sinks = []
        outs = []
        for j in range(n_sinks):
            out = open(tmp_path / f"part{j}", "wb")
            proc = subprocess.Popen(["cat"], stdin=subprocess.PIPE,
                                    stdout=out, bufsize=0)
            sinks.append(proc)
            outs.append(out)
        writer = subprocess.Popen(["cat"], stdin=subprocess.PIPE,
                                  stdout=src_w, bufsize=0)
        os.close(src_w)
        writer.stdin.write(fasta)
        writer.stdin.close()
        n_bytes, n_records = _stream_records_to_pipes(
            src_r, [s.stdin.fileno() for s in sinks], block_bytes=block)
        for s in sinks:
            s.stdin.close()
            s.wait()
        for out in outs:
            out.close()
        writer.wait()
        os.close(src_r)
        assert (n_bytes, n_records) == (len(fasta), 2000)
        got = []
        for j in range(n_sinks):
            text = (tmp_path / f"part{j}").read_bytes().decode()
            lines = text.splitlines()
            assert len(lines) % 2 == 0
            for h, s in zip(lines[::2], lines[1::2]):
                assert h.startswith(">") and len(s) == 31
                got.append((h[1:], s))
        assert sorted(got) == sorted(records)
        if n_sinks > 1 and block == 64:
            assert all((tmp_path / f"part{j}").stat().st_size
                       for j in range(n_sinks)), "every sink got records"

    def test_end_to_end_against_the_one_query_way(self, tmp_path):
        """The parallel streamed subtraction keeps exactly the k-mers a
        single ``jellyfish query`` over the dumped FASTA keeps."""
        rng = random.Random(8)
        ref = "".join(rng.choice("ACGT") for _ in range(20000))
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(f">chr1\n{ref}\n")
        ref_jf = str(tmp_path / "ref.jf")
        subprocess.run(["jellyfish", "count", "-m", "11", "-s", "1M", "-t", "1",
                        "-C", "-o", ref_jf, str(ref_fa)], check=True)
        # the child: the reference three times, plus novel sequence twice
        novel = "".join(rng.choice("ACGT") for _ in range(3000))
        child_fa = tmp_path / "child.fa"
        child_fa.write_text(f">a\n{ref}\n>b\n{ref}\n>c\n{ref}\n>n\n{novel}\n"
                            f">m\n{novel}\n>once\n{novel[:500]}ACGT{novel[500:]}\n")
        child_jf = str(tmp_path / "child.jf")
        subprocess.run(["jellyfish", "count", "-m", "11", "-s", "1M", "-t", "1",
                        "-C", "-o", child_jf, str(child_fa)], check=True)

        dumped = subprocess.run(["jellyfish", "dump", "-c", "-L", "2", child_jf],
                                capture_output=True, text=True, check=True).stdout
        candidates = [line.split()[0] for line in dumped.splitlines() if line]
        cand_fa = tmp_path / "cand.fa"
        cand_fa.write_text("".join(f">{i}\n{k}\n" for i, k in enumerate(candidates)))
        queried = subprocess.run(["jellyfish", "query", ref_jf, "-s", str(cand_fa)],
                                 capture_output=True, text=True, check=True).stdout
        expected = sorted(line.split()[0] for line in queried.splitlines()
                          if line.split()[1] == "0")
        assert 2500 < len(expected) < len(candidates)

        out_dir = tmp_path / "run"
        out_dir.mkdir()
        non_ref_fa, n_candidates, n_non_ref = _subtract_reference_kmers(
            ref_jf, child_jf, 2, str(out_dir), threads=3, min_child_count=2,
        )
        assert n_candidates == len(candidates)
        assert n_non_ref == len(expected)
        with open(non_ref_fa) as fh:
            lines = fh.read().splitlines()
        headers, kmers = lines[::2], lines[1::2]
        assert sorted(kmers) == expected
        assert len(set(headers)) == len(headers) and all(
            h.startswith(">") for h in headers)
        assert not os.path.exists(child_jf), "the child index is removed"
        assert not os.path.exists(out_dir / "non_ref_parts")

    def test_a_failing_query_fails_the_run(self, tmp_path, monkeypatch):
        rng = random.Random(9)
        ref = "".join(rng.choice("ACGT") for _ in range(2000))
        (tmp_path / "ref.fa").write_text(f">chr1\n{ref}\n")
        jf = str(tmp_path / "x.jf")
        subprocess.run(["jellyfish", "count", "-m", "11", "-s", "1M", "-t", "1",
                        "-C", "-o", jf, str(tmp_path / "ref.fa")], check=True)
        with pytest.raises(RuntimeError, match="jellyfish query .*failed"):
            _subtract_reference_kmers(str(tmp_path / "missing.jf"), jf, 1,
                                      str(tmp_path), threads=2)

    def test_an_empty_dump(self, tmp_path):
        rng = random.Random(10)
        ref = "".join(rng.choice("ACGT") for _ in range(2000))
        (tmp_path / "ref.fa").write_text(f">chr1\n{ref}\n")
        jf = str(tmp_path / "x.jf")
        subprocess.run(["jellyfish", "count", "-m", "11", "-s", "1M", "-t", "1",
                        "-C", "-o", jf, str(tmp_path / "ref.fa")], check=True)
        ref_jf = str(tmp_path / "ref.jf")
        subprocess.run(["jellyfish", "count", "-m", "11", "-s", "1M", "-t", "1",
                        "-C", "-o", ref_jf, str(tmp_path / "ref.fa")], check=True)
        # nothing has count >= 50
        non_ref_fa, n_candidates, n_non_ref = _subtract_reference_kmers(
            ref_jf, jf, 50, str(tmp_path), threads=2)
        assert (n_candidates, n_non_ref) == (0, 0)
        assert os.path.getsize(non_ref_fa) == 0
