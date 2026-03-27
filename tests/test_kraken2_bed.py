"""Tests for per-read Kraken2 classification detail BED output."""

import gzip
import os
import tempfile

import pysam
import pytest

from kmer_denovo_filter.kmer_utils import Kraken2Runner
from kmer_denovo_filter.pipeline import (
    _extract_softclips,
    _parse_kmer_votes,
    _write_kraken2_read_detail_bed,
    _write_kraken2_span_bed,
    _write_kraken2_expanded_span_bed,
)


class TestParseKmerVotes:
    """Tests for _parse_kmer_votes helper."""

    def test_empty_string(self):
        votes, named, total, human = _parse_kmer_votes("")
        assert votes == ""
        assert named == ""
        assert total == 0
        assert human == 0

    def test_simple_votes(self):
        votes, named, total, human = _parse_kmer_votes("562:10 0:5")
        assert "562:10" in votes
        assert "0:5" in votes
        assert total == 15
        assert human == 0

    def test_human_kmer_count(self):
        votes, named, total, human = _parse_kmer_votes("562:10 9606:3 0:5")
        assert human == 3
        assert total == 18

    def test_with_name_map(self):
        name_map = {562: "Escherichia_coli", 9606: "Homo_sapiens"}
        votes, named, total, human = _parse_kmer_votes(
            "562:10 9606:3 0:5", name_map=name_map,
        )
        assert "Escherichia_coli:10" in named
        assert "Homo_sapiens:3" in named
        assert "unclassified:5" in named

    def test_paired_end_delimiter(self):
        """Paired-end reads use |:| delimiter between mates."""
        votes, named, total, human = _parse_kmer_votes(
            "562:5 0:2 |:| 562:3 0:1",
        )
        assert total == 11
        # 562:5 + 562:3 = 562:8
        assert "562:8" in votes

    def test_ambiguous_tokens_excluded(self):
        """'A' tokens (ambiguous) should be excluded."""
        votes, named, total, human = _parse_kmer_votes("562:10 A:5 0:3")
        # A:5 is excluded (not a valid int taxid)
        assert total == 13
        assert "A" not in votes

    def test_top_n_truncation(self):
        """Truncated to top_n entries."""
        kmer_str = " ".join(f"{i}:{100-i}" for i in range(20))
        votes, named, total, human = _parse_kmer_votes(
            kmer_str, top_n=5,
        )
        assert len(votes.split(";")) == 5

    def test_descending_sort_order(self):
        """Votes are sorted descending by count."""
        votes, named, total, human = _parse_kmer_votes("2:5 562:20 0:3")
        parts = votes.split(";")
        assert parts[0] == "562:20"
        assert parts[1] == "2:5"
        assert parts[2] == "0:3"


class TestWriteKraken2ReadDetailBed:
    """Tests for _write_kraken2_read_detail_bed."""

    def _make_result_with_detail(self, details):
        """Create a Kraken2Runner.Result with given per_read_detail."""
        result = Kraken2Runner.Result()
        result.per_read_detail = details
        for rname, d in details.items():
            if d["is_nonhuman"]:
                result.nonhuman_read_names.add(rname)
                result.nonhuman_count += 1
        return result

    def test_basic_bed_output(self):
        """BED file has correct header and data rows."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.kraken2_reads.bed.gz")

            informative = {
                "chr1:1000:A:T": {"read001", "read002"},
            }
            alt_informative = {
                "chr1:1000:A:T": {"read001"},
            }
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C",
                    "taxid": 562,
                    "domain": "Bacteria",
                    "guard_status": "PASS",
                    "is_nonhuman": True,
                    "kmer_string": "562:25 2:8 0:3",
                },
                "read002": {
                    "status": "C",
                    "taxid": 9606,
                    "domain": "Human",
                    "guard_status": "HUMAN",
                    "is_nonhuman": False,
                    "kmer_string": "9606:40 0:2",
                },
            })
            name_map = {562: "Escherichia_coli", 9606: "Homo_sapiens", 2: "Bacteria"}

            _write_kraken2_read_detail_bed(
                out_path, informative, alt_informative, result, name_map,
            )

            # File exists and is bgzipped
            assert os.path.isfile(out_path)
            # Tabix index exists
            assert os.path.isfile(out_path + ".tbi")

            # Read contents
            with gzip.open(out_path, "rt") as fh:
                lines = fh.readlines()

            # Header line
            assert lines[0].startswith("#chrom\t")
            header_cols = lines[0].strip().split("\t")
            assert "variant" in header_cols
            assert "read_name" in header_cols
            assert "kmer_votes" in header_cols

            # Data rows (2 reads)
            data_lines = [l for l in lines if not l.startswith("#")]
            assert len(data_lines) == 2

            # Verify first data row
            fields = data_lines[0].strip().split("\t")
            assert fields[0] == "chr1"  # chrom
            assert fields[1] == "1000"  # chromStart
            assert fields[2] == "1001"  # chromEnd (pos + len("A"))
            assert fields[3] == "chr1:1000:A:T"  # variant
            assert fields[4] == "read001"  # read_name (sorted first)
            assert fields[5] == "DKA"  # read_set

    def test_tabix_queryable(self):
        """BED file can be queried with tabix."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            informative = {
                "chr1:1000:A:T": {"read001"},
                "chr2:2000:G:C": {"read002"},
            }
            alt_informative = {}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562, "domain": "Bacteria",
                    "guard_status": "PASS", "is_nonhuman": True,
                    "kmer_string": "562:10",
                },
                "read002": {
                    "status": "C", "taxid": 562, "domain": "Bacteria",
                    "guard_status": "PASS", "is_nonhuman": True,
                    "kmer_string": "562:10",
                },
            })

            _write_kraken2_read_detail_bed(
                out_path, informative, alt_informative, result, None,
            )

            # Query with pysam tabix
            tbx = pysam.TabixFile(out_path)
            assert "chr1" in tbx.contigs
            assert "chr2" in tbx.contigs

            rows = list(tbx.fetch("chr1", 999, 1002))
            assert len(rows) == 1
            assert "read001" in rows[0]

            rows2 = list(tbx.fetch("chr2", 1999, 2002))
            assert len(rows2) == 1
            assert "read002" in rows2[0]
            tbx.close()

    def test_sort_order(self):
        """Rows are sorted by chrom, pos, then read_name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            informative = {
                "chr2:500:A:T": {"readB"},
                "chr1:200:G:C": {"readA", "readC"},
                "chr1:100:T:A": {"readD"},
            }
            detail = {}
            for var_reads in informative.values():
                for rname in var_reads:
                    detail[rname] = {
                        "status": "C", "taxid": 562, "domain": "Bacteria",
                        "guard_status": "PASS", "is_nonhuman": True,
                        "kmer_string": "562:10",
                    }
            result = self._make_result_with_detail(detail)

            _write_kraken2_read_detail_bed(
                out_path, informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                lines = [l for l in fh if not l.startswith("#")]

            # Expected order: chr1:100 readD, chr1:200 readA, chr1:200 readC, chr2:500 readB
            assert "chr1\t100\t" in lines[0]
            assert "readD" in lines[0]
            assert "chr1\t200\t" in lines[1]
            assert "readA" in lines[1]
            assert "chr1\t200\t" in lines[2]
            assert "readC" in lines[2]
            assert "chr2\t500\t" in lines[3]
            assert "readB" in lines[3]

    def test_unclassified_reads(self):
        """Unclassified reads have correct empty fields."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            informative = {"chr1:1000:A:T": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "U", "taxid": 0, "domain": "Unclassified",
                    "guard_status": "UNCLASSIFIED", "is_nonhuman": False,
                    "kmer_string": "",
                },
            })

            _write_kraken2_read_detail_bed(
                out_path, informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                lines = [l for l in fh if not l.startswith("#")]

            fields = lines[0].strip().split("\t")
            assert fields[6] == "U"  # kraken2_status
            assert fields[7] == "0"  # assigned_taxid
            assert fields[8] == "."  # assigned_taxon
            assert fields[9] == "Unclassified"  # domain
            assert fields[10] == "UNCLASSIFIED"  # guard_status
            assert fields[11] == "false"  # is_nonhuman

    def test_multi_variant_same_read(self):
        """A read spanning multiple variants appears once per variant."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            informative = {
                "chr1:1000:A:T": {"shared_read"},
                "chr1:1005:G:C": {"shared_read"},
            }
            result = self._make_result_with_detail({
                "shared_read": {
                    "status": "C", "taxid": 562, "domain": "Bacteria",
                    "guard_status": "PASS", "is_nonhuman": True,
                    "kmer_string": "562:10",
                },
            })

            _write_kraken2_read_detail_bed(
                out_path, informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                lines = [l for l in fh if not l.startswith("#")]

            assert len(lines) == 2  # one row per variant
            assert "chr1:1000:A:T" in lines[0]
            assert "chr1:1005:G:C" in lines[1]

    def test_indel_chrom_end(self):
        """Indel variants have correct chromEnd (pos + len(ref))."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            informative = {"chr1:1000:ATG:A": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562, "domain": "Bacteria",
                    "guard_status": "PASS", "is_nonhuman": True,
                    "kmer_string": "562:10",
                },
            })

            _write_kraken2_read_detail_bed(
                out_path, informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                lines = [l for l in fh if not l.startswith("#")]

            fields = lines[0].strip().split("\t")
            assert fields[1] == "1000"  # chromStart
            assert fields[2] == "1003"  # chromEnd = 1000 + len("ATG")


class TestExtractSoftclips:
    """Tests for _extract_softclips helper."""

    def test_no_cigar(self):
        assert _extract_softclips(None) == (0, 0)

    def test_empty_cigar(self):
        assert _extract_softclips([]) == (0, 0)

    def test_no_clips(self):
        # 100M
        assert _extract_softclips([(0, 100)]) == (0, 0)

    def test_left_clip_only(self):
        # 42S58M
        assert _extract_softclips([(4, 42), (0, 58)]) == (42, 0)

    def test_right_clip_only(self):
        # 58M5S
        assert _extract_softclips([(0, 58), (4, 5)]) == (0, 5)

    def test_both_clips(self):
        # 10S80M10S
        assert _extract_softclips([(4, 10), (0, 80), (4, 10)]) == (10, 10)

    def test_single_softclip_op(self):
        # Edge case: entire read is soft-clipped (single op)
        assert _extract_softclips([(4, 100)]) == (100, 0)

    def test_hard_clips_ignored(self):
        # 5H10S80M5S3H — hard clips (op 5) do not mask soft clips
        assert _extract_softclips([(5, 5), (4, 10), (0, 80), (4, 5), (5, 3)]) == (10, 5)


class TestWriteKraken2SpanBed:
    """Tests for _write_kraken2_span_bed."""

    def _make_result_with_detail(self, details):
        """Create a Kraken2Runner.Result with given per_read_detail."""
        result = Kraken2Runner.Result()
        result.per_read_detail = details
        for rname, d in details.items():
            if d["is_nonhuman"]:
                result.nonhuman_read_names.add(rname)
                result.nonhuman_count += 1
        return result

    def test_basic_span_bed_output(self):
        """Span BED file has correct header and data rows."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.kraken2_spans.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 100000, "end": 100150,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 5,
                    "has_sa": False, "is_supplementary": False,
                }],
                "read002": [{
                    "chrom": "chr1", "start": 100020, "end": 100170,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {
                "chr1:100050:A:T": {"read001", "read002"},
            }
            alt_informative = {
                "chr1:100050:A:T": {"read001"},
            }
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562,
                    "domain": "Bacteria", "guard_status": "PASS",
                    "is_nonhuman": True, "kmer_string": "562:25",
                },
                "read002": {
                    "status": "C", "taxid": 9606,
                    "domain": "Human", "guard_status": "HUMAN",
                    "is_nonhuman": False, "kmer_string": "9606:40",
                },
            })
            name_map = {562: "Escherichia_coli", 9606: "Homo_sapiens"}

            _write_kraken2_span_bed(
                out_path, alignment_meta,
                informative, alt_informative, result, name_map,
            )

            assert os.path.isfile(out_path)
            assert os.path.isfile(out_path + ".tbi")

            with gzip.open(out_path, "rt") as fh:
                lines = fh.readlines()

            # Header
            header = lines[0].strip().split("\t")
            assert header[0] == "#chrom"
            assert "taxon_name" in header
            assert "softclip_left" in header
            assert "is_split" in header
            assert "is_supplementary" in header

            data = [l for l in lines if not l.startswith("#")]
            assert len(data) == 2

            # First row (sorted by start → read001 at 100000)
            f = data[0].strip().split("\t")
            assert f[0] == "chr1"
            assert f[1] == "100000"
            assert f[2] == "100150"
            assert f[3] == "Escherichia_coli"
            assert f[4] == "Bacteria"
            assert f[5] == "PASS"
            assert f[6] == "true"
            assert f[7] == "read001"
            assert f[8] == "chr1:100050:A:T"
            assert f[9] == "DKA"
            assert f[10] == "60"
            assert f[11] == "0"
            assert f[12] == "5"
            assert f[13] == "false"
            assert f[14] == "false"

    def test_split_read_produces_two_rows(self):
        """A split read with SA tag produces two BED rows."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            alignment_meta = {
                "read003": [
                    {
                        "chrom": "chr1", "start": 100030, "end": 100100,
                        "mapq": 25, "softclip_left": 42, "softclip_right": 0,
                        "has_sa": True, "is_supplementary": False,
                    },
                    {
                        "chrom": "chr7", "start": 50000, "end": 50070,
                        "mapq": 0, "softclip_left": 0, "softclip_right": 42,
                        "has_sa": True, "is_supplementary": True,
                    },
                ],
            }
            informative = {"chr1:100050:A:T": {"read003"}}
            result = self._make_result_with_detail({
                "read003": {
                    "status": "C", "taxid": 2,
                    "domain": "Bacteria", "guard_status": "HHG",
                    "is_nonhuman": False, "kmer_string": "2:10 9606:2",
                },
            })
            name_map = {2: "Bacteria"}

            _write_kraken2_span_bed(
                out_path, alignment_meta,
                informative, {}, result, name_map,
            )

            with gzip.open(out_path, "rt") as fh:
                data = [l for l in fh if not l.startswith("#")]

            assert len(data) == 2

            # Primary (chr1)
            f1 = data[0].strip().split("\t")
            assert f1[0] == "chr1"
            assert f1[1] == "100030"
            assert f1[2] == "100100"
            assert f1[13] == "true"   # is_split
            assert f1[14] == "false"  # is_supplementary (primary)
            assert f1[11] == "42"     # softclip_left

            # Supplementary (chr7)
            f2 = data[1].strip().split("\t")
            assert f2[0] == "chr7"
            assert f2[1] == "50000"
            assert f2[2] == "50070"
            assert f2[13] == "true"   # is_split
            assert f2[14] == "true"   # is_supplementary
            assert f2[12] == "42"     # softclip_right

            # Both have same read_name and classification
            assert f1[7] == f2[7] == "read003"
            assert f1[3] == f2[3] == "Bacteria"

    def test_unclassified_read(self):
        """Unclassified reads have taxon_name=Unclassified."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 1000, "end": 1100,
                    "mapq": 30, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {"chr1:1050:A:T": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "U", "taxid": 0,
                    "domain": "Unclassified", "guard_status": "UNCLASSIFIED",
                    "is_nonhuman": False, "kmer_string": "",
                },
            })

            _write_kraken2_span_bed(
                out_path, alignment_meta,
                informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                data = [l for l in fh if not l.startswith("#")]

            f = data[0].strip().split("\t")
            assert f[3] == "Unclassified"
            assert f[4] == "Unclassified"
            assert f[5] == "UNCLASSIFIED"
            assert f[6] == "false"

    def test_unknown_taxid_fallback(self):
        """Without name_map, taxon uses Unknown_taxid_NNN format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 1000, "end": 1100,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {"chr1:1050:A:T": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 999999,
                    "domain": "Bacteria", "guard_status": "PASS",
                    "is_nonhuman": True, "kmer_string": "999999:10",
                },
            })

            _write_kraken2_span_bed(
                out_path, alignment_meta,
                informative, {}, result, None,  # no name_map
            )

            with gzip.open(out_path, "rt") as fh:
                data = [l for l in fh if not l.startswith("#")]

            f = data[0].strip().split("\t")
            assert f[3] == "Unknown_taxid_999999"

    def test_multi_variant_comma_separated(self):
        """Read spanning multiple variants has comma-separated variant keys."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            alignment_meta = {
                "shared_read": [{
                    "chrom": "chr1", "start": 1000, "end": 1200,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {
                "chr1:1050:A:T": {"shared_read"},
                "chr1:1100:G:C": {"shared_read"},
            }
            result = self._make_result_with_detail({
                "shared_read": {
                    "status": "C", "taxid": 562,
                    "domain": "Bacteria", "guard_status": "PASS",
                    "is_nonhuman": True, "kmer_string": "562:10",
                },
            })

            _write_kraken2_span_bed(
                out_path, alignment_meta,
                informative, {}, result, {562: "Escherichia_coli"},
            )

            with gzip.open(out_path, "rt") as fh:
                data = [l for l in fh if not l.startswith("#")]

            # Single row (one alignment record)
            assert len(data) == 1
            f = data[0].strip().split("\t")
            # Variant column has comma-separated sorted keys
            variants = f[8].split(",")
            assert len(variants) == 2
            assert "chr1:1050:A:T" in variants
            assert "chr1:1100:G:C" in variants

    def test_tabix_queryable(self):
        """Span BED file can be queried with tabix."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 1000, "end": 1100,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
                "read002": [{
                    "chrom": "chr2", "start": 2000, "end": 2100,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {
                "chr1:1050:A:T": {"read001"},
                "chr2:2050:G:C": {"read002"},
            }
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562, "domain": "Bacteria",
                    "guard_status": "PASS", "is_nonhuman": True,
                    "kmer_string": "562:10",
                },
                "read002": {
                    "status": "C", "taxid": 562, "domain": "Bacteria",
                    "guard_status": "PASS", "is_nonhuman": True,
                    "kmer_string": "562:10",
                },
            })

            _write_kraken2_span_bed(
                out_path, alignment_meta,
                informative, {}, result, None,
            )

            tbx = pysam.TabixFile(out_path)
            assert "chr1" in tbx.contigs
            assert "chr2" in tbx.contigs

            rows = list(tbx.fetch("chr1", 999, 1101))
            assert len(rows) == 1
            assert "read001" in rows[0]

            rows2 = list(tbx.fetch("chr2", 1999, 2101))
            assert len(rows2) == 1
            assert "read002" in rows2[0]
            tbx.close()

    def test_sort_order(self):
        """Rows sorted by chrom, start, read_name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.bed.gz")

            alignment_meta = {
                "readB": [{
                    "chrom": "chr2", "start": 500, "end": 600,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
                "readA": [{
                    "chrom": "chr1", "start": 200, "end": 300,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
                "readD": [{
                    "chrom": "chr1", "start": 100, "end": 200,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {
                "chr2:550:A:T": {"readB"},
                "chr1:250:G:C": {"readA"},
                "chr1:150:T:A": {"readD"},
            }
            detail = {}
            for rname in alignment_meta:
                detail[rname] = {
                    "status": "C", "taxid": 562, "domain": "Bacteria",
                    "guard_status": "PASS", "is_nonhuman": True,
                    "kmer_string": "562:10",
                }
            result = self._make_result_with_detail(detail)

            _write_kraken2_span_bed(
                out_path, alignment_meta,
                informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                data = [l for l in fh if not l.startswith("#")]

            assert "chr1\t100\t" in data[0]
            assert "readD" in data[0]
            assert "chr1\t200\t" in data[1]
            assert "readA" in data[1]
            assert "chr2\t500\t" in data[2]
            assert "readB" in data[2]


class TestWriteKraken2ExpandedSpanBed:
    """Tests for _write_kraken2_expanded_span_bed."""

    def _make_result_with_detail(self, details):
        """Create a Kraken2Runner.Result with given per_read_detail."""
        result = Kraken2Runner.Result()
        result.per_read_detail = details
        for rname, d in details.items():
            if d["is_nonhuman"]:
                result.nonhuman_read_names.add(rname)
                result.nonhuman_count += 1
        return result

    def test_basic_expanded_coordinates(self):
        """Expanded BED extends start/end by soft-clip lengths."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.expanded.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 100000, "end": 100150,
                    "mapq": 60, "softclip_left": 20, "softclip_right": 5,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {"chr1:100050:A:T": {"read001"}}
            alt_informative = {"chr1:100050:A:T": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562,
                    "domain": "Bacteria", "guard_status": "PASS",
                    "is_nonhuman": True, "kmer_string": "562:25",
                },
            })
            name_map = {562: "Escherichia_coli"}

            _write_kraken2_expanded_span_bed(
                out_path, alignment_meta,
                informative, alt_informative, result, name_map,
            )

            assert os.path.isfile(out_path)
            assert os.path.isfile(out_path + ".tbi")

            with gzip.open(out_path, "rt") as fh:
                lines = fh.readlines()

            # Header
            header = lines[0].strip().split("\t")
            assert header[0] == "#chrom"
            assert "aligned_start" in header
            assert "aligned_end" in header

            data = [l for l in lines if not l.startswith("#")]
            assert len(data) == 1

            f = data[0].strip().split("\t")
            # Expanded coordinates: start=100000-20=99980, end=100150+5=100155
            assert f[0] == "chr1"
            assert f[1] == "99980"
            assert f[2] == "100155"
            assert f[3] == "Escherichia_coli"
            assert f[4] == "Bacteria"
            assert f[5] == "PASS"
            assert f[6] == "true"
            assert f[7] == "read001"
            assert f[8] == "chr1:100050:A:T"
            assert f[9] == "DKA"
            assert f[10] == "60"
            assert f[11] == "20"  # softclip_left
            assert f[12] == "5"   # softclip_right
            assert f[13] == "false"
            assert f[14] == "false"
            # aligned_start and aligned_end (original mapped coordinates)
            assert f[15] == "100000"
            assert f[16] == "100150"

    def test_expanded_start_clamped_at_zero(self):
        """Expanded start is clamped to 0 when softclip exceeds start."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.expanded.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 10, "end": 160,
                    "mapq": 60, "softclip_left": 50, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {"chr1:50:A:T": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562,
                    "domain": "Bacteria", "guard_status": "PASS",
                    "is_nonhuman": True, "kmer_string": "562:10",
                },
            })

            _write_kraken2_expanded_span_bed(
                out_path, alignment_meta,
                informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                data = [l for l in fh if not l.startswith("#")]

            f = data[0].strip().split("\t")
            # max(0, 10 - 50) = 0
            assert f[1] == "0"
            assert f[2] == "160"
            assert f[15] == "10"   # aligned_start
            assert f[16] == "160"  # aligned_end

    def test_no_clips_matches_standard(self):
        """Without soft clips, expanded coords equal aligned coords."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.expanded.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 1000, "end": 1100,
                    "mapq": 60, "softclip_left": 0, "softclip_right": 0,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {"chr1:1050:A:T": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562,
                    "domain": "Bacteria", "guard_status": "PASS",
                    "is_nonhuman": True, "kmer_string": "562:10",
                },
            })

            _write_kraken2_expanded_span_bed(
                out_path, alignment_meta,
                informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                data = [l for l in fh if not l.startswith("#")]

            f = data[0].strip().split("\t")
            # No clips: expanded = aligned
            assert f[1] == "1000"
            assert f[2] == "1100"
            assert f[15] == "1000"
            assert f[16] == "1100"

    def test_split_read_expanded(self):
        """Split read expanded BED produces two rows with correct expansion."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.expanded.bed.gz")

            alignment_meta = {
                "read003": [
                    {
                        "chrom": "chr1", "start": 100030, "end": 100100,
                        "mapq": 25, "softclip_left": 42, "softclip_right": 0,
                        "has_sa": True, "is_supplementary": False,
                    },
                    {
                        "chrom": "chr7", "start": 50000, "end": 50070,
                        "mapq": 0, "softclip_left": 0, "softclip_right": 42,
                        "has_sa": True, "is_supplementary": True,
                    },
                ],
            }
            informative = {"chr1:100050:A:T": {"read003"}}
            result = self._make_result_with_detail({
                "read003": {
                    "status": "C", "taxid": 2,
                    "domain": "Bacteria", "guard_status": "HHG",
                    "is_nonhuman": False, "kmer_string": "2:10 9606:2",
                },
            })

            _write_kraken2_expanded_span_bed(
                out_path, alignment_meta,
                informative, {}, result, {2: "Bacteria"},
            )

            with gzip.open(out_path, "rt") as fh:
                data = [l for l in fh if not l.startswith("#")]

            assert len(data) == 2

            # Primary (chr1): start=100030-42=99988, end=100100+0=100100
            f1 = data[0].strip().split("\t")
            assert f1[0] == "chr1"
            assert f1[1] == "99988"
            assert f1[2] == "100100"
            assert f1[15] == "100030"  # aligned_start
            assert f1[16] == "100100"  # aligned_end
            assert f1[13] == "true"    # is_split

            # Supplementary (chr7): start=50000-0=50000, end=50070+42=50112
            f2 = data[1].strip().split("\t")
            assert f2[0] == "chr7"
            assert f2[1] == "50000"
            assert f2[2] == "50112"
            assert f2[15] == "50000"   # aligned_start
            assert f2[16] == "50070"   # aligned_end
            assert f2[14] == "true"    # is_supplementary

    def test_tabix_queryable(self):
        """Expanded BED file can be queried with tabix."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.expanded.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 1000, "end": 1100,
                    "mapq": 60, "softclip_left": 50, "softclip_right": 30,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {"chr1:1050:A:T": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562, "domain": "Bacteria",
                    "guard_status": "PASS", "is_nonhuman": True,
                    "kmer_string": "562:10",
                },
            })

            _write_kraken2_expanded_span_bed(
                out_path, alignment_meta,
                informative, {}, result, None,
            )

            tbx = pysam.TabixFile(out_path)
            assert "chr1" in tbx.contigs
            # expanded: start=950, end=1130
            rows = list(tbx.fetch("chr1", 949, 1131))
            assert len(rows) == 1
            assert "read001" in rows[0]
            tbx.close()

    def test_column_count(self):
        """Expanded BED has 17 columns (15 standard + 2 extra)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test.expanded.bed.gz")

            alignment_meta = {
                "read001": [{
                    "chrom": "chr1", "start": 5000, "end": 5150,
                    "mapq": 60, "softclip_left": 10, "softclip_right": 20,
                    "has_sa": False, "is_supplementary": False,
                }],
            }
            informative = {"chr1:5050:A:T": {"read001"}}
            result = self._make_result_with_detail({
                "read001": {
                    "status": "C", "taxid": 562,
                    "domain": "Bacteria", "guard_status": "PASS",
                    "is_nonhuman": True, "kmer_string": "562:10",
                },
            })

            _write_kraken2_expanded_span_bed(
                out_path, alignment_meta,
                informative, {}, result, None,
            )

            with gzip.open(out_path, "rt") as fh:
                lines = fh.readlines()

            header = lines[0].strip().split("\t")
            assert len(header) == 17

            data = [l for l in lines if not l.startswith("#")]
            f = data[0].strip().split("\t")
            assert len(f) == 17
