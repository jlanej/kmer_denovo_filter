"""Tests for per-read Kraken2 classification detail BED output."""

import gzip
import os
import tempfile

import pysam
import pytest

from kmer_denovo_filter.kmer_utils import Kraken2Runner
from kmer_denovo_filter.pipeline import (
    _parse_kmer_votes,
    _write_kraken2_read_detail_bed,
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
