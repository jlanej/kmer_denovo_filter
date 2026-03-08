"""Tests for Kraken2Runner bacterial content classification."""

import os
import subprocess
import tempfile
from unittest import mock

import pytest

from kmer_denovo_filter.kmer_utils import Kraken2Runner, _BACTERIA_TAXID


class TestKraken2Result:
    """Tests for Kraken2Runner.Result container."""

    def test_empty_result(self):
        r = Kraken2Runner.Result()
        assert r.total == 0
        assert r.classified == 0
        assert r.unclassified == 0
        assert r.bacterial_count == 0
        assert r.human_count == 0
        assert r.root_count == 0
        assert r.bacterial_read_names == set()

    def test_summary_empty(self):
        r = Kraken2Runner.Result()
        s = r.summary()
        assert "0 reads" in s
        assert "0 bacterial" in s

    def test_summary_with_counts(self):
        r = Kraken2Runner.Result()
        r.total = 100
        r.classified = 80
        r.unclassified = 20
        r.bacterial_count = 10
        r.human_count = 60
        r.root_count = 10
        s = r.summary()
        assert "100 reads" in s
        assert "80 classified" in s
        assert "10 bacterial" in s
        assert "10.0%" in s
        assert "60 human" in s
        assert "10 root" in s


class TestKraken2RunnerInit:
    """Tests for Kraken2Runner constructor."""

    def test_default_params(self):
        kr = Kraken2Runner("/fake/db")
        assert kr.db_path == "/fake/db"
        assert kr.confidence == 0.0
        assert kr.threads == 1

    def test_custom_params(self):
        kr = Kraken2Runner("/db", confidence=0.2, threads=8)
        assert kr.confidence == 0.2
        assert kr.threads == 8


class TestKraken2RunnerClassify:
    """Tests for classify_sequences with mocked subprocess."""

    def test_empty_sequences_returns_empty_result(self):
        kr = Kraken2Runner("/fake/db")
        result = kr.classify_sequences({})
        assert result.total == 0

    def test_empty_list_returns_empty_result(self):
        kr = Kraken2Runner("/fake/db")
        result = kr.classify_sequences([])
        assert result.total == 0

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_classify_parses_output(self, mock_popen):
        """Test that kraken2 output is correctly parsed."""
        # Simulate kraken2 per-read output:
        # C = classified, U = unclassified
        # Format: status\tread_name\ttaxid\tlength\tkmers
        kraken2_output = (
            "C\tread1\t2\t100\t2:10 0:5\n"    # bacterial (taxid 2)
            "C\tread2\t9606\t100\t9606:20\n"   # human
            "U\tread3\t0\t100\t0:15\n"          # unclassified
            "C\tread4\t1\t100\t1:8\n"           # root
            "C\tread5\t562\t100\t562:12\n"      # E. coli (bacterial)
        )

        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (
            kraken2_output.encode(), b""
        )
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")

        # Without taxonomy files, only exact taxid 2 is bacterial
        with mock.patch.object(
            Kraken2Runner, '_load_bacterial_taxids', return_value=None,
        ):
            result = kr.classify_sequences({
                "read1": "ACGTACGTACGT",
                "read2": "ACGTACGTACGT",
                "read3": "ACGTACGTACGT",
                "read4": "ACGTACGTACGT",
                "read5": "ACGTACGTACGT",
            })

        assert result.total == 5
        assert result.classified == 4
        assert result.unclassified == 1
        assert result.bacterial_count == 1  # only exact taxid 2
        assert result.human_count == 1
        assert result.root_count == 1
        assert "read1" in result.bacterial_read_names

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_classify_with_bacterial_taxids(self, mock_popen):
        """Test bacterial lineage-aware matching."""
        kraken2_output = (
            "C\tread1\t2\t100\t2:10\n"      # Bacteria
            "C\tread2\t562\t100\t562:12\n"   # E. coli (descendant)
            "C\tread3\t9606\t100\t9606:20\n" # human
        )

        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (
            kraken2_output.encode(), b""
        )
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")

        # Simulate bacterial taxid set including descendants
        bacterial_set = {2, 562, 1234}
        with mock.patch.object(
            Kraken2Runner, '_load_bacterial_taxids',
            return_value=bacterial_set,
        ):
            result = kr.classify_sequences({
                "read1": "ACGTACGTACGT",
                "read2": "ACGTACGTACGT",
                "read3": "ACGTACGTACGT",
            })

        assert result.bacterial_count == 2
        assert "read1" in result.bacterial_read_names
        assert "read2" in result.bacterial_read_names
        assert result.human_count == 1

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_classify_nonzero_exit(self, mock_popen):
        """Non-zero exit code returns empty result with correct total."""
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (b"", b"error msg")
        mock_proc.returncode = 1
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        result = kr.classify_sequences({"r1": "ACGT"})
        assert result.total == 1
        assert result.classified == 0

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_classify_list_input(self, mock_popen):
        """Accepts list of (name, seq) tuples."""
        kraken2_output = "C\tread1\t9606\t100\t9606:20\n"
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (
            kraken2_output.encode(), b""
        )
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        with mock.patch.object(
            Kraken2Runner, '_load_bacterial_taxids', return_value=None,
        ):
            result = kr.classify_sequences([("read1", "ACGTACGT")])

        assert result.total == 1
        assert result.human_count == 1

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_temp_fastq_cleaned_up(self, mock_popen):
        """Temporary FASTQ file is cleaned up after classification."""
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (b"", b"")
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        with mock.patch.object(
            Kraken2Runner, '_load_bacterial_taxids', return_value=None,
        ):
            with tempfile.TemporaryDirectory() as tmpdir:
                kr.classify_sequences({"r1": "ACGT"}, tmpdir=tmpdir)
                # After classification, no .fq files should remain
                remaining = [
                    f for f in os.listdir(tmpdir) if f.endswith(".fq")
                ]
                assert remaining == []

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_fastq_content(self, mock_popen):
        """Verify FASTQ written to kraken2 has correct format."""
        captured_args = {}

        def capture_popen(cmd, **kwargs):
            # Read the fastq file that was just written
            fq_path = cmd[-1]  # last argument is the fastq path
            with open(fq_path) as fh:
                captured_args["fastq_content"] = fh.read()
            proc = mock.MagicMock()
            proc.communicate.return_value = (b"", b"")
            proc.returncode = 0
            return proc

        mock_popen.side_effect = capture_popen

        kr = Kraken2Runner("/fake/db")
        with mock.patch.object(
            Kraken2Runner, '_load_bacterial_taxids', return_value=None,
        ):
            kr.classify_sequences({"myread": "ACGTACGT"})

        content = captured_args["fastq_content"]
        assert "@myread\n" in content
        assert "ACGTACGT\n" in content
        assert "+\n" in content

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_command_includes_confidence(self, mock_popen):
        """Verify confidence flag is passed to kraken2."""
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (b"", b"")
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/my/db", confidence=0.3, threads=4)
        with mock.patch.object(
            Kraken2Runner, '_load_bacterial_taxids', return_value=None,
        ):
            kr.classify_sequences({"r1": "ACGT"})

        cmd = mock_popen.call_args[0][0]
        assert "kraken2" in cmd
        assert "--db" in cmd
        assert "/my/db" in cmd
        assert "--confidence" in cmd
        assert "0.3" in cmd
        assert "--threads" in cmd
        assert "4" in cmd


class TestLoadBacterialTaxids:
    """Tests for _load_bacterial_taxids taxonomy parsing."""

    def test_missing_file_returns_none(self):
        result = Kraken2Runner._load_bacterial_taxids("/nonexistent/db")
        assert result is None

    def test_parses_nodes_dmp(self):
        """Build a simple taxonomy and verify bacteria detection."""
        with tempfile.TemporaryDirectory() as db:
            tax_dir = os.path.join(db, "taxonomy")
            os.makedirs(tax_dir)
            # Write a minimal nodes.dmp:
            # taxid | parent_taxid | rank | ...
            # 1 = root (parent of self)
            # 2 = Bacteria (parent = 1)
            # 562 = E. coli (parent = 2)
            # 9606 = Homo sapiens (parent = 1)
            nodes = os.path.join(tax_dir, "nodes.dmp")
            with open(nodes, "w") as fh:
                fh.write("1\t|\t1\t|\tno rank\t|\n")
                fh.write("2\t|\t1\t|\tsuperkingdom\t|\n")
                fh.write("562\t|\t2\t|\tspecies\t|\n")
                fh.write("9606\t|\t1\t|\tspecies\t|\n")

            bacterial = Kraken2Runner._load_bacterial_taxids(db)
            assert bacterial is not None
            assert 2 in bacterial       # Bacteria itself
            assert 562 in bacterial     # E. coli (descendant)
            assert 9606 not in bacterial  # Human
            assert 1 not in bacterial     # Root


class TestBacterialTaxidConstant:
    """Verify the bacterial taxid constant."""

    def test_bacteria_taxid(self):
        assert _BACTERIA_TAXID == 2
