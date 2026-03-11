"""Tests for Kraken2Runner non-human content classification."""

import os
import subprocess
import tempfile
from unittest import mock

import pytest

from kmer_denovo_filter.kmer_utils import (
    Kraken2Runner,
    _ARCHAEA_TAXID,
    _BACTERIA_TAXID,
    _EUKARYOTA_TAXID,
    _FUNGI_TAXID,
    _HUMAN_TAXID,
    _METAZOA_TAXID,
    _VIRIDIPLANTAE_TAXID,
)


class TestKraken2Result:
    """Tests for Kraken2Runner.Result container."""

    def test_empty_result(self):
        r = Kraken2Runner.Result()
        assert r.total == 0
        assert r.classified == 0
        assert r.unclassified == 0
        assert r.bacterial_count == 0
        assert r.archaeal_count == 0
        assert r.fungal_count == 0
        assert r.protist_count == 0
        assert r.nonhuman_count == 0
        assert r.human_count == 0
        assert r.root_count == 0
        assert r.bacterial_read_names == set()
        assert r.archaeal_read_names == set()
        assert r.fungal_read_names == set()
        assert r.protist_read_names == set()
        assert r.nonhuman_read_names == set()

    def test_summary_empty(self):
        r = Kraken2Runner.Result()
        s = r.summary()
        assert "0 reads" in s
        assert "0 bacterial" in s
        assert "0 non-human" in s

    def test_summary_with_counts(self):
        r = Kraken2Runner.Result()
        r.total = 100
        r.classified = 80
        r.unclassified = 20
        r.bacterial_count = 10
        r.archaeal_count = 2
        r.fungal_count = 3
        r.protist_count = 1
        r.nonhuman_count = 16
        r.human_count = 60
        r.root_count = 4
        s = r.summary()
        assert "100 reads" in s
        assert "80 classified" in s
        assert "10 bacterial" in s
        assert "10.0%" in s
        assert "2 archaeal" in s
        assert "3 fungal" in s
        assert "1 protist" in s
        assert "16 non-human" in s
        assert "16.0%" in s
        assert "60 human" in s
        assert "4 root" in s

    def test_bacterial_fraction_empty(self):
        r = Kraken2Runner.Result()
        assert r.bacterial_fraction == 0.0

    def test_bacterial_fraction_with_counts(self):
        r = Kraken2Runner.Result()
        r.total = 200
        r.bacterial_count = 50
        assert r.bacterial_fraction == 0.25


class TestKraken2RunnerInit:
    """Tests for Kraken2Runner constructor."""

    def test_default_params(self):
        kr = Kraken2Runner("/fake/db")
        assert kr.db_path == "/fake/db"
        assert kr.confidence == 0.0
        assert kr.threads == 1
        assert kr.memory_mapping is False

    def test_custom_params(self):
        kr = Kraken2Runner(
            "/db", confidence=0.2, threads=8, memory_mapping=True,
        )
        assert kr.confidence == 0.2
        assert kr.threads == 8
        assert kr.memory_mapping is True


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
    def test_classify_parses_output(self, mock_popen, caplog):
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

        # Without taxonomy files, only exact taxid matching
        caplog.set_level("WARNING", logger="kmer_denovo_filter.kmer_utils")
        with mock.patch.object(
            Kraken2Runner, '_load_all_taxid_sets', return_value=None,
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
        assert "lineage matching is unavailable" in caplog.text
        assert "exact taxid matching only" in caplog.text

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

        # Simulate taxonomy sets including bacterial descendants
        taxid_sets = {
            "bacterial": {2, 562, 1234},
            "archaeal": set(),
            "fungal": set(),
            "protist": set(),
            "human_lineage": {9606, 1},
            "human_clade": {9606},
        }
        with mock.patch.object(
            Kraken2Runner, '_load_all_taxid_sets',
            return_value=taxid_sets,
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
        # Both bacterial reads should be non-human too
        assert result.nonhuman_count == 2
        assert "read1" in result.nonhuman_read_names
        assert "read2" in result.nonhuman_read_names

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
    def test_classify_nonzero_exit_logs_elapsed(self, mock_popen, caplog):
        """Non-zero exit code warning should include elapsed time."""
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (b"", b"kraken2 OOM")
        mock_proc.returncode = -9  # SIGKILL (OOM)
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        caplog.set_level("WARNING", logger="kmer_denovo_filter.kmer_utils")
        result = kr.classify_sequences({"r1": "ACGT"})
        assert result.total == 1
        assert result.classified == 0
        # Warning should include exit code and elapsed time
        assert "-9" in caplog.text
        assert "after" in caplog.text

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_classify_success_logs_completion(self, mock_popen, caplog):
        """Successful classification should log a completion message."""
        kraken2_output = "C\tread1\t9606\t100\t9606:20\n"
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (kraken2_output.encode(), b"")
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        caplog.set_level("INFO", logger="kmer_denovo_filter.kmer_utils")
        with mock.patch.object(
            Kraken2Runner, "_load_all_taxid_sets", return_value=None,
        ):
            result = kr.classify_sequences({"read1": "ACGTACGT"})

        assert result.total == 1
        # Completion message should appear in logs
        assert "classification complete" in caplog.text

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_classify_malformed_kmer_field(self, mock_popen):
        """Malformed kmer information field is ignored during parsing."""
        kraken2_output = (
            "C\tread1\t2\t100\tGARBAGE_KMER_DATA!!!\n"
            "C\tread2\t9606\t100\t\n"
            "\n"
            "short\n"
            "C\tread3\tbadtaxid\t100\tfoo\n"
        )
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (
            kraken2_output.encode(), b""
        )
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        with mock.patch.object(
            Kraken2Runner, '_load_all_taxid_sets', return_value=None,
        ):
            result = kr.classify_sequences({
                "read1": "ACGT", "read2": "ACGT", "read3": "ACGT",
            })
        # read1 classified bacterial, read2 classified human,
        # read3 has bad taxid so skipped, short line skipped
        assert result.classified == 2
        assert result.bacterial_count == 1
        assert result.human_count == 1

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
            Kraken2Runner, '_load_all_taxid_sets', return_value=None,
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
            Kraken2Runner, '_load_all_taxid_sets', return_value=None,
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
            Kraken2Runner, '_load_all_taxid_sets', return_value=None,
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
            Kraken2Runner, '_load_all_taxid_sets', return_value=None,
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
        assert "--memory-mapping" not in cmd

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_command_includes_memory_mapping_when_enabled(self, mock_popen):
        """Verify memory-mapping flag is passed to kraken2 when enabled."""
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (b"", b"")
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/my/db", memory_mapping=True)
        with mock.patch.object(
            Kraken2Runner, '_load_all_taxid_sets', return_value=None,
        ):
            kr.classify_sequences({"r1": "ACGT"})

        cmd = mock_popen.call_args[0][0]
        assert "--memory-mapping" in cmd


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

    def test_parses_root_level_nodes_dmp(self):
        """nodes.dmp at the DB root (PrackenDB layout) is found."""
        with tempfile.TemporaryDirectory() as db:
            # PrackenDB places nodes.dmp directly in the DB directory,
            # not under a taxonomy/ subdirectory.
            nodes = os.path.join(db, "nodes.dmp")
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

    def test_taxonomy_subdir_preferred_over_root(self):
        """taxonomy/nodes.dmp is preferred when both locations exist."""
        with tempfile.TemporaryDirectory() as db:
            # Create nodes.dmp in both locations with different content.
            # taxonomy/ version has taxid 562 (E. coli).
            tax_dir = os.path.join(db, "taxonomy")
            os.makedirs(tax_dir)
            with open(os.path.join(tax_dir, "nodes.dmp"), "w") as fh:
                fh.write("1\t|\t1\t|\tno rank\t|\n")
                fh.write("2\t|\t1\t|\tsuperkingdom\t|\n")
                fh.write("562\t|\t2\t|\tspecies\t|\n")

            # Root version has taxid 1280 (S. aureus) instead.
            with open(os.path.join(db, "nodes.dmp"), "w") as fh:
                fh.write("1\t|\t1\t|\tno rank\t|\n")
                fh.write("2\t|\t1\t|\tsuperkingdom\t|\n")
                fh.write("1280\t|\t2\t|\tspecies\t|\n")

            bacterial = Kraken2Runner._load_bacterial_taxids(db)
            assert bacterial is not None
            # Should use taxonomy/ version (562 present, 1280 absent)
            assert 562 in bacterial
            assert 1280 not in bacterial


class TestLoadParentMap:
    """Tests for _load_parent_map taxonomy parsing."""

    def test_missing_file_returns_none(self):
        assert Kraken2Runner._load_parent_map("/nonexistent/db") is None

    def test_parses_nodes_dmp(self):
        with tempfile.TemporaryDirectory() as db:
            tax_dir = os.path.join(db, "taxonomy")
            os.makedirs(tax_dir)
            nodes = os.path.join(tax_dir, "nodes.dmp")
            with open(nodes, "w") as fh:
                fh.write("1\t|\t1\t|\tno rank\t|\n")
                fh.write("2\t|\t1\t|\tsuperkingdom\t|\n")
                fh.write("562\t|\t2\t|\tspecies\t|\n")

            pm = Kraken2Runner._load_parent_map(db)
            assert pm is not None
            assert pm[1] == 1  # root
            assert pm[2] == 1  # bacteria -> root
            assert pm[562] == 2  # E. coli -> bacteria


class TestDescendantsOf:
    """Tests for _descendants_of helper."""

    def test_simple_tree(self):
        parent_map = {1: 1, 2: 1, 562: 2, 9606: 1}
        desc = Kraken2Runner._descendants_of(parent_map, 2)
        assert 2 in desc
        assert 562 in desc
        assert 9606 not in desc
        assert 1 not in desc


class TestAncestorsOf:
    """Tests for _ancestors_of helper."""

    def test_human_lineage(self):
        parent_map = {1: 1, 131567: 1, 2759: 131567, 33208: 2759, 9606: 33208}
        anc = Kraken2Runner._ancestors_of(parent_map, 9606)
        assert 9606 in anc
        assert 33208 in anc
        assert 2759 in anc
        assert 131567 in anc
        assert 1 in anc


class TestLoadAllTaxidSets:
    """Tests for _load_all_taxid_sets taxonomy computation."""

    def _write_mini_taxonomy(self, db_path):
        """Write a minimal nodes.dmp with all major domains."""
        nodes = os.path.join(db_path, "nodes.dmp")
        with open(nodes, "w") as fh:
            # Root and cellular organisms
            fh.write("1\t|\t1\t|\tno rank\t|\n")
            fh.write("131567\t|\t1\t|\tno rank\t|\n")
            # Bacteria
            fh.write("2\t|\t131567\t|\tsuperkingdom\t|\n")
            fh.write("562\t|\t2\t|\tspecies\t|\n")
            # Archaea
            fh.write("2157\t|\t131567\t|\tsuperkingdom\t|\n")
            fh.write("2287\t|\t2157\t|\tspecies\t|\n")
            # Eukaryota
            fh.write("2759\t|\t131567\t|\tsuperkingdom\t|\n")
            # Fungi
            fh.write("4751\t|\t2759\t|\tkingdom\t|\n")
            fh.write("5061\t|\t4751\t|\tspecies\t|\n")
            # Metazoa
            fh.write("33208\t|\t2759\t|\tkingdom\t|\n")
            fh.write("9606\t|\t33208\t|\tspecies\t|\n")
            # Viridiplantae
            fh.write("33090\t|\t2759\t|\tkingdom\t|\n")
            fh.write("3702\t|\t33090\t|\tspecies\t|\n")
            # Protist (eukaryotic, not metazoan/fungal/plant)
            fh.write("5794\t|\t2759\t|\tphylum\t|\n")
            fh.write("5820\t|\t5794\t|\tspecies\t|\n")
            # Virus (not under cellular organisms)
            fh.write("10239\t|\t1\t|\tsuperkingdom\t|\n")
            fh.write("11676\t|\t10239\t|\tspecies\t|\n")

    def test_missing_db_returns_none(self):
        assert Kraken2Runner._load_all_taxid_sets("/nonexistent") is None

    def test_all_domains_classified(self):
        with tempfile.TemporaryDirectory() as db:
            self._write_mini_taxonomy(db)
            sets = Kraken2Runner._load_all_taxid_sets(db)
            assert sets is not None

            # Bacteria
            assert 2 in sets["bacterial"]
            assert 562 in sets["bacterial"]
            assert 9606 not in sets["bacterial"]

            # Archaea
            assert 2157 in sets["archaeal"]
            assert 2287 in sets["archaeal"]
            assert 562 not in sets["archaeal"]

            # Fungi
            assert 4751 in sets["fungal"]
            assert 5061 in sets["fungal"]
            assert 562 not in sets["fungal"]

            # Protist: eukaryotic but not metazoan/fungal/plant
            assert 5794 in sets["protist"]
            assert 5820 in sets["protist"]
            assert 9606 not in sets["protist"]
            assert 4751 not in sets["protist"]
            assert 33208 not in sets["protist"]

            # Human lineage
            assert 9606 in sets["human_lineage"]
            assert 33208 in sets["human_lineage"]
            assert 2759 in sets["human_lineage"]
            assert 1 in sets["human_lineage"]
            assert 562 not in sets["human_lineage"]

            # Human clade
            assert 9606 in sets["human_clade"]
            assert 33208 not in sets["human_clade"]


class TestKrakenKmerTaxidParsing:
    """Tests for parsing Kraken2 k-mer detail field."""

    def test_extract_taxids_from_kmer_string(self):
        parsed = Kraken2Runner._extract_taxids_from_kmer_string(
            "562:12 A:6 0:3 |:| 9606:4 10239:2",
        )
        assert parsed == {0, 562, 9606, 10239}

    def test_extract_taxids_from_empty_kmer_string(self):
        assert Kraken2Runner._extract_taxids_from_kmer_string("") == set()


class TestTaxidConstants:
    """Verify the taxonomy ID constants."""

    def test_bacteria_taxid(self):
        assert _BACTERIA_TAXID == 2

    def test_archaea_taxid(self):
        assert _ARCHAEA_TAXID == 2157

    def test_fungi_taxid(self):
        assert _FUNGI_TAXID == 4751

    def test_eukaryota_taxid(self):
        assert _EUKARYOTA_TAXID == 2759

    def test_metazoa_taxid(self):
        assert _METAZOA_TAXID == 33208

    def test_viridiplantae_taxid(self):
        assert _VIRIDIPLANTAE_TAXID == 33090

    def test_human_taxid(self):
        assert _HUMAN_TAXID == 9606


class TestKrakenHomologyGuard:
    """Ensure human-homologous reads are not over-flagged as non-human."""

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_bacterial_assignment_with_human_kmers_not_flagged(
        self, mock_popen,
    ):
        kraken2_output = (
            "C\tread1\t562\t100\t562:8 9606:4\n"   # mixed; skip bacterial
            "C\tread2\t562\t100\t562:10 0:2\n"     # bacterial only
        )
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (kraken2_output.encode(), b"")
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        taxid_sets = {
            "bacterial": {2, 562},
            "archaeal": set(),
            "fungal": set(),
            "protist": set(),
            "human_lineage": {9606, 1},
            "human_clade": {9606},
        }
        with mock.patch.object(
            Kraken2Runner, "_load_all_taxid_sets", return_value=taxid_sets,
        ):
            result = kr.classify_sequences({
                "read1": "ACGTACGTACGT",
                "read2": "ACGTACGTACGT",
            })

        assert result.bacterial_count == 1
        assert "read1" not in result.bacterial_read_names
        assert "read2" in result.bacterial_read_names
        # Human homology guard also excludes from non-human
        assert "read1" not in result.nonhuman_read_names
        assert "read2" in result.nonhuman_read_names
        assert result.nonhuman_count == 1

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_archaeal_assignment_with_human_kmers_not_flagged(
        self, mock_popen,
    ):
        """Archaeal reads with human k-mer evidence are excluded."""
        kraken2_output = (
            "C\tread1\t2157\t100\t2157:8 9606:4\n"  # mixed; skip
            "C\tread2\t2157\t100\t2157:10 0:2\n"    # archaeal only
        )
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (kraken2_output.encode(), b"")
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        taxid_sets = {
            "bacterial": set(),
            "archaeal": {2157},
            "fungal": set(),
            "protist": set(),
            "human_lineage": {9606, 1},
            "human_clade": {9606},
        }
        with mock.patch.object(
            Kraken2Runner, "_load_all_taxid_sets", return_value=taxid_sets,
        ):
            result = kr.classify_sequences({
                "read1": "ACGTACGTACGT",
                "read2": "ACGTACGTACGT",
            })

        assert result.archaeal_count == 1
        assert "read1" not in result.archaeal_read_names
        assert "read2" in result.archaeal_read_names
        assert result.nonhuman_count == 1


class TestMultiDomainClassification:
    """Test classification across bacteria, archaea, fungi, protist."""

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_multi_domain_classification(self, mock_popen):
        """Reads from different domains are classified correctly."""
        kraken2_output = (
            "C\tread_bact\t562\t100\t562:10\n"     # E. coli (bacterial)
            "C\tread_arch\t2157\t100\t2157:10\n"   # Archaea
            "C\tread_fung\t4751\t100\t4751:10\n"   # Fungi
            "C\tread_prot\t5794\t100\t5794:10\n"   # Protist (e.g. Apicomplexa)
            "C\tread_hum\t9606\t100\t9606:20\n"    # Human
            "U\tread_unk\t0\t100\t0:15\n"           # Unclassified
        )
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (kraken2_output.encode(), b"")
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        taxid_sets = {
            "bacterial": {2, 562},
            "archaeal": {2157},
            "fungal": {4751},
            "protist": {5794},
            "human_lineage": {9606, 9605, 33208, 2759, 131567, 1},
            "human_clade": {9606},
        }
        with mock.patch.object(
            Kraken2Runner, "_load_all_taxid_sets", return_value=taxid_sets,
        ):
            result = kr.classify_sequences({
                "read_bact": "ACGTACGT",
                "read_arch": "ACGTACGT",
                "read_fung": "ACGTACGT",
                "read_prot": "ACGTACGT",
                "read_hum": "ACGTACGT",
                "read_unk": "ACGTACGT",
            })

        assert result.total == 6
        assert result.classified == 5
        assert result.unclassified == 1
        assert result.bacterial_count == 1
        assert result.archaeal_count == 1
        assert result.fungal_count == 1
        assert result.protist_count == 1
        assert result.human_count == 1
        assert result.nonhuman_count == 4

        assert "read_bact" in result.bacterial_read_names
        assert "read_arch" in result.archaeal_read_names
        assert "read_fung" in result.fungal_read_names
        assert "read_prot" in result.protist_read_names
        assert "read_hum" not in result.nonhuman_read_names
        for name in ("read_bact", "read_arch", "read_fung", "read_prot"):
            assert name in result.nonhuman_read_names

    @mock.patch("kmer_denovo_filter.kmer_utils.subprocess.Popen")
    def test_ambiguous_eukaryota_not_counted_nonhuman(self, mock_popen):
        """Reads classified at Eukaryota level are NOT counted as non-human."""
        kraken2_output = (
            "C\tread1\t2759\t100\t2759:10\n"   # Eukaryota (ambiguous)
            "C\tread2\t33208\t100\t33208:10\n"  # Metazoa (ancestor of human)
        )
        mock_proc = mock.MagicMock()
        mock_proc.communicate.return_value = (kraken2_output.encode(), b"")
        mock_proc.returncode = 0
        mock_popen.return_value = mock_proc

        kr = Kraken2Runner("/fake/db")
        taxid_sets = {
            "bacterial": set(),
            "archaeal": set(),
            "fungal": set(),
            "protist": set(),
            "human_lineage": {9606, 9605, 33208, 2759, 131567, 1},
            "human_clade": {9606},
        }
        with mock.patch.object(
            Kraken2Runner, "_load_all_taxid_sets", return_value=taxid_sets,
        ):
            result = kr.classify_sequences({
                "read1": "ACGTACGT",
                "read2": "ACGTACGT",
            })

        # Eukaryota and Metazoa are ancestors of human → not non-human
        assert result.nonhuman_count == 0
        assert result.bacterial_count == 0
        assert result.archaeal_count == 0
