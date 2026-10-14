"""Module 1's child k-mer count stays inside its memory budget.

The jobs that motivated this died at the third in-memory doubling of a 4G-entry jellyfish table (23, 68, 45, 133,
88 GB of RSS, then a 264-GB allocation against a 248-GB cgroup).  The count now runs with ``--disk`` (a full table
is dumped and merged, never doubled), a hash sized from the memory budget and the input, and a Bloom filter that
keeps k-mers seen once out of the hash.  The Bloom filter admits a k-mer from its second occurrence, so the hash's
counts are one less than the truth and the dump bound is lowered by one.  The jellyfish facts are pinned here
against the real binary (skipped where jellyfish is not installed; CI installs it)."""
import os
import shutil
import subprocess

import pytest

from kmer_denovo_filter.core.jellyfish_wrappers import (
    BLOOM_MIN_ELEMENTS,
    HASH_MIN_ENTRIES,
    _estimate_jf_hash_size,
    _jf_bytes_per_entry,
    _parse_size,
    _pow2_at_least,
    size_child_count,
)
from kmer_denovo_filter.discovery.pipeline import _child_count_command

HAVE_JELLYFISH = shutil.which("jellyfish") is not None


def _fake_file(path, size):
    with open(path, "wb") as fh:
        fh.seek(size - 1)
        fh.write(b"\0")
    return path


class TestSizing:
    def test_bytes_per_entry_matches_the_jobs_that_died(self):
        # 4G entries at k=31 took 23.1 GB, 8G 45.2, 16G 88.0: about 5.8, 5.6 and 5.5 bytes per entry
        assert abs(_jf_bytes_per_entry(31, 2**32) * 2**32 / 1024**3 - 23.1) < 2.0
        assert abs(_jf_bytes_per_entry(31, 2**33) * 2**33 / 1024**3 - 45.2) < 2.0
        assert abs(_jf_bytes_per_entry(31, 2**34) * 2**34 / 1024**3 - 88.0) < 2.0

    def test_parse_size(self):
        assert _parse_size("4G") == 4 * 2**30 and _parse_size("500M") == 500 * 2**20 and _parse_size("1234") == 1234

    def test_a_wgs_child_fits_the_budget_with_the_filter(self, tmp_path):
        # the child of the killed job: a 41.5-GB CRAM, k=31, a 244-GB budget
        cram = _fake_file(str(tmp_path / "child.cram"), 1)
        os.truncate(cram, int(41.5 * 1024**3))
        z = size_child_count(cram, 31, memory_gb=244.0, bloom_filter=True)
        assert 1.5e10 < z["est_distinct"] < 4e10                       # errors dominate the distinct k-mers
        assert z["bloom_elements"] >= z["est_distinct"] and 20 < z["bloom_gb"] < 60
        assert z["hash_entries"] >= 2**32 and z["hash_gb"] + z["bloom_gb"] <= 244 * 0.6 + 1e-6
        assert not z["overflow"], z["note"]                            # the table holds what passes the filter
        assert "Bloom filter" in z["note"] and "budget 146 GB of 244" in z["note"]

    def test_without_the_filter_the_table_is_capped_by_the_budget_and_overflows_to_disk(self, tmp_path):
        cram = _fake_file(str(tmp_path / "child.cram"), 1)
        os.truncate(cram, int(41.5 * 1024**3))
        z = size_child_count(cram, 31, memory_gb=244.0, bloom_filter=False)
        assert z["bloom_elements"] == 0 and z["hash_gb"] <= 244 * 0.6 + 1e-6
        assert z["overflow"] and "dump to disk" in z["note"]           # what the old code could never say: it doubled instead

    def test_small_inputs_get_the_floors(self, tmp_path):
        bam = _fake_file(str(tmp_path / "tiny.bam"), 50_000)
        z = size_child_count(bam, 31, memory_gb=16.0, bloom_filter=True)
        assert z["hash_entries"] == _pow2_at_least(HASH_MIN_ENTRIES) and z["bloom_elements"] == BLOOM_MIN_ELEMENTS
        assert z["hash_gb"] < 1.0 and z["bloom_gb"] < 0.05

    def test_a_given_hash_size_is_honoured(self, tmp_path):
        bam = _fake_file(str(tmp_path / "tiny.bam"), 50_000)
        z = size_child_count(bam, 31, memory_gb=16.0, bloom_filter=True, jf_hash_size="2G")
        assert z["hash_entries"] == 2 * 2**30 and "(given)" in z["note"]
        assert _estimate_jf_hash_size(bam, 31, memory_gb=16.0) == str(_pow2_at_least(HASH_MIN_ENTRIES))

    def test_a_reference_fasta_is_error_free(self, tmp_path):
        # the reference index (_ensure_ref_jf) sizes its count with the same estimator: GRCh38's 3.2 GB of sequence hold
        # the genome's distinct k-mers and no sequencing errors - one 2^32 table (about 25 GB), as the driver's 64-GB
        # reference job allows; a read file of that size would be sized for its errors
        fa = _fake_file(str(tmp_path / "ref.fa"), 1)
        os.truncate(fa, int(3.2 * 1e9))
        z = size_child_count(fa, 31, memory_gb=64.0, bloom_filter=False)
        assert 3.0e9 <= z["est_distinct"] <= 3.1e9 and z["hash_entries"] == 2**32 and z["hash_gb"] < 30 and not z["overflow"]
        assert _estimate_jf_hash_size(fa, 31, memory_gb=64.0) == str(2**32)   # without a budget the system's RAM caps it (a 16-GB CI runner: 2^30)

    def test_a_missing_file_falls_back(self):
        assert _estimate_jf_hash_size("/nonexistent/child.bam", 31, memory_gb=16.0) == str(_pow2_at_least(HASH_MIN_ENTRIES))


class TestCommand:
    def test_disk_filter_and_lower_bound(self):
        cmd = _child_count_command("/t/child.jf", 31, 22, 8589934592, 30000000000, 2)
        assert cmd[:2] == ["jellyfish", "count"] and "--disk" in cmd and cmd[-1] == "/dev/fd/0"
        assert cmd[cmd.index("-s") + 1] == "8589934592" and cmd[cmd.index("--bf-size") + 1] == "30000000000"
        assert cmd[cmd.index("-L") + 1] == "2" and cmd[cmd.index("-o") + 1] == "/t/child.jf" and "-C" in cmd

    def test_without_the_filter(self):
        cmd = _child_count_command("/t/child.jf", 31, 4, 100000000, 0, 3)
        assert "--bf-size" not in cmd and "--disk" in cmd and cmd[cmd.index("-L") + 1] == "3"


@pytest.mark.skipif(not HAVE_JELLYFISH, reason="jellyfish not installed")
class TestJellyfishFacts:
    """The behaviours the fix rests on, against the installed jellyfish."""

    @pytest.fixture()
    def reads(self, tmp_path):
        # four distinct 7-mers at multiplicities 1, 2, 3 and 5
        seqs = {"CCGTAAT": 1, "GCCTTTC": 2, "CCTAACA": 3, "GAGTTTT": 5}
        fa = tmp_path / "reads.fa"
        with open(fa, "w") as fh:
            i = 0
            for s, n in seqs.items():
                for _ in range(n):
                    fh.write(">r%d\n%s\n" % (i, s))
                    i += 1
        return str(fa), seqs

    def _dump(self, jf):
        out = subprocess.run(["jellyfish", "dump", "-c", jf], capture_output=True, text=True, check=True).stdout
        return {line.split()[0]: int(line.split()[1]) for line in out.splitlines() if line.strip()}

    def _canon(self, s):
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        rc = "".join(comp[c] for c in reversed(s))
        return min(s, rc)

    def test_bloom_filter_drops_singletons_and_counts_one_low(self, tmp_path, reads):
        fa, seqs = reads
        plain, bf = str(tmp_path / "plain.jf"), str(tmp_path / "bf.jf")
        subprocess.run(["jellyfish", "count", "-m", "7", "-s", "1000", "-t", "1", "-C", "-o", plain, fa], check=True)
        subprocess.run(["jellyfish", "count", "-m", "7", "-s", "1000", "-t", "1", "-C", "--bf-size", "1000", "-o", bf, fa], check=True)
        truth = {self._canon(s): n for s, n in seqs.items()}
        assert self._dump(plain) == truth
        assert self._dump(bf) == {k: n - 1 for k, n in truth.items() if n >= 2}

    def test_disk_mode_dumps_merges_and_applies_the_lower_bound(self, tmp_path, reads):
        fa, seqs = reads
        out = str(tmp_path / "disk.jf")
        cmd = _child_count_command(out, 7, 1, 2, 0, 2)                  # a two-entry table: fills on every read
        cmd[-1] = fa                                                    # the file instead of stdin
        subprocess.run(cmd, check=True)
        assert sorted(os.listdir(tmp_path)) == ["disk.jf", "reads.fa"]  # the overflow files were merged and removed
        truth = {self._canon(s): n for s, n in seqs.items() if n >= 2}
        assert self._dump(out) == truth

    def test_filter_and_lowered_bound_recover_the_true_threshold(self, tmp_path, reads):
        fa, seqs = reads
        out = str(tmp_path / "bf3.jf")
        cmd = _child_count_command(out, 7, 1, 1000, 1000, 3 - 1)        # min count 3 through the filter: -L 2 on hash counts
        cmd[-1] = fa
        subprocess.run(cmd, check=True)
        got = self._dump(out)
        assert set(got) == {self._canon(s) for s, n in seqs.items() if n >= 3}
