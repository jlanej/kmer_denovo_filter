"""Module 1b's scratch and its failures.

The reference subtraction writes the non-reference parts beside the child's
index on the tmp filesystem.  Before the dump starts, the run reads the
k-mer count from the index's header, estimates the parts (a tenth of the
dump's text) and stops with the numbers when the free space is short; while
the dump streams, a ``[Disk]`` line samples the filesystem.  When a part
runs out of space anyway, the error names the process that failed first and
on its own - awk, with its "No space left on device" - and never the dump
or a query the cleanup stopped.
"""

import logging
import os
import random
import stat
import subprocess
import time

import pytest

import kmer_denovo_filter.discovery.pipeline as pipeline_mod
from kmer_denovo_filter.discovery.pipeline import (
    _NON_REF_PARTS_FRACTION,
    _ProcessEnd,
    _SUBTRACT_STAGE_DUMP,
    _SUBTRACT_STAGE_FILTER,
    _SUBTRACT_STAGE_QUERY,
    _attribute_subtraction_failure,
    _ended_on_their_own,
    _estimate_subtraction_scratch,
    _jf_index_records,
    _subtract_reference_kmers,
    _subtraction_failure_message,
)

LOGGER = "kmer_denovo_filter.discovery.pipeline"


def _random_sequence(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _count(fasta, jf, k):
    subprocess.run(["jellyfish", "count", "-m", str(k), "-s", "4M", "-t", "1",
                    "-C", "-o", jf, str(fasta)], check=True)


def _distinct(jf):
    stats = subprocess.run(["jellyfish", "stats", jf], capture_output=True,
                           text=True, check=True).stdout
    return int(next(line.split()[1] for line in stats.splitlines()
                    if line.startswith("Distinct")))


def _trio_indexes(tmp_path, k=21, ref_len=20_000, novel_len=30_000, seed=3):
    """A reference index and a child index (the reference plus novel
    sequence, so most of the child's k-mers are non-reference)."""
    ref = _random_sequence(ref_len, seed)
    novel = _random_sequence(novel_len, seed + 1)
    (tmp_path / "ref.fa").write_text(f">chr1\n{ref}\n")
    (tmp_path / "child.fa").write_text(f">a\n{ref}\n>b\n{novel}\n")
    ref_jf = str(tmp_path / "ref.jf")
    child_jf = str(tmp_path / "child.jf")
    _count(tmp_path / "ref.fa", ref_jf, k)
    _count(tmp_path / "child.fa", child_jf, k)
    return ref_jf, child_jf


class TestIndexHeader:
    """The k-mer count read from a jellyfish index's header and size."""

    @pytest.mark.parametrize("k,record_bytes", [(31, 12), (21, 10), (33, 13)])
    def test_the_header_gives_the_kmer_count(self, tmp_path, k, record_bytes):
        (tmp_path / "ref.fa").write_text(f">chr1\n{_random_sequence(50_000, k)}\n")
        jf = str(tmp_path / "ref.jf")
        _count(tmp_path / "ref.fa", jf, k)
        assert _jf_index_records(jf) == (_distinct(jf), record_bytes, k)

    def test_what_is_not_an_index(self, tmp_path):
        text = tmp_path / "dump.fa"
        text.write_text(">1\nACGTACGTACG\n" * 100)
        assert _jf_index_records(str(text)) is None
        empty = tmp_path / "empty.jf"
        empty.write_bytes(b"")
        assert _jf_index_records(str(empty)) is None
        other_format = tmp_path / "text.jf"
        header = b'{"format":"text/sorted","key_len":62,"counter_len":4}'
        other_format.write_bytes(f"{len(header):09d}".encode() + header + b"x" * 120)
        assert _jf_index_records(str(other_format)) is None
        assert _jf_index_records(str(tmp_path / "missing.jf")) is None


class TestEstimate:
    """The scratch the subtraction will need, from the index."""

    def test_a_30x_genome(self, tmp_path, monkeypatch):
        """2.6 G k-mers at k = 31 in a 29 GB index: parts of about 9 GB
        (the measured 8 GB with the margin), and that is what must be
        free, the index being far bigger than the parts."""
        child_jf = str(tmp_path / "child.jf")
        n_kmers = 2_600_000_000
        index_bytes = n_kmers * 12
        monkeypatch.setattr(pipeline_mod, "_jf_index_records",
                            lambda path: (n_kmers, 12, 31))
        monkeypatch.setattr(pipeline_mod, "_child_index_files",
                            lambda tmpdir, jf: [jf])
        monkeypatch.setattr(pipeline_mod, "_file_bytes", lambda path: index_bytes)
        est = _estimate_subtraction_scratch(child_jf, str(tmp_path))
        assert est["from_header"] and est["n_kmers"] == n_kmers and est["k"] == 31
        assert est["index_bytes"] == index_bytes
        assert est["dump_text_bytes"] == n_kmers * 36
        assert est["parts_bytes"] == int(n_kmers * 36 * _NON_REF_PARTS_FRACTION)
        assert 8.0e9 < est["parts_bytes"] < 10.0e9
        assert est["needed_bytes"] == est["parts_bytes"]
        assert 26.0 < index_bytes / 1024**3 < 30.0, "the index of the measured runs"

    def test_the_size_alone_when_the_header_cannot_be_read(self, tmp_path, monkeypatch):
        _ref_jf, child_jf = _trio_indexes(tmp_path)
        from_header = _estimate_subtraction_scratch(child_jf, str(tmp_path))
        assert from_header["from_header"] and from_header["k"] == 21
        monkeypatch.setattr(pipeline_mod, "_jf_index_records", lambda path: None)
        assert _estimate_subtraction_scratch(child_jf, str(tmp_path)) is None
        fallback = _estimate_subtraction_scratch(child_jf, str(tmp_path), kmer_size=21)
        assert not fallback["from_header"] and fallback["k"] == 21
        # 10 bytes per record at k = 21; the header's bytes are the only error
        assert 0 <= fallback["n_kmers"] - from_header["n_kmers"] <= 400

    def test_chunk_files_count_as_index(self, tmp_path):
        _ref_jf, child_jf = _trio_indexes(tmp_path)
        alone = _estimate_subtraction_scratch(child_jf, str(tmp_path))
        chunk = tmp_path / "child.jf_0"
        chunk.write_bytes(b"\0" * 5000)
        with_chunk = _estimate_subtraction_scratch(child_jf, str(tmp_path))
        assert with_chunk["index_bytes"] == alone["index_bytes"] + 5000
        assert with_chunk["n_kmers"] == alone["n_kmers"]


class TestSpaceCheck:
    """Before the dump starts."""

    def test_too_little_free_space_stops_before_anything_starts(self, tmp_path, monkeypatch):
        ref_jf, child_jf = _trio_indexes(tmp_path)
        est = _estimate_subtraction_scratch(child_jf, str(tmp_path))
        assert est["needed_bytes"] > 0
        monkeypatch.setattr(pipeline_mod, "_tmp_space",
                            lambda path: (est["needed_bytes"] // 2, 10 * 1024**3))
        run_dir = tmp_path / "run"
        run_dir.mkdir()
        with pytest.raises(RuntimeError) as exc:
            _subtract_reference_kmers(ref_jf, child_jf, 1, str(run_dir), threads=2,
                                      kmer_size=21)
        msg = str(exc.value)
        assert msg.startswith("Module 1b cannot start: the tmp filesystem holding")
        for phrase in (
            "free of 10.0 GB",
            "needs about",
            "non-reference parts, estimated at",
            f"{est['n_kmers']:,} dumped k-mers make at k=21",
            "beside the", "child index",
            "concatenation takes as much again once the index is removed",
            "Nothing of the subtraction has started",
            "scratch of at least",
            "about 40 GB for a 30x genome at k=31",
            "--tmp-dir at a larger filesystem",
        ):
            assert phrase in msg, phrase
        assert os.path.exists(child_jf), "the index is left as it was"
        assert not os.path.exists(run_dir / "non_ref_parts")
        assert not os.path.exists(run_dir / "child_non_ref_kmers.fa")

    def test_enough_space_logs_the_estimate_and_samples_the_disk(self, tmp_path, caplog):
        ref_jf, child_jf = _trio_indexes(tmp_path)
        caplog.set_level(logging.INFO)  # the [Disk] lines come from core.memory_utils
        _non_ref_fa, n_candidates, n_non_ref = _subtract_reference_kmers(
            ref_jf, child_jf, 1, str(tmp_path), threads=2, kmer_size=21,
            poll_interval=0,
        )
        assert n_candidates > n_non_ref > 20_000
        scratch = [r.getMessage() for r in caplog.records
                   if r.getMessage().startswith("Module 1b scratch:")]
        assert len(scratch) == 1
        assert "k-mers at k=21) stays on the tmp filesystem" in scratch[0]
        assert "10% of the" in scratch[0] and "free on" in scratch[0]
        assert "[Disk] tmpdir during reference subtraction" in caplog.text
        assert "[TmpDir] non-reference parts so far" in caplog.text
        assert not any(r.levelno >= logging.WARNING for r in caplog.records)

    def test_tight_space_warns_and_runs(self, tmp_path, monkeypatch, caplog):
        ref_jf, child_jf = _trio_indexes(tmp_path)
        est = _estimate_subtraction_scratch(child_jf, str(tmp_path))
        monkeypatch.setattr(pipeline_mod, "_tmp_space",
                            lambda path: (int(est["needed_bytes"] * 1.2), 10 * 1024**3))
        caplog.set_level(logging.INFO, logger=LOGGER)
        _non_ref_fa, n_candidates, _n_non_ref = _subtract_reference_kmers(
            ref_jf, child_jf, 1, str(tmp_path), threads=2, kmer_size=21)
        assert n_candidates > 0
        warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1 and warnings[0].startswith("Module 1b scratch is tight:")

    def test_no_header_and_no_k_skips_the_check(self, tmp_path, monkeypatch, caplog):
        ref_jf, child_jf = _trio_indexes(tmp_path)
        monkeypatch.setattr(pipeline_mod, "_jf_index_records", lambda path: None)
        monkeypatch.setattr(pipeline_mod, "_tmp_space", lambda path: (0, 10 * 1024**3))
        caplog.set_level(logging.INFO, logger=LOGGER)
        _non_ref_fa, n_candidates, _n_non_ref = _subtract_reference_kmers(
            ref_jf, child_jf, 1, str(tmp_path), threads=2)
        assert n_candidates > 0
        assert "could not size the scratch" in caplog.text


def _fake_awk(bin_dir, failing_pipe=0, after_records=40):
    """An ``awk`` whose part runs out of space after a few records on one
    pipe; on the others it reads to the end and exits 0."""
    bin_dir.mkdir()
    script = bin_dir / "awk"
    script.write_text(
        "#!/bin/sh\n"
        f"case \"$*\" in *p={failing_pipe}\\ *|*p={failing_pipe}) fail=1;; *) fail=0;; esac\n"
        "if [ \"$fail\" = 1 ]; then\n"
        "  n=0\n"
        "  while IFS= read -r line; do\n"
        "    n=$((n+1))\n"
        f"    if [ \"$n\" -ge {after_records} ]; then\n"
        "      echo 'awk: write failure (No space left on device)' >&2\n"
        "      exit 2\n"
        "    fi\n"
        "  done\n"
        "  exit 0\n"
        "fi\n"
        "cat >/dev/null\n"
        "exit 0\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return bin_dir


class TestFailureAttribution:
    """Who is named when the subtraction fails."""

    def test_a_part_out_of_space_names_awk_not_the_dump(self, tmp_path, monkeypatch):
        """A child whose dump is far bigger than the pipes hold, so the
        splitter is still streaming when pipe 0's awk hits ENOSPC: the
        query behind it dies of SIGPIPE, the splitter's write breaks, and
        the cleanup stops the dump and the other pipes.  The message
        names awk's write error, never the dump."""
        ref_jf, child_jf = _trio_indexes(tmp_path, novel_len=600_000)
        monkeypatch.setenv(
            "PATH", str(_fake_awk(tmp_path / "bin")) + os.pathsep + os.environ["PATH"])
        run_dir = tmp_path / "run"
        run_dir.mkdir()
        started = time.monotonic()
        with pytest.raises(RuntimeError) as exc:
            _subtract_reference_kmers(ref_jf, child_jf, 1, str(run_dir), threads=3,
                                      kmer_size=21)
        assert time.monotonic() - started < 30
        msg = str(exc.value)
        lines = msg.splitlines()
        assert lines[0] == (
            "reference subtraction failed: awk (ref subtraction, pipe 0) failed with "
            "exit code 2 - the first failure, on its own:")
        assert lines[1] == "awk: write failure (No space left on device)"
        assert "jellyfish dump (child) failed" not in msg
        stopped = next(line for line in lines
                       if line.startswith("Stopped by the cleanup, not the cause: "))
        assert "jellyfish dump (child) (-" in stopped
        # the query of pipe 0, left without a reader, died of SIGPIPE before
        # the splitter noticed: it is listed, not named
        assert [line for line in lines if line.startswith("Also ended on their own: ")] == [
            "Also ended on their own: jellyfish query (ref subtraction, pipe 0) (-13 (SIGPIPE))"]
        assert "The scratch ran out:" in msg
        assert "the subtraction needs about" in msg and "child index" in msg
        assert "point --tmp-dir at a larger filesystem" in msg
        assert "non-reference parts in" in msg and "non_ref_parts hold" in msg
        assert isinstance(exc.value.__cause__, BrokenPipeError)

    def test_a_failing_query_is_still_named(self, tmp_path):
        """The reference index is missing: every query exits 1 with its
        message, the dump is stopped by the cleanup; a query is named."""
        _ref_jf, child_jf = _trio_indexes(tmp_path)
        with pytest.raises(RuntimeError) as exc:
            _subtract_reference_kmers(str(tmp_path / "missing.jf"), child_jf, 1,
                                      str(tmp_path), threads=2, kmer_size=21)
        msg = str(exc.value)
        assert msg.startswith(
            "reference subtraction failed: jellyfish query (ref subtraction, pipe 0) "
            "failed with exit code 1 - the first failure, on its own:")
        assert "jellyfish dump (child) failed" not in msg
        assert "[jellyfish query (ref subtraction, pipe 1)] stderr:" in msg

    def test_the_ranking(self):
        D, Q, F = _SUBTRACT_STAGE_DUMP, _SUBTRACT_STAGE_QUERY, _SUBTRACT_STAGE_FILTER
        enospc = "awk: write failure (No space left on device)"
        ends = [
            _ProcessEnd("jellyfish dump (child)", -13, "", False, D),
            _ProcessEnd("jellyfish query (ref subtraction, pipe 0)", -13, "", True, Q),
            _ProcessEnd("awk (ref subtraction, pipe 0)", 2, enospc, True, F),
            _ProcessEnd("jellyfish query (ref subtraction, pipe 1)", -9, "", False, Q),
            _ProcessEnd("awk (ref subtraction, pipe 1)", 0, "", False, F),
        ]
        cause, also, stopped = _attribute_subtraction_failure(ends)
        assert cause.label == "awk (ref subtraction, pipe 0)"
        assert [e.label for e in also] == ["jellyfish query (ref subtraction, pipe 0)"]
        assert [e.label for e in stopped] == [
            "jellyfish dump (child)", "jellyfish query (ref subtraction, pipe 1)"]

    def test_a_process_the_cleanup_stopped_is_never_the_cause(self):
        D, Q, F = _SUBTRACT_STAGE_DUMP, _SUBTRACT_STAGE_QUERY, _SUBTRACT_STAGE_FILTER
        ends = [
            _ProcessEnd("jellyfish dump (child)", -9, "", False, D),
            _ProcessEnd("jellyfish query (ref subtraction, pipe 0)", -9, "", False, Q),
            _ProcessEnd("awk (ref subtraction, pipe 0)", -9, "some noise", False, F),
        ]
        cause, also, stopped = _attribute_subtraction_failure(ends)
        assert cause is None and also == [] and len(stopped) == 3
        msg = _subtraction_failure_message(ends, "/nonexistent/tmp", "/nonexistent/tmp/parts")
        assert msg.startswith("reference subtraction failed: a query pipe closed early "
                              "although no process failed on its own")
        assert "Stopped by the cleanup, not the cause: jellyfish dump (child) (-9 (SIGKILL))" in msg
        assert "[awk (ref subtraction, pipe 0)] stderr:\nsome noise" in msg

    @pytest.mark.parametrize("ends,named", [
        # an error exit before a death by signal
        ([_ProcessEnd("jellyfish dump (child)", 1, "dump: bad index", True, 0),
          _ProcessEnd("jellyfish query (ref subtraction, pipe 0)", -13, "", True, 1)],
         "jellyfish dump (child)"),
        # stderr before silence
        ([_ProcessEnd("jellyfish query (ref subtraction, pipe 0)", 1, "", True, 1),
          _ProcessEnd("jellyfish query (ref subtraction, pipe 1)", 1, "query: boom", True, 1)],
         "jellyfish query (ref subtraction, pipe 1)"),
        # the reader of a pipe before its writer
        ([_ProcessEnd("jellyfish query (ref subtraction, pipe 0)", 1, "query: write error", True, 1),
          _ProcessEnd("awk (ref subtraction, pipe 0)", 2, "awk: write failure", True, 2)],
         "awk (ref subtraction, pipe 0)"),
        # then pipeline order
        ([_ProcessEnd("awk (ref subtraction, pipe 0)", 2, "awk: write failure", True, 2),
          _ProcessEnd("awk (ref subtraction, pipe 1)", 2, "awk: write failure", True, 2)],
         "awk (ref subtraction, pipe 0)"),
    ])
    def test_who_failed_first(self, ends, named):
        cause, _also, _stopped = _attribute_subtraction_failure(ends)
        assert cause.label == named

    def test_the_message_carries_every_stderr_tail_and_the_disk(self, tmp_path, monkeypatch):
        D, Q, F = _SUBTRACT_STAGE_DUMP, _SUBTRACT_STAGE_QUERY, _SUBTRACT_STAGE_FILTER
        parts_dir = tmp_path / "non_ref_parts"
        parts_dir.mkdir()
        (parts_dir / "part_000.fa").write_bytes(b"x" * 3000)
        (parts_dir / "part_001.fa").write_bytes(b"x" * 1000)
        ends = [
            _ProcessEnd("jellyfish dump (child)", -13, "", False, D),
            _ProcessEnd("jellyfish query (ref subtraction, pipe 0)", -13, "", True, Q),
            _ProcessEnd("awk (ref subtraction, pipe 0)", 2,
                        "awk: write failure (No space left on device)", True, F),
            _ProcessEnd("jellyfish query (ref subtraction, pipe 1)", -9, "", False, Q),
            _ProcessEnd("awk (ref subtraction, pipe 1)", 2,
                        "awk: write failure (No space left on device)", True, F),
        ]
        monkeypatch.setattr(pipeline_mod, "_tmp_space", lambda path: (2 * 1024**3, 31 * 1024**3))
        estimate = {"parts_bytes": 9 * 1024**3, "index_bytes": 29 * 1024**3}
        msg = _subtraction_failure_message(ends, str(tmp_path), str(parts_dir), estimate)
        assert msg.splitlines()[0].startswith(
            "reference subtraction failed: awk (ref subtraction, pipe 0) failed with exit code 2")
        assert ("Also ended on their own: jellyfish query (ref subtraction, pipe 0) "
                "(-13 (SIGPIPE)), awk (ref subtraction, pipe 1) (2)") in msg
        assert ("Stopped by the cleanup, not the cause: jellyfish dump (child) "
                "(-13 (SIGPIPE)), jellyfish query (ref subtraction, pipe 1) (-9 (SIGKILL))") in msg
        assert "[awk (ref subtraction, pipe 1)] stderr:\nawk: write failure" in msg
        assert msg.count("No space left on device") == 2
        assert (f"The tmp filesystem holding {tmp_path} has 2.0 GB free of 31.0 GB; the "
                f"non-reference parts in {parts_dir} hold 0.0 MB.") in msg
        assert ("The scratch ran out: the subtraction needs about 9.0 GB for the parts "
                "beside the 29.0 GB child index, then as much again for their "
                "concatenation. Give the job a bigger scratch (about 40 GB for a 30x "
                "genome at k=31), or point --tmp-dir at a larger filesystem.") in msg

    def test_nearly_nothing_free_gives_the_remedy_without_the_words(self, tmp_path, monkeypatch):
        ends = [_ProcessEnd("awk (ref subtraction, pipe 0)", 2, "awk: write failure", True,
                            _SUBTRACT_STAGE_FILTER)]
        monkeypatch.setattr(pipeline_mod, "_tmp_space", lambda path: (10 * 1024**2, 31 * 1024**3))
        msg = _subtraction_failure_message(ends, str(tmp_path), str(tmp_path / "parts"))
        assert "The scratch ran out: Give the job a bigger scratch" in msg
        monkeypatch.setattr(pipeline_mod, "_tmp_space", lambda path: (20 * 1024**3, 31 * 1024**3))
        msg = _subtraction_failure_message(ends, str(tmp_path), str(tmp_path / "parts"))
        assert "has 20.0 GB free of 31.0 GB" in msg and "The scratch ran out" not in msg
        monkeypatch.setattr(pipeline_mod, "_tmp_space", lambda path: None)
        msg = _subtraction_failure_message(ends, str(tmp_path), str(tmp_path / "parts"))
        assert "tmp filesystem" not in msg

    def test_ended_on_their_own(self):
        gone = subprocess.Popen(["sh", "-c", "exit 3"])
        staying = subprocess.Popen(["sleep", "30"])
        try:
            started = time.monotonic()
            ended = _ended_on_their_own([gone, staying], grace_s=0.3)
            assert 0.25 < time.monotonic() - started < 5
            assert ended == {id(gone)} and gone.returncode == 3
            assert staying.poll() is None
        finally:
            staying.kill()
            staying.wait()
        both = subprocess.Popen(["sh", "-c", "exit 0"])
        both.wait()
        started = time.monotonic()
        assert _ended_on_their_own([both, gone], grace_s=5) == {id(both), id(gone)}
        assert time.monotonic() - started < 1, "every process ended: no waiting"
