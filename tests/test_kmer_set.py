"""The k-mer set behind discovery mode's read scan.

:mod:`kmer_denovo_filter.core.kmer_set` encodes k-mers as 2-bit codes and
scans batches of reads with numpy.  Its result must be exactly that of the
per-k-mer reference (:func:`kmer_utils._extract_read_kmers` and a set
lookup): the same reads, the same k-mer positions, the same distinct
k-mers, with N, lower case, short reads and read boundaries handled the
same way.
"""

import random

import numpy as np
import pytest

from kmer_denovo_filter.core.kmer_set import (
    MAX_CODED_K,
    KmerSet,
    _window_codes,
    _window_mask,
    _BASE_CODE,
    decode_kmer,
    encode_kmer,
)
from kmer_denovo_filter.kmer_utils import (
    _extract_read_kmers,
    canonicalize,
    reverse_complement,
)


def _random_seq(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def _reference_scan(reads, kmers, k):
    """What the per-k-mer reference finds: {read index: (positions, k-mers)}."""
    found = {}
    for index, read in enumerate(reads):
        canon_at_pos, _ = _extract_read_kmers(read, k)
        hits = {(pos, canon) for pos, canon in canon_at_pos.items()
                if canon in kmers}
        if hits:
            found[index] = ({pos for pos, _ in hits},
                            {canon for _, canon in hits})
    return found


def _as_strings(result, k):
    """The scan's result with the k-mers decoded to strings."""
    out = {}
    for index, positions, codes in result:
        if hasattr(codes, "tolist"):
            kmers = {decode_kmer(code, k) for code in codes.tolist()}
        else:
            kmers = set(codes)
        out[index] = (set(positions.tolist()), kmers)
    return out


class TestEncoding:
    """The 2-bit code orders k-mers as their strings do."""

    @pytest.mark.parametrize("k", [1, 5, 11, 21, 31, 32])
    def test_code_order_is_string_order(self, k):
        rng = random.Random(k)
        kmers = [_random_seq(rng, k) for _ in range(300)]
        codes = [encode_kmer(s) for s in kmers]
        assert [decode_kmer(c, k) for c in sorted(codes)] == sorted(kmers)

    @pytest.mark.parametrize("k", [5, 11, 31, 32])
    def test_canonical_code_is_min_of_both_strands(self, k):
        rng = random.Random(k)
        kmers = [_random_seq(rng, k) for _ in range(200)]
        buf = "".join(kmers).encode()
        fwd, rc = _window_codes(_BASE_CODE[np.frombuffer(buf, dtype=np.uint8)], k)
        firsts = np.arange(len(kmers)) * k
        for kmer, f, r in zip(kmers, fwd[firsts], rc[firsts]):
            assert int(f) == encode_kmer(kmer)
            assert int(r) == encode_kmer(reverse_complement(kmer))
            assert min(int(f), int(r)) == encode_kmer(canonicalize(kmer))

    def test_lower_case_and_other_letters(self):
        codes = _BASE_CODE[np.frombuffer(b"acgtACGTNnRX.-", dtype=np.uint8)]
        assert codes.tolist() == [0, 1, 2, 3, 0, 1, 2, 3, 4, 4, 4, 4, 4, 4]

    def test_window_mask_keeps_whole_kmers_of_one_read(self):
        reads = [b"ACGTACGT", b"ACGNACGTAC", b"AC"]
        buf = b"".join(reads)
        codes = _BASE_CODE[np.frombuffer(buf, dtype=np.uint8)]
        keep, read_of_window, starts = _window_mask(
            codes, 4, np.array([len(r) for r in reads]))
        kept = {(int(read_of_window[i]), int(i - starts[read_of_window[i]]))
                for i in np.flatnonzero(keep)}
        # read 0: every window; read 1: those after the N; read 2: none
        assert kept == ({(0, p) for p in range(5)}
                        | {(1, p) for p in range(4, 7)})


class TestKmerSet:
    def test_from_kmers_canonicalizes_and_deduplicates(self):
        kmers = ["AAACC", "GGTTT", "aaacc", "ACGTA"]
        s = KmerSet.from_kmers(kmers, 5)
        assert len(s) == 2
        assert s.nbytes == 16
        assert sorted(decode_kmer(int(c), 5) for c in s.codes) == [
            "AAACC", "ACGTA"]

    def test_wrong_length_is_rejected(self):
        with pytest.raises(ValueError, match="length 4, not 5"):
            KmerSet.from_kmers(["ACGT"], 5)
        with pytest.raises(ValueError, match="other than A, C, G and T"):
            KmerSet.from_kmers(["ACGTN"], 5)

    def test_empty(self):
        s = KmerSet.from_kmers([], 5)
        assert len(s) == 0
        assert s.scan([b"ACGTACGTACGT"]) == []

    def test_from_fasta_matches_from_kmers(self, tmp_path):
        rng = random.Random(5)
        kmers = [_random_seq(rng, 31) for _ in range(5000)]
        path = tmp_path / "kmers.fa"
        with open(path, "w") as fh:
            for i, kmer in enumerate(kmers):
                fh.write(f">{i}\n{kmer}\n")
            # a blank line, a repeat and a lower-case copy change nothing
            fh.write(f"\n>dup\n{kmers[0]}\n>lower\n{kmers[1].lower()}\n")
        # small chunks so that records straddle chunk boundaries
        loaded = KmerSet.from_fasta(str(path), 31, chunk_bytes=777)
        built = KmerSet.from_kmers(kmers, 31)
        assert np.array_equal(loaded.codes, built.codes)
        assert len(loaded) == len({canonicalize(k) for k in kmers})

    def test_from_fasta_rejects_another_k(self, tmp_path):
        path = tmp_path / "kmers.fa"
        path.write_text(">0\nACGTACGTACG\n")
        with pytest.raises(ValueError, match="length 11, not 31"):
            KmerSet.from_fasta(str(path), 31)

    def test_contains(self):
        s = KmerSet.from_kmers(["AAACC", "ACGTA"], 5)
        codes = np.array([encode_kmer("AAACC"), encode_kmer("ACGTA"),
                          encode_kmer("TTTTT"), (1 << 10) - 1],
                         dtype=np.uint64)
        assert s.contains(codes).tolist() == [True, True, False, False]

    def test_pickles(self):
        import pickle
        s = KmerSet.from_kmers(["AAACC", "ACGTA"], 5)
        t = pickle.loads(pickle.dumps(s))
        assert np.array_equal(s.codes, t.codes) and t.k == 5


class TestScanMatchesTheReference:
    """On random reads with Ns, lower case, every length around k and
    k-mers of the set planted in them."""

    @pytest.mark.parametrize("k", [5, 11, 21, 31, 32, 33])
    def test_random_reads(self, k):
        rng = random.Random(k)
        universe = [_random_seq(rng, k) for _ in range(1000)]
        chosen = {canonicalize(s) for s in rng.sample(universe, 200)}
        kmer_set = KmerSet.from_kmers(chosen, k)
        assert (kmer_set.codes is None) == (k > MAX_CODED_K)
        reads = []
        for _ in range(1500):
            length = rng.choice([k - 1, k, k + 1, 50, 100, 150, 151, 250])
            parts = []
            while sum(len(p) for p in parts) < length:
                if rng.random() < 0.7:
                    parts.append(rng.choice(universe))
                else:
                    parts.append(_random_seq(rng, rng.randint(1, 10), "ACGTN"))
            read = "".join(parts)[:length]
            if rng.random() < 0.2:
                read = read.lower()
            if rng.random() < 0.3:
                i = rng.randrange(len(read))
                read = read[:i] + "N" + read[i + 1:]
            reads.append(read)
        expected = _reference_scan(reads, chosen, k)
        assert expected, "the fixture plants hits"
        got = _as_strings(kmer_set.scan([r.encode() for r in reads]), k)
        assert got == expected

    def test_no_window_crosses_two_reads(self):
        # The last 4 bases of read 1 and the first of read 2 form a k-mer
        # of the set; laid end to end they must not match.
        k = 5
        kmer_set = KmerSet.from_kmers(["ACGTT"], k)
        assert kmer_set.scan([b"GGGGACGT", b"TGGGG"]) == []
        assert len(kmer_set.scan([b"GGGGACGTT", b"TGGGG"])) == 1

    def test_positions_and_distinct_kmers(self):
        k = 5
        kmer_set = KmerSet.from_kmers(["AAAAA", "ACGTA"], k)
        # AAAAAAA holds AAAAA at 0, 1, 2 (one distinct k-mer); TACGT at 3
        # is the reverse complement of ACGTA
        [(index, positions, codes)] = kmer_set.scan([b"AAAAAAATACGT"])
        assert index == 0
        assert positions.tolist() == [0, 1, 2, 7]
        assert [decode_kmer(c, k) for c in codes.tolist()] == [
            "AAAAA", "AAAAA", "AAAAA", "ACGTA"]
        assert len(set(codes.tolist())) == 2

    def test_reads_shorter_than_k_and_empty_batch(self):
        kmer_set = KmerSet.from_kmers(["AAAAA"], 5)
        assert kmer_set.scan([]) == []
        assert kmer_set.scan([b"AAAA", b""]) == []
        assert kmer_set.scan([b"AAAA", b"AAAAA"])[0][0] == 1
