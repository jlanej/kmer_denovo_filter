"""A set of canonical k-mers as sorted 2-bit codes, and the vectorized scan
of read batches against it.

This is the engine behind Module 3 of discovery mode.  The proband-unique
k-mers (a few million for a 30x genome) are loaded once, in the parent,
as a sorted ``numpy`` array of ``uint64`` codes: two bits per base, so a
31-mer takes 62 bits and the set about 8 bytes per k-mer.  Workers forked
from the parent share that array copy-on-write.  Each batch of reads is
encoded with ``numpy`` too: every k-mer window of every read is derived
by shifts, canonicalized as ``min(forward, reverse complement)``, masked
where it holds an N (or any base other than A, C, G and T), and tested
for membership with ``numpy.searchsorted`` against the sorted set.  The
transient arrays of a batch are a few tens of bytes per base, so a batch
of 1.5 million bases costs about 100 MB.

The encoding orders k-mers exactly as their strings do (A < C < G < T),
so the canonical code is the code of :func:`kmer_utils.canonicalize`.
The result of a scan is the same as the per-k-mer reference
implementation :func:`kmer_utils._extract_read_kmers` followed by a set
lookup, which the tests check against random reads.

Codes fit a ``uint64`` for k <= 32.  A larger k falls back to the Python
reference implementation (:class:`KmerSet` keeps the strings then), which
is correct but slow; discovery mode's default and the real runs use 31.
"""

import numpy as np

from kmer_denovo_filter.kmer_utils import _extract_read_kmers, canonicalize

#: The largest k whose 2-bit code fits a uint64.
MAX_CODED_K = 32

#: Byte value -> 2-bit base code; 4 marks N and anything else.
_BASE_CODE = np.full(256, 4, dtype=np.uint8)
for _base, _code in ((b"A", 0), (b"C", 1), (b"G", 2), (b"T", 3)):
    _BASE_CODE[_base[0]] = _code
    _BASE_CODE[_base.lower()[0]] = _code

_U64 = np.uint64


def encode_kmer(kmer):
    """The 2-bit code of one k-mer string (not canonicalized)."""
    code = 0
    for base in kmer.upper():
        code = (code << 2) | "ACGT".index(base)
    return code


def decode_kmer(code, k):
    """The k-mer string of a 2-bit code."""
    bases = []
    for _ in range(k):
        bases.append("ACGT"[code & 3])
        code >>= 2
    return "".join(reversed(bases))


def _powers_of_two(k):
    """The powers of two that sum to *k*, largest first."""
    return [1 << bit for bit in range(k.bit_length() - 1, -1, -1)
            if k & (1 << bit)]


def _window_codes(codes, k):
    """Forward and reverse-complement 2-bit codes of every k-window.

    *codes* is a ``uint8`` array of base codes (0-3, or 4 for a base that
    is not A, C, G or T) over a buffer of ``T`` bases; the result holds
    ``T - k + 1`` windows, window ``i`` covering ``codes[i:i + k]``.  A
    window holding a 4 gets an arbitrary code: callers mask those windows
    out (see :func:`_window_mask`).

    Built by doubling: codes of 2-windows from 1-windows, 4-windows from
    2-windows, and so on, then the parts of *k* (its powers of two) are
    combined - about ``2 * log2(k)`` passes over the buffer rather than
    *k*.
    """
    n = codes.shape[0] - k + 1
    if n <= 0:
        return np.empty(0, dtype=_U64), np.empty(0, dtype=_U64)
    parts = _powers_of_two(k)

    # Forward: the first base of a window is the most significant.
    by_width = {1: codes.astype(_U64)}
    width = 1
    while width * 2 <= parts[0]:
        prev = by_width[width]
        by_width[width * 2] = (prev[:-width] << _U64(2 * width)) | prev[width:]
        width *= 2
    fwd = np.zeros(n, dtype=_U64)
    offset = 0
    for part in parts:
        fwd |= by_width[part][offset:offset + n] << _U64(2 * (k - offset - part))
        offset += part
    del by_width

    # Reverse complement: complement every base, and the first base of a
    # window is the least significant.
    comp = (3 - codes.astype(np.int16)).clip(0, 3).astype(_U64)
    by_width = {1: comp}
    width = 1
    while width * 2 <= parts[0]:
        prev = by_width[width]
        by_width[width * 2] = prev[:-width] | (prev[width:] << _U64(2 * width))
        width *= 2
    rc = np.zeros(n, dtype=_U64)
    offset = 0
    for part in parts:
        rc |= by_width[part][offset:offset + n] << _U64(2 * offset)
        offset += part
    return fwd, rc


def _window_mask(codes, k, lengths):
    """Which k-windows of a concatenated buffer are real k-mers of one read.

    *codes* is the ``uint8`` base-code buffer of the reads laid end to
    end, *lengths* their lengths in order.  A window is kept when it lies
    within a single read and holds no base outside A, C, G and T.
    Returns the boolean mask over the ``T - k + 1`` windows, and the
    read index and read start of every window.
    """
    n = codes.shape[0] - k + 1
    starts = np.zeros(len(lengths), dtype=np.int64)
    if len(lengths) > 1:
        np.cumsum(lengths[:-1], out=starts[1:])
    ends = starts + lengths
    read_of_base = np.repeat(np.arange(len(lengths), dtype=np.int64), lengths)
    read_of_window = read_of_base[:n]
    within_read = (np.arange(n, dtype=np.int64) + k) <= ends[read_of_window]
    bad = np.cumsum(codes >= 4, dtype=np.int32)
    with_n = bad[k - 1:] - np.concatenate((np.zeros(1, dtype=np.int32), bad[:n - 1]))
    return within_read & (with_n == 0), read_of_window, starts


class KmerSet:
    """Canonical k-mers as a sorted array of 2-bit codes.

    ``KmerSet.from_fasta`` reads a FASTA of k-mers (one per record, as
    ``jellyfish dump`` and Module 2 write them); ``KmerSet.from_kmers``
    takes any iterable of k-mer strings.  Both canonicalize, so a
    non-canonical input k-mer matches reads carrying either strand.
    ``scan`` tests a batch of reads.

    For ``k > MAX_CODED_K`` the set holds the canonical strings instead
    and ``scan`` runs the per-k-mer reference implementation.
    """

    def __init__(self, k, codes=None, strings=None):
        self.k = int(k)
        if self.k > MAX_CODED_K:
            self.codes = None
            self.strings = frozenset(strings or ())
        else:
            if codes is None:
                codes = np.empty(0, dtype=_U64)
            self.codes = np.unique(np.asarray(codes, dtype=_U64))
            self.strings = None

    # ── construction ──────────────────────────────────────────────

    @classmethod
    def from_kmers(cls, kmers, k):
        """From an iterable of k-mer strings (each of length *k*)."""
        kmers = list(kmers)
        for kmer in kmers:
            if len(kmer) != k:
                raise ValueError(
                    f"k-mer {kmer!r} has length {len(kmer)}, not {k}")
        if k > MAX_CODED_K:
            return cls(k, strings=(canonicalize(s.upper()) for s in kmers))
        if not kmers:
            return cls(k)
        buf = "".join(kmers).encode("ascii")
        return cls(k, codes=_canonical_codes_of_records(buf, k))

    @classmethod
    def from_fasta(cls, path, k, chunk_bytes=32 << 20):
        """From a FASTA of k-mers: every sequence line is one k-mer of
        length *k* (header lines and blank lines are ignored).  Read in
        chunks, so the file is never held in memory at once."""
        if k > MAX_CODED_K:
            strings = []
            with open(path) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line and not line.startswith(">"):
                        if len(line) != k:
                            raise ValueError(
                                f"{path}: k-mer {line!r} has length "
                                f"{len(line)}, not {k}")
                        strings.append(canonicalize(line.upper()))
            return cls(k, strings=strings)
        pieces = []
        carry = b""
        with open(path, "rb") as fh:
            while True:
                block = fh.read(chunk_bytes)
                if not block:
                    break
                block = carry + block
                cut = block.rfind(b"\n")
                if cut < 0:
                    carry = block
                    continue
                carry = block[cut + 1:]
                pieces.append(_codes_of_fasta_block(block[:cut + 1], k, path))
        if carry:
            pieces.append(_codes_of_fasta_block(carry + b"\n", k, path))
        if not pieces:
            return cls(k)
        return cls(k, codes=np.concatenate(pieces))

    # ── size ──────────────────────────────────────────────────────

    def __len__(self):
        if self.codes is None:
            return len(self.strings)
        return int(self.codes.shape[0])

    @property
    def nbytes(self):
        """Bytes of the array (0 for the string fallback)."""
        return int(self.codes.nbytes) if self.codes is not None else 0

    # ── membership ────────────────────────────────────────────────

    def contains(self, canonical_codes):
        """Boolean array: which of *canonical_codes* are in the set."""
        canonical_codes = np.asarray(canonical_codes, dtype=_U64)
        if self.codes is None or self.codes.shape[0] == 0:
            return np.zeros(canonical_codes.shape[0], dtype=bool)
        idx = np.searchsorted(self.codes, canonical_codes)
        idx[idx == self.codes.shape[0]] = 0
        return self.codes[idx] == canonical_codes

    def scan(self, seqs):
        """Scan a batch of read sequences.

        *seqs* is a list of ``bytes`` (ASCII; lower case is read as upper
        case).  Returns a list of ``(index, positions, kmers)`` for the
        reads holding at least one k-mer of the set, in batch order:
        *positions* is the sorted array of 0-based start positions in the
        read of the k-mer windows found, and *kmers* the canonical code of
        each (strings for the fallback; equal codes mean one k-mer).  A
        window holding an N (or any base other than A, C, G and T) is
        never a hit; a read shorter than k has no window.
        """
        k = self.k
        if self.codes is None:
            return self._scan_strings(seqs)
        lengths = np.fromiter((len(s) for s in seqs), dtype=np.int64,
                              count=len(seqs))
        if lengths.sum() < k or len(self) == 0:
            return []
        buf = b"".join(seqs)
        codes = _BASE_CODE[np.frombuffer(buf, dtype=np.uint8)]
        keep, read_of_window, starts = _window_mask(codes, k, lengths)
        selected = np.flatnonzero(keep)
        if selected.shape[0] == 0:
            return []
        fwd, rc = _window_codes(codes, k)
        canon = np.minimum(fwd[selected], rc[selected])
        del fwd, rc
        hit = self.contains(canon)
        hit_windows = selected[hit]
        if hit_windows.shape[0] == 0:
            return []
        hit_codes = canon[hit]
        hit_reads = read_of_window[hit_windows]
        positions = hit_windows - starts[hit_reads]
        # hit_windows ascend, so hit_reads is non-decreasing: cut at the
        # changes of read index.
        bounds = np.flatnonzero(np.diff(hit_reads)) + 1
        bounds = np.concatenate(([0], bounds, [hit_reads.shape[0]]))
        return [
            (int(hit_reads[lo]), positions[lo:hi], hit_codes[lo:hi])
            for lo, hi in zip(bounds[:-1], bounds[1:])
        ]

    def _scan_strings(self, seqs):
        out = []
        for index, seq in enumerate(seqs):
            canon_at_pos, _ = _extract_read_kmers(
                seq.decode("ascii"), self.k)
            found = [(pos, canon) for pos, canon in canon_at_pos.items()
                     if canon in self.strings]
            if found:
                out.append((index, np.array([p for p, _ in found]),
                            [c for _, c in found]))
        return out


def _canonical_codes_of_columns(column, n, k):
    """Canonical codes of *n* k-mers given column by column: ``column(j)``
    returns the ``uint8`` base codes of base *j* of every k-mer.  Works
    on one column at a time, so the transient memory is a few arrays of
    *n* elements, never an ``(n, k)`` table."""
    fwd = np.zeros(n, dtype=_U64)
    rc = np.zeros(n, dtype=_U64)
    bad = np.zeros(n, dtype=bool)
    for j in range(k):
        head = column(j)
        bad |= head >= 4
        fwd <<= _U64(2)
        fwd |= head
        tail = column(k - 1 - j)
        rc <<= _U64(2)
        rc |= (tail ^ 3) & 3          # the complement of a base code
    if bad.any():
        row = int(np.flatnonzero(bad)[0])
        letters = "ACGTN"
        kmer = "".join(letters[min(int(column(j)[row]), 4)] for j in range(k))
        raise ValueError(
            f"k-mer {kmer!r} holds a base other than A, C, G and T")
    return np.minimum(fwd, rc)


def _canonical_codes_of_records(buf, k):
    """Canonical codes of the k-mers laid end to end in *buf* (ASCII)."""
    codes = _BASE_CODE[np.frombuffer(buf, dtype=np.uint8)]
    n_records = codes.shape[0] // k
    if codes.shape[0] != n_records * k:
        raise ValueError("k-mer buffer is not a multiple of k")
    rows = codes.reshape(n_records, k)
    return _canonical_codes_of_columns(lambda j: rows[:, j], n_records, k)


def _codes_of_fasta_block(block, k, path):
    """Canonical codes of the sequence lines of a FASTA block that ends
    with a newline: every line not starting with ``>`` is one k-mer."""
    arr = np.frombuffer(block, dtype=np.uint8)
    newlines = np.flatnonzero(arr == 10)
    if newlines.shape[0] == 0:
        return np.empty(0, dtype=_U64)
    starts = np.concatenate((np.zeros(1, dtype=np.int64), newlines[:-1] + 1))
    lengths = newlines - starts
    is_seq = (lengths > 0) & (arr[np.minimum(starts, arr.shape[0] - 1)] != ord(">"))
    starts, lengths = starts[is_seq], lengths[is_seq]
    if starts.shape[0] == 0:
        return np.empty(0, dtype=_U64)
    if np.any(lengths != k):
        bad = int(np.flatnonzero(lengths != k)[0])
        line = bytes(arr[starts[bad]:starts[bad] + lengths[bad]])
        raise ValueError(
            f"{path}: k-mer {line.decode('ascii', 'replace')!r} has "
            f"length {lengths[bad]}, not {k}")
    try:
        return _canonical_codes_of_columns(
            lambda j: _BASE_CODE[arr[starts + j]], starts.shape[0], k)
    except ValueError as exc:
        raise ValueError(f"{path}: {exc}") from None
