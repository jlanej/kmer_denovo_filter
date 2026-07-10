"""BAM-aware entry points: classify reads, and classify per-variant allele support.

This module is the convenient, high-level surface most callers want.  It pulls
read sequences out of a BAM/CRAM and hands them to the pysam-free
:class:`~nonhuman_screen.engine.Kraken2Runner`, so no consumer has to
re-implement the "which reads support this allele, and how non-human are they"
machinery.

Requires the ``[bam]`` extra (pysam).  Import errors here are deliberately
actionable.

Coordinate convention: every ``pos`` is **0-based** (pysam reference
coordinates), and variant keys are ``"{chrom}:{pos}:{ref}:{alt}"`` to match the
internal conventions of consuming pipelines.
"""

from __future__ import annotations

try:
    import pysam
except ImportError as exc:  # pragma: no cover - import-guard
    raise ImportError(
        "nonhuman_screen.bam requires pysam. Install with: "
        "pip install 'nonhuman-screen[bam]'"
    ) from exc

from nonhuman_screen.alleles import read_supports_alt
from nonhuman_screen.engine import Kraken2Runner
from nonhuman_screen.result import TaxonomicFractions, VariantNHF


def _as_variant(v):
    """Normalize a variant into a ``(chrom, pos, ref, alt)`` tuple (pos 0-based)."""
    if isinstance(v, (tuple, list)):
        chrom, pos, ref, alt = v
        return (chrom, int(pos), ref, alt)
    if isinstance(v, dict):
        return (v["chrom"], int(v["pos"]), v["ref"], v["alt"])
    return (v.chrom, int(v.pos), v.ref, v.alt)


def _open(bam_path, ref_fasta):
    return pysam.AlignmentFile(
        bam_path, reference_filename=ref_fasta if ref_fasta else None
    )


def _alt_support_seqs(bam, chrom, pos, ref, alt, *, min_baseq=0):
    """Return ``{read_name: sequence}`` for reads supporting *alt* at 0-based *pos*.

    Fragments are de-duplicated by read name (first supporting mate wins), so a
    paired-end fragment contributes at most one sequence — consistent with how
    fragment-level fractions are computed downstream.
    """
    out = {}
    span = pos + max(len(ref) if ref else 1, 1)
    for read in bam.fetch(chrom, pos, span):
        name = read.query_name
        if name in out:
            continue
        seq = read.query_sequence
        if not seq:
            continue
        if read_supports_alt(read, pos, ref, alt, min_baseq):
            out[name] = seq
    return out


def reads_supporting_alt(
    bam_path, chrom, pos, ref, alt, *, ref_fasta=None, min_baseq=0
):
    """Return the set of read names supporting *alt* at 0-based *pos*.

    A thin, classification-free helper for callers that only need the
    allele-support read set (e.g. to feed their own logic).
    """
    with _open(bam_path, ref_fasta) as bam:
        return set(_alt_support_seqs(bam, chrom, pos, ref, alt, min_baseq=min_baseq))


def classify_reads_from_bam(
    bam_path,
    db_path,
    *,
    read_names=None,
    loci=None,
    ref_fasta=None,
    confidence=0.0,
    threads=1,
    memory_mapping=False,
    tmpdir=None,
):
    """Classify a set of reads from a BAM/CRAM and return a ``ClassificationResult``.

    Args:
        read_names: Set of read names to classify.  When *loci* is not given,
            the whole file is scanned for these names.
        loci: Optional ``{(chrom, pos): {read_name, ...}}`` mapping (0-based
            ``pos``) to restrict fetching to specific loci — far cheaper than a
            whole-file scan.  When both are given, *read_names* acts as an
            additional filter on the fetched reads.
    """
    if not read_names and not loci:
        return Kraken2Runner.Result()

    sequences = {}
    with _open(bam_path, ref_fasta) as bam:
        if loci:
            name_filter = set(read_names) if read_names else None
            for (chrom, pos), names in sorted(loci.items()):
                targets = set(names)
                if name_filter is not None:
                    targets &= name_filter
                if not targets:
                    continue
                for read in bam.fetch(chrom, pos, pos + 1):
                    n = read.query_name
                    if n in targets and read.query_sequence and n not in sequences:
                        sequences[n] = read.query_sequence
        else:
            wanted = set(read_names)
            for read in bam.fetch(until_eof=True):
                n = read.query_name
                if n in wanted and read.query_sequence and n not in sequences:
                    sequences[n] = read.query_sequence

    if not sequences:
        return Kraken2Runner.Result()

    runner = Kraken2Runner(
        db_path, confidence=confidence, threads=threads,
        memory_mapping=memory_mapping,
    )
    return runner.classify_sequences(sequences, tmpdir=tmpdir)


def classify_variants_alt_reads(
    bam_path,
    db_path,
    variants,
    *,
    ref_fasta=None,
    min_baseq=0,
    confidence=0.0,
    threads=1,
    memory_mapping=False,
    tmpdir=None,
):
    """Compute the allele-based non-human fraction for many variants at once.

    For each variant, the reads supporting its ALT allele are gathered; the
    union of those reads is classified with **a single kraken2 invocation**
    (kraken2's database load dominates runtime, so batching matters), then the
    result is split back per variant by read-name intersection.

    Args:
        variants: Iterable of variants, each a ``(chrom, pos, ref, alt)`` tuple
            (0-based ``pos``), a mapping with those keys, or any object with
            ``chrom``/``pos``/``ref``/``alt`` attributes.

    Returns:
        A list of :class:`~nonhuman_screen.result.VariantNHF`, in input order.
    """
    normalized = [_as_variant(v) for v in variants]

    sequences = {}
    per_variant_names = []
    with _open(bam_path, ref_fasta) as bam:
        for chrom, pos, ref, alt in normalized:
            vs = _alt_support_seqs(bam, chrom, pos, ref, alt, min_baseq=min_baseq)
            per_variant_names.append(set(vs))
            for name, seq in vs.items():
                sequences.setdefault(name, seq)

    runner = Kraken2Runner(
        db_path, confidence=confidence, threads=threads,
        memory_mapping=memory_mapping,
    )
    result = runner.classify_sequences(sequences, tmpdir=tmpdir) if sequences \
        else Kraken2Runner.Result()

    out = []
    for (chrom, pos, ref, alt), names in zip(normalized, per_variant_names):
        out.append(
            VariantNHF(
                chrom=chrom, pos=pos, ref=ref, alt=alt,
                supporting_read_names=frozenset(names),
                fractions=TaxonomicFractions.over_reads(result, names),
            )
        )
    return out


def classify_variant_alt_reads(
    bam_path, db_path, chrom, pos, ref, alt, *,
    ref_fasta=None, min_baseq=0, confidence=0.0, threads=1,
    memory_mapping=False, tmpdir=None,
):
    """Allele-based non-human fraction for a single variant.

    Convenience wrapper over :func:`classify_variants_alt_reads`.

    Example::

        vnhf = classify_variant_alt_reads(
            "sample.bam", "kraken2_db", "chr1", 12345, "A", "T",
            ref_fasta="ref.fa",
        )
        vnhf.nonhuman_fraction   # 0.0-1.0 over ALT-supporting reads
        vnhf.fractions.bacterial # per-domain breakdown
    """
    return classify_variants_alt_reads(
        bam_path, db_path, [(chrom, pos, ref, alt)],
        ref_fasta=ref_fasta, min_baseq=min_baseq, confidence=confidence,
        threads=threads, memory_mapping=memory_mapping, tmpdir=tmpdir,
    )[0]
