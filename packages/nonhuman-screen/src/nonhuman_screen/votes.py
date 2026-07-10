"""Parse a Kraken2 per-read k-mer detail string into taxid vote summaries."""

from nonhuman_screen.engine import _HUMAN_TAXID


def parse_kmer_votes(kmer_string, name_map=None, top_n=10):
    """Parse a Kraken2 kmer_detail_string into vote summaries.

    Args:
        kmer_string: Raw Kraken2 kmer detail string (space-separated
            ``taxid:count`` tokens; ``|:|`` separates paired-end mates).
        name_map: Optional ``{taxid: name}`` dict for named output.
        top_n: Maximum number of entries to include (sorted descending
            by count).

    Returns:
        Tuple ``(kmer_votes, kmer_votes_named, total_kmers,
        human_kmer_count)`` where ``kmer_votes`` is
        ``taxid1:count1;taxid2:count2;...`` and ``kmer_votes_named``
        replaces taxids with names.  Taxid ``0`` is rendered as
        ``unclassified`` in the named column.  ``A`` (ambiguous) tokens
        are excluded.
    """
    if not kmer_string:
        return ("", "", 0, 0)

    counts = {}  # taxid (int) -> total count
    for token in kmer_string.replace("|:|", " ").split():
        taxid_str, _, count_str = token.partition(":")
        if not taxid_str or not count_str:
            continue
        try:
            tid = int(taxid_str)
            cnt = int(count_str)
        except ValueError:
            continue
        counts[tid] = counts.get(tid, 0) + cnt

    total_kmers = sum(counts.values())
    human_kmer_count = counts.get(_HUMAN_TAXID, 0)

    # Sort descending by count, truncate
    sorted_votes = sorted(counts.items(), key=lambda x: (-x[1], x[0]))
    top_votes = sorted_votes[:top_n]

    kmer_votes = ";".join(f"{tid}:{cnt}" for tid, cnt in top_votes)

    def _name_for(tid):
        if tid == 0:
            return "unclassified"
        if name_map and tid in name_map:
            return name_map[tid]
        return str(tid)

    kmer_votes_named = ";".join(
        f"{_name_for(tid)}:{cnt}" for tid, cnt in top_votes
    )

    return (kmer_votes, kmer_votes_named, total_kmers, human_kmer_count)
