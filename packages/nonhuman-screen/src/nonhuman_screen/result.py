"""Result value objects for non-human read screening.

These dataclasses turn the low-level :class:`~nonhuman_screen.engine.Kraken2Runner`
``Result`` (a bag of per-domain read-name sets) into the two reductions
callers actually want:

* :class:`TaxonomicFractions` — per-domain fractions over an arbitrary set of
  reads, either the whole batch or the reads supporting one variant allele.
* :class:`VariantNHF` — the non-human fraction of the reads supporting a single
  variant's ALT allele, i.e. the "allele-based NHF" primitive.

``result.py`` intentionally imports nothing from :mod:`nonhuman_screen.engine`
so that ``engine`` can depend on it without a cycle.  It reads domain read-name
sets off any object exposing the ``*_read_names`` attributes.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

# Ordered so ``fractions()`` output and BED/JSON columns are stable.
DOMAINS = (
    "nonhuman",
    "bacterial",
    "archaeal",
    "fungal",
    "protist",
    "viral",
    "univec_core",
    "human_lineage",
    "unclassified",
)

# Map each public domain name to the read-name set attribute and count
# attribute on an engine ``Result``.  Every domain uses a ``{d}_count``
# attribute except ``unclassified``, which the engine stores as ``unclassified``.
_SET_ATTR = {d: f"{d}_read_names" for d in DOMAINS}
_COUNT_ATTR = {d: f"{d}_count" for d in DOMAINS}
_COUNT_ATTR["unclassified"] = "unclassified"


def _round(numerator: int, denominator: int) -> float:
    return round(numerator / denominator, 4) if denominator else 0.0


@dataclass(frozen=True)
class TaxonomicFractions:
    """Per-domain fractions in ``[0, 1]`` over some denominator of reads.

    The domains partition classified reads such that
    ``nonhuman + univec_core + human_lineage + unclassified == 1.0`` (modulo
    rounding); ``bacterial``/``archaeal``/``fungal``/``protist``/``viral`` are
    sub-breakdowns of ``nonhuman``.
    """

    total: int
    nonhuman: float
    bacterial: float
    archaeal: float
    fungal: float
    protist: float
    viral: float
    univec_core: float
    human_lineage: float
    unclassified: float

    @classmethod
    def from_result(cls, result) -> "TaxonomicFractions":
        """Whole-batch fractions: each domain count divided by ``result.total``."""
        total = result.total
        return cls(
            total=total,
            **{d: _round(getattr(result, _COUNT_ATTR[d]), total) for d in DOMAINS},
        )

    @classmethod
    def over_reads(cls, result, read_names) -> "TaxonomicFractions":
        """Fractions of *read_names* falling in each domain of *result*.

        This is the per-variant reduction: *read_names* is the set of reads
        supporting a variant allele, and the denominator is ``len(read_names)``.
        Names absent from *result* (e.g. never classified) count toward the
        denominator but no numerator.
        """
        names = set(read_names)
        denom = len(names)
        return cls(
            total=denom,
            **{
                d: _round(len(names & getattr(result, _SET_ATTR[d])), denom)
                for d in DOMAINS
            },
        )

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class VariantNHF:
    """Non-human content of the reads supporting one variant's ALT allele.

    ``pos`` is 0-based, matching pysam's reference coordinates and the internal
    key convention of consuming pipelines.  ``variant_key`` is
    ``"{chrom}:{pos}:{ref}:{alt}"``.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    supporting_read_names: frozenset
    fractions: TaxonomicFractions

    @property
    def variant_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def supporting_reads(self) -> int:
        return len(self.supporting_read_names)

    @property
    def nonhuman_fraction(self) -> float:
        """Fraction of ALT-supporting reads classified as non-human."""
        return self.fractions.nonhuman

    def to_dict(self) -> dict:
        return {
            "variant_key": self.variant_key,
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "supporting_reads": self.supporting_reads,
            "nonhuman_fraction": self.nonhuman_fraction,
            "fractions": self.fractions.to_dict(),
        }
