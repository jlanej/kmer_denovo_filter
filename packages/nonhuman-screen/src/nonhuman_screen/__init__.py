"""nonhuman-screen: classify sequencing reads by non-human taxonomic content.

Wraps kraken2's LCA classification and reduces per-read verdicts to a
non-human fraction (NHF) plus per-domain breakdowns, with a human-homology
guard and UniVec-Core exclusion.  Sample-agnostic: it classifies whatever
named reads you hand it, so it works equally on a proband, a parent, a tumour,
or a plain FASTQ.

Two levels of API:

* Engine (stdlib only) — ``Kraken2Runner.classify_sequences({name: seq})``.
* BAM/allele helpers (needs the ``[bam]`` extra) — ``classify_variant_alt_reads``
  and ``classify_variants_alt_reads`` compute the non-human fraction of the
  reads supporting a variant's ALT allele.
"""

from nonhuman_screen.engine import (
    ClassificationResult,
    Kraken2Runner,
    _ARCHAEA_TAXID,
    _BACTERIA_TAXID,
    _EUKARYOTA_TAXID,
    _FUNGI_TAXID,
    _HUMAN_TAXID,
    _METAZOA_TAXID,
    _UNIVEC_CORE_TAXID,
    _VIRIDIPLANTAE_TAXID,
    _VIRUSES_TAXID,
)
from nonhuman_screen.alleles import read_supports_alt
from nonhuman_screen.result import TaxonomicFractions, VariantNHF
from nonhuman_screen.votes import parse_kmer_votes

__version__ = "0.1.0"

__all__ = [
    "Kraken2Runner",
    "ClassificationResult",
    "TaxonomicFractions",
    "VariantNHF",
    "read_supports_alt",
    "parse_kmer_votes",
    "__version__",
]

# BAM/allele helpers require pysam ([bam] extra). Expose them at the top level
# when available, but never make importing the core engine depend on pysam.
try:
    from nonhuman_screen.bam import (  # noqa: F401
        classify_reads_from_bam,
        classify_variant_alt_reads,
        classify_variants_alt_reads,
        reads_supporting_alt,
    )

    __all__ += [
        "classify_reads_from_bam",
        "classify_variant_alt_reads",
        "classify_variants_alt_reads",
        "reads_supporting_alt",
    ]
except ImportError:
    pass
