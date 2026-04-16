"""Backward-compatibility shim for the monolithic pipeline module.

The pipeline logic has been split into two sub-packages:

* :mod:`kmer_denovo_filter.vcf.pipeline` — VCF-mode pipeline
* :mod:`kmer_denovo_filter.discovery.pipeline` — Discovery-mode pipeline

This module re-exports all public and private names so that existing
``from kmer_denovo_filter.pipeline import ...`` statements continue to
work without modification.
"""

# Re-export everything from the VCF sub-package
from kmer_denovo_filter.vcf.pipeline import (  # noqa: F401
    _FRACTION_PRECISION,
    _build_span_bed_rows,
    _collect_child_kmers,
    _format_expanded_span_row,
    _format_span_row,
    _parse_kmer_votes,
    _parse_vcf_variants,
    _run_kraken2_on_reads,
    _write_annotated_vcf,
    _write_bed_from_rows,
    _write_informative_reads,
    _write_kraken2_expanded_span_bed,
    _write_kraken2_read_detail_bed,
    _write_kraken2_span_bed,
    _write_summary,
    run_pipeline,
)

# Re-export everything from the Discovery sub-package
from kmer_denovo_filter.discovery.pipeline import (  # noqa: F401
    SULOVARI_DNM_REGIONS,
    _anchor_and_cluster,
    _annotate_and_link_from_metadata,
    _classify_regions,
    _compare_candidates_to_regions,
    _count_parent_jellyfish,
    _evaluate_dnm_regions,
    _extract_child_kmers_discovery,
    _filter_parents_discovery,
    _parse_candidate_summary,
    _subtract_reference_kmers,
    _write_bed,
    _write_bedgraph,
    _write_bedpe,
    _write_discovery_summary,
    _write_empty_discovery_outputs,
    _write_informative_reads_discovery,
    _write_read_coverage_bed,
    run_discovery_pipeline,
)

# Re-export shared validation and utilities accessed via this module
from kmer_denovo_filter.utils import _check_tool, _validate_inputs  # noqa: F401

# Re-export names from kmer_utils that were previously in this
# module's namespace (tests use ``pipeline_mod.Kraken2Runner``, etc.)
from kmer_denovo_filter.kmer_utils import Kraken2Runner  # noqa: F401

# Re-export names from core.bam_scanner that were previously importable
# from this module (e.g. _extract_softclips used in test_kraken2_bed.py)
from kmer_denovo_filter.core.bam_scanner import _extract_softclips  # noqa: F401

# Re-export the report generator
from kmer_denovo_filter.report import generate_report  # noqa: F401
