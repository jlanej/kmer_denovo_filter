#!/usr/bin/env python3
"""List a trio's child-only SV genotypes from the 1000 Genomes SV call set.

The NYGC 3,202-sample Illumina SV call set (GATK-SV, svtools and Absinthe
merged; freeze V3 with allele frequencies) genotypes every sample at every
SV.  For one trio this script keeps the SVs genotyped non-reference in the
child and homozygous reference in both parents, and writes them in the
``--dnm-regions`` format of ``kmer-discovery``::

    chrom  pos  size  event_type

so a discovery run on the trio reports, per SV, whether it found reads at
the locus.  A second, wider table keeps the call set's own evidence for
review.

A child-only genotype in a population call set is more often a parent
genotyped wrongly than a de novo SV, so by default a call must be PASS,
have a population allele frequency of at most ``--max-af`` and genotype
qualities of at least ``--min-gq`` in all three members.  The counts before
and after these filters are printed, and both are worth reading: a locus
discovery supports despite failing them is as interesting as one it does
not support despite passing.

Sizes: deletions, duplications, inversions and complex SVs get their
reference span (|SVLEN|, else END - POS); insertions and breakends get
``.``, so discovery counts regions within ``--cluster-distance`` of the
position.  Multiallelic copy-number variants (``<CNV>``) carry copy states
rather than genotypes and are skipped unless ``--include-cnv`` is given.

Needs pysam.  Subsetting the 1.7 GB call set to a trio takes a few minutes;
00_setup.sh subsets it once to the cohort's samples, and the per-trio job
reads that smaller file.
"""

import argparse
import sys

import pysam

DEFAULT_TYPES = "DEL,DUP,INS,INV,BND,CTX,CPX"
REGIONS_HEADER = ("chrom", "pos", "size", "event_type")
DETAIL_HEADER = ("ID", "CHROM", "POS", "END", "SVTYPE", "ALT", "SVLEN", "AF",
                 "FILTER", "CHILD_GT", "CHILD_GQ", "FATHER_GT", "FATHER_GQ",
                 "MOTHER_GT", "MOTHER_GQ", "ALGORITHMS", "EVIDENCE")


def _gt(sample):
    alleles = sample.get("GT")
    if not alleles or any(a is None for a in alleles):
        return None
    return tuple(alleles)


def _fmt_gt(sample):
    alleles = sample.get("GT")
    if not alleles:
        return "."
    sep = "|" if getattr(sample, "phased", False) else "/"
    return sep.join("." if a is None else str(a) for a in alleles)


def _gq(sample):
    try:
        value = sample.get("GQ")
    except KeyError:
        return None
    return None if value is None else int(value)


def _first(value):
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def _join(value):
    if value is None:
        return "."
    if isinstance(value, tuple):
        return ",".join(str(v) for v in value)
    return str(value)


def event_size(record, svtype):
    """Reference span in bp for the --dnm-regions file, or None."""
    if svtype in ("INS", "BND", "CTX"):
        return None
    svlen = _first(record.info.get("SVLEN"))
    if svlen:
        return abs(int(svlen))
    # pysam folds INFO/END into record.stop
    if record.stop is not None and record.stop > record.pos:
        return record.stop - record.pos
    return None


def child_only_svs(vcf_path, child, father, mother, sv_types, pass_only=True,
                   max_af=None, min_gq=0, include_cnv=False):
    """Yield (record, child_only, passes_filters) for every SV record.

    *child_only* is True when the child carries an alternate allele and both
    parents are homozygous reference; *passes_filters* applies the quality
    filters on top of that.
    """
    with pysam.VariantFile(vcf_path) as vcf:
        missing = [s for s in (child, father, mother)
                   if s not in vcf.header.samples]
        if missing:
            raise ValueError("samples not in the VCF: " + ", ".join(missing))
        vcf.subset_samples([child, father, mother])
        for record in vcf:
            svtype = record.info.get("SVTYPE")
            if svtype is None:
                continue
            if svtype == "CNV" and not include_cnv:
                continue
            if sv_types and svtype not in sv_types:
                continue
            c, f, m = (_gt(record.samples[s]) for s in (child, father, mother))
            is_child_only = (
                c is not None and f is not None and m is not None
                and any(a > 0 for a in c)
                and all(a == 0 for a in f) and all(a == 0 for a in m))
            if not is_child_only:
                yield record, False, False
                continue
            ok = True
            if pass_only:
                filters = list(record.filter.keys())
                ok = filters in ([], ["PASS"])
            if ok and max_af is not None:
                af = record.info.get("AF")
                if af is not None:
                    af = af if isinstance(af, tuple) else (af,)
                    ok = max(a for a in af if a is not None) <= max_af
            if ok and min_gq > 0 and "GQ" in record.format:
                gqs = [_gq(record.samples[s]) for s in (child, father, mother)]
                ok = all(g is not None and g >= min_gq for g in gqs)
            yield record, True, ok


def write_outputs(records, out_regions, out_table, child, father, mother,
                  params):
    """Write the --dnm-regions file and the detail table; return the counts."""
    n_scanned = n_child_only = n_pass = 0
    with open(out_regions, "w") as regions, \
            (open(out_table, "w") if out_table else _null()) as table:
        regions.write(
            "# child-only SVs for %s (father %s, mother %s) from %s\n"
            % (child, father, mother, params["vcf"]))
        regions.write("# filters: %s\n" % params["filters"])
        regions.write("# columns: " + "\t".join(REGIONS_HEADER) + "\n")
        if table:
            table.write("\t".join(DETAIL_HEADER) + "\n")
        for record, is_child_only, ok in records:
            n_scanned += 1
            if not is_child_only:
                continue
            n_child_only += 1
            svtype = record.info["SVTYPE"]
            alt = (record.alts[0] if record.alts else ".")
            label = alt.strip("<>") if alt.startswith("<") else svtype
            if table:
                table.write("\t".join((
                    record.id or ".", record.chrom, str(record.pos),
                    str(record.stop), svtype, alt,
                    _join(record.info.get("SVLEN")),
                    _join(record.info.get("AF")),
                    ";".join(record.filter.keys()) or "PASS",
                    _fmt_gt(record.samples[child]),
                    _join(_gq(record.samples[child])),
                    _fmt_gt(record.samples[father]),
                    _join(_gq(record.samples[father])),
                    _fmt_gt(record.samples[mother]),
                    _join(_gq(record.samples[mother])),
                    _join(record.info.get("ALGORITHMS")),
                    _join(record.info.get("EVIDENCE")),
                )) + ("\tPASS_FILTERS\n" if ok else "\tFILTERED\n"))
            if not ok:
                continue
            n_pass += 1
            size = event_size(record, svtype)
            regions.write("%s\t%d\t%s\t%s\n" % (
                record.chrom, record.pos, "." if size is None else size, label))
    return n_scanned, n_child_only, n_pass


class _null:
    """A stand-in for an unwanted output file."""

    def __enter__(self):
        return None

    def __exit__(self, *exc):
        return False


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", required=True, help="SV call set (bgzipped)")
    parser.add_argument("--child", required=True)
    parser.add_argument("--father", required=True)
    parser.add_argument("--mother", required=True)
    parser.add_argument("--out-regions", required=True,
                        help="--dnm-regions file for kmer-discovery")
    parser.add_argument("--out-table", default=None,
                        help="detail table of every child-only call, "
                             "filtered or not")
    parser.add_argument("--sv-types", default=DEFAULT_TYPES,
                        help="SVTYPEs to keep (default: %(default)s)")
    parser.add_argument("--include-cnv", action="store_true",
                        help="also keep multiallelic <CNV> records")
    parser.add_argument("--max-af", type=float, default=0.01,
                        help="maximum population AF (default: %(default)s; "
                             "negative = no limit)")
    parser.add_argument("--min-gq", type=int, default=20,
                        help="minimum GQ in all three members "
                             "(default: %(default)s)")
    parser.add_argument("--no-pass-only", action="store_true",
                        help="keep calls that failed the call set's filters")
    args = parser.parse_args(argv)

    sv_types = {t.strip() for t in args.sv_types.split(",") if t.strip()}
    max_af = None if args.max_af < 0 else args.max_af
    filters = "pass_only=%s max_af=%s min_gq=%d sv_types=%s cnv=%s" % (
        not args.no_pass_only, max_af, args.min_gq, ",".join(sorted(sv_types)),
        args.include_cnv)
    try:
        records = child_only_svs(
            args.vcf, args.child, args.father, args.mother, sv_types,
            pass_only=not args.no_pass_only, max_af=max_af,
            min_gq=args.min_gq, include_cnv=args.include_cnv)
        counts = write_outputs(
            records, args.out_regions, args.out_table, args.child,
            args.father, args.mother, {"vcf": args.vcf, "filters": filters})
    except (ValueError, OSError) as exc:
        raise SystemExit(f"error: {exc}")
    n_scanned, n_child_only, n_pass = counts
    print(f"{args.child}: {n_scanned} SV records scanned, {n_child_only} "
          f"child-only, {n_pass} after filters ({filters}) -> "
          f"{args.out_regions}", file=sys.stderr)


if __name__ == "__main__":
    main()
