#!/usr/bin/env python3
"""Summarise kmer-discovery runs across a cohort of trios.

Reads every trio listed in ``trios.tsv`` whose run finished (the per-trio
job leaves ``<trio>.done`` beside its outputs) and writes, under
``--out-dir``:

``trio_summary.tsv``
    One row per trio: region counts by class (``SV``, ``AMBIGUOUS``,
    ``SMALL``), SV-class regions by type, BEDPE junctions, k-mer and read
    totals, the child-only SVs of the published call set and how many
    discovery supported, wall time and per-module times.
``recurrent_regions.bed``
    Loci where regions from at least ``--recurrent-min-trios`` different
    children fall within ``--slack`` bp of one another.  Unrelated children
    do not share de novo events, so these are reference gaps, systematic
    alignment artefacts or common variation missing from the parents'
    calls, not de novos; they are dropped from the candidate list below.
``sv_candidates.tsv`` / ``sv_candidates.bed``
    Every ``SV``-class region that is not recurrent, with its evidence
    columns and the BEDPE junctions that start or end in it: the de novo
    SV candidates to review (IGV with the trio's informative BAM, or a
    dedicated SV caller on the staged CRAMs).
``sv_candidates.bedpe``
    The junctions of those candidates, with a trio column appended.
``cohort_summary.md``
    The headline numbers and the per-trio table, in Markdown.

Standard library only, so it runs on any login node.  Re-run it while trios
are still finishing; unfinished trios are listed as missing.
"""

import argparse
import collections
import csv
import json
import os
import re
import statistics
import sys

CLASSES = ("SV", "AMBIGUOUS", "SMALL")
SV_TYPES = ("DEL", "DUP", "INV", "INS", "BND")
EVIDENCE_COLUMNS = ("reads", "unique_kmers", "split_reads", "discordant_pairs",
                    "unmapped_mates", "breakpoint_reads", "large_indel_reads",
                    "max_clip_len")
MODULES = ("1", "2", "2b", "3")
MODULE_RE = re.compile(r"\[Module (\w+)\] Complete \(([^)]*)\)")
FINISHED_RE = re.compile(r"Pipeline finished in (.+)$")


# -- parsing ------------------------------------------------------------------

def parse_elapsed(text):
    """Seconds from kmer-discovery's elapsed strings ('1h 2m 3s', '5m 1.2s', '3.4s')."""
    total = 0.0
    for value, unit in re.findall(r"([\d.]+)\s*([hms])", text):
        total += float(value) * {"h": 3600, "m": 60, "s": 1}[unit]
    return total


def read_trios(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_bedpe(path):
    """Junctions of a trio's BEDPE as dicts; [] when the file is missing."""
    links = []
    if not os.path.isfile(path):
        return links
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            links.append({
                "chrom1": f[0], "start1": int(f[1]), "end1": int(f[2]),
                "chrom2": f[3], "start2": int(f[4]), "end2": int(f[5]),
                "sv_id": f[6], "reads": int(f[7]), "strand1": f[8],
                "strand2": f[9], "sv_type": f[10] if len(f) > 10 else ".",
            })
    return links


def read_timing(path):
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[-1] if rows else {}


def read_module_times(path):
    """{'1': seconds, ..., 'total': seconds} from the saved kmer-discovery log."""
    times = {}
    if not os.path.isfile(path):
        return times
    with open(path, errors="replace") as fh:
        for line in fh:
            m = MODULE_RE.search(line)
            if m:
                times[m.group(1)] = parse_elapsed(m.group(2))
                continue
            m = FINISHED_RE.search(line.rstrip())
            if m:
                times["total"] = parse_elapsed(m.group(1))
    return times


def read_child_only_counts(detail_path):
    """(child-only calls, calls passing the filters) from the detail table."""
    if not os.path.isfile(detail_path):
        return None, None
    n_all = n_pass = 0
    with open(detail_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return 0, 0
        for row in reader:
            if not row:
                continue
            n_all += 1
            n_pass += row[-1] == "PASS_FILTERS"
    return n_all, n_pass


def load_trio(results_dir, trio_id):
    """Everything the summary needs about one trio, or None when unfinished."""
    prefix = os.path.join(results_dir, trio_id, trio_id)
    if not (os.path.isfile(prefix + ".done")
            and os.path.isfile(prefix + ".metrics.json")):
        return None
    with open(prefix + ".metrics.json") as fh:
        metrics = json.load(fh)
    regions = metrics.get("regions", [])
    for r in regions:
        r["trio"] = trio_id
    n_all, n_pass = read_child_only_counts(prefix + ".child_only_svs.detail.tsv")
    return {
        "trio": trio_id,
        "metrics": metrics,
        "regions": regions,
        "links": read_bedpe(prefix + ".sv.bedpe"),
        "timing": read_timing(prefix + ".timing.tsv"),
        "module_times": read_module_times(prefix + ".log"),
        "dnm": metrics.get("dnm_evaluation"),
        "child_only_all": n_all,
        "child_only_pass": n_pass,
    }


# -- analysis -----------------------------------------------------------------

def cluster_regions(regions, slack):
    """Group regions from all trios that lie within *slack* bp of each other.

    Returns a list of clusters; each is a dict with the merged interval, the
    member regions and the set of trios.  Every region gets a ``cluster``
    index and an ``n_trios`` count.
    """
    clusters = []
    ordered = sorted(regions, key=lambda r: (r["chrom"], r["start"], r["end"]))
    current = None
    for r in ordered:
        if (current is not None and r["chrom"] == current["chrom"]
                and r["start"] <= current["end"] + slack):
            current["end"] = max(current["end"], r["end"])
            current["members"].append(r)
        else:
            current = {"chrom": r["chrom"], "start": r["start"],
                       "end": r["end"], "members": [r]}
            clusters.append(current)
    for i, c in enumerate(clusters):
        c["trios"] = sorted({m["trio"] for m in c["members"]})
        for m in c["members"]:
            m["cluster"] = i
            m["n_trios"] = len(c["trios"])
    return clusters


def links_for_region(region, links):
    """BEDPE lines of the region's trio that start or end in the region."""
    out = []
    key = (region["chrom"], region["start"], region["end"])
    for l in links:
        if (l["chrom1"], l["start1"], l["end1"]) == key:
            out.append(l)
        elif (l["chrom2"], l["start2"], l["end2"]) == key:
            out.append(l)
    return out


def summarise_trio(info):
    """The per-trio row of trio_summary.tsv as an ordered dict."""
    m = info["metrics"]
    regions = info["regions"]
    by_class = collections.Counter(r.get("class", ".") for r in regions)
    sv_types = collections.Counter(
        r.get("sv_type", ".") for r in regions if r.get("class") == "SV")
    links = info["links"]
    row = collections.OrderedDict()
    row["trio"] = info["trio"]
    row["status"] = "done"
    row["regions"] = m.get("candidate_regions", len(regions))
    for cls in CLASSES:
        row[cls] = by_class.get(cls, 0)
    for t in SV_TYPES:
        row["SV_" + t] = sv_types.get(t, 0)
    row["SV_untyped"] = sv_types.get(".", 0)
    row["SV_nonrecurrent"] = sum(
        1 for r in regions
        if r.get("class") == "SV" and r.get("n_trios", 1) < 2)
    row["junctions"] = len(links)
    row["child_candidate_kmers"] = m.get("child_candidate_kmers", "")
    row["non_ref_kmers"] = m.get("non_ref_kmers", "")
    row["proband_unique_kmers"] = m.get("proband_unique_kmers", "")
    row["informative_reads"] = m.get("informative_reads", "")
    row["child_only_svs"] = ("" if info["child_only_all"] is None
                             else info["child_only_all"])
    row["child_only_svs_filtered"] = ("" if info["child_only_pass"] is None
                                      else info["child_only_pass"])
    dnm = info["dnm"] or {}
    row["known_svs_tested"] = dnm.get("total_loci", 0 if info["dnm"] is None
                                      else dnm.get("total_loci", ""))
    row["known_svs_detected"] = dnm.get("detected", 0)
    timing = info["timing"]
    wall = timing.get("wall_s")
    row["wall_hours"] = f"{int(wall) / 3600:.2f}" if wall else ""
    row["threads"] = timing.get("threads", "")
    row["memory_gb"] = timing.get("memory_gb", "")
    row["node"] = timing.get("node", "")
    mt = info["module_times"]
    for mod in MODULES:
        row["module%s_hours" % mod] = (f"{mt[mod] / 3600:.2f}" if mod in mt
                                      else "")
    return row


def missing_row(trio_id):
    row = collections.OrderedDict()
    row["trio"] = trio_id
    row["status"] = "missing"
    return row


def _median(values):
    values = [v for v in values if v is not None]
    return statistics.median(values) if values else None


def _fmt(value, digits=1):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


# -- output -------------------------------------------------------------------

def write_tsv(path, rows, header):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(h, "")) for h in header) + "\n")


def write_recurrent(path, clusters, min_trios):
    recurrent = [c for c in clusters if len(c["trios"]) >= min_trios]
    with open(path, "w") as fh:
        fh.write("#chrom\tstart\tend\tn_trios\tn_regions\tclasses\tsv_types\ttrios\n")
        for c in recurrent:
            classes = collections.Counter(m.get("class", ".") for m in c["members"])
            types = collections.Counter(m.get("sv_type", ".") for m in c["members"])
            fh.write("\t".join((
                c["chrom"], str(c["start"]), str(c["end"]),
                str(len(c["trios"])), str(len(c["members"])),
                ",".join(f"{k}:{v}" for k, v in sorted(classes.items())),
                ",".join(f"{k}:{v}" for k, v in sorted(types.items())),
                ",".join(c["trios"]),
            )) + "\n")
    return recurrent


def write_candidates(out_dir, trios, min_trios):
    """The non-recurrent SV-class regions, with their junctions."""
    header = ["chrom", "start", "end", "trio", "class", "sv_type"] + \
        list(EVIDENCE_COLUMNS) + ["n_trios_at_locus", "junctions"]
    candidates = []
    bedpe_lines = []
    for info in trios:
        for r in info["regions"]:
            if r.get("class") != "SV" or r.get("n_trios", 1) >= min_trios:
                continue
            links = links_for_region(r, info["links"])
            row = collections.OrderedDict()
            row["chrom"], row["start"], row["end"] = r["chrom"], r["start"], r["end"]
            row["trio"] = info["trio"]
            row["class"], row["sv_type"] = r.get("class", "."), r.get("sv_type", ".")
            for col in EVIDENCE_COLUMNS:
                row[col] = r.get(col, "")
            row["n_trios_at_locus"] = r.get("n_trios", 1)
            row["junctions"] = ";".join(
                "%s:%d-%d|%s%s|%s|%d" % (
                    l["chrom2"] if (l["chrom1"], l["start1"], l["end1"])
                    == (r["chrom"], r["start"], r["end"]) else l["chrom1"],
                    l["start2"] if (l["chrom1"], l["start1"], l["end1"])
                    == (r["chrom"], r["start"], r["end"]) else l["start1"],
                    l["end2"] if (l["chrom1"], l["start1"], l["end1"])
                    == (r["chrom"], r["start"], r["end"]) else l["end1"],
                    l["strand1"], l["strand2"], l["sv_type"], l["reads"])
                for l in links) or "."
            candidates.append(row)
            for l in links:
                line = (l["chrom1"], l["start1"], l["end1"], l["chrom2"],
                        l["start2"], l["end2"], l["sv_id"], l["reads"],
                        l["strand1"], l["strand2"], l["sv_type"], info["trio"])
                if line not in bedpe_lines:
                    bedpe_lines.append(line)
    candidates.sort(key=lambda r: (r["chrom"], r["start"], r["trio"]))
    write_tsv(os.path.join(out_dir, "sv_candidates.tsv"), candidates, header)
    with open(os.path.join(out_dir, "sv_candidates.bed"), "w") as fh:
        fh.write("#chrom\tstart\tend\tname\tscore\tstrand\n")
        for r in candidates:
            name = "%s|%s|%s" % (r["trio"], r["sv_type"], r["class"])
            fh.write("\t".join((r["chrom"], str(r["start"]), str(r["end"]),
                                name, str(r["reads"]), ".")) + "\n")
    with open(os.path.join(out_dir, "sv_candidates.bedpe"), "w") as fh:
        fh.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\t"
                 "supporting_reads\tstrand1\tstrand2\tsv_type\ttrio\n")
        for line in sorted(bedpe_lines, key=lambda l: (l[0], l[1], l[11])):
            fh.write("\t".join(str(v) for v in line) + "\n")
    return candidates


def write_markdown(path, rows, trio_meta, candidates, recurrent, clusters,
                   min_trios, slack, n_selected):
    done = [r for r in rows if r["status"] == "done"]
    lines = ["# 1000 Genomes trio screen: cohort summary", ""]
    lines.append(f"* Trios selected: {n_selected}; finished: {len(done)}; "
                 f"missing: {n_selected - len(done)}")
    if not done:
        lines.append("")
        lines.append("No finished trio yet.")
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return
    regions = [r["regions"] for r in done]
    sv = [r["SV"] for r in done]
    sv_nr = [r["SV_nonrecurrent"] for r in done]
    lines.append(f"* Candidate regions per trio: median {_fmt(_median(regions))}, "
                 f"range {min(regions)}-{max(regions)}")
    lines.append(f"* `SV`-class regions per trio: median {_fmt(_median(sv))}, "
                 f"range {min(sv)}-{max(sv)}")
    lines.append(f"* Recurrent loci (regions from >= {min_trios} children within "
                 f"{slack} bp): {len(recurrent)} of {len(clusters)} clusters")
    lines.append(f"* De novo SV candidates after dropping recurrent loci: "
                 f"{len(candidates)} (median {_fmt(_median(sv_nr))} per trio)")
    types = collections.Counter(c["sv_type"] for c in candidates)
    if candidates:
        lines.append("* Candidate types: " + ", ".join(
            f"{t}: {n}" for t, n in sorted(types.items(),
                                           key=lambda kv: -kv[1])))
    tested = sum(r["known_svs_tested"] or 0 for r in done)
    detected = sum(r["known_svs_detected"] or 0 for r in done)
    if tested:
        lines.append(f"* Published child-only SVs passing the filters: {tested}; "
                     f"supported by discovery: {detected} "
                     f"({100.0 * detected / tested:.0f}%)")
    walls = [float(r["wall_hours"]) for r in done if r["wall_hours"]]
    if walls:
        lines.append(f"* Wall time per trio: median {_fmt(_median(walls), 2)} h, "
                     f"max {max(walls):.2f} h")
        for mod in MODULES:
            vals = [float(r["module%s_hours" % mod]) for r in done
                    if r.get("module%s_hours" % mod)]
            if vals:
                lines.append(f"  * Module {mod}: median {_fmt(_median(vals), 2)} h")
    lines.append("")
    lines.append("| Trio | Pop | Regions | SV | AMB | SMALL | DEL | DUP | INV | INS | BND | SV non-rec. | Junctions | Known SVs found | Wall h |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        meta = trio_meta.get(r["trio"], {})
        pop = meta.get("POPULATION", "")
        if r["status"] != "done":
            lines.append(f"| {r['trio']} | {pop} | (missing) |" + " |" * 12)
            continue
        known = (f"{r['known_svs_detected']}/{r['known_svs_tested']}"
                 if r["known_svs_tested"] else "-")
        lines.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            r["trio"], pop, r["regions"], r["SV"], r["AMBIGUOUS"], r["SMALL"],
            r["SV_DEL"], r["SV_DUP"], r["SV_INV"], r["SV_INS"], r["SV_BND"],
            r["SV_nonrecurrent"], r["junctions"], known, r["wall_hours"] or "-"))
    lines.append("")
    lines.append("Files: `trio_summary.tsv`, `recurrent_regions.bed`, "
                 "`sv_candidates.tsv` (with evidence and junctions), "
                 "`sv_candidates.bed`, `sv_candidates.bedpe`.")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--trios", required=True, help="trios.tsv manifest")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--recurrent-min-trios", type=int, default=2,
                        help="children sharing a locus for it to count as "
                             "recurrent (default: %(default)s)")
    parser.add_argument("--slack", type=int, default=500,
                        help="bp between regions of different trios for "
                             "them to share a locus (default: %(default)s)")
    args = parser.parse_args(argv)

    trio_rows = read_trios(args.trios)
    trio_meta = {t["TRIO_ID"]: t for t in trio_rows}
    os.makedirs(args.out_dir, exist_ok=True)

    loaded, rows = [], []
    for t in trio_rows:
        info = load_trio(args.results_dir, t["TRIO_ID"])
        if info is None:
            rows.append(missing_row(t["TRIO_ID"]))
            continue
        loaded.append(info)
    all_regions = [r for info in loaded for r in info["regions"]]
    clusters = cluster_regions(all_regions, args.slack)
    recurrent = write_recurrent(
        os.path.join(args.out_dir, "recurrent_regions.bed"), clusters,
        args.recurrent_min_trios)
    for info in loaded:
        rows.append(summarise_trio(info))
    rows.sort(key=lambda r: [t["TRIO_ID"] for t in trio_rows].index(r["trio"]))
    header = list(summarise_trio(loaded[0]).keys()) if loaded else ["trio", "status"]
    write_tsv(os.path.join(args.out_dir, "trio_summary.tsv"), rows, header)
    candidates = write_candidates(args.out_dir, loaded, args.recurrent_min_trios)
    write_markdown(os.path.join(args.out_dir, "cohort_summary.md"), rows,
                   trio_meta, candidates, recurrent, clusters,
                   args.recurrent_min_trios, args.slack, len(trio_rows))
    print(f"{len(loaded)} of {len(trio_rows)} trios finished; "
          f"{len(recurrent)} recurrent loci; {len(candidates)} SV candidates "
          f"-> {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
