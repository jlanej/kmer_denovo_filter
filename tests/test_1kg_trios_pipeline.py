"""Hermetic end-to-end test of the ``examples/1kg_trios`` SLURM scripts.

Every external dependency is shimmed on PATH: ``sbatch`` runs the submitted
script inline with a SLURM-like environment, ``squeue`` reports nothing,
``curl`` and ``aria2c`` fail so that ``wget`` (serving files from a local
store) carries every download, and ``apptainer`` fakes ``pull``,
``jellyfish count``, ``samtools faidx`` and ``kmer-discovery`` while running
``bcftools`` and this example's Python helpers for real.  Three trios, two of
them siblings sharing both parents and one with a corrupt CRAM, cover the
full cycle: setup, manifest, reference index, staging with MD5 checks, the
shared-parent bookkeeping, a failed trio healed on a resweep, cleanup, and
the cohort summary.  Nothing touches the network or a scheduler.
"""

import hashlib
import os
import shutil
import stat
import subprocess
import sys
import textwrap

import pysam
import pytest

from tests.test_1kg_trios_helpers import (
    INDEX_HEADER, INDEX_SAMPLES, PED, SV_VCF_HEADER, index_row, sv_record,
)

EXAMPLE_DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir, "examples", "1kg_trios"))
BCFTOOLS = shutil.which("bcftools")


def _executable(path, content):
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)


def make_shims(bin_dir, remote_dir, log_dir):
    """Fake sbatch, squeue, curl, wget, aria2c and apptainer in *bin_dir*."""
    bin_dir.mkdir()
    os.symlink(sys.executable, bin_dir / "python3")
    _executable(bin_dir / "squeue", "#!/usr/bin/env bash\nexit 0\n")
    _executable(bin_dir / "aria2c", "#!/usr/bin/env bash\nexit 1\n")
    # curl: serves the store for setup's fetch (-o DEST URL); every other
    # use (HEAD requests, range downloads) fails so the manager falls
    # through to wget.
    _executable(bin_dir / "curl", textwrap.dedent(f"""\
        #!/usr/bin/env bash
        dest=""; url=""; head=0
        while [[ $# -gt 0 ]]; do
          case "$1" in
            -o) dest="$2"; shift 2 ;;
            -I) head=1; shift ;;
            -r|--retry|--retry-delay|--max-time) shift 2 ;;
            -C) shift 2 ;;
            -*) shift ;;
            *) url="$1"; shift ;;
          esac
        done
        [[ -n "$dest" && $head -eq 0 ]] || exit 1
        src="{remote_dir}/$(basename "$url")"
        [[ -f "$src" ]] || exit 22
        cp "$src" "$dest"
        """))
    _executable(bin_dir / "wget", textwrap.dedent(f"""\
        #!/usr/bin/env bash
        dest=""; url=""
        while [[ $# -gt 0 ]]; do
          case "$1" in
            -O) dest="$2"; shift 2 ;;
            -*) shift ;;
            *) url="$1"; shift ;;
          esac
        done
        src="{remote_dir}/$(basename "$url")"
        [[ -f "$src" ]] || exit 1
        cp "$src" "$dest"
        """))
    _executable(bin_dir / "sbatch", textwrap.dedent(f"""\
        #!/usr/bin/env bash
        # Runs the submitted script inline under a SLURM-like environment and
        # prints a job id, as --parsable would.
        script=""; args=(); cpus=2; mem=4096
        while [[ $# -gt 0 ]]; do
          if [[ -n "$script" ]]; then args+=("$1"); shift; continue; fi
          case "$1" in
            --cpus-per-task=*) cpus="${{1#*=}}"; shift ;;
            --mem=*) m="${{1#*=}}"; m="${{m%G}}"; mem=$(( m * 1024 )); shift ;;
            --*) shift ;;
            *) script="$1"; shift ;;
          esac
        done
        jobid=$(( RANDOM + 1000 ))
        echo "$jobid $(basename "$script") ${{args[*]:-}}" >> "{log_dir}/sbatch.log"
        # Like a site whose sbatch exports nothing: the job sees none of the
        # submitter's configuration variables, only what --config restores.
        scrub=()
        for v in $(env | cut -d= -f1 | grep -E '^(WORK_DIR|CRAM_DIR|RESULTS_DIR|LOG_DIR|STATE_DIR|TMP_ROOT|N_TRIOS|TRIO_|DISCOVERY_|USE_SV_CALLSET|SBATCH_EXTRA_ARGS|EXPECTED_COMPLETE_TRIOS|KMER_TRIOS_CONFIG|MODULES_TO_LOAD)'); do scrub+=("-u" "$v"); done
        env ${{scrub[@]+"${{scrub[@]}}"}} SLURM_JOB_ID=$jobid SLURM_CPUS_PER_TASK=$cpus SLURM_MEM_PER_NODE=$mem \\
          bash "$script" ${{args[@]+"${{args[@]}}"}} >> "{log_dir}/jobs.log" 2>&1 \\
          || echo "job $jobid failed" >> "{log_dir}/jobs.log"
        echo "$jobid"
        """))
    _executable(bin_dir / "fake-kmer-discovery", textwrap.dedent("""\
        #!/usr/bin/env bash
        # Writes the outputs kmer-discovery would, records its arguments,
        # and fails for the child named in FAKE_DISCOVERY_FAIL_CHILD.
        prefix=""; dnm=""; child=""; prev=""
        for a in "$@"; do
          case "$prev" in
            --out-prefix) prefix="$a" ;;
            --dnm-regions) dnm="$a" ;;
            --child) child="$a" ;;
          esac
          prev="$a"
        done
        sample="$(basename "$child" .cram)"
        echo "$*" > "${prefix}.args"
        if [[ "$sample" == "${FAKE_DISCOVERY_FAIL_CHILD:-}" ]]; then
          echo "fake failure for $sample"; exit 3
        fi
        pos=$(( $(printf '%s' "$sample" | cksum | cut -d' ' -f1) % 100000 + 10000 ))
        printf '#filters: min_distinct_kmers=1\\n#chrom\\tstart\\tend\\treads\\tunique_kmers\\tsplit_reads\\tdiscordant_pairs\\tmax_clip_len\\tunmapped_mates\\tclass\\tbreakpoint_reads\\tlarge_indel_reads\\tsv_type\\n' > "${prefix}.bed"
        printf 'chr1\\t1000\\t1300\\t12\\t40\\t5\\t3\\t60\\t0\\tSV\\t6\\t0\\tDEL\\n' >> "${prefix}.bed"
        printf 'chr2\\t%d\\t%d\\t9\\t30\\t4\\t2\\t55\\t0\\tSV\\t5\\t0\\tDUP\\n' "$pos" $(( pos + 250 )) >> "${prefix}.bed"
        printf '#chrom1\\tstart1\\tend1\\tchrom2\\tstart2\\tend2\\tsv_id\\tsupporting_reads\\tstrand1\\tstrand2\\tsv_type\\n' > "${prefix}.sv.bedpe"
        printf 'chr2\\t%d\\t%d\\tchr2\\t%d\\t%d\\tSV_1\\t7\\t-\\t+\\tDUP\\n' "$pos" $(( pos + 250 )) $(( pos + 900 )) $(( pos + 1200 )) >> "${prefix}.sv.bedpe"
        dnm_json=""
        if [[ -n "$dnm" ]]; then
          n=$(grep -cv '^#' "$dnm")
          dnm_json=", \\"dnm_evaluation\\": {\\"source\\": \\"$(basename "$dnm")\\", \\"slack_bp\\": 500, \\"total_loci\\": $n, \\"detected\\": 1, \\"detection_rate\\": 0.5, \\"loci\\": []}"
        fi
        cat > "${prefix}.metrics.json" <<JSON
        {"mode": "discovery", "child_candidate_kmers": 5000, "non_ref_kmers": 900,
         "proband_unique_kmers": 120, "informative_reads": 60, "unmapped_informative_reads": 2,
         "candidate_regions": 2,
         "regions": [
          {"chrom": "chr1", "start": 1000, "end": 1300, "size": 300, "reads": 12, "unique_kmers": 40,
           "split_reads": 5, "discordant_pairs": 3, "max_clip_len": 60, "unmapped_mates": 0,
           "breakpoint_reads": 6, "large_indel_reads": 0, "sv_type": "DEL", "class": "SV"},
          {"chrom": "chr2", "start": $pos, "end": $(( pos + 250 )), "size": 250, "reads": 9, "unique_kmers": 30,
           "split_reads": 4, "discordant_pairs": 2, "max_clip_len": 55, "unmapped_mates": 0,
           "breakpoint_reads": 5, "large_indel_reads": 0, "sv_type": "DUP", "class": "SV"}
         ]$dnm_json}
        JSON
        echo "fake summary" > "${prefix}.summary.txt"
        : > "${prefix}.informative.bam"
        echo "[Module 1] Complete (2m 0.0s)"
        echo "[Module 2] Complete (1m 0.0s)"
        echo "[Module 2b] Complete (3.0s, index: 1 MB)"
        echo "[Module 3] Complete (30.0s)"
        echo "Pipeline finished in 3m 33.0s"
        """))
    _executable(bin_dir / "apptainer", textwrap.dedent(f"""\
        #!/usr/bin/env bash
        # pull: writes the .sif. exec: strips --bind and the image, then
        # dispatches on the tool.
        cmd="$1"; shift
        if [[ "$cmd" == "pull" ]]; then echo "fake sif" > "$1"; exit 0; fi
        [[ "$cmd" == "exec" ]] || exit 1
        while [[ $# -gt 0 ]]; do
          case "$1" in --bind) shift 2 ;; -*) shift ;; *) break ;; esac
        done
        shift   # the image
        tool="$1"; shift
        case "$tool" in
          aria2c) exit 1 ;;
          samtools) [[ "$1" == faidx ]] && echo "chr1\\t8\\t6\\t8\\t9" > "$2.fai"; exit 0 ;;
          jellyfish)
            out=""; prev=""
            for a in "$@"; do [[ "$prev" == "-o" ]] && out="$a"; prev="$a"; done
            echo "fake jellyfish index" > "$out"; exit 0 ;;
          bcftools) exec "{BCFTOOLS or 'bcftools'}" "$@" ;;
          python3)
            script="$1"; shift
            script="${{script/\\/app\\/examples\\/1kg_trios/{EXAMPLE_DIR}}}"
            exec "{sys.executable}" "$script" "$@" ;;
          kmer-discovery) exec "{bin_dir}/fake-kmer-discovery" "$@" ;;
          *) echo "fake apptainer: unexpected tool $tool" >&2; exit 1 ;;
        esac
        """))


def write_remote_store(remote, cram_bytes):
    """The public files the setup and the manager fetch, by basename."""
    remote.mkdir()
    (remote / "GRCh38_full_analysis_set_plus_decoy_hla.fa").write_text(">chr1\nACGTACGT\n")
    (remote / "GRCh38_full_analysis_set_plus_decoy_hla.fa.fai").write_text("chr1\t8\t6\t8\t9\n")
    (remote / "20130606_g1k_3202_samples_ped_population.txt").write_text(PED)
    for name, batch in (("1000G_2504_high_coverage.sequence.index", "2504"),
                        ("1000G_698_related_high_coverage.sequence.index", "698")):
        rows = [index_row(s, p, b, md5=hashlib.md5(cram_bytes[s.strip()]).hexdigest())
                for s, p, b in INDEX_SAMPLES if b == batch]
        (remote / name).write_text("##FileDate=20200820\n" + INDEX_HEADER + "\n"
                                   + "\n".join(rows) + "\n")
    for sample, data in cram_bytes.items():
        (remote / f"{sample}.final.cram").write_bytes(data)
        (remote / f"{sample}.final.cram.crai").write_bytes(b"crai " + sample.encode())
    # The SV call set: one child-only deletion for HG00405, one common
    # (filtered) insertion for HG00406, nothing for NA19240.
    samples = [s.strip() for s, _, _ in INDEX_SAMPLES]
    def gts(**per_sample):
        return [per_sample.get(s, "0/0:60") for s in samples]
    vcf = remote / "sv.vcf"
    vcf.write_text(
        SV_VCF_HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples) + "\n"
        + sv_record("chr1", 1000, "DEL1", "<DEL>", "PASS",
                    "SVTYPE=DEL;SVLEN=-300;END=1300;AF=0.001", gts(HG00405="0/1:70"))
        + sv_record("chr1", 5000, "INS1", "<INS:ME:ALU>", "PASS",
                    "SVTYPE=INS;SVLEN=300;END=5001;AF=0.2", gts(HG00406="0/1:70")))
    gz = str(remote / "1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz")
    pysam.tabix_compress(str(vcf), gz, force=True)
    pysam.tabix_index(gz, preset="vcf", force=True)


@pytest.fixture
def workspace(tmp_path):
    cram_bytes = {s.strip(): (s.strip() * 100).encode() for s, _, _ in INDEX_SAMPLES}
    remote = tmp_path / "remote"
    write_remote_store(remote, cram_bytes)
    # NA19240's CRAM on the "server" does not match the index MD5
    (remote / "NA19240.final.cram").write_bytes(b"corrupt")
    logs = tmp_path / "shimlogs"
    logs.mkdir()
    bin_dir = tmp_path / "bin"
    make_shims(bin_dir, remote, logs)
    work = tmp_path / "work"
    env = {
        "PATH": os.pathsep.join([str(bin_dir), "/usr/bin", "/bin", "/usr/sbin", "/sbin"]
                                + ([os.path.dirname(BCFTOOLS)] if BCFTOOLS else [])),
        "HOME": str(tmp_path),
        "WORK_DIR": str(work),
        "N_TRIOS": "0",
        "TRIO_SEED": "",
        "TRIO_ALL_SIBLINGS": "1",
        "EXPECTED_COMPLETE_TRIOS": "3",
        "TRIO_SLOTS": "1",
        "DISCOVERY_MEM": "4G",
        "DISCOVERY_CPUS": "2",
        "USE_SV_CALLSET": "1" if BCFTOOLS else "0",
        "SBATCH_EXTRA_ARGS": "--partition=fake --account=fake",
    }
    return {"env": env, "work": work, "remote": remote, "logs": logs,
            "cram_bytes": cram_bytes}


def run(script, env, *args, check=True):
    result = subprocess.run(
        ["bash", os.path.join(EXAMPLE_DIR, script)] + list(args),
        env=env, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise AssertionError(
            f"{script} exited {result.returncode}\n--- stdout ---\n"
            f"{result.stdout}\n--- stderr ---\n{result.stderr}")
    return result


def test_setup_stage_discover_and_summarise(workspace):
    env, work, remote, logs = (workspace[k] for k in ("env", "work", "remote", "logs"))

    # ── 00_setup.sh ───────────────────────────────────────────────────────
    out = run("00_setup.sh", env).stdout
    assert (work / "kmer_denovo_filter.sif").exists()
    ref = work / "reference" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
    assert ref.exists() and (ref.parent / (ref.name + ".fai")).exists()
    trios = (work / "trios.tsv").read_text().splitlines()
    assert [t.split("\t")[0] for t in trios[1:]] == ["HG00405", "HG00406", "NA19240"]
    assert len((work / "trio_files.tsv").read_text().splitlines()) == 10
    sbatch_log = (logs / "sbatch.log").read_text()
    assert "00b_build_ref_index.sh" in sbatch_log
    assert (ref.parent / (ref.name + ".k31.jf")).exists()
    assert "3 trios selected" in out
    if BCFTOOLS:
        cohort_vcf = work / "sv_callset" / "cohort_trios.sv.vcf.gz"
        assert cohort_vcf.exists() and (cohort_vcf.parent / (cohort_vcf.name + ".tbi")).exists()
        with pysam.VariantFile(str(cohort_vcf)) as vf:
            assert sorted(vf.header.samples) == sorted(
                ["HG00403", "HG00404", "HG00405", "HG00406", "NA19238", "NA19239",
                 "NA19240"])
    # idempotent: a second run changes nothing and resubmits nothing
    run("00_setup.sh", env)
    assert (logs / "sbatch.log").read_text() == sbatch_log

    # ── globus batch files for the whole cohort ───────────────────────────
    run("globus_batch_from_manifest.sh", env, str(work / "gb"))
    crai_lines = (work / "gb_crai.txt").read_text().splitlines()
    cram_lines = (work / "gb_cram.txt").read_text().splitlines()
    assert len(crai_lines) == len(cram_lines) == 7     # shared parents once
    assert cram_lines[0].startswith('"/vol1/run/') and cram_lines[0].endswith(
        f'"{work}/crams/HG00405.cram"')

    # ── 01b by hand with "verify" on a corrupt pre-staged CRAM ────────────
    crams = work / "crams"
    for sample in ("NA19240", "NA19239", "NA19238"):
        (crams / f"{sample}.cram").write_bytes(
            b"corrupt" if sample == "NA19240" else workspace["cram_bytes"][sample])
        (crams / f"{sample}.cram.crai").write_bytes(b"crai")
    result = run("01b_discover_trio.sh", env, "NA19240", "verify", check=False)
    assert result.returncode == 1
    assert "failed MD5" in result.stdout
    assert (crams / "NA19240.cram").exists()           # left in place
    assert not (work / "results" / "NA19240" / "NA19240.done").exists()
    for sample in ("NA19240", "NA19239", "NA19238"):
        (crams / f"{sample}.cram").unlink()
        (crams / f"{sample}.cram.crai").unlink()

    # ── sweep 1: siblings succeed, NA19240 fails its MD5 twice ────────────
    env1 = dict(env, MANAGER_LOCAL="1")
    result = run("01_stage_and_discover.sh", env1, check=False)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "2 trios dispatched to discovery, 1 failed" in result.stdout
    for trio in ("HG00405", "HG00406"):
        prefix = work / "results" / trio / trio
        assert (work / "results" / trio / f"{trio}.done").exists()
        for ext in ("bed", "sv.bedpe", "metrics.json", "summary.txt", "log",
                    "timing.tsv", "args"):
            assert (work / "results" / trio / f"{trio}.{ext}").exists(), ext
        args = (work / "results" / trio / f"{trio}.args").read_text()
        assert f"--child {crams}/{trio}.cram" in args
        assert "--father " + str(crams / "HG00403.cram") in args
        assert "--ref-jf " + str(ref) + ".k31.jf" in args
        assert "--threads 2 --memory 4" in args
        assert "--report " + str(prefix) + ".report.html" in args
    args_405 = (work / "results" / "HG00405" / "HG00405.args").read_text()
    args_406 = (work / "results" / "HG00406" / "HG00406.args").read_text()
    if BCFTOOLS:
        assert "--dnm-regions" in args_405        # one rare child-only DEL
        assert "--dnm-regions" not in args_406    # its only call is common
        detail = (work / "results" / "HG00406" / "HG00406.child_only_svs.detail.tsv")
        assert detail.read_text().splitlines()[1].startswith("INS1\tchr1\t5000")
    else:
        assert "--dnm-regions" not in args_405
    stage_log = (work / "logs" / "stage_NA19240.log").read_text()
    assert stage_log.count("MD5 mismatch") == 2
    assert "failed MD5 after a fresh download" in stage_log
    assert (work / "state" / "results" / "NA19240").read_text().strip() == "fail"
    jobs_log = (logs / "jobs.log").read_text()
    # the first sibling kept the shared parents for the second
    assert "Keeping HG00403.cram: trio HG00406 still needs it" in jobs_log
    assert "Removed HG00405.cram" in jobs_log
    assert "Removed HG00403.cram" in jobs_log       # by HG00406's job
    assert sorted(os.listdir(crams)) == []
    timing = (work / "results" / "HG00405" / "HG00405.timing.tsv").read_text().splitlines()
    assert timing[0].startswith("trio\tchild\twall_s")
    assert timing[1].split("\t")[:2] == ["HG00405", "HG00405"]
    assert timing[1].split("\t")[6:8] == ["2", "4"]
    assert timing[1].split("\t")[10] == "2"           # regions

    # ── sweep 2 after the server's copy is fixed ──────────────────────────
    (remote / "NA19240.final.cram").write_bytes(workspace["cram_bytes"]["NA19240"])
    result = run("01_stage_and_discover.sh", env1)
    assert "Finished:            2" in result.stdout
    assert "1 trios dispatched to discovery, 0 failed" in result.stdout
    assert (work / "results" / "NA19240" / "NA19240.done").exists()
    assert sorted(os.listdir(crams)) == []
    submissions = [l for l in (logs / "sbatch.log").read_text().splitlines()
                   if "01b_discover_trio.sh" in l]
    assert [l.split()[-1] for l in submissions] == ["HG00405", "HG00406", "NA19240"]
    assert all(f"--config {work}/state/config.env" in l for l in submissions)

    # a third sweep finds nothing to do
    result = run("01_stage_and_discover.sh", env1)
    assert "Nothing to stage." in result.stdout

    # ── 02_summarize.sh ───────────────────────────────────────────────────
    run("02_summarize.sh", env)
    cohort = work / "results" / "cohort"
    rows = [l.split("\t") for l in (cohort / "trio_summary.tsv").read_text().splitlines()]
    header = rows[0]
    by_trio = {r[0]: dict(zip(header, r)) for r in rows[1:]}
    assert all(by_trio[t]["status"] == "done" for t in ("HG00405", "HG00406", "NA19240"))
    assert by_trio["HG00405"]["SV"] == "2"
    assert by_trio["HG00405"]["SV_nonrecurrent"] == "1"
    assert by_trio["HG00405"]["module1_hours"] == "0.03"
    recurrent = [l for l in (cohort / "recurrent_regions.bed").read_text().splitlines()
                 if not l.startswith("#")]
    assert len(recurrent) == 1
    assert recurrent[0].startswith("chr1\t1000\t1300\t3\t3\tSV:3\tDEL:3\t")
    candidates = [l for l in (cohort / "sv_candidates.tsv").read_text().splitlines()[1:]]
    assert len(candidates) == 3
    assert all("\tDUP\t" in c and c.endswith("|-+|DUP|7") for c in candidates)
    md = (cohort / "cohort_summary.md").read_text()
    assert "Trios selected: 3; finished: 3; missing: 0" in md
    assert "Recurrent loci (regions from >= 2 children within 500 bp): 1 of" in md
    if BCFTOOLS:
        assert by_trio["HG00405"]["known_svs_tested"] == "1"
        assert "Published child-only SVs passing the filters: 1; supported by discovery: 1" in md


def test_failed_discovery_keeps_the_crams(workspace):
    env, work, remote, logs = (workspace[k] for k in ("env", "work", "remote", "logs"))
    (remote / "NA19240.final.cram").write_bytes(workspace["cram_bytes"]["NA19240"])
    run("00_setup.sh", dict(env, TRIO_CHILDREN="NA19240"))
    env1 = dict(env, MANAGER_LOCAL="1", FAKE_DISCOVERY_FAIL_CHILD="NA19240")
    result = run("01_stage_and_discover.sh", env1)
    assert "1 trios dispatched to discovery, 0 failed" in result.stdout
    jobs_log = (logs / "jobs.log").read_text()
    assert "kmer-discovery exited 3" in jobs_log
    assert not (work / "results" / "NA19240" / "NA19240.done").exists()
    assert sorted(os.listdir(work / "crams")) == [
        "NA19238.cram", "NA19238.cram.crai", "NA19239.cram", "NA19239.cram.crai",
        "NA19240.cram", "NA19240.cram.crai"]
    # the resweep reuses the verified CRAMs without downloading
    result = run("01_stage_and_discover.sh", dict(env, MANAGER_LOCAL="1"))
    stage_log = (work / "logs" / "stage_NA19240.log").read_text()
    assert stage_log.count("CRAM already present") == 3
    assert "downloading CRAM" not in stage_log
    assert (work / "results" / "NA19240" / "NA19240.done").exists()
    assert sorted(os.listdir(work / "crams")) == []


def test_flags_site_file_and_snapshot(workspace):
    """Paths come from flags, a site file or the environment, and jobs get
    them through the snapshot even when nothing is exported to them."""
    env, logs = workspace["env"], workspace["logs"]
    work2 = workspace["work"].parent / "work2"
    site = workspace["work"].parent / "site.sh"
    site.write_text(
        'N_TRIOS="${N_TRIOS:-1}"\n'                 # a flag below overrides it
        'TRIO_ALL_SIBLINGS="${TRIO_ALL_SIBLINGS:-1}"\n'
        'TRIO_SEED="${TRIO_SEED-}"\n'
        'EXPECTED_COMPLETE_TRIOS="${EXPECTED_COMPLETE_TRIOS:-3}"\n'
        'USE_SV_CALLSET="${USE_SV_CALLSET:-0}"\n'
        'MODULES_TO_LOAD="${MODULES_TO_LOAD:-}"\n')
    bare = {k: v for k, v in env.items() if k in ("PATH", "HOME")}

    # --help prints the header and the effective values, and runs nothing
    result = run("00_setup.sh", bare, "--help", "--work-dir", str(work2))
    assert "one-time setup" in result.stdout
    assert f"--work-dir                 {work2}" in result.stdout
    assert not work2.exists()

    # flags beat the site file, which beats the defaults
    run("00_setup.sh", bare, "--config", str(site), "--work-dir", str(work2),
        "--n-trios", "2", "--set", "DISCOVERY_MEM=4G")
    trios = (work2 / "trios.tsv").read_text().splitlines()
    assert [t.split("\t")[0] for t in trios[1:]] == ["HG00405", "HG00406"]
    snapshot = (work2 / "state" / "config.env").read_text()
    assert f"export WORK_DIR={work2}" in snapshot
    assert "export DISCOVERY_MEM=4G" in snapshot
    assert "export N_TRIOS=2" in snapshot
    assert "export USE_SV_CALLSET=0" in snapshot
    assert "00b_build_ref_index.sh --config " + str(work2 / "state" / "config.env") \
        in (logs / "sbatch.log").read_text()
    assert (work2 / "reference" / "GRCh38_full_analysis_set_plus_decoy_hla.fa.k31.jf").exists()

    # the site file through the environment, flags in --name=value form, and
    # a manager whose jobs see nothing but the snapshot
    env2 = dict(bare, KMER_TRIOS_CONFIG=str(site), MANAGER_LOCAL="1",
                TRIO_SLOTS="1", DISCOVERY_CPUS="2")
    result = run("01_stage_and_discover.sh", env2, f"--work-dir={work2}",
                 "--trio-limit=1", "--set=DISCOVERY_MEM=4G")
    assert "1 trios dispatched to discovery, 0 failed" in result.stdout
    submissions = [l for l in (logs / "sbatch.log").read_text().splitlines()
                   if "01b_discover_trio.sh" in l]
    assert submissions[-1].endswith(
        f"01b_discover_trio.sh --config {work2}/state/config.env HG00405")
    assert (work2 / "results" / "HG00405" / "HG00405.done").exists()
    args = (work2 / "results" / "HG00405" / "HG00405.args").read_text()
    assert f"--child {work2}/crams/HG00405.cram" in args
    assert "--threads 2 --memory 4" in args

    # a missing site file, and an unknown --flag, are errors not silent defaults
    result = run("02_summarize.sh", bare, "--config", str(work2 / "nope.sh"), check=False)
    assert result.returncode == 2 and "config file not found" in result.stderr
    result = run("02_summarize.sh", dict(bare, WORK_DIR=str(work2)), "--bogus", check=False)
    assert result.returncode != 0


def test_manager_refuses_without_the_reference_index(workspace):
    env, work = workspace["env"], workspace["work"]
    run("00_setup.sh", env)
    os.remove(work / "reference" / "GRCh38_full_analysis_set_plus_decoy_hla.fa.k31.jf")
    result = run("01_stage_and_discover.sh", dict(env, MANAGER_LOCAL="1"), check=False)
    assert result.returncode == 1
    assert "reference Jellyfish index missing" in result.stdout
