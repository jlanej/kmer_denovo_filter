"""Tests for scripts/download_kraken2_db.sh behavior."""

import os
import stat
import subprocess
from pathlib import Path


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_download_kraken2_db_uses_ftp_mode(tmp_path):
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "scripts" / "download_kraken2_db.sh"

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    db_path = tmp_path / "db"
    args_log = tmp_path / "kraken2-build-args.log"
    _write_executable(
        fake_bin / "kraken2-build",
        f"""#!/usr/bin/env bash
set -euo pipefail
echo "$*" >> "{args_log}"
db=""
while [[ $# -gt 0 ]]; do
  if [[ "$1" == "--db" ]]; then
    db="$2"
    shift 2
  else
    shift
  fi
done
mkdir -p "$db/taxonomy"
touch "$db/hash.k2d" "$db/opts.k2d" "$db/taxo.k2d" "$db/taxonomy/nodes.dmp"
""",
    )

    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}:{env['PATH']}"

    result = subprocess.run(
        [str(script_path), "--db", str(db_path), "--threads", "2"],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    args = args_log.read_text(encoding="utf-8").splitlines()
    assert len(args) == 1
    assert "--use-ftp" in args[0]


def test_download_kraken2_db_does_not_retry_when_pub_module_error_is_logged(tmp_path):
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "scripts" / "download_kraken2_db.sh"

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    db_path = tmp_path / "db"
    args_log = tmp_path / "kraken2-build-args.log"
    _write_executable(
        fake_bin / "kraken2-build",
        f"""#!/usr/bin/env bash
set -euo pipefail
echo "$*" >> "{args_log}"
echo "Downloading nucleotide gb accession to taxon map...rsync: @ERROR: Unknown module 'pub'" >&2
echo "rsync error: error starting client-server protocol (code 5) at main.c(1850) [Receiver=3.4.1]" >&2
db=""
while [[ $# -gt 0 ]]; do
  if [[ "$1" == "--db" ]]; then
    db="$2"
    shift 2
  else
    shift
  fi
done
mkdir -p "$db/taxonomy"
touch "$db/hash.k2d" "$db/opts.k2d" "$db/taxo.k2d" "$db/taxonomy/nodes.dmp"
""",
    )

    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}:{env['PATH']}"

    result = subprocess.run(
        [str(script_path), "--db", str(db_path), "--threads", "2"],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    args = args_log.read_text(encoding="utf-8").splitlines()
    assert len(args) == 1
    assert "--use-ftp" in args[0]
