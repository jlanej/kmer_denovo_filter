"""Tests for scripts/download_kraken2_db.sh behavior."""

import os
import stat
import subprocess
from pathlib import Path

PRACKENDB_URL = (
    "https://genome-idx.s3.amazonaws.com/kraken/"
    "k2_NCBI_reference_20251007.tar.gz"
)


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_download_kraken2_db_uses_wget(tmp_path):
    """Script uses wget to download the PrackenDB tarball."""
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "scripts" / "download_kraken2_db.sh"

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    db_path = tmp_path / "db"
    wget_args_log = tmp_path / "wget-args.log"
    # Fake wget: log args, create DB files to satisfy validation.
    _write_executable(
        fake_bin / "wget",
        f"""#!/usr/bin/env bash
set -euo pipefail
echo "$*" >> "{wget_args_log}"
# Parse -O <dest> to find db directory.
dest=""
for arg in "$@"; do
  if [[ -n "${{prev:-}}" && "$prev" == "-O" ]]; then
    dest="$arg"
  fi
  prev="$arg"
done
db_dir="$(dirname "$dest")"
mkdir -p "$db_dir/taxonomy"
touch "$db_dir/hash.k2d" "$db_dir/opts.k2d" "$db_dir/taxo.k2d" "$db_dir/taxonomy/nodes.dmp"
# Create a dummy tarball so tar succeeds.
tar -czf "$dest" -C "$db_dir" hash.k2d opts.k2d taxo.k2d taxonomy
""",
    )

    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}:{env['PATH']}"

    result = subprocess.run(
        [str(script_path), "--db", str(db_path)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    args = wget_args_log.read_text(encoding="utf-8").strip()
    assert PRACKENDB_URL in args
    assert "PrackenDB" in result.stdout


def test_download_kraken2_db_custom_url(tmp_path):
    """--url overrides the default PrackenDB URL."""
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "scripts" / "download_kraken2_db.sh"

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    db_path = tmp_path / "db"
    wget_args_log = tmp_path / "wget-args.log"
    custom_url = "https://example.com/custom_kraken2.tar.gz"
    _write_executable(
        fake_bin / "wget",
        f"""#!/usr/bin/env bash
set -euo pipefail
echo "$*" >> "{wget_args_log}"
dest=""
for arg in "$@"; do
  if [[ -n "${{prev:-}}" && "$prev" == "-O" ]]; then
    dest="$arg"
  fi
  prev="$arg"
done
db_dir="$(dirname "$dest")"
mkdir -p "$db_dir/taxonomy"
touch "$db_dir/hash.k2d" "$db_dir/opts.k2d" "$db_dir/taxo.k2d" "$db_dir/taxonomy/nodes.dmp"
tar -czf "$dest" -C "$db_dir" hash.k2d opts.k2d taxo.k2d taxonomy
""",
    )

    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}:{env['PATH']}"

    result = subprocess.run(
        [str(script_path), "--db", str(db_path), "--url", custom_url],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    args = wget_args_log.read_text(encoding="utf-8").strip()
    assert custom_url in args
    assert PRACKENDB_URL not in args


def test_download_kraken2_db_fails_without_wget(tmp_path):
    """Script exits with an error when wget is not available."""
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "scripts" / "download_kraken2_db.sh"

    import shutil

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    db_path = tmp_path / "db"

    # Symlink only bash so the script can run but wget is absent.
    bash_path = shutil.which("bash")
    if bash_path:
        (fake_bin / "bash").symlink_to(bash_path)

    env = os.environ.copy()
    env["PATH"] = str(fake_bin)

    result = subprocess.run(
        [str(script_path), "--db", str(db_path)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode != 0
    assert "wget not found" in result.stderr


def test_download_kraken2_db_fails_on_missing_db_files(tmp_path):
    """Script fails when required database files are missing after extract."""
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "scripts" / "download_kraken2_db.sh"

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    db_path = tmp_path / "db"
    # Fake wget that creates a tarball with only an empty file.
    _write_executable(
        fake_bin / "wget",
        f"""#!/usr/bin/env bash
set -euo pipefail
dest=""
for arg in "$@"; do
  if [[ -n "${{prev:-}}" && "$prev" == "-O" ]]; then
    dest="$arg"
  fi
  prev="$arg"
done
db_dir="$(dirname "$dest")"
mkdir -p "$db_dir"
touch "$db_dir/dummy.txt"
tar -czf "$dest" -C "$db_dir" dummy.txt
""",
    )

    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}:{env['PATH']}"

    result = subprocess.run(
        [str(script_path), "--db", str(db_path)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode != 0
    assert "missing required database file" in result.stderr
