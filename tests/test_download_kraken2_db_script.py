"""Tests for scripts/download_kraken2_db.sh behavior."""

import os
import stat
import subprocess
from pathlib import Path


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_download_kraken2_db_prefers_k2(tmp_path):
    """When k2 is on PATH the script uses 'k2 build' (HTTP mode)."""
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "scripts" / "download_kraken2_db.sh"

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    db_path = tmp_path / "db"
    k2_args_log = tmp_path / "k2-args.log"
    _write_executable(
        fake_bin / "k2",
        f"""#!/usr/bin/env bash
set -euo pipefail
echo "$*" >> "{k2_args_log}"
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
    # Also place a kraken2-build shim to verify it is NOT called.
    kb_args_log = tmp_path / "kraken2-build-args.log"
    _write_executable(
        fake_bin / "kraken2-build",
        f"""#!/usr/bin/env bash
echo "$*" >> "{kb_args_log}"
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
    # k2 should have been invoked
    args = k2_args_log.read_text(encoding="utf-8").splitlines()
    assert len(args) == 1
    assert "--standard" in args[0]
    assert str(db_path) in args[0]
    # kraken2-build should not have been invoked
    assert not kb_args_log.exists()


def test_download_kraken2_db_fails_without_k2(tmp_path):
    """Script exits with an error when k2 is not available."""
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "scripts" / "download_kraken2_db.sh"

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    db_path = tmp_path / "db"

    env = os.environ.copy()
    # Keep /usr/bin (for bash, env, etc.) but remove everything else so
    # k2 is not found.
    env["PATH"] = f"{fake_bin}:/usr/bin:/bin"

    result = subprocess.run(
        [str(script_path), "--db", str(db_path)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode != 0
    assert "k2 not found on PATH" in result.stderr
    assert "Kraken2 >= 2.17" in result.stderr
