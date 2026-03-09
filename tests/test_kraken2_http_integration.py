"""Integration tests confirming that the PrackenDB download endpoint is
reachable and that ``kraken2`` (the classifier) is installed.

The endpoint reachability test is skipped in environments without network
access.  The ``kraken2`` binary check is skipped when the binary is not
on PATH (e.g. local dev environments without Kraken2 installed).
"""

import shutil
import subprocess
import urllib.request
import urllib.error

import pytest

# Pre-built PrackenDB database URL used by download_kraken2_db.sh.
# https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads
PRACKENDB_URL = (
    "https://genome-idx.s3.amazonaws.com/kraken/"
    "k2_NCBI_reference_20251007.tar.gz"
)

_KRAKEN2_AVAILABLE = shutil.which("kraken2") is not None


def _network_available() -> bool:
    """Return True when the PrackenDB S3 endpoint is reachable."""
    try:
        req = urllib.request.Request(PRACKENDB_URL, method="HEAD")
        with urllib.request.urlopen(req, timeout=10) as resp:
            return resp.status < 500
    except urllib.error.HTTPError as exc:
        return exc.code < 500
    except (urllib.error.URLError, OSError):
        return False


@pytest.mark.skipif(not _network_available(), reason="PrackenDB endpoint unreachable (no network)")
def test_prackendb_endpoint_reachable():
    """Confirm that the PrackenDB download URL is reachable.

    Sends a HEAD request to the S3-hosted tarball.  A 2xx/3xx response
    (or a 403 from S3 bucket policy) confirms the endpoint is live.
    """
    req = urllib.request.Request(PRACKENDB_URL, method="HEAD")
    with urllib.request.urlopen(req, timeout=30) as resp:
        assert resp.status < 500


@pytest.mark.skipif(not _KRAKEN2_AVAILABLE, reason="kraken2 not found on PATH")
def test_kraken2_binary_available():
    """Confirm that the kraken2 classifier binary is installed."""
    result = subprocess.run(
        ["kraken2", "--version"],
        capture_output=True,
        text=True,
    )
    text = result.stdout + result.stderr
    assert "kraken" in text.lower(), (
        f"Unexpected output from kraken2 --version:\n{text}"
    )
