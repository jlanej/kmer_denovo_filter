"""Integration tests confirming k2 (Kraken2 >= 2.17) is installed and that
the HTTP download endpoint used by ``k2 build`` is reachable.

These tests are skipped when ``k2`` is not found on PATH (e.g. in local dev
environments without Kraken2 installed).  They run automatically in CI where
the container image is built with the required Kraken2 version.
"""

import re
import shutil
import subprocess
import urllib.request
import urllib.error

import pytest

# URL used by k2 to list available pre-built standard databases.
# Verifying this endpoint is reachable confirms HTTP downloads will work.
K2_INDEX_URL = "https://genome-idx.s3.amazonaws.com/kraken/"

# Minimum Kraken2 version that ships the k2 wrapper with HTTP download support.
MIN_VERSION = (2, 17)

_K2_AVAILABLE = shutil.which("k2") is not None


@pytest.fixture(scope="module")
def k2_version() -> tuple[int, ...]:
    """Return the installed k2/Kraken2 version as a tuple of ints."""
    result = subprocess.run(
        ["k2", "--version"],
        capture_output=True,
        text=True,
    )
    # k2 / kraken2 typically prints something like "Kraken version 2.1.6" or
    # "k2 2.17.0".  Extract the first dotted-number sequence found.
    text = result.stdout + result.stderr
    match = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text)
    if not match:
        pytest.fail(f"Could not parse k2 --version output:\n{text}")
    parts = tuple(int(g) for g in match.groups() if g is not None)
    return parts


@pytest.mark.skipif(not _K2_AVAILABLE, reason="k2 not found on PATH")
def test_k2_version_meets_minimum(k2_version):
    """Confirm that the installed k2 version is >= 2.17 (HTTP-download capable)."""
    major_minor = k2_version[:2]
    assert major_minor >= MIN_VERSION, (
        f"k2 version {'.'.join(str(v) for v in k2_version)} is below the "
        f"minimum required version {'.'.join(str(v) for v in MIN_VERSION)}. "
        "Please upgrade Kraken2: https://github.com/DerrickWood/kraken2"
    )


@pytest.mark.skipif(not _K2_AVAILABLE, reason="k2 not found on PATH")
def test_k2_http_download_endpoint_reachable():
    """Confirm that the k2 HTTP download index is reachable.

    This test sends a HEAD request to the S3 bucket index used by ``k2 build``
    to download pre-built standard databases.  A successful response (any 2xx
    or 3xx status) confirms that HTTP downloads are functional in this
    environment.
    """
    try:
        req = urllib.request.Request(K2_INDEX_URL, method="HEAD")
        with urllib.request.urlopen(req, timeout=30) as resp:
            status = resp.status
    except urllib.error.HTTPError as exc:
        # HTTPError still has a status code; anything below 500 is acceptable
        # (e.g. 403 Forbidden from S3 for a bucket listing is still reachable).
        status = exc.code
        if status >= 500:
            pytest.fail(
                f"k2 HTTP download endpoint returned server error {status}: "
                f"{K2_INDEX_URL}"
            )
    except urllib.error.URLError as exc:
        pytest.fail(
            f"k2 HTTP download endpoint is unreachable: {K2_INDEX_URL}\n"
            f"Reason: {exc.reason}"
        )
