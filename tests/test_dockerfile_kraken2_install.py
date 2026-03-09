from pathlib import Path


def test_dockerfile_builds_kraken2_from_source() -> None:
    dockerfile = Path(__file__).resolve().parent.parent / "Dockerfile"
    content = dockerfile.read_text(encoding="utf-8")

    assert "\n        kraken2 \\" not in content
    assert "install_kraken2.sh /usr/local/bin" in content
    assert "https://github.com/DerrickWood/kraken2.git" in content
