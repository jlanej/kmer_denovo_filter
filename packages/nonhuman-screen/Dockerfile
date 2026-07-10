FROM python:3.12-slim

# kraken2 version this package is developed and tested against.
ARG KRAKEN2_VERSION=v2.17.1

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        perl \
        zlib1g-dev \
        wget \
    && git clone --depth 1 --branch "${KRAKEN2_VERSION}" \
        https://github.com/DerrickWood/kraken2.git /tmp/kraken2 \
    && /tmp/kraken2/install_kraken2.sh /usr/local/bin \
    && rm -rf /tmp/kraken2 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml README.md LICENSE ./
COPY src/ src/
COPY scripts/ scripts/

# Install with the [bam] extra so the allele-based NHF helpers are available.
RUN pip install --no-cache-dir '.[bam]'

ENTRYPOINT ["nonhuman-screen"]
