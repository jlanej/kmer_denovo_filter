FROM python:3.12-slim

ARG KRAKEN2_VERSION=v2.1.6

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        perl \
        zlib1g-dev \
        samtools \
        jellyfish \
        wget \
    && git clone --depth 1 --branch "${KRAKEN2_VERSION}" \
        https://github.com/DerrickWood/kraken2.git /tmp/kraken2 \
    && /tmp/kraken2/install_kraken2.sh /usr/local/bin \
    && rm -rf /tmp/kraken2 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml README.md ./
COPY src/ src/
COPY scripts/ scripts/

RUN pip install --no-cache-dir .

ENTRYPOINT ["kmer-denovo"]
