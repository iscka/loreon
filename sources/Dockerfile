FROM ubuntu:22.04 AS tools

ENV DEBIAN_FRONTEND=noninteractive

ARG MINIMAP2_VERSION=2.28
ARG MINIMAP2_SHA256=ffa5712735d229119f8c05722a0638ae0cc15aeb8938e29a3e52d5da5c92a0b4
ARG SAMTOOLS_VERSION=1.21
ARG SAMTOOLS_SHA256=05724b083a6b6f0305fcae5243a056cc36cf826309c3cb9347a6b89ee3fc5ada

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    ca-certificates \
    curl \
    libbz2-dev \
    libdeflate-dev \
    liblzma-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

RUN curl -fsSL -o samtools.tar.bz2 \
        "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    && echo "${SAMTOOLS_SHA256}  samtools.tar.bz2" | sha256sum -c - \
    && tar -xjf samtools.tar.bz2 \
    && cd "samtools-${SAMTOOLS_VERSION}" \
    && ./configure --prefix=/opt/tools --without-curses --disable-libcurl \
    && make -j"$(nproc)" \
    && make install

RUN curl -fsSL -o minimap2.tar.bz2 \
        "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2" \
    && echo "${MINIMAP2_SHA256}  minimap2.tar.bz2" | sha256sum -c - \
    && tar -xjf minimap2.tar.bz2 \
    && cd "minimap2-${MINIMAP2_VERSION}" \
    && if [ "$(uname -m)" = "aarch64" ]; then \
           make -j"$(nproc)" arm_neon=1 aarch64=1; \
       else \
           make -j"$(nproc)"; \
       fi \
    && install -m 0755 minimap2 /opt/tools/bin/minimap2


FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    gawk \
    libdeflate0 \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY --from=tools /opt/tools/bin/minimap2 /opt/tools/bin/samtools /usr/local/bin/

WORKDIR /app

COPY requirements.txt .
RUN grep -v -i '^pyqt5' requirements.txt > requirements-pipeline.txt \
    && pip3 install --no-cache-dir -r requirements-pipeline.txt

COPY metaGenomics_new.py \
     OtuUtils.py \
     ResultsReader.py \
     PipelineLogger.py \
     report_generator.py \
     pipeline_profiler.py \
     template.html \
     ./
COPY example_data ./example_data

RUN mkdir -p /data/input /data/output /data/db

CMD ["python3", "/app/metaGenomics_new.py", "--help"]
