FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    minimap2 \
    samtools \
    gawk \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

COPY metaGenomics_new.py \
     OtuUtils.py \
     ResultsReader.py \
     PipelineLogger.py \
     report_generator.py \
     pipeline_profiler.py \
     template.html \
     ./

RUN mkdir -p /data/input /data/output /data/db

CMD ["python3", "/app/metaGenomics_new.py", "--help"]
