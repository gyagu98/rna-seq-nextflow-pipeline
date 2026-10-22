# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.10.0-1

LABEL maintainer="Gretta Y. Rozenberg <https://github.com/gyagu98>"
LABEL description="All-in-one container for rna-seq-nextflow-pipeline"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    curl \
    default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools via conda / bioconda
COPY environment.yml /tmp/environment.yml
RUN conda env update -n base -f /tmp/environment.yml && conda clean -afy

WORKDIR /pipeline
ENTRYPOINT ["/bin/bash"]
