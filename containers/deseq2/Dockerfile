FROM mambaorg/micromamba:1.5.8

LABEL software="deseq2"
LABEL version="1.16"
LABEL description="DESeq2 1.16 environment (R 3.4.1)"

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        r-base=3.4.1 \
        bioconductor-deseq2=1.16.* \
        r-ggplot2=3.0.* \
        libgfortran=3.0.0 \
        libgcc-ng=7.* \
        libstdcxx-ng=7.* \
    && micromamba clean -a -y

SHELL ["/bin/bash", "-lc"]

USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["micromamba", "run", "-n", "base"]

WORKDIR /data

