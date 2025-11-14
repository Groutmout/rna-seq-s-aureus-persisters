FROM mambaorg/micromamba:1.5.8

USER root
ENV DEBIAN_FRONTEND=noninteractive
ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN apt-get update && apt-get install -y procps && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /tmp/mamba && chmod -R 777 /tmp/mamba
ENV MAMBA_ROOT_PREFIX=/tmp/mamba
ENV HOME=/tmp

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        sra-tools=3.0.* \
    && micromamba clean -a -y

SHELL ["/bin/bash", "-lc"]

WORKDIR /data

ENV PATH="/tmp/mamba/bin:${PATH}"
ENTRYPOINT ["fasterq-dump"]
