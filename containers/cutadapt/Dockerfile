FROM mambaorg/micromamba:1.5.8

LABEL software="cutadapt"
LABEL version="1.11"

RUN micromamba install -y -n base -c bioconda -c conda-forge cutadapt=1.11 python=3.5 && \
    micromamba clean -a -y

WORKDIR /data

ENTRYPOINT ["micromamba", "run", "-n", "base", "cutadapt"]