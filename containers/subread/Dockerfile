FROM ubuntu:24.04

LABEL software="subread"
LABEL version="1.4.6-p3"
LABEL description="Subread 1.4.6-p3 (featureCounts) for paper reproduction"

RUN apt-get update && apt-get install -y \
    wget make gzip gcc zlib1g-dev libpthread-stubs0-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local

RUN wget -O subread.tar.gz https://sourceforge.net/projects/subread/files/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz/download && \
    tar -xzf subread.tar.gz && \
    rm subread.tar.gz

ENV PATH="/usr/local/subread-1.4.6-p3-Linux-x86_64/bin:${PATH}"

WORKDIR /data
