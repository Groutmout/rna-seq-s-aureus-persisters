FROM ubuntu:16.04

LABEL software="bowtie"
LABEL version="0.12.7"
LABEL description="Bowtie 0.12.7 + Samtools 0.1.19 for paper reproduction"

RUN apt-get update && apt-get install -y \
    build-essential zlib1g-dev wget unzip python \
    libncurses5-dev libbz2-dev liblzma-dev git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# Bowtie
RUN wget -q https://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.7/bowtie-0.12.7-src.zip && \
    unzip bowtie-0.12.7-src.zip && \
    cd bowtie-0.12.7 && \
    sed -i '1i #include <algorithm>' alphabet.h && \
    sed -i 's/-O3/-O3 -fpermissive/' Makefile && \
    make && \
    ln -s /opt/bowtie-0.12.7/bowtie* /usr/local/bin/

# Samtools
RUN wget -q https://github.com/samtools/samtools/archive/refs/tags/0.1.19.tar.gz -O samtools.tar.gz && \
    tar -xzf samtools.tar.gz && \
    cd samtools-0.1.19 && \
    sed -i 's/-lcurses/-lncurses/' Makefile && \
    make CFLAGS="-O2 -fcommon -D_CURSES_LIB=1" && \
    cp samtools /usr/local/bin/ && \
    cp bcftools/bcftools /usr/local/bin/ && \
    cp misc/md5sum-lite /usr/local/bin/

WORKDIR /data

USER root

CMD ["bash"]
