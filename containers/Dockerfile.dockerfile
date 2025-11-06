FROM ubuntu:24.04

RUN apt-get update && apt-get install -y wget make gzip gcc zlib1g zlib1g-dev libpthread-stubs0-dev

WORKDIR /usr/local
RUN wget https://sourceforge.net/projects/subread/files/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz/download
RUN tar -zxvf download

ENTRYPOINT ["/usr/local/subread-1.4.6-p3-Linux-x86_64/bin"]
