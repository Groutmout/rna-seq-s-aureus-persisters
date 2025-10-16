# =========================================================
#  Dockerfile for RNA-seq reproducibility hackathon
#  Common environment for all team members
# =========================================================
FROM ubuntu:22.04

LABEL maintainer="François Faramus <francois@exemple.com>"
LABEL description="Common Docker environment for RNA-seq workflow (Nextflow)"

# --------------------------
# Install system dependencies
# --------------------------
RUN apt-get update && apt-get install -y \
    wget curl git gzip tar bzip2 build-essential \
    python3 python3-pip \
    samtools bowtie subread fastqc \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# --------------------------
# Install BioConda tools via Miniconda
# --------------------------
RUN curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    /opt/conda/bin/conda install -y -c bioconda -c conda-forge \
        sra-tools trim-galore multiqc fastqc bowtie samtools subread \
        && /opt/conda/bin/conda clean -a -y

# Add Conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# --------------------------
# Create working directory
# --------------------------
WORKDIR /data

# --------------------------
# Default command
# --------------------------
CMD ["bash"]
