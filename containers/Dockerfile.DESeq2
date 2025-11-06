FROM mambaorg/micromamba:1.5.8 

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    r-base=3.4.1 \
    bioconductor-deseq2=1.16.* \
    r-ggplot2=3.0.* \
    libgfortran=3.0.0 \
    libgcc-ng=7.* \
    libstdcxx-ng=7.* \
 && micromamba clean -a -y


SHELL ["/bin/bash", "-lc"]
