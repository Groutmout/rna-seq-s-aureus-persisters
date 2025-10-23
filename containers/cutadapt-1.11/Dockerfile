# =========================================================
# Cutadapt (v1.11, paper-faithful) in a reproducible container
# =========================================================
FROM mambaorg/micromamba:1.5.8

LABEL maintainer="Ton Nom <ton.email@exemple.com>"
LABEL description="Cutadapt 1.11 (trimming) with min length 25, as per paper"
LABEL org.opencontainers.image.source="https://github.com/Groutmout/rna-seq-s-aureus-persisters"

# Active l'env base automatiquement dans les RUN
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Dépendances conda (versions figées dans environment.yml)
COPY --chown=mambauser:mambauser environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y

# EntryPoint (copié en root, exécuté en user non-root)
USER root
COPY --chmod=0755 entrypoint.sh /usr/local/bin/entrypoint.sh
USER mambauser

# Répertoire de travail par défaut
WORKDIR /data

# Entrée : wrapper qui appelle "micromamba run -n base cutadapt"
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

