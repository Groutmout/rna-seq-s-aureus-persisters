FROM mambaorg/micromamba:1.5.8
ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=mambauser:mambauser environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y

USER root
COPY --chmod=0755 entrypoint.sh /usr/local/bin/entrypoint.sh
USER mambauser

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

