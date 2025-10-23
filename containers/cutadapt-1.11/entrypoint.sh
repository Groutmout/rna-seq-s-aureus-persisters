#!/usr/bin/env bash
set -euo pipefail

ARGS=("$@")

# Usage minimal:
#   docker run --rm -v $PWD:/data image \
#     -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     -o out.fastq.gz in.fastq.gz
#
# Tout ce que vous passez ici est répercuté à cutadapt,
# avec -m 25 ajouté systématiquement pour coller au papier.

# Ajouter -m 25 si absent
if ! printf '%s\0' "${ARGS[@]}" | grep -q -- '-m\0\|--minimum-length\0' ; then
  ARGS+=("-m" "25")
fi

# Log des versions pour la traçabilité
echo "== Cutadapt version =="
cutadapt --version
echo "======================"

exec micromamba run -n base cutadapt "${ARGS[@]}"
