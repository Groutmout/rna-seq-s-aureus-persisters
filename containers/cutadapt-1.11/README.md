# Cutadapt 1.11 (trimming)
Versions figées pour reproduire l’article :
- cutadapt **1.11** ; seules les séquences **≥ 25 nt** sont conservées.

Référence : « Reads were cleaned of adapter sequences and low-quality sequences using **cutadapt version 1.11**. **Only sequences at least 25 nt** were considered. » :contentReference[oaicite:0]{index=0}

Build:
  docker build -t persisters/cutadapt:1.11 containers/cutadapt-1.11

Run (exemple):
  docker run --rm -v "$PWD":/work persisters/cutadapt:1.11 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -o /work/out.trimmed.fastq.gz /work/in.fastq.gz
