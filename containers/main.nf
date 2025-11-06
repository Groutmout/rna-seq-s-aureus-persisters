#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// =====================
// PARAMÈTRES GLOBAUX
// =====================
params.reference_index = params.reference_index ?: "data/reference/S_aureus"
params.annotation_gtf  = params.annotation_gtf  ?: "data/reference/S_aureus.gtf"
params.outdir          = params.outdir ?: "results"

// Debug info
println "Pipeline started"
println "Output directory: ${params.outdir}"

// =====================
// PROCESS 1: Téléchargement des données SRA
// =====================
process downloadFASTQ {
    container 'alantrbt/sratoolkit:latest'
    publishDir "${params.outdir}/raw_fastq", mode: 'copy'

    output:
    path "*.fastq.gz"

    script:
    """
    echo "=== Téléchargement des fichiers FASTQ depuis SRA ==="
    for id in SRR10379726 SRR10379725 SRR10379724 SRR10379723 SRR10379722 SRR10379721; do
        echo "Téléchargement de \$id..."
        fasterq-dump --split-files --threads 4 \$id -O .
        gzip *.fastq
    done
    echo "=== Téléchargement terminé ==="
    """
}

// =====================
// PROCESS 2: Trimming (Cutadapt)
// =====================
process TRIM {
    container 'alantrbt/cutadapt:1.11'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path reads

    output:
    path "*.trimmed.fastq.gz"

    script:
    """
    for f in ${reads}; do
        base=\$(basename \$f .fastq.gz)
        cutadapt -q 20 -m 20 -o \${base}.trimmed.fastq.gz \$f
    done
    """
}

// =====================
// PROCESS 3: Téléchargement du génome + index (Bowtie2)
// =====================
process DOWNLOAD_REFERENCE {
    container 'alantrbt/bowtie2:latest'
    publishDir "data/reference", mode: 'copy'

    output:
    path "S_aureus.fasta"
    path "S_aureus.gtf"
    path "S_aureus.*.bt2"

    script:
    """
    mkdir -p data/reference
    cd data/reference

    # Télécharger le génome et l’annotation
    wget -O S_aureus.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
    wget -O S_aureus.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz

    gunzip -f S_aureus.fasta.gz
    gunzip -f S_aureus.gtf.gz

    # Construire l’index Bowtie2
    bowtie2-build S_aureus.fasta S_aureus
    """
}


// =====================
// PROCESS 4: Alignement (Bowtie2 + Samtools)
// =====================
process ALIGN {
    container 'alantrbt/bowtie2:latest'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    path reads
    path ref

    output:
    path "*.bam"

    script:
    """
    for f in ${reads}; do
        base=\$(basename \$f .trimmed.fastq.gz)
        bowtie2 -x S_aureus -U \$f | samtools sort -o \${base}.bam
    done
    """
}

// =====================
// PROCESS 5: Comptage (FeatureCounts)
// =====================
process COUNT {
    container 'alantrbt/subread:latest'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path bam_files

    output:
    path "counts.txt"

    script:
    """
    featureCounts -a ${params.annotation_gtf} -o counts.txt ${bam_files}
    """
}

// =====================
// PROCESS 6: Analyse différentielle (DESeq2)
// =====================
process DESEQ2 {
    container 'alantrbt/deseq2:latest'
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path counts

    output:
    path "deseq2_results.csv"

    script:
    """
    echo "=== Analyse différentielle avec DESeq2 ==="
    Rscript /scripts/run_deseq2.R counts.txt deseq2_results.csv
    echo "=== Analyse DESeq2 terminée ==="
    """
}



// =====================
// WORKFLOW PRINCIPAL
// =====================
workflow {
    reads = downloadFASTQ()
    trimmed = TRIM(reads.out)
    ref = DOWNLOAD_REFERENCE()
    aligned = ALIGN(trimmed.out, ref.out)
    counts = COUNT(aligned.out)
    DESEQ2(counts.out)
}
