#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================
// PARAMÈTRES
// =====================
params.outdir  = params.outdir  ?: "results"
params.sra_ids = params.sra_ids ?: ["SRR10379726","SRR10379725"]

println "Pipeline started"
println "Output directory: ${params.outdir}"

// =====================
// PROCESS 1 – Téléchargement FASTQ (SRA Toolkit)
// Sortie: tuple (id, R1.fastq.gz, R2.fastq.gz)
// =====================

process DOWNLOAD_FASTQ {
    container 'alantrbt/sratoolkit:latest'
    publishDir "${params.outdir}/raw_fastq", mode: 'copy'

    input:
        val id

    output:
        tuple val(id), path("${id}_1.fastq.gz"), path("${id}_2.fastq.gz"), optional: true

    script:
        """
        echo "=== Téléchargement du FASTQ pour $id ==="
        fasterq-dump --split-files --threads 4 $id -O .
        for f in ${id}_*.fastq; do gzip "\$f"; done
        """
}

// =====================
// PROCESS 2 – Trimming (Cutadapt)
// Entrée:  tuple (id, R1, R2)
// Sortie:  tuple (id, R1.trimmed, R2.trimmed)
// =====================
process TRIM {
    container 'alantrbt/cutadapt:1.11'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
        tuple val(id), path(read1), path(read2)

    output:
        tuple val(id), path("${id}_1.trimmed.fastq.gz"), path("${id}_2.trimmed.fastq.gz")

    script:
        """
        echo "=== Trimming pour $id ==="
        cutadapt -q 20 -m 20 -o ${id}_1.trimmed.fastq.gz ${read1}
        cutadapt -q 20 -m 20 -o ${id}_2.trimmed.fastq.gz ${read2}
        """
}

// =====================
// PROCESS 3 – Référence + Index Bowtie1
// Sorties nommées (emit):
//   index -> index.*  (pour ALIGN)
//   gtf   -> S_aureus.gtf (pour COUNT)
// =====================

process DOWNLOAD_REFERENCE {
    container 'alantrbt/bowtie:latest'
    publishDir "data/reference", mode: 'copy'

    output:
        path "*.ebwt",      emit: index
        path "S_aureus.gtf", emit: gtf

    script:
        """
        set -e
        wget -O S_aureus.fasta.gz \
          https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
        gunzip -f S_aureus.fasta.gz

        wget -O S_aureus.gtf.gz \
          https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz
        gunzip -f S_aureus.gtf.gz

        bowtie-build S_aureus.fasta index
        """
}

// =====================
// PROCESS 4 – Alignement (Bowtie1 + Samtools)
// Entrées : tuple(id,R1trim,R2trim) + path(index.*)
// Sortie  : tuple(id, BAM)
// =====================
process ALIGN {
    container 'alantrbt/bowtie:latest'
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
        tuple val(id), path(read1), path(read2)
        path index_files

    output:
        tuple val(id), path("${id}.bam")

    script:
        """
        echo "=== Alignement de $id ==="
        # Bowtie1 : l'index est référencé par son prefix ("index"), pas par chaque fichier
        bowtie -S -p 4 index -1 ${read1} -2 ${read2} \
          | samtools view -Sb - \
          | samtools sort -o ${id}.bam
        samtools index ${id}.bam
        """
}

// =====================
// PROCESS 5 – Comptage (featureCounts)
// Entrées : tuple(id,BAM) + path(GTF)
// Sortie  : path counts_<id>.txt
// =====================
process COUNT {
    container 'alantrbt/subread:latest'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
        tuple val(id), path(bam_file)
        path gtf_file

    output:
        path "counts_${id}.txt"

    script:
        """
        echo "=== Comptage pour $id ==="
        featureCounts -a ${gtf_file} -o counts_${id}.txt ${bam_file}
        """
}

// =====================
// PROCESS 6 – DESeq2
// Entrée : liste collectée des fichiers de counts
// Sortie : deseq2_results.csv
// =====================
process DESEQ2 {
    container 'alantrbt/deseq2:latest'
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
        path counts_files

    output:
        path "deseq2_results.csv"

    script:
        """
        echo "=== Analyse différentielle avec DESeq2 ==="
        Rscript /scripts/run_deseq2.R ${counts_files} deseq2_results.csv
        """
}

// =====================
// WORKFLOW PRINCIPAL
// =====================
workflow {
    Channel.fromList(params.sra_ids).set { sra_ch }

    reads_ch   = DOWNLOAD_FASTQ(sra_ch)
    trimmed_ch = TRIM(reads_ch)

    ref      = DOWNLOAD_REFERENCE()
    index_ch = ref.index                 // -> tous les *.ebwt
    gtf_one  = ref.gtf.first()           // -> 1 seul fichier, répliqué tout seul par Nextflow

    aligned_ch = ALIGN(trimmed_ch, index_ch)
    counts_ch  = COUNT(aligned_ch, gtf_one)  // on passe le gtf “single-value channel”

    counts_files = counts_ch.collect()
    DESEQ2(counts_files)
}
