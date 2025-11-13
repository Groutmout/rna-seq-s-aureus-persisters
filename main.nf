#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//
// Paramètres
//
params.outdir  = params.outdir  ?: "results"


/*
 * 1) Téléchargement des lectures
 */
process DOWNLOAD_FASTQ {
    tag "$id"
    publishDir "${params.outdir}/raw_fastq", mode: 'copy'

    input:
        val id

    output:
        tuple val(id), path("${id}.fastq.gz")

    script:
        """
        fasterq-dump ${id} --threads 4 -O .
        gzip ${id}.fastq
        """

}

/*
 * 2) Trimming (Cutadapt)
 */
process TRIM {
    tag "$id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
        tuple val(id), path(read)

    output:
        tuple val(id), path("${id}.trimmed.fastq.gz")

    script:
        """
        cutadapt -q 20 --phred33 --length 25 \
            -o ${id}.trimmed.fastq.gz \
            ${read}
        """
}

/*
 * 3) Téléchargement génome
 */
process DOWNLOAD_REFERENCE {
    publishDir "data/reference", mode: 'copy'

    output:
        path "reference.fasta", emit: fasta
        path "reference.gff3",  emit: gff

    script:
        """
        wget -q -O reference.fasta \
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"

        wget -q -O reference.gff3 \
          "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
        """
}

/*
 * 4) Index Bowtie1
 */
process INDEX {
    tag "bowtie-index"
    publishDir "${params.outdir}/index", mode: 'copy'

    input:
        path fasta

    output:
        path "index.*", emit: index

    script:
        """
        bowtie-build ${fasta} index
        """
}

/*
 * 5) Alignement
 */
process ALIGN {
    tag "$id"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
        tuple val(id), path(read)
        path index_files

    output:
        tuple val(id), path("${id}.bam")

    script:
        """
        bowtie -S -p 4 index ${read} \
            | samtools sort -o ${id}.bam

        samtools index ${id}.bam
        """
}

/*
 * 6) Comptage (featureCounts, GFF3)
 */
process COUNT {
    tag "$id"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
        tuple val(id), path(bam)
        path gff

    output:
        path "counts_${id}.txt"

    script:
        """
        featureCounts \
            -F GFF \
            -t gene \
            -g ID \
            -T 4 \
            -a ${gff} \
            -o counts_${id}.txt \
            ${bam}
        """
}

/*
 * 7) DESeq2
 */
process DESEQ2 {
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
        path counts_files
        val samples_metadata

    output:
        path "deseq2_results.csv"

    script:
        """
        Rscript /scripts/run_deseq2.R ${counts_files.join(" ")} deseq2_results.csv
        """
}

/*
 * Workflow principal
 */
workflow {

    /*
     * Charger le sample sheet dans le workflow (obligatoire en DSL2)
     */
    def sample_list = file("$baseDir/samples.tsv").splitCsv(header:true, sep:'\t')

    samples_ch = Channel.from(sample_list)
    // samples_ch.view()

    sra_ids_ch = samples_ch.map { row -> row.sra }
    // sra_ids_ch.view()


    reads_ch = DOWNLOAD_FASTQ(sra_ids_ch)
    trimmed_ch = TRIM(reads_ch)

    ref        = DOWNLOAD_REFERENCE()
    fasta_ch   = ref.fasta
    gff_ch     = ref.gff

    index_ch   = INDEX(fasta_ch)
    aligned_ch = ALIGN(trimmed_ch, index_ch)
    counts_ch  = COUNT(aligned_ch, gff_ch)

    all_counts = counts_ch.collect()

    samples_ch.collect().set { samples_metadata }
    DESEQ2(all_counts, samples_metadata)
}

