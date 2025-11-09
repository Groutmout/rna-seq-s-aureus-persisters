#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process hello {
    container 'ubuntu:22.04'
    output:
        path "hello.txt"
    script:
        """
        echo '✅ Nextflow fonctionne avec Docker local !' > hello.txt
        """
}

process bowtie_build {
    container 'bouty:latest'          // <-- ton image Docker locale
    publishDir "results/bowtie_index", mode: 'copy'

    input:
        path reference_fasta
    output:
        path "index.*"

    script:
        """
        echo "🔹 Construction de l'index Bowtie pour $reference_fasta"
        which bowtie-build
        bowtie-build $reference_fasta index
        """
}

workflow {
    hello()
    def reference_ch = Channel.fromPath('data/reference.fasta', checkIfExists: true)
    bowtie_build(reference_ch)
}
