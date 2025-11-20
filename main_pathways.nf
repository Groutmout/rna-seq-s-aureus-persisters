#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.deseq2_results = "results/deseq2_results.csv"
params.gff_file       = "data/reference/reference.gff3"

// ----------------------------
// CHECK INPUT FILES
// ----------------------------

if( !file(params.deseq2_results).exists() ) {
    error "Fichier DESeq2 manquant : ${params.deseq2_results}"
}

if( !file(params.gff_file).exists() ) {
    error "Fichier GFF manquant : ${params.gff_file}"
}

// ----------------------------
// CHANNELS
// ----------------------------

ch_deseq = Channel.fromPath(params.deseq2_results)
ch_gff   = Channel.fromPath(params.gff_file)

// ----------------------------
// PROCESS ANNOTATE_GENES
// ----------------------------

process ANNOTATE_GENES {

    tag "annotate"
    container "bioconductor/bioconductor_docker:RELEASE_3_17"

    input:
        file deseq from ch_deseq
        file gff   from ch_gff

    output:
        file "deseq2_results_annotated.csv" into ch_annotated

    script:
    """
    Rscript scripts/annotate_genes.R \
        $deseq \
        $gff \
        deseq2_results_annotated.csv
    """
}

// ----------------------------
// PROCESS PATHWAYS
// ----------------------------

process PATHWAYS {

    tag "pathways"
    container "bioconductor/bioconductor_docker:RELEASE_3_17"

    input:
        file annotated from ch_annotated

    output:
        file "kegg_pathways_up.csv"
        file "kegg_pathways_down.csv"
        file "kegg_up_barplot.png"
        file "kegg_down_barplot.png"

    script:
    """
    Rscript scripts/run_pathways.R \
        $annotated \
        pathways_results

    cp pathways_results/* .
    """
}

// ----------------------------
// WORKFLOW
// ----------------------------

workflow {
    ANNOTATE_GENES()
    PATHWAYS()
}
