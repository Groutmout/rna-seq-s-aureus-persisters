process downloadFASTQ {

    output:
    path "*.fastq.gz" 

    container = 'alantrbt/sratoolkit:latest'

    script:
    """
    fasterq-dump SRR10379726
    fasterq-dump SRR10379725
    fasterq-dump SRR10379724
    fasterq-dump SRR10379723
    fasterq-dump SRR10379722
    fasterq-dump SRR10379721
    """
}

workflow {
    fastq_files = downloadFASTQ()
}