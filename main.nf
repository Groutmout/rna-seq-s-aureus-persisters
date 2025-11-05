#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process hello {
    container 'ubuntu:22.04'
    output:
        path "hello.txt"
    script:
        """
        echo 'Youpiyoup Nextflow fonctionne !' > hello.txt
        """
}

workflow {
    hello()
}
