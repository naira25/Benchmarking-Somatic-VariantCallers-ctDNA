#!/bin/bash

/*
 * Fastq Pre Processing Subworkflow
*/

nextflow.enable.dsl=2

workflow {
    nextflow run nf-core/sarek -r 3.8.1 <VERSION> \ 
        --profile <docker/singulatiry/.../institute> \
        --input samplesheet.csv
        --step "mapping" # Starting step of the workflow, can be "mapping"
        --outdir
}


# Version is laterst version
# Started from FastQ files gzip.compressed, paired-end sequencing data of ctDNA samples B, D and E
# No GPU accelerated alingment
# No initial duplicate marking