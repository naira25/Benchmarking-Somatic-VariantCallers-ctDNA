#!/bin/bash

BASE_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA"
OUTPUT_DIR="${BASE_DIR}/Results/Test_Results"

nextflow run nf-core/sarek -r 3.8.1 \
    -c ${BASE_DIR}/Conf/nextflow.config \
    -profile docker \
    --input ${BASE_DIR}/Metadata/samplesheet.csv \
    --outdir ${OUTPUT_DIR} \
    --genome GATK.GRCh38 \
    --wes true \
    --intervals ${BASE_DIR}/Metadata/ROC2.bed \
    --save_output_as_bam true
    -resume

# Version is laterst version
# Started from FastQ files gzip.compressed, paired-end sequencing data of ctDNA samples B, D and E
# No GPU accelerated alingment
# No initial duplicate marking