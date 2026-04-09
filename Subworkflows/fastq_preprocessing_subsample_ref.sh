#!/bin/bash

BASE_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA"
OUTPUT_DIR="${BASE_DIR}/Results/Subsample_1"

nextflow run nf-core/sarek -r 3.8.1 \
    -c ${BASE_DIR}/Conf/nextflow.config \
    -profile docker \
    -resume \
    --input ${BASE_DIR}/Metadata/SampleSheets_Subsamples/samplesheet_1.csv \
    --outdir ${OUTPUT_DIR} \
    --genome hg38 \
    --igenomes_ignore true \
    --fasta https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
    --wes true \
    --intervals ${BASE_DIR}/Metadata/ROC2.bed \
    --save_output_as_bam true \
    --save_reference true \
    --split_fastq 0 \
    --skip_tools baserecalibrator \
    --resource_limits true